"""Annotate the top-K candidates from a scored JSONL with regulatory context.

Turns a ranked scorer output into a target shortlist by attaching, per
candidate:

    * nearest gene symbol + distance to TSS (+ve = upstream of TSS,
      i.e. in the promoter direction; -ve = downstream / inside gene)
    * feature class (promoter / 5'UTR / gene_body / 3'UTR / intergenic)
    * CpG-island context (island / shore / shelf / open_sea) + bp
      distance to nearest CGI boundary

Stdlib only. Uses two UCSC hg19 tables: refGene.txt.gz (transcripts →
gene symbols) and cpgIslandExt.txt.gz (CGI intervals). These are not
committed (reproducible artifacts); fetch once:

    mkdir -p data/raw/ucsc
    curl -sSo data/raw/ucsc/refGene.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
    curl -sSo data/raw/ucsc/cpgIslandExt.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz

The script is intentionally small. It is the smallest thing that turns
"good AUC" into "here is a real target shortlist, here is what each
top hit actually *is*". If the v1 output changes decisions, we extend
with repeat / mappability / cCRE axes in v2.

Output is a tab-separated table; the first row carries column names.
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import heapq
import json
from pathlib import Path

#: Promoter = ±this many bp from the TSS. Standard cutoff.
PROMOTER_BP = 1000

#: Shore = ±this many bp from the CGI edge. Shelf = next band of same
#: width. Standard Illumina "Relation_to_UCSC_CpG_Island" convention.
SHORE_BP = 2000
SHELF_BP = 2000


# ---------- loaders ----------


def load_refgene(path: Path, chroms: set[str]) -> dict[str, list[tuple]]:
    """chrom → list of (txStart, txEnd, strand, gene_symbol) sorted by txStart.

    Ranges are 0-indexed half-open, matching UCSC's genePred convention.
    """
    out: dict[str, list[tuple]] = {c: [] for c in chroms}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[2]
            if chrom not in chroms:
                continue
            strand = parts[3]
            tx_start = int(parts[4])
            tx_end = int(parts[5])
            gene = parts[12]
            out[chrom].append((tx_start, tx_end, strand, gene))
    for c in out:
        out[c].sort()
    return out


def load_cgi(path: Path, chroms: set[str]) -> dict[str, list[tuple[int, int]]]:
    """chrom → list of (start, end) sorted by start."""
    out: dict[str, list[tuple[int, int]]] = {c: [] for c in chroms}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[1]
            if chrom not in chroms:
                continue
            out[chrom].append((int(parts[2]), int(parts[3])))
    for c in out:
        out[c].sort()
    return out


# ---------- per-candidate annotation ----------


def _tss(strand: str, tx_start: int, tx_end: int) -> int:
    return tx_start if strand == "+" else tx_end - 1


def annotate_gene(
    chrom: str, pos: int, tx_list: list[tuple]
) -> tuple[str, int, str]:
    """Return (nearest_gene, tss_distance_bp, feature_class).

    `tss_distance_bp` is positive when the candidate lies on the upstream
    (promoter) side of the TSS with respect to the gene's strand, i.e.
    standard regulatory-biology convention.
    """

    if not tx_list:
        return ("", 0, "intergenic")

    # Collect overlapping transcripts. UCSC genePred intervals are 0-indexed
    # half-open: a tx at [tx_start, tx_end) covers `pos` iff
    # tx_start <= pos < tx_end. Earlier failures of this function treated
    # pos == tx_end as still inside, pulling boundary loci into gene_body.
    #
    # Do NOT early-break on tx_end. The list is sorted by tx_start, not tx_end,
    # so a short later-starting transcript can end before `pos` while an
    # earlier-starting long transcript still overlaps. We have to scan every
    # transcript whose start is <= pos.
    starts = [t[0] for t in tx_list]
    idx = bisect.bisect_right(starts, pos) - 1
    overlapping: list[tuple] = []
    j = idx
    while j >= 0:
        if tx_list[j][0] <= pos < tx_list[j][1]:
            overlapping.append(tx_list[j])
        j -= 1

    # Compute the candidate's best match: prefer overlapping, then nearest
    # by TSS distance. Ties broken by tx length (prefer compact isoform)
    # then by gene symbol for determinism.

    def _tss_distance(tx: tuple) -> tuple[int, int]:
        tx_start, tx_end, strand, _ = tx
        tss = _tss(strand, tx_start, tx_end)
        # upstream-positive convention:
        if strand == "+":
            signed = tss - pos
        else:
            signed = pos - tss
        return abs(signed), signed

    if overlapping:
        chosen = min(
            overlapping,
            key=lambda t: (_tss_distance(t)[0], t[1] - t[0], t[3]),
        )
    else:
        # nearest by TSS across all transcripts on the chrom
        chosen = min(tx_list, key=lambda t: (_tss_distance(t)[0], t[1] - t[0], t[3]))

    tx_start, tx_end, strand, gene = chosen
    abs_d, signed_d = _tss_distance(chosen)

    # Feature class
    if chosen in overlapping:
        if abs_d <= PROMOTER_BP:
            feat = "promoter"  # TSS-adjacent
        else:
            feat = "gene_body"
    else:
        if abs_d <= PROMOTER_BP:
            feat = "promoter"
        else:
            feat = "intergenic"

    return (gene, signed_d, feat)


def annotate_cgi(
    chrom: str, pos: int, cgi_list: list[tuple[int, int]]
) -> tuple[str, int]:
    """Return (cpg_island_context, abs_distance_to_nearest_CGI_boundary)."""

    if not cgi_list:
        return ("open_sea", -1)

    # Is pos inside any CGI? UCSC CGI intervals are 0-indexed half-open;
    # a CGI at [start, end) contains pos iff start <= pos < end.
    starts = [c[0] for c in cgi_list]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx >= 0 and cgi_list[idx][0] <= pos < cgi_list[idx][1]:
        return ("island", 0)

    # Distance to nearest CGI boundary — check the two nearest intervals.
    # `end` is exclusive, so positions at `>= end` are measured from `end`.
    cands: list[int] = []
    for j in (idx, idx + 1):
        if 0 <= j < len(cgi_list):
            s, e = cgi_list[j]
            if pos < s:
                cands.append(s - pos)
            elif pos >= e:
                cands.append(pos - e + 1)
    if not cands:
        return ("open_sea", -1)
    d = min(cands)
    if d <= SHORE_BP:
        ctx = "shore"
    elif d <= SHORE_BP + SHELF_BP:
        ctx = "shelf"
    else:
        ctx = "open_sea"
    return (ctx, d)


# ---------- top-K stream ----------


def top_k_by(score_field: str, scored_path: Path, k: int) -> list[dict]:
    """Return top-K ScoredCandidate records by the selected score field.

    Tie-break is deterministic and matches `evaluate_ranking`: score
    descending, candidate_id **ascending** within tied score bands.

    Implementation note: use `heapq.nsmallest` with the exact benchmark key
    `(-score, candidate_id)`. This preserves the same ordering contract as
    `evaluate_ranking`, including prefix-related candidate IDs like `a` < `ab`.
    """
    if score_field not in ("final_score", "p_therapeutic_selectivity"):
        raise ValueError(f"unsupported score_field: {score_field!r}")

    def _iter_records():
        with scored_path.open() as f:
            for line in f:
                r = json.loads(line)
                if score_field == "final_score":
                    score = r["final_score"]
                else:
                    prob = r.get("probabilistic")
                    if not prob:
                        continue
                    score = prob["p_therapeutic_selectivity"]
                yield {
                    "record": r,
                    "score": score,
                    "candidate_id": r["candidate"]["candidate_id"],
                }

    top = heapq.nsmallest(
        k,
        _iter_records(),
        key=lambda x: (-x["score"], x["candidate_id"]),
    )
    return [x["record"] for x in top]


# ---------- CLI ----------


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scored", required=True, type=Path)
    p.add_argument("--score-field", default="p_therapeutic_selectivity",
                   choices=("final_score", "p_therapeutic_selectivity"))
    p.add_argument("--top-k", type=int, default=20)
    p.add_argument("--refgene", type=Path,
                   default=Path("data/raw/ucsc/refGene.txt.gz"))
    p.add_argument("--cpg-islands", type=Path,
                   default=Path("data/raw/ucsc/cpgIslandExt.txt.gz"))
    p.add_argument("--positives", type=Path, default=None,
                   help="Optional file of candidate_ids to flag in the output")
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    # Stream once to know which chroms appear — small catalog, just read the
    # top-K first and learn chroms from them.
    top = top_k_by(args.score_field, args.scored, args.top_k)
    chroms = sorted({r["candidate"]["chrom"] for r in top})
    print(f"top-{args.top_k} spans chroms: {chroms}", flush=True)

    print(f"loading refGene for {chroms}...", flush=True)
    tx_by_chrom = load_refgene(args.refgene, set(chroms))
    for c in chroms:
        print(f"  {c}: {len(tx_by_chrom[c]):,} transcripts", flush=True)

    print(f"loading CpG islands for {chroms}...", flush=True)
    cgi_by_chrom = load_cgi(args.cpg_islands, set(chroms))
    for c in chroms:
        print(f"  {c}: {len(cgi_by_chrom[c]):,} islands", flush=True)

    positives: set[str] = set()
    if args.positives:
        positives = set(args.positives.read_text().split())

    cols = [
        "rank", "candidate_id", "chrom", "critical_c_pos", "strand",
        "pam_family", "pam", "score",
        "beta_tumor_mean", "beta_normal_mean", "delta_beta",
        "p_targ", "p_diff", "p_prot", "p_trust",
        "nearest_gene", "tss_distance_bp", "feature_class",
        "cpg_island_context", "cpg_island_distance_bp",
        "is_positive",
    ]
    out_rows: list[list[str]] = [cols]

    for rank, r in enumerate(top, 1):
        cand = r["candidate"]
        obs = r["observation"]
        prob = r.get("probabilistic") or {}
        cid = cand["candidate_id"]
        chrom = cand["chrom"]
        pos = cand["critical_c_pos"]

        gene, tss_d, feat = annotate_gene(chrom, pos, tx_by_chrom[chrom])
        cgi_ctx, cgi_d = annotate_cgi(chrom, pos, cgi_by_chrom[chrom])

        if args.score_field == "final_score":
            score = r["final_score"]
        else:
            score = prob.get("p_therapeutic_selectivity")

        bt = obs.get("beta_tumor_mean")
        bn = obs.get("beta_normal_mean")
        db = (bn - bt) if (bt is not None and bn is not None) else None

        def f(x, fmt=".4f"):
            if x is None:
                return "NA"
            return format(x, fmt)

        out_rows.append([
            str(rank),
            cid, chrom, str(pos), cand["strand"],
            cand["pam_family"], cand["pam"],
            f(score, ".6g"),
            f(bt, ".3f"), f(bn, ".3f"), f(db, "+.3f"),
            f(prob.get("p_targetable_tumor"), ".3f"),
            f(prob.get("p_differential_protection"), ".3f"),
            f(prob.get("p_protected_normal"), ".3f"),
            f(prob.get("p_observation_trustworthy"), ".3f"),
            gene or "-", str(tss_d) if gene else "NA", feat,
            cgi_ctx, str(cgi_d) if cgi_d >= 0 else "NA",
            "YES" if cid in positives else "-",
        ])

    with args.output.open("w") as fh:
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    print(f"→ {args.output}  ({len(out_rows) - 1} candidates annotated)", flush=True)

    # Also print to stdout for quick inspection — but compact: drop some columns
    compact_cols = ["rank", "candidate_id", "score",
                    "beta_tumor_mean", "beta_normal_mean", "delta_beta",
                    "nearest_gene", "tss_distance_bp", "feature_class",
                    "cpg_island_context", "is_positive"]
    compact_idxs = [cols.index(c) for c in compact_cols]
    widths = [max(len(row[i]) for row in out_rows) for i in compact_idxs]
    for row in out_rows:
        print("  ".join(row[i].ljust(widths[j]) for j, i in enumerate(compact_idxs)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
