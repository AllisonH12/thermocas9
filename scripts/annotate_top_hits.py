"""Annotate the top-K candidates from a scored JSONL with regulatory context.

Turns a ranked scorer output into a target shortlist by attaching, per
candidate:

    * nearest gene symbol + distance to TSS (+ve = upstream of TSS,
      i.e. in the promoter direction; -ve = downstream / inside gene)
    * feature class (promoter / 5'UTR / gene_body / 3'UTR / intergenic)
    * CpG-island context (island / shore / shelf / open_sea) + bp
      distance to nearest CGI boundary
    * repeat overlap (RepeatMasker class/family/name, or "-")
    * DNase-hypersensitivity cluster overlap (an ENCODE-wide regulatory-
      activity proxy; v1 stand-in for the formal cCRE Registry which
      lives behind a JS challenge on screen.wenglab.org)

Stdlib only. UCSC hg19 table sources — fetch once:

    mkdir -p data/raw/ucsc
    curl -sSo data/raw/ucsc/refGene.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz
    curl -sSo data/raw/ucsc/cpgIslandExt.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cpgIslandExt.txt.gz
    curl -sSo data/raw/ucsc/rmsk.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz
    curl -sSo data/raw/ucsc/wgEncodeRegDnaseClusteredV3.txt.gz \\
        https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeRegDnaseClusteredV3.txt.gz

The --rmsk and --dnase arguments are optional; if omitted, the repeat
and DNase columns are emitted as "NA" so downstream consumers can tell
"annotation not available" from "no overlap".

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


def load_repeats(
    path: Path, chroms: set[str]
) -> dict[str, list[tuple[int, int, str, str, str]]]:
    """chrom → list of (start, end, repName, repClass, repFamily) sorted by start.

    UCSC rmsk.txt.gz schema (0-indexed half-open genomic coords):
        bin swScore milliDiv milliDel milliIns
        genoName genoStart genoEnd genoLeft strand
        repName repClass repFamily repStart repEnd repLeft id
    """
    out: dict[str, list[tuple[int, int, str, str, str]]] = {c: [] for c in chroms}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[5]
            if chrom not in chroms:
                continue
            out[chrom].append(
                (int(parts[6]), int(parts[7]), parts[10], parts[11], parts[12])
            )
    for c in out:
        out[c].sort()
    return out


def load_regulatory(
    path: Path, chroms: set[str]
) -> dict[str, list[tuple[int, int, int]]]:
    """chrom → list of (start, end, sourceCount) sorted by start.

    UCSC wgEncodeRegDnaseClusteredV3 schema (0-indexed half-open):
        bin chrom chromStart chromEnd name score sourceCount sourceIds sourceScores

    `sourceCount` is the number of ENCODE cell types in which the DNase
    hypersensitivity peak was observed — a rough proxy for how ubiquitous
    the open-chromatin region is (higher = more constitutively open).
    """
    out: dict[str, list[tuple[int, int, int]]] = {c: [] for c in chroms}
    with gzip.open(path, "rt") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[1]
            if chrom not in chroms:
                continue
            out[chrom].append((int(parts[2]), int(parts[3]), int(parts[6])))
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


def annotate_repeat(
    chrom: str, pos: int, rep_list: list[tuple[int, int, str, str, str]]
) -> tuple[bool, str, str, str]:
    """Return (in_repeat, repeat_class, repeat_family, repeat_name).

    UCSC rmsk intervals are 0-indexed half-open: a repeat at [start, end)
    contains pos iff start <= pos < end. rmsk intervals can nest or
    overlap — RepeatMasker typically fragments an outer repeat at each
    inner insertion, so deep nesting is uncommon in practice, but it is
    not structurally forbidden. A defensive full backward scan costs
    nothing at top-K and guarantees correctness regardless of track-
    specific invariants. Same shape as `annotate_gene`. First match in
    backward-scan order wins, so the innermost (largest-start)
    overlapping repeat is preferred — typically the most specific
    annotation a reader wants.
    """
    if not rep_list:
        return (False, "-", "-", "-")
    starts = [r[0] for r in rep_list]
    idx = bisect.bisect_right(starts, pos) - 1
    j = idx
    while j >= 0:
        s, e, name, cls, fam = rep_list[j]
        if s <= pos < e:
            return (True, cls, fam, name)
        j -= 1
    return (False, "-", "-", "-")


def annotate_regulatory(
    chrom: str, pos: int, reg_list: list[tuple[int, int, int]]
) -> tuple[bool, int]:
    """Return (in_dnase_cluster, source_count).

    UCSC ENCODE DNase-clustered intervals are 0-indexed half-open. This
    is the v1 proxy for the SCREEN/ENCODE cCRE Registry (which we would
    prefer but currently cannot fetch directly from screen.wenglab.org).
    """
    if not reg_list:
        return (False, 0)
    starts = [r[0] for r in reg_list]
    idx = bisect.bisect_right(starts, pos) - 1
    if idx >= 0:
        s, e, count = reg_list[idx]
        if s <= pos < e:
            return (True, count)
    return (False, 0)


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


# ---------- flag rules & markdown emitter ----------

#: Thresholds the flag rules depend on. Pulled out here so a reader can see
#: "promoter" vs "small Δβ" vs "sparse evidence" as one coherent filter set.
FLAG_LOW_PTRUST = 0.10
FLAG_SMALL_DELTA_BETA = 0.30


def compute_flags(card: dict) -> list[str]:
    """Return a list of short tagged strings summarizing each candidate.

    The rules are intentionally conservative and deterministic from the
    annotated fields — they are a triage aid for the person picking
    guides, not a re-ranking. A card can carry multiple flags; a card
    with no flags is neither good nor bad, it just didn't trip any rule.
    """
    flags: list[str] = []

    feat = card.get("feature_class")
    cgi = card.get("cpg_island_context")
    in_rep = card.get("in_repeat")  # True / False / None
    rep_fam = card.get("repeat_family")
    in_dnase = card.get("in_dnase_cluster")  # True / False / None
    dnase_count = card.get("dnase_source_count")  # int / None
    p_trust = card.get("p_trust")
    delta_beta = card.get("delta_beta")

    # Positives — traits a bench reader cares to see front-loaded.
    if feat == "promoter" and cgi == "island":
        flags.append("STRONG: island-localized promoter")
    if feat == "promoter" and in_dnase is True and (dnase_count or 0) >= 10:
        flags.append(
            f"STRONG: active promoter (DNase support in {dnase_count} cell types)"
        )

    # Cautions — things the reader should eyeball before ordering.
    if in_rep is True:
        flags.append(f"CAUTION: overlaps {rep_fam or 'a'} repeat")
    if p_trust is not None and p_trust < FLAG_LOW_PTRUST:
        flags.append(f"CAUTION: sparse evidence (p_trust={p_trust:.3f})")
    if delta_beta is not None and delta_beta < FLAG_SMALL_DELTA_BETA:
        flags.append(f"CAUTION: small differential (Delta_beta={delta_beta:+.3f})")
    if feat == "gene_body" and in_dnase is False:
        flags.append("NOTE: gene body with no DNase support")

    return flags


def emit_markdown(cards: list[dict], path: Path, *, cohort_label: str) -> None:
    """Write a per-candidate Markdown shortlist.

    Structure: a short summary header, followed by one card per
    candidate. Each card shows the decision-critical fields (coords,
    gene, feature class, methylation numbers, score components) plus
    any flags from `compute_flags(card)`.
    """
    n = len(cards)
    n_strong = sum(1 for c in cards if any(f.startswith("STRONG:") for f in c["_flags"]))
    n_caution = sum(1 for c in cards if any(f.startswith("CAUTION:") for f in c["_flags"]))
    n_clean = sum(1 for c in cards if not any(f.startswith("CAUTION:") for f in c["_flags"]))

    lines: list[str] = []
    lines.append(f"# Top-{n} ThermoCas9 target shortlist — {cohort_label}\n")
    lines.append(
        "Human-readable companion to the TSV. For each candidate: "
        "coordinates, nearest gene, feature class, methylation numbers, "
        "score components, and rule-based flags. "
        "STRONG flags are traits a bench reader wants to see; "
        "CAUTION flags are worth eyeballing before ordering guides. "
        "Flags are a triage aid, not a re-ranking.\n"
    )
    lines.append(
        f"- **{n_strong}** candidate(s) with at least one STRONG flag.\n"
        f"- **{n_caution}** candidate(s) with at least one CAUTION flag.\n"
        f"- **{n_clean}** candidate(s) with no CAUTION flags.\n"
    )
    lines.append("---\n")

    for c in cards:
        pos_str = c.get("is_positive") or "-"
        pos_tag = " — **labeled positive**" if pos_str == "YES" else ""
        lines.append(
            f"## Rank {c['rank']} — `{c['candidate_id']}`{pos_tag}\n"
        )

        gene = c.get("nearest_gene") or "-"
        tss_d = c.get("tss_distance_bp")
        tss_d_str = f"{tss_d:+d} bp to TSS" if isinstance(tss_d, int) else "TSS distance: NA"
        feat = c.get("feature_class") or "-"
        cgi = c.get("cpg_island_context") or "-"
        cgi_d = c.get("cpg_island_distance_bp")
        cgi_str = f"{cgi}" + (f" ({cgi_d} bp to boundary)" if isinstance(cgi_d, int) and cgi_d > 0 else "")

        lines.append(
            f"**Locus.** {c['chrom']}:{c['critical_c_pos']} ({c['strand']}), "
            f"PAM `{c['pam']}` (family `{c['pam_family']}`).  \n"
            f"**Gene / feature.** {gene}, {feat} ({tss_d_str}); CpG context: {cgi_str}.\n"
        )

        bt = c.get("beta_tumor_mean")
        bn = c.get("beta_normal_mean")
        db = c.get("delta_beta")
        beta_parts = []
        if isinstance(bt, float):
            beta_parts.append(f"beta_tumor = {bt:.3f}")
        if isinstance(bn, float):
            beta_parts.append(f"beta_normal = {bn:.3f}")
        if isinstance(db, float):
            beta_parts.append(f"Delta_beta = {db:+.3f}")
        lines.append("**Methylation.** " + ", ".join(beta_parts) + ".\n")

        prob_parts = []
        for key, label in (
            ("p_targ", "p_targ"),
            ("p_diff", "p_diff"),
            ("p_trust", "p_trust"),
        ):
            v = c.get(key)
            if isinstance(v, float):
                prob_parts.append(f"{label} = {v:.3f}")
        score = c.get("score")
        score_str = f"{score:.4g}" if isinstance(score, float) else "NA"
        lines.append(
            f"**Score.** {score_str}  ({', '.join(prob_parts)}).\n"
        )

        rep_overlay = []
        if c.get("in_repeat") is True:
            rep_overlay.append(
                f"repeat: {c.get('repeat_class')}/{c.get('repeat_family')}/{c.get('repeat_name')}"
            )
        elif c.get("in_repeat") is False:
            rep_overlay.append("no repeat overlap")
        if c.get("in_dnase_cluster") is True:
            rep_overlay.append(
                f"DNase cluster in {c.get('dnase_source_count')} cell type(s)"
            )
        elif c.get("in_dnase_cluster") is False:
            rep_overlay.append("no DNase-cluster support")
        if rep_overlay:
            lines.append("**Structural context.** " + "; ".join(rep_overlay) + ".\n")

        flags = c.get("_flags") or []
        if flags:
            lines.append("**Flags.**")
            for f in flags:
                lines.append(f"- {f}")
            lines.append("")  # blank line after list
        lines.append("")  # blank line between cards

    path.write_text("\n".join(lines))
    print(f"-> {path}  ({n} candidate cards)", flush=True)


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
    p.add_argument("--rmsk", type=Path, default=None,
                   help="Optional UCSC rmsk.txt.gz (RepeatMasker track, hg19). "
                        "If omitted, repeat columns are 'NA'.")
    p.add_argument("--dnase", type=Path, default=None,
                   help="Optional UCSC wgEncodeRegDnaseClusteredV3.txt.gz "
                        "(ENCODE DNase-HS clusters, hg19). Stand-in for the "
                        "formal SCREEN cCRE Registry. If omitted, DNase "
                        "columns are 'NA'.")
    p.add_argument("--positives", type=Path, default=None,
                   help="Optional file of candidate_ids to flag in the output")
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--markdown", type=Path, default=None,
                   help="Optional human-readable Markdown companion. Emits a "
                        "per-candidate card with decision-critical fields and "
                        "rule-based flags (island-localized promoter / in-repeat "
                        "caution / sparse-evidence caution / etc.). Intended for "
                        "the person who will actually order guides.")
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

    rep_by_chrom: dict[str, list[tuple[int, int, str, str, str]]] | None = None
    if args.rmsk:
        print(f"loading RepeatMasker for {chroms}...", flush=True)
        rep_by_chrom = load_repeats(args.rmsk, set(chroms))
        for c in chroms:
            print(f"  {c}: {len(rep_by_chrom[c]):,} repeats", flush=True)

    reg_by_chrom: dict[str, list[tuple[int, int, int]]] | None = None
    if args.dnase:
        print(f"loading DNase clusters for {chroms}...", flush=True)
        reg_by_chrom = load_regulatory(args.dnase, set(chroms))
        for c in chroms:
            print(f"  {c}: {len(reg_by_chrom[c]):,} DNase clusters", flush=True)

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
        "in_repeat", "repeat_class", "repeat_family", "repeat_name",
        "in_dnase_cluster", "dnase_source_count",
        "flags",
        "is_positive",
    ]

    # Build typed cards first; derive TSV rows and optional Markdown from
    # the same source of truth so the two emissions can't drift.
    cards: list[dict] = []

    for rank, r in enumerate(top, 1):
        cand = r["candidate"]
        obs = r["observation"]
        prob = r.get("probabilistic") or {}
        cid = cand["candidate_id"]
        chrom = cand["chrom"]
        pos = cand["critical_c_pos"]

        gene, tss_d, feat = annotate_gene(chrom, pos, tx_by_chrom[chrom])
        cgi_ctx, cgi_d = annotate_cgi(chrom, pos, cgi_by_chrom[chrom])

        if rep_by_chrom is not None:
            in_rep_bool, rep_cls, rep_fam, rep_name = annotate_repeat(
                chrom, pos, rep_by_chrom[chrom]
            )
            in_rep_val: bool | None = in_rep_bool
        else:
            in_rep_val, rep_cls, rep_fam, rep_name = None, "NA", "NA", "NA"

        if reg_by_chrom is not None:
            in_dnase_bool, dnase_count = annotate_regulatory(
                chrom, pos, reg_by_chrom[chrom]
            )
            in_dnase_val: bool | None = in_dnase_bool
        else:
            in_dnase_val, dnase_count = None, 0

        if args.score_field == "final_score":
            score = r["final_score"]
        else:
            score = prob.get("p_therapeutic_selectivity")

        bt = obs.get("beta_tumor_mean")
        bn = obs.get("beta_normal_mean")
        db = (bn - bt) if (bt is not None and bn is not None) else None

        card = {
            "rank": rank,
            "candidate_id": cid,
            "chrom": chrom,
            "critical_c_pos": pos,
            "strand": cand["strand"],
            "pam_family": cand["pam_family"],
            "pam": cand["pam"],
            "score": score,
            "beta_tumor_mean": bt,
            "beta_normal_mean": bn,
            "delta_beta": db,
            "p_targ": prob.get("p_targetable_tumor"),
            "p_diff": prob.get("p_differential_protection"),
            "p_prot": prob.get("p_protected_normal"),
            "p_trust": prob.get("p_observation_trustworthy"),
            "nearest_gene": gene or None,
            "tss_distance_bp": tss_d if gene else None,
            "feature_class": feat,
            "cpg_island_context": cgi_ctx,
            "cpg_island_distance_bp": cgi_d if cgi_d >= 0 else None,
            "in_repeat": in_rep_val,
            "repeat_class": rep_cls,
            "repeat_family": rep_fam,
            "repeat_name": rep_name,
            "in_dnase_cluster": in_dnase_val,
            "dnase_source_count": dnase_count if in_dnase_val else 0,
            "is_positive": "YES" if cid in positives else "-",
        }
        card["_flags"] = compute_flags(card)
        cards.append(card)

    def f(x, fmt=".4f"):
        if x is None:
            return "NA"
        return format(x, fmt)

    def tri(v: bool | None) -> str:
        if v is None:
            return "NA"
        return "YES" if v else "-"

    out_rows: list[list[str]] = [cols]
    for c in cards:
        out_rows.append([
            str(c["rank"]),
            c["candidate_id"], c["chrom"], str(c["critical_c_pos"]), c["strand"],
            c["pam_family"], c["pam"],
            f(c["score"], ".6g"),
            f(c["beta_tumor_mean"], ".3f"), f(c["beta_normal_mean"], ".3f"),
            f(c["delta_beta"], "+.3f"),
            f(c["p_targ"], ".3f"), f(c["p_diff"], ".3f"),
            f(c["p_prot"], ".3f"), f(c["p_trust"], ".3f"),
            c["nearest_gene"] or "-",
            str(c["tss_distance_bp"]) if c["tss_distance_bp"] is not None else "NA",
            c["feature_class"],
            c["cpg_island_context"],
            str(c["cpg_island_distance_bp"]) if c["cpg_island_distance_bp"] is not None else "NA",
            tri(c["in_repeat"]), c["repeat_class"], c["repeat_family"], c["repeat_name"],
            tri(c["in_dnase_cluster"]),
            str(c["dnase_source_count"]) if c["in_dnase_cluster"] is True else
                ("0" if c["in_dnase_cluster"] is False else "NA"),
            "; ".join(c["_flags"]) if c["_flags"] else "-",
            c["is_positive"],
        ])

    with args.output.open("w") as fh:
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    print(f"-> {args.output}  ({len(out_rows) - 1} candidates annotated)", flush=True)

    if args.markdown is not None:
        # Cohort label = scored JSONL stem minus the leading "scored_" prefix
        # if present. Purely cosmetic.
        stem = args.scored.stem
        cohort_label = stem[len("scored_"):] if stem.startswith("scored_") else stem
        emit_markdown(cards, args.markdown, cohort_label=cohort_label)

    # Also print to stdout for quick inspection — but compact: drop some columns
    compact_cols = ["rank", "candidate_id", "score",
                    "beta_tumor_mean", "beta_normal_mean", "delta_beta",
                    "nearest_gene", "tss_distance_bp", "feature_class",
                    "cpg_island_context",
                    "in_repeat", "repeat_family",
                    "in_dnase_cluster",
                    "is_positive"]
    compact_idxs = [cols.index(c) for c in compact_cols]
    widths = [max(len(row[i]) for row in out_rows) for i in compact_idxs]
    for row in out_rows:
        print("  ".join(row[i].ljust(widths[j]) for j, i in enumerate(compact_idxs)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
