"""Phase 4 diagnostics — four follow-up investigations on the Phase 4 results:

  1. V2 probabilistic underperformance — inspect the distribution of each
     factor (P(targetable), P(protected), P(trust)) for positives vs negatives.
  2. P@100 = 0 — show the top-20 by final_score and what they actually hit.
  3. Tighter positives set — restrict to HM450 probes in TSS200/TSS1500
     (promoter) groups, not anywhere in the gene, and re-benchmark.
  4. Roth coordinate system — hg19 vs hg38? Check ESR1 + GATA3 genomic
     positions against the paper's stated coordinates.

No framework changes; this is purely analysis glue. Writes a plain-text
report to data/derived/phase4_diagnostics.txt and per-step supporting files.
"""

from __future__ import annotations

import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DERIVED = REPO / "data" / "derived"
RAW = REPO / "data" / "raw"
REPORT = DERIVED / "phase4_diagnostics.txt"


def _out(fh, *args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=fh)


def step1_probabilistic_distribution(fh):
    """Distribution of V2 probabilistic factors for positives vs negatives."""
    _out(fh, "")
    _out(fh, "=" * 72)
    _out(fh, "STEP 1 — V2 probabilistic: why does it underperform naive?")
    _out(fh, "=" * 72)

    positives = set((DERIVED / "positives.txt").read_text().split())

    buckets = {
        "positives": {"t": [], "p": [], "r": [], "sel": []},
        "negatives": {"t": [], "p": [], "r": [], "sel": []},
    }
    with (DERIVED / "scored_brca_full.jsonl").open() as f:
        for ln in f:
            r = json.loads(ln)
            prob = r.get("probabilistic")
            if not prob:
                continue
            key = "positives" if r["candidate"]["candidate_id"] in positives else "negatives"
            b = buckets[key]
            b["t"].append(prob["p_targetable_tumor"])
            b["p"].append(prob["p_protected_normal"])
            b["r"].append(prob["p_observation_trustworthy"])
            b["sel"].append(
                prob["p_targetable_tumor"]
                * prob["p_protected_normal"]
                * prob["p_observation_trustworthy"]
            )

    def _summary(xs):
        if not xs:
            return "n=0"
        s = sorted(xs)
        return (f"n={len(xs):>7}  mean={statistics.fmean(xs):.4f}  "
                f"med={s[len(s)//2]:.4f}  q90={s[int(len(s)*0.9)]:.4f}  "
                f"max={max(xs):.4f}")

    _out(fh, "")
    _out(fh, "V2 probabilistic factor distributions across 3M-candidate full-BRCA run:")
    for label in ("positives", "negatives"):
        _out(fh, f"\n{label.upper()} (n={len(buckets[label]['t'])})")
        _out(fh, f"  P(targetable_tumor) :  {_summary(buckets[label]['t'])}")
        _out(fh, f"  P(protected_normal):  {_summary(buckets[label]['p'])}")
        _out(fh, f"  P(observation_trust):  {_summary(buckets[label]['r'])}")
        _out(fh, f"  P(therapeutic_sel):   {_summary(buckets[label]['sel'])}")

    # Key diagnostic: how often is P(trust) the zeroing factor?
    n_trust_low = sum(1 for r in buckets["negatives"]["r"] if r < 0.1)
    n_trust_zero = sum(1 for r in buckets["negatives"]["r"] if r == 0.0)
    _out(fh, "")
    _out(fh, "Trust-factor zero-out:")
    _out(fh, f"  negatives with P(trust) == 0.0: {n_trust_zero:>7}  "
          f"({n_trust_zero/len(buckets['negatives']['r']):.1%})")
    _out(fh, f"  negatives with P(trust) < 0.1 : {n_trust_low:>7}  "
          f"({n_trust_low/len(buckets['negatives']['r']):.1%})")

    # If huge fraction is zeroed by trust, the product is degenerate.
    pos_trust_zero = sum(1 for r in buckets["positives"]["r"] if r == 0.0)
    _out(fh, f"  positives with P(trust) == 0.0: {pos_trust_zero:>7}  "
          f"({pos_trust_zero/len(buckets['positives']['r']):.1%})")


def step2_top20_inspection(fh):
    """Top-20 by final_score on LumA-full. What are they, and why do they win?"""
    _out(fh, "")
    _out(fh, "=" * 72)
    _out(fh, "STEP 2 — top-20 by final_score on LumA-full. Why is P@100 = 0?")
    _out(fh, "=" * 72)

    # Load probe annotations (probe_id → gene)
    probe_to_gene: dict[str, str] = {}
    with (DERIVED / "probes_at_roth_genes.tsv").open() as f:
        next(f)
        for ln in f:
            pid, chrom, pos, gene = ln.rstrip().split("\t")
            probe_to_gene[pid] = gene

    positives = set((DERIVED / "positives.txt").read_text().split())

    # Keep top-100 by final_score via a bounded heap. Use candidate_id as
    # tiebreaker so heapq can compare tuples deterministically.
    import heapq
    heap: list[tuple[float, str, dict]] = []
    K = 100
    with (DERIVED / "scored_brca_luma_full.jsonl").open() as f:
        for ln in f:
            r = json.loads(ln)
            key = (r["final_score"], r["candidate"]["candidate_id"], r)
            if len(heap) < K:
                heapq.heappush(heap, key)
            elif key > heap[0]:
                heapq.heapreplace(heap, key)
    top = [(score, r) for score, _cid, r in sorted(heap, reverse=True)]

    # Show top 20 with context
    _out(fh, "")
    _out(fh, f"{'rank':>4}  {'candidate_id':<36}  {'β_tumor':>7}  {'β_normal':>8}  "
          f"{'Δ':>6}  {'evidence':<14}  {'probe':<12}  {'gene':<10}  {'pos?':<6}  {'final':>7}")
    _out(fh, "-" * 125)
    for i, (score, r) in enumerate(top[:20], 1):
        cand = r["candidate"]
        obs = r["observation"]
        bt = obs.get("beta_tumor_mean")
        bn = obs.get("beta_normal_mean")
        pid = obs.get("probe_id") or "-"
        gene = probe_to_gene.get(pid, "")
        delta = (bn - bt) if (bn is not None and bt is not None) else None
        is_pos = "YES" if cand["candidate_id"] in positives else "no"
        _out(fh,
            f"{i:>4}  {cand['candidate_id']:<36}  "
            f"{(f'{bt:.3f}' if bt is not None else '   --'):>7}  "
            f"{(f'{bn:.3f}' if bn is not None else '    --'):>8}  "
            f"{(f'{delta:+.3f}' if delta is not None else '  --'):>6}  "
            f"{obs['evidence_class']:<14}  {pid:<12}  {gene:<10}  {is_pos:<6}  "
            f"{score:>7.3f}")

    # How many of top-100 are in positives?
    top100_pos = sum(1 for _, r in top if r["candidate"]["candidate_id"] in positives)
    _out(fh, "")
    _out(fh, f"Top-100 hit counts:  positives={top100_pos}/100  P@100={top100_pos/100:.3f}")

    # Gene distribution in top-100
    gene_counts: dict[str, int] = defaultdict(int)
    for _score, r in top:
        pid = r["observation"].get("probe_id") or "-"
        gene_counts[probe_to_gene.get(pid, "(no Roth gene)")] += 1
    _out(fh, f"Top-100 gene breakdown: {dict(sorted(gene_counts.items(), key=lambda t: -t[1]))}")


def step3_tighter_positives(fh):
    """Re-build positives.txt restricted to TSS200/TSS1500 probes, rerun benchmark."""
    _out(fh, "")
    _out(fh, "=" * 72)
    _out(fh, "STEP 3 — tighter positives (promoter-only HM450 probes)")
    _out(fh, "=" * 72)

    gene_targets = {"ESR1", "GATA3", "EGFLAM", "VEGFA", "EMX1", "PRDX4"}
    promoter_groups = {"TSS200", "TSS1500", "5'UTR"}

    promoter_probes: list[tuple[str, str, int, str, str]] = []  # probe_id, chrom, pos, gene, group
    with (RAW / "HM450_manifest.csv").open() as f:
        line = ""
        while not line.startswith("IlmnID"):
            line = f.readline()
        header = line.rstrip("\n").split(",")
        idx_id = header.index("IlmnID")
        idx_chr = header.index("CHR")
        idx_pos = header.index("MAPINFO")
        idx_gene = header.index("UCSC_RefGene_Name")
        idx_group = header.index("UCSC_RefGene_Group")
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith("[Controls]"):
                break
            if len(row) <= max(idx_id, idx_chr, idx_pos, idx_gene, idx_group):
                continue
            names = row[idx_gene].split(";")
            groups = row[idx_group].split(";")
            # Each gene/group pairs up position-wise
            for name, grp in zip(names, groups, strict=False):
                if name.strip() in gene_targets and grp.strip() in promoter_groups:
                    try:
                        pos = int(row[idx_pos])
                    except ValueError:
                        continue
                    promoter_probes.append(
                        (row[idx_id], f"chr{row[idx_chr].strip()}", pos, name.strip(), grp.strip())
                    )
                    break

    _out(fh, f"\nPromoter-region probes in {len(gene_targets)} target genes: {len(promoter_probes)}")
    by_gene: dict[str, int] = defaultdict(int)
    for _, _, _, gene, _ in promoter_probes:
        by_gene[gene] += 1
    _out(fh, f"  breakdown: {dict(sorted(by_gene.items(), key=lambda t: -t[1]))}")

    # Scan LumA-full scored candidates for hits within ±50 bp of any promoter probe
    from bisect import bisect_left
    probe_positions: dict[str, list[tuple[int, str]]] = defaultdict(list)
    for _, chrom, pos, gene, _ in promoter_probes:
        probe_positions[chrom].append((pos, gene))
    for chrom in probe_positions:
        probe_positions[chrom].sort()

    WINDOW = 50
    tight_positives: list[str] = []
    tight_by_gene: dict[str, int] = defaultdict(int)
    with (DERIVED / "scored_brca_luma_full.jsonl").open() as f:
        for ln in f:
            r = json.loads(ln)
            cand = r["candidate"]
            chrom = cand["chrom"]
            if chrom not in probe_positions:
                continue
            pos = cand["critical_c_pos"]
            lst = probe_positions[chrom]
            positions = [p for p, _ in lst]
            idx = bisect_left(positions, pos)
            for j in (idx - 1, idx):
                if 0 <= j < len(lst):
                    p, gene = lst[j]
                    if abs(p - pos) <= WINDOW:
                        tight_positives.append(cand["candidate_id"])
                        tight_by_gene[gene] += 1
                        break

    tight_path = DERIVED / "positives_tight.txt"
    tight_path.write_text("\n".join(tight_positives) + "\n")
    _out(fh, f"\nTight positives (within ±{WINDOW} bp of promoter probe): {len(tight_positives)}")
    _out(fh, f"  by gene: {dict(sorted(tight_by_gene.items(), key=lambda t: -t[1]))}")
    _out(fh, f"  (loose positives from Phase 3 were 1687 — a {len(tight_positives)/1687:.1%} shrink)")
    _out(fh, f"  written → {tight_path}")


def step4_coordinate_system(fh):
    """Check whether Roth's coords are hg19 or hg38 by looking up gene spans."""
    _out(fh, "")
    _out(fh, "=" * 72)
    _out(fh, "STEP 4 — Roth coordinate system: hg19 or hg38?")
    _out(fh, "=" * 72)

    # The Roth paper reports these target coordinates (from Fig. 5d caption):
    #   ESR1  T17  Chr 6 : 151,690,043
    #   GATA3 T18  Chr 10 :  8,045,463
    #   EGFLAM T11 Chr 5 : 38,258,782
    roth = {
        "ESR1":   ("chr6",  151_690_043),
        "GATA3":  ("chr10",  8_045_463),
        "EGFLAM": ("chr5",  38_258_782),
    }
    # Known gene spans. Source: UCSC genome browser (Ensembl canonical transcripts).
    gene_spans = {
        "ESR1": {
            "hg19": ("chr6",  152_126_182, 152_420_094),
            "hg38": ("chr6",  151_656_691, 151_845_984),
        },
        "GATA3": {
            "hg19": ("chr10",  8_096_668,   8_117_164),
            "hg38": ("chr10",  8_054_471,   8_074_967),
        },
        "EGFLAM": {
            "hg19": ("chr5",  38_258_583,  38_513_500),
            "hg38": ("chr5",  38_258_470,  38_513_387),
        },
    }

    _out(fh, "")
    _out(fh, f"{'gene':<8}  {'roth coord':>14}  {'hg19 span':>25}  {'hg38 span':>25}  verdict")
    _out(fh, "-" * 100)
    verdict_counts = {"hg19": 0, "hg38": 0, "both": 0, "neither": 0}
    for gene, (_, pos) in roth.items():
        hg19_chrom, hg19_s, hg19_e = gene_spans[gene]["hg19"]
        hg38_chrom, hg38_s, hg38_e = gene_spans[gene]["hg38"]
        # Accept matches within gene body ±50 kb (covers promoter + nearby regulatory)
        SLOP = 50_000
        in_hg19 = (hg19_s - SLOP) <= pos <= (hg19_e + SLOP)
        in_hg38 = (hg38_s - SLOP) <= pos <= (hg38_e + SLOP)
        if in_hg19 and in_hg38:
            v = "both (ambiguous)"; verdict_counts["both"] += 1
        elif in_hg38:
            v = "hg38"; verdict_counts["hg38"] += 1
        elif in_hg19:
            v = "hg19"; verdict_counts["hg19"] += 1
        else:
            v = "neither"; verdict_counts["neither"] += 1
        _out(fh, f"{gene:<8}  {pos:>14,}  "
             f"{hg19_s:>12,}–{hg19_e:<12,}  "
             f"{hg38_s:>12,}–{hg38_e:<12,}  {v}")

    _out(fh, "")
    _out(fh, f"Verdict counts: {verdict_counts}")
    if verdict_counts["hg38"] > verdict_counts["hg19"]:
        _out(fh, "→ Roth paper coords are hg38.")
    elif verdict_counts["hg19"] > verdict_counts["hg38"]:
        _out(fh, "→ Roth paper coords are hg19.")
    else:
        _out(fh, "→ Inconclusive.")
    _out(fh, "")
    _out(fh, "Impact on Phase 3/4 analysis:")
    _out(fh, "  * Catalog + HM450 manifest both hg19.")
    _out(fh, "  * GDC methylation data is hg38, but matched by probe_id → build-agnostic.")
    _out(fh, "  * Positives selected by gene-symbol intersection → robust to build.")
    _out(fh, "  * Only the raw Roth (chrom, pos) coordinates from the paper captions")
    _out(fh, "    would need liftOver before being used directly; we never did.")


def main() -> int:
    DERIVED.mkdir(parents=True, exist_ok=True)
    with REPORT.open("w") as fh:
        _out(fh, "Phase 4 diagnostics")
        _out(fh, "Generated by scripts/phase4_diagnostics.py")
        step1_probabilistic_distribution(fh)
        step2_top20_inspection(fh)
        step3_tighter_positives(fh)
        step4_coordinate_system(fh)
    print(f"\nreport → {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
