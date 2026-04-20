"""Build the Roth et al. (Nature 2026) paper-comparable cohort from GSE322563.

GSE322563 is the actual Infinium MethylationEPIC v2.0 data Roth et al.
deposited (the main paper's data-availability statement prints `GSE32256`,
which is wrong — the supplement's reporting summary and the public GEO
record agree that the correct accession is `GSE322563`).

Layout:
    MCF-7   R1 / R2   → tumor arm   (n=2)
    MCF-10A R1 / R2   → normal arm  (n=2)

Platform is EPIC v2.0 (GPL33022), not HM450. Fast probe-intersect path:
EPIC v2 probe IDs are `cgNNNNNNNN_BCXX`; we strip the `_BCXX` suffix to
recover the canonical HM450-compatible ID, and keep only probes that
exist in our HM450 annotation (`data/raw/probes_hg19.tsv`). Probes that
appear multiple times under different `_BCXX` suffixes are collapsed by
per-sample mean across the non-NA variants — same underlying CpG, multiple
bead designs.

Caveats worth citing in any downstream report:
  * n=2 per side means `p_trust` saturates at 0.95 * 2/30 = 0.0633 for
    every EXACT-evidence candidate — tie-band is expected to be larger
    than the n=3 GSE77348 surrogate (299) and much larger than the
    n=305/50 GSE69914 cohort (2).
  * EPIC v2 drops some HM450 probes. Retention for the Roth positives
    (ESR1, GATA3, EGFLAM, VEGFA) is reported in PROVENANCE.md.
  * This is the paper-comparable biological test, NOT a high-n
    rank-metric validation.

Outputs (LocalSummaryBackend-compatible):
    <out>/tumor_summary.tsv
    <out>/normal_summary.tsv
    <out>/probes.tsv
    <out>/PROVENANCE.md
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_surrogate_cohort import summarize, write_summary_tsv  # noqa: E402

#: EPIC v2 probe IDs append a two-letter beadchip/design code plus a numeric
#: index (`_BC11`, `_TC21`, `_TO22`, `_BO11`, etc.) to the canonical
#: HM450-style identifier. Strip any such suffix to intersect with HM450.
_EPIC_V2_SUFFIX = re.compile(r"_[A-Z]{2}\d+$")


def _canonical(probe_id: str) -> str:
    return _EPIC_V2_SUFFIX.sub("", probe_id)


def _beta(s: str) -> float | None:
    s = s.strip()
    if s in ("", "NA", "null", "nan", "NaN"):
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if 0.0 <= v <= 1.0:
        return v
    return None


def read_beta_matrix(
    path: Path,
    tumor_cols: list[str],
    normal_cols: list[str],
) -> tuple[int, dict[str, list[float | None]], dict[str, list[float | None]]]:
    """Return (n_epic_probes_read, tumor_by_probe, normal_by_probe).

    Each output dict maps canonical (HM450-compatible) probe_id →
    [β per sample], with `None` for missing values. Multi-variant probes
    (same canonical ID under different `_BCxx` suffixes) are collapsed by
    per-sample mean across non-NA variants.
    """

    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    tumor_idx = [header.index(c) for c in tumor_cols]
    normal_idx = [header.index(c) for c in normal_cols]

    # Accumulators: per canonical probe, per sample index → list of β values
    # across beadchip variants. We mean them at the end.
    tumor_acc: dict[str, list[list[float]]] = defaultdict(
        lambda: [[] for _ in tumor_cols]
    )
    normal_acc: dict[str, list[list[float]]] = defaultdict(
        lambda: [[] for _ in normal_cols]
    )
    n_read = 0

    with gzip.open(path, "rt") as f:
        f.readline()  # header consumed above
        for line in f:
            parts = line.rstrip("\n").split("\t")
            canon = _canonical(parts[0])
            for j, col_idx in enumerate(tumor_idx):
                v = _beta(parts[col_idx])
                if v is not None:
                    tumor_acc[canon][j].append(v)
            for j, col_idx in enumerate(normal_idx):
                v = _beta(parts[col_idx])
                if v is not None:
                    normal_acc[canon][j].append(v)
            n_read += 1

    # Collapse per-sample lists to a scalar β or None.
    def _collapse(acc: dict[str, list[list[float]]]) -> dict[str, list[float | None]]:
        out: dict[str, list[float | None]] = {}
        for canon, per_sample in acc.items():
            out[canon] = [
                (statistics.fmean(xs) if xs else None) for xs in per_sample
            ]
        return out

    return n_read, _collapse(tumor_acc), _collapse(normal_acc)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--beta-matrix", required=True, type=Path,
                   help="GSE322563_beta_matrix_EPIC_v2.txt.gz")
    p.add_argument("--tumor-cols", nargs="+", default=["MCF71", "MCF72"])
    p.add_argument("--normal-cols", nargs="+", default=["MCF10A1", "MCF10A2"])
    p.add_argument("--hm450-probes", required=True, type=Path,
                   help="data/raw/probes_hg19.tsv — canonical HM450 annotation")
    p.add_argument("--roth-gene-probes", type=Path, default=None,
                   help="Optional data/derived/probes_at_roth_genes.tsv for "
                        "per-gene retention reporting")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"streaming {args.beta_matrix.name}...", flush=True)
    n_epic, tumor_by, normal_by = read_beta_matrix(
        args.beta_matrix, args.tumor_cols, args.normal_cols,
    )
    print(f"  EPIC v2 probe rows read:             {n_epic:,}", flush=True)
    print(f"  canonical probes (tumor side):       {len(tumor_by):,}", flush=True)

    # Load HM450 probe universe for the fast intersect.
    hm450 = set()
    with args.hm450_probes.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        pid_col = header.index("probe_id")
        for line in f:
            parts = line.rstrip("\n").split("\t")
            hm450.add(parts[pid_col])
    print(f"  HM450 probe universe:                {len(hm450):,}", flush=True)

    # Probes that appear on EITHER side (avoid losing records where one
    # cohort has all-NA for a probe covered on the other).
    all_probes = (set(tumor_by) | set(normal_by))
    overlap = all_probes & hm450
    print(f"  EPIC-v2 ∩ HM450 (canonical match):    {len(overlap):,} "
          f"({100 * len(overlap) / len(hm450):.1f}% of HM450)", flush=True)

    # Filter summaries to the HM450-resolvable probes — that's what
    # LocalSummaryBackend + our positives list consumes. Fill missing
    # sides with all-None arrays so the summarize() helper emits a
    # well-defined zero-n row.
    n_tumor_cols = len(args.tumor_cols)
    n_normal_cols = len(args.normal_cols)
    tumor_by = {
        pid: tumor_by.get(pid, [None] * n_tumor_cols) for pid in overlap
    }
    normal_by = {
        pid: normal_by.get(pid, [None] * n_normal_cols) for pid in overlap
    }

    # Roth gene retention report
    retention_by_gene: dict[str, tuple[int, int]] = {}
    if args.roth_gene_probes and args.roth_gene_probes.exists():
        gene_probes: dict[str, set[str]] = defaultdict(set)
        with args.roth_gene_probes.open() as f:
            f.readline()  # header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                pid, _chrom, _pos, gene = parts[0], parts[1], parts[2], parts[3]
                gene_probes[gene].add(pid)
        for gene, probes in gene_probes.items():
            retained = len(probes & overlap)
            retention_by_gene[gene] = (retained, len(probes))
        print("  Roth-gene probe retention:", flush=True)
        for gene, (kept, total) in sorted(retention_by_gene.items()):
            print(f"    {gene:<10}  {kept:>3} / {total:<3}  ({100 * kept / total:.0f}%)", flush=True)

    print("summarizing tumor arm (MCF-7 R1/R2)...", flush=True)
    t_summary = summarize(tumor_by)
    n_tumor = write_summary_tsv(args.output_dir / "tumor_summary.tsv", t_summary)

    print("summarizing normal arm (MCF-10A R1/R2)...", flush=True)
    n_summary = summarize(normal_by)
    n_normal = write_summary_tsv(args.output_dir / "normal_summary.tsv", n_summary)

    shutil.copy(args.hm450_probes, args.output_dir / "probes.tsv")

    prov_lines = [
        "# GSE322563 cohort build — Roth et al. Nature 2026 actual samples",
        "",
        f"- β matrix: {args.beta_matrix}",
        "- tumor arm (MCF-7):   "
        f"{', '.join(args.tumor_cols)}   (n={len(args.tumor_cols)})",
        "- normal arm (MCF-10A): "
        f"{', '.join(args.normal_cols)}  (n={len(args.normal_cols)})",
        "- platform: Infinium MethylationEPIC v2.0 (GPL33022)",
        "- intersect policy: canonical probe_id match after stripping `_BCxx` suffix",
        f"- EPIC-v2 probe rows read:           {n_epic:,}",
        f"- HM450 probe universe:              {len(hm450):,}",
        f"- probes retained (EPIC-v2 ∩ HM450): "
        f"{len(overlap):,}  ({100 * len(overlap) / len(hm450):.1f}% of HM450)",
    ]
    if retention_by_gene:
        prov_lines.append("- Roth-gene probe retention:")
        for gene, (kept, total) in sorted(retention_by_gene.items()):
            prov_lines.append(
                f"    - {gene}: {kept}/{total} ({100 * kept / total:.0f}%)"
            )
    prov_lines += [
        "",
        "## Caveats",
        "",
        "- **n=2 per side** means `p_trust` saturates at 0.0633 for EXACT-evidence",
        "  candidates. Tie-band at the top of any p_trust-multiplied composite is",
        "  expected to be large; treat P@K as secondary to AUC and top-list biology.",
        "- EPIC v2 drops some HM450 probes; retention reported above.",
        "- This is the paper-comparable biological test, NOT a high-n rank-metric",
        "  validation.",
        "",
        "## Citation",
        "",
        "Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.",
        "Molecular basis for methylation-sensitive editing by Cas9.",
        "Nature (2026). DOI 10.1038/s41586-026-10384-z.",
        "GEO accession: GSE322563 (per supplementary reporting summary).",
        "Main paper's data-availability statement prints 'GSE32256' (typo).",
    ]
    (args.output_dir / "PROVENANCE.md").write_text("\n".join(prov_lines) + "\n")

    print(f"→ {args.output_dir}/  "
          f"(tumor: {n_tumor:,} probes, normal: {n_normal:,} probes)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
