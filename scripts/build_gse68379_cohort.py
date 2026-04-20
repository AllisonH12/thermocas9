"""Build a GSE68379 breast cell-line tumor cohort on HM450.

GSE68379 (Iorio et al., Sanger GDSC1000) profiles ~1028 cancer cell lines
on HM450. No in-study normal. We use the 52 breast-primary-site cell lines
as the tumor arm and pair with an external healthy-normal cohort
(typically GSE69914 status=0, n=50) as the normal arm.

This is a cross-series composition — biologically a "breast cancer cell
panel vs primary healthy breast epithelium" comparison. Batch effects
between Sanger (GSE68379) and the external normal study are a real
caveat and are flagged in PROVENANCE.md. The test purpose is whether
V2.5's ranking behavior on matched cell-line cohorts survives when
the tumor side has 52 heterogeneous breast lines instead of just MCF-7.

Inputs:
    --beta-matrix       GSE68379_Matrix.processed.txt.gz
    --breast-samples    one GSM per line (52 breast lines, from metadata)
    --hm450-probes      probes_hg19.tsv (copied through as probes.tsv)

Outputs:
    <out>/tumor_summary.tsv
    <out>/probes.tsv
    <out>/PROVENANCE.md

(The normal_summary.tsv is taken from the external normal cohort — copy
it in after this script runs.)
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_surrogate_cohort import summarize, write_summary_tsv  # noqa: E402


def _beta(s: str) -> float | None:
    s = s.strip().strip('"')
    if s in ("", "NA", "null", "nan", "NaN"):
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    if 0.0 <= v <= 1.0:
        return v
    return None


def _resolve_columns(header: list[str], wanted: set[str]) -> dict[str, int]:
    """Map each wanted cell-line name → column index in the matrix header.

    GSE68379's processed matrix names columns `{cell_line}_AVG.Beta`
    (e.g. `MCF7_AVG.Beta`). Strip the suffix to match against the
    cell-line-name set; also accept exact matches if the caller supplies
    full column names.
    """
    out: dict[str, int] = {}
    suffix = "_AVG.Beta"
    for i, h in enumerate(header):
        if h in wanted:
            out[h] = i
            continue
        if h.endswith(suffix):
            core = h[: -len(suffix)]
            if core in wanted:
                out[core] = i
    return out


def read_beta_matrix(
    path: Path, sample_ids: set[str]
) -> tuple[int, dict[str, list[float | None]]]:
    """Stream probe × sample β matrix, extracting β values for the requested
    sample IDs only. Returns (n_probe_rows, probe → [β per sample])."""

    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        # First column is a probe identifier; strip surrounding quotes on every cell.
        header = [h.strip('"') for h in header]
    # Pass 1: resolve columns by exact match.
    col_map = _resolve_columns(header, sample_ids)
    missing = sample_ids - set(col_map)
    if missing:
        raise ValueError(
            f"{len(missing)} requested samples not found in matrix header: "
            f"{sorted(missing)[:5]}..."
        )
    ordered = sorted(col_map, key=lambda s: col_map[s])
    idxs = [col_map[s] for s in ordered]
    print(f"  resolved {len(ordered)} samples → cols {idxs[0]}..{idxs[-1]}", flush=True)

    by_probe: dict[str, list[float | None]] = {}
    n_rows = 0
    with gzip.open(path, "rt") as f:
        f.readline()  # consumed
        for line in f:
            parts = line.rstrip("\n").split("\t")
            probe_id = parts[0].strip('"')
            by_probe[probe_id] = [_beta(parts[i]) for i in idxs]
            n_rows += 1
            if n_rows % 100_000 == 0:
                print(f"  streamed {n_rows:,} probes", flush=True)

    return n_rows, by_probe


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--beta-matrix", required=True, type=Path)
    p.add_argument("--breast-samples", required=True, type=Path)
    p.add_argument("--hm450-probes", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    breast_ids = set(args.breast_samples.read_text().split())
    print(f"read {len(breast_ids)} breast sample IDs", flush=True)

    n_rows, tumor_by = read_beta_matrix(args.beta_matrix, breast_ids)
    print(f"collected β arrays for {len(tumor_by):,} probes", flush=True)

    print("summarizing tumor arm (52 breast cell lines)...", flush=True)
    t_summary = summarize(tumor_by)
    n_tumor = write_summary_tsv(args.output_dir / "tumor_summary.tsv", t_summary)

    shutil.copy(args.hm450_probes, args.output_dir / "probes.tsv")

    prov = args.output_dir / "PROVENANCE.md"
    prov.write_text(
        "# GSE68379 breast cell-line tumor cohort\n\n"
        f"- β matrix: {args.beta_matrix}\n"
        f"- tumor samples: {len(breast_ids)} breast cell lines "
        "(primary_site=breast per GSE68379 metadata)\n"
        f"- probes with summaries: {n_tumor:,}\n"
        f"- platform: HM450 (GPL13534)\n\n"
        "## Normal arm note\n\n"
        "GSE68379 has no in-study normal. Pair the emitted `tumor_summary.tsv`\n"
        "with an external healthy-normal `normal_summary.tsv` (typically\n"
        "`data/derived/gse69914_cohort/normal_summary.tsv`, n=50 healthy\n"
        "donor breast on the same HM450 platform).\n\n"
        "Cross-series composition is a known caveat — Sanger (GSE68379) and\n"
        "the external normal cohort were processed by different labs with\n"
        "different scanners. Treat AUC and tie-band behavior as the primary\n"
        "signal; interpret absolute β values with caution.\n"
    )

    print(f"→ {args.output_dir}/tumor_summary.tsv  ({n_tumor:,} probes)", flush=True)
    print(f"  next: copy a normal_summary.tsv from an external HM450 normal cohort.", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
