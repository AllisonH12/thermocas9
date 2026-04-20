"""Build a V2.5 higher-n validation cohort from GEO series GSE69914.

GSE69914 (Teschendorff et al.) is a HM450 methylation series with:
    n=305  breast cancer (non-BRCA1)       — status code 2
    n=50   healthy-donor breast normal      — status code 0
    n=42   tumor-adjacent normal            — status code 1   (excluded)
    n=7    BRCA1-carrier normal             — status code 3   (excluded)
    n=3    BRCA1-carrier cancer             — status code 4   (excluded)

We keep (status=2) as the tumor arm and (status=0) as the normal arm.
tumor-adjacent-normal is excluded because TCGA-BRCA already showed that
adjacent-normal bulk produces anti-predictive p_prot; we want clean
healthy-vs-tumor. BRCA1 arms are excluded because their methylation
biology is germline-driven, not representative of sporadic breast cancer.

With n=305 and n=50 on the two sides, p_trust should NOT saturate at the
low-n ceiling that afflicts the n=3 GSE77348 surrogate, so V2.5's
differential factor gets a real chance to drive ranking rather than
being trapped inside a tied band at the top score.

Inputs:
    --series-matrix   GSE69914_series_matrix.txt.gz  (from GEO FTP)
    --tumor-samples   one GSM per line (status=2 filter)
    --normal-samples  one GSM per line (status=0 filter)
    --probe-annotation  HM450 probes.tsv (copied through)

Outputs (LocalSummaryBackend-compatible):
    <out>/tumor_summary.tsv
    <out>/normal_summary.tsv
    <out>/probes.tsv
    <out>/PROVENANCE.md
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# Reuse helpers from the existing surrogate builder — identical semantics.
from build_surrogate_cohort import summarize, write_summary_tsv  # noqa: E402


def read_series_matrix(
    series_matrix: Path,
    tumor_ids: set[str],
    normal_ids: set[str],
) -> tuple[dict[str, list[float | None]], dict[str, list[float | None]]]:
    """Stream the GEO SeriesMatrix and collect β values for the two arms.

    The SeriesMatrix is gzipped; the probe table block sits between
    `!series_matrix_table_begin` and `!series_matrix_table_end`. The table
    header is `"ID_REF"` + per-sample GSM IDs (in quotes); each data row is
    `probe_id` + β values in the same column order.
    """

    tumor_by_probe: dict[str, list[float | None]] = {}
    normal_by_probe: dict[str, list[float | None]] = {}

    with gzip.open(series_matrix, "rt") as f:
        # 1. Fast-forward to the data table.
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
        # 2. Read the header row — strip surrounding quotes on each cell.
        header = f.readline().rstrip("\n").split("\t")
        header = [h.strip('"') for h in header]
        assert header[0] == "ID_REF", f"unexpected first header cell: {header[0]!r}"
        sample_ids = header[1:]

        tumor_idxs = [i for i, g in enumerate(sample_ids) if g in tumor_ids]
        normal_idxs = [i for i, g in enumerate(sample_ids) if g in normal_ids]
        if len(tumor_idxs) != len(tumor_ids):
            missing = tumor_ids - set(sample_ids)
            raise ValueError(f"{len(missing)} tumor samples not found in SeriesMatrix: {sorted(missing)[:5]}...")
        if len(normal_idxs) != len(normal_ids):
            missing = normal_ids - set(sample_ids)
            raise ValueError(f"{len(missing)} normal samples not found in SeriesMatrix: {sorted(missing)[:5]}...")
        print(f"  resolved {len(tumor_idxs)} tumor + {len(normal_idxs)} normal columns", flush=True)

        # 3. Stream rows, accumulate β arrays per probe.
        n_rows = 0
        for line in f:
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            probe_id = parts[0].strip('"')
            # β values are unquoted decimals; GEO uses "NA" / "null" / empty.
            row = parts[1:]
            def _beta(s: str) -> float | None:
                s = s.strip()
                if s in ("", "NA", "null", "nan"):
                    return None
                try:
                    v = float(s)
                except ValueError:
                    return None
                if 0.0 <= v <= 1.0:
                    return v
                return None

            tumor_by_probe[probe_id] = [_beta(row[i]) for i in tumor_idxs]
            normal_by_probe[probe_id] = [_beta(row[i]) for i in normal_idxs]
            n_rows += 1
            if n_rows % 100_000 == 0:
                print(f"  streamed {n_rows:,} probes", flush=True)

    return tumor_by_probe, normal_by_probe


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--series-matrix", required=True, type=Path)
    p.add_argument("--tumor-samples", required=True, type=Path)
    p.add_argument("--normal-samples", required=True, type=Path)
    p.add_argument("--probe-annotation", required=True, type=Path,
                   help="probes_hg19.tsv — copied through as probes.tsv")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    tumor_ids = set(args.tumor_samples.read_text().split())
    normal_ids = set(args.normal_samples.read_text().split())
    print(f"read {len(tumor_ids)} tumor + {len(normal_ids)} normal sample IDs", flush=True)

    print("streaming SeriesMatrix...", flush=True)
    tumor_by, normal_by = read_series_matrix(args.series_matrix, tumor_ids, normal_ids)
    print(f"collected β arrays for {len(tumor_by):,} probes", flush=True)

    print("summarizing tumor arm...", flush=True)
    tumor_summary = summarize(tumor_by)
    n_tumor = write_summary_tsv(args.output_dir / "tumor_summary.tsv", tumor_summary)

    print("summarizing normal arm...", flush=True)
    normal_summary = summarize(normal_by)
    n_normal = write_summary_tsv(args.output_dir / "normal_summary.tsv", normal_summary)

    shutil.copy(args.probe_annotation, args.output_dir / "probes.tsv")

    prov = args.output_dir / "PROVENANCE.md"
    prov.write_text(
        "# GSE69914 cohort build\n\n"
        f"- SeriesMatrix: {args.series_matrix}\n"
        f"- tumor samples (status=2, non-BRCA1): {len(tumor_ids)}\n"
        f"- normal samples (status=0, healthy donor): {len(normal_ids)}\n"
        f"- probes with summaries: tumor={n_tumor:,}, normal={n_normal:,}\n"
        f"- excluded: tumor-adjacent-normal (status=1), BRCA1 carrier arms (status=3,4)\n"
        "- purpose: higher-n matched tissue validation of V2.5 "
        "differential-protection mode (n-saturation test, not a Roth-comparable cohort)\n"
    )

    print(f"→ {args.output_dir}/  ({n_tumor:,} tumor probes, {n_normal:,} normal probes)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
