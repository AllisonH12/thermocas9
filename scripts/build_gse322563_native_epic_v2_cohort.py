"""Build GSE322563 against the NATIVE EPIC v2 probe set — no HM450 intersect.

Sibling to `scripts/build_gse322563_cohort.py`. The sibling stripped
each probe's `_BC##`/`_TC##`/`_TO##`/`_BO##` beadchip-design suffix and
intersected with the HM450 universe (~80.7 % retention). This script
keeps every probe's native EPIC v2 identifier (suffix intact) and emits
summaries for every probe in the EPIC v2 hg19 annotation
(`data/raw/probes_hg19_epic_v2.tsv`, produced by
`scripts/build_epic_v2_probes.py`).

Outputs (LocalSummaryBackend-compatible):
    <out>/tumor_summary.tsv         (MCF-7 replicates, n=2, native probe IDs)
    <out>/normal_summary.tsv        (MCF-10A replicates, n=2, native probe IDs)
    <out>/probes.tsv                (copy of probes_hg19_epic_v2.tsv)
    <out>/PROVENANCE.md

Same caveats as the sibling (n=2 per side, p_trust saturation, platform-
drift MCF-7 subline considerations); differs only in the probe-ID
preservation path.
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_surrogate_cohort import summarize, write_summary_tsv  # noqa: E402


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


def read_native_beta_matrix(
    path: Path, tumor_cols: list[str], normal_cols: list[str],
    probe_universe: set[str],
) -> tuple[int, dict[str, list[float | None]], dict[str, list[float | None]]]:
    """Stream the GSE322563 β matrix with native EPIC v2 probe IDs.

    Unlike the sibling HM450-intersect ingest, this function does NOT
    strip beadchip suffixes. Each emitted row keeps the full native
    probe_id (e.g. `cg00381604_BC11`). We retain only rows whose
    probe_id is in `probe_universe` (the hg19-lifted EPIC v2 probes).
    """
    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    tumor_idx = [header.index(c) for c in tumor_cols]
    normal_idx = [header.index(c) for c in normal_cols]

    tumor_by: dict[str, list[float | None]] = {}
    normal_by: dict[str, list[float | None]] = {}
    n_rows = 0
    n_in_universe = 0

    with gzip.open(path, "rt") as f:
        f.readline()  # consume header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            probe_id = parts[0]
            n_rows += 1
            if probe_id not in probe_universe:
                continue
            tumor_by[probe_id] = [_beta(parts[i]) for i in tumor_idx]
            normal_by[probe_id] = [_beta(parts[i]) for i in normal_idx]
            n_in_universe += 1
    return n_rows, tumor_by, normal_by, n_in_universe


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--beta-matrix", required=True, type=Path,
                   help="GSE322563_beta_matrix_EPIC_v2.txt.gz")
    p.add_argument("--tumor-cols", nargs="+", default=["MCF71", "MCF72"])
    p.add_argument("--normal-cols", nargs="+", default=["MCF10A1", "MCF10A2"])
    p.add_argument("--epic-v2-probes", required=True, type=Path,
                   help="probes_hg19_epic_v2.tsv (native EPIC v2 IDs + hg19 pos)")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load the native EPIC v2 probe universe.
    universe: set[str] = set()
    with args.epic_v2_probes.open() as f:
        next(f)
        for ln in f:
            universe.add(ln.rstrip("\n").split("\t", 1)[0])
    print(f"EPIC v2 native probe universe: {len(universe):,}", flush=True)

    print(f"streaming {args.beta_matrix.name}...", flush=True)
    n_rows, tumor_by, normal_by, n_in_u = read_native_beta_matrix(
        args.beta_matrix, args.tumor_cols, args.normal_cols, universe,
    )
    print(f"  rows streamed:          {n_rows:,}", flush=True)
    print(f"  rows in native universe: {n_in_u:,}  "
          f"({100 * n_in_u / len(universe):.1f}% of universe)", flush=True)

    print("summarizing tumor arm...", flush=True)
    t_summary = summarize(tumor_by)
    n_tumor = write_summary_tsv(args.output_dir / "tumor_summary.tsv", t_summary)

    print("summarizing normal arm...", flush=True)
    n_summary = summarize(normal_by)
    n_normal = write_summary_tsv(args.output_dir / "normal_summary.tsv", n_summary)

    shutil.copy(args.epic_v2_probes, args.output_dir / "probes.tsv")

    (args.output_dir / "PROVENANCE.md").write_text(
        "# GSE322563 native EPIC v2 cohort build\n\n"
        f"- β matrix: {args.beta_matrix}\n"
        f"- tumor arm:  {', '.join(args.tumor_cols)}  (n={len(args.tumor_cols)})\n"
        f"- normal arm: {', '.join(args.normal_cols)}  (n={len(args.normal_cols)})\n"
        f"- probe universe: {args.epic_v2_probes} (EPIC v2, hg19-lifted)\n"
        f"- rows streamed from β matrix: {n_rows:,}\n"
        f"- rows kept (probe in native universe): {n_in_u:,}\n"
        "- probe IDs are kept in canonical EPIC v2 form with `_BC##` / `_TC##` /\n"
        "  `_TO##` / `_BO##` suffixes preserved (no HM450 intersect).\n"
        f"- probes with summaries: tumor={n_tumor:,}, normal={n_normal:,}\n"
    )
    print(f"→ {args.output_dir}/", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
