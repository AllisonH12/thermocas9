#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2.0", "scipy>=1.13"]
# ///
"""Cohort-aware driver: extract probe × sample β matrix + run limma-eBayes.

Per-cohort adapters (GSE322563 beta matrix; GSE77348 U/M intensities;
GSE69914 series matrix) → common (probes, samples, β-matrix, group
vector) form → `limma_ebayes.limma_ebayes(...)` → per-probe moderated-t
TSV.

Output: `data/derived/limma_{cohort_slug}_probes.tsv` with columns
`probe_id, delta_beta, s_sq, s_tilde_sq, t_mod, df_tot, p_value`. Plus
a `.meta.json` companion with the fitted prior parameters.

Usage:
  uv run scripts/run_limma_per_cohort.py --cohort gse322563
  uv run scripts/run_limma_per_cohort.py --cohort gse77348
  uv run scripts/run_limma_per_cohort.py --cohort gse69914
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from limma_ebayes import limma_ebayes  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
RAW = REPO / "data" / "raw"
DERIVED = REPO / "data" / "derived"


# ---------- adapters ----------


def _open_gz(path: Path):
    return gzip.open(path, "rt")


def extract_gse322563() -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
    """Roth MCF-7/MCF-10A (n=2/2 EPIC v2). Output groups: MCF7 = 1 (tumor), MCF10A = 0.

    Source: GSE322563_beta_matrix_EPIC_v2.txt.gz (probe_id + 4 sample columns).
    """
    matrix = RAW / "gse322563" / "GSE322563_beta_matrix_EPIC_v2.txt.gz"
    with _open_gz(matrix) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        probes: list[str] = []
        rows: list[list[float]] = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            probes.append(parts[0])
            vals = [np.nan if v in ("", "NA", "nan", "NaN") else float(v) for v in parts[1:]]
            rows.append(vals)
    betas = np.asarray(rows, dtype=np.float64)
    # Samples expected: MCF71, MCF72, MCF10A1, MCF10A2
    group = np.array([1 if s.upper().startswith("MCF7") and not s.upper().startswith("MCF10A") else 0 for s in samples])
    return probes, samples, betas, group


def extract_gse77348() -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
    """MCF-7/MCF-10A surrogate, restricted to untreated replicates only
    (3 MCF7 + 3 MCF10A) — matches `build_surrogate_cohort.py`'s filter.

    Beta = M / (M + U + 100). Two signal-intensity files, one per cell line.
    """
    import re

    mcf7_path = RAW / "surrogate" / "GSE77348_MCF7_SignalIntensities.txt.gz"
    mcf10a_path = RAW / "surrogate" / "GSE77348_MCF10A_SignalIntensities.txt.gz"

    def _read_signal(path: Path, cell_prefix: str) -> tuple[list[str], list[str], np.ndarray]:
        """Keep only untreated replicates — column names `{prefix}-<digit>_…`
        (case-insensitive), matching the cohort-builder heuristic."""
        pattern = re.compile(rf"^{re.escape(cell_prefix)}-\d+_\S+_Unmethylated signal$", re.IGNORECASE)
        with _open_gz(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            # Column triplets: (Unmethylated, Methylated, Detection Pval); keep indices of
            # untreated Unmethylated-signal columns.
            untreated_u_cols = [i for i, h in enumerate(header) if pattern.match(h)]
            # Sample name: strip "_Unmethylated signal" suffix.
            samples: list[str] = []
            for i in untreated_u_cols:
                samples.append(re.sub(r"_Unmethylated signal$", "", header[i], flags=re.IGNORECASE))
            probes: list[str] = []
            betas: list[list[float]] = []
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                probes.append(parts[0])
                row_betas = []
                for u_idx in untreated_u_cols:
                    m_idx = u_idx + 1
                    try:
                        u = float(parts[u_idx])
                        m = float(parts[m_idx])
                    except (ValueError, IndexError):
                        row_betas.append(np.nan); continue
                    denom = m + u + 100.0
                    row_betas.append(m / denom if denom > 0 else np.nan)
                betas.append(row_betas)
        return probes, samples, np.asarray(betas, dtype=np.float64)

    probes7, samples7, betas7 = _read_signal(mcf7_path, "MCF7")
    probes10, samples10, betas10 = _read_signal(mcf10a_path, "MCF10A")
    assert probes7 == probes10, "probe order mismatch between MCF7 and MCF10A intensity files"
    samples = samples7 + samples10
    betas = np.concatenate([betas7, betas10], axis=1)
    group = np.array([1] * len(samples7) + [0] * len(samples10))
    return probes7, samples, betas, group


def extract_gse69914() -> tuple[list[str], list[str], np.ndarray, np.ndarray]:
    """GSE69914: parse the full series matrix, keep only tumor (status=2) + normal (status=0) samples."""
    series = RAW / "gse69914" / "GSE69914_series_matrix.txt.gz"
    tumors = set((RAW / "gse69914" / "samples_tumor.txt").read_text().split())
    normals = set((RAW / "gse69914" / "samples_normal.txt").read_text().split())
    keep = tumors | normals
    header: list[str] | None = None
    samples_kept: list[str] = []
    keep_cols: list[int] | None = None
    probes: list[str] = []
    rows: list[list[float]] = []
    with _open_gz(series) as fh:
        in_data = False
        for line in fh:
            if not in_data:
                if line.startswith("!series_matrix_table_begin"):
                    in_data = True
                    header_line = fh.readline().rstrip("\n").split("\t")
                    header = [h.strip('"') for h in header_line]
                    keep_cols = [i for i, h in enumerate(header) if h in keep]
                    samples_kept = [header[i] for i in keep_cols]
                    print(f"  GSE69914: header has {len(header)-1} samples; kept {len(samples_kept)}", flush=True)
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            probes.append(parts[0].strip('"'))
            vals = [np.nan if parts[i] in ("", "NA", '""', "null") else float(parts[i].strip('"')) for i in keep_cols]
            rows.append(vals)
            if len(rows) % 100_000 == 0:
                print(f"  GSE69914: read {len(rows):,} probes...", flush=True)
    betas = np.asarray(rows, dtype=np.float64)
    group = np.array([1 if s in tumors else 0 for s in samples_kept])
    return probes, samples_kept, betas, group


_ADAPTERS = {
    "gse322563": extract_gse322563,
    "gse77348": extract_gse77348,
    "gse69914": extract_gse69914,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(_ADAPTERS.keys()))
    ap.add_argument("--output", type=Path,
                    help="Per-probe TSV (default: data/derived/limma_{cohort}_probes.tsv)")
    args = ap.parse_args()

    adapter = _ADAPTERS[args.cohort]
    out = args.output or (DERIVED / f"limma_{args.cohort}_probes.tsv")

    print(f"extracting {args.cohort} ...", flush=True)
    probes, samples, betas, group = adapter()
    print(f"  shape: {betas.shape}; group=1: {int((group==1).sum())}, group=0: {int((group==0).sum())}", flush=True)

    print(f"fitting limma-eBayes ...", flush=True)
    r = limma_ebayes(betas, group)
    print(f"  d = {r['d']}; prior: s0² = {r['s0_sq']:.4g}, d0 = {r['d0']!r}", flush=True)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        fh.write("probe_id\tdelta_beta\ts_sq\ts_tilde_sq\tt_mod\tdf_tot\tp_value\n")
        for i, pid in enumerate(probes):
            fh.write(
                f"{pid}\t"
                f"{r['delta'][i]:.6g}\t"
                f"{r['s_sq'][i]:.6g}\t"
                f"{r['s_tilde_sq'][i]:.6g}\t"
                f"{r['t_mod'][i]:.6g}\t"
                f"{r['df_tot'][i]:.6g}\t"
                f"{r['p_value'][i]:.6g}\n"
            )
    print(f"wrote {out}", flush=True)

    meta = out.with_suffix(".meta.json")
    meta.write_text(json.dumps({
        "cohort": args.cohort,
        "n_probes": int(len(probes)),
        "n_samples_used": int(len(samples)),
        "n0": int((group == 0).sum()),
        "n1": int((group == 1).sum()),
        "d": r["d"],
        "s0_sq": r["s0_sq"],
        "d0": None if np.isinf(r["d0"]) else r["d0"],
        "d0_is_infinite": bool(np.isinf(r["d0"])),
    }, indent=2))
    print(f"wrote {meta}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
