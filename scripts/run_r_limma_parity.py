#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2.0", "scipy>=1.13"]
# ///
"""Driver: run canonical R `limma::lmFit + eBayes` on the same (β,
group) inputs the Python `scripts/run_limma_per_cohort.py` consumed,
then compare R t-statistics vs the committed Python `data/derived/
limma_<cohort>_probes.tsv` outputs.

Pipeline:
  1. Re-extract (probes, samples, β, group) for the requested cohort
     using the same adapters in `scripts/run_limma_per_cohort.py`.
  2. Write β matrix + group vector to a temp TSV pair.
  3. Invoke `scripts/r_limma_parity.R` to produce per-probe R t-stats.
  4. Join R vs Python by probe_id, compute Spearman/Pearson,
     top-K overlap (K ∈ {100, 1000}), and produce both a TSV and
     a markdown summary under `examples/r_limma_parity_<cohort>.{tsv,md}`.

Used to address the third-party reviewer's R-parity ask
(memo-2026-04-22-az / Bioinformatics-revision request: confirm the
pure-Python Smyth (2004) implementation matches the canonical
limma::lmFit + eBayes math on real cohort inputs).

Usage:
  uv run python scripts/run_r_limma_parity.py --cohort gse69914
  uv run python scripts/run_r_limma_parity.py --cohort gse322563
  uv run python scripts/run_r_limma_parity.py --cohort gse77348
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"
DERIVED = REPO / "data" / "derived"

sys.path.insert(0, str(Path(__file__).parent))
from run_limma_per_cohort import _ADAPTERS  # noqa: E402


def _spearman_pearson(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    """Return (spearman, pearson) on finite-paired entries.

    No SciPy dependency — small per-probe arrays, NumPy-only. The
    finite-mask is critical for n = 2/2 cell-line cohorts where R emits
    NA for ~12,500 probes the eBayes fit cannot resolve; without it,
    NumPy's corrcoef returns NaN.
    """
    mask = np.isfinite(x) & np.isfinite(y)
    xs, ys = x[mask], y[mask]
    if xs.size < 3:
        return float("nan"), float("nan")
    pearson = float(np.corrcoef(xs, ys)[0, 1])
    rx = np.argsort(np.argsort(xs))
    ry = np.argsort(np.argsort(ys))
    spearman = float(np.corrcoef(rx, ry)[0, 1])
    return spearman, pearson


def _topk_overlap(x: np.ndarray, y: np.ndarray, k: int) -> float:
    """Jaccard overlap of top-K |t|-scoring probes between two arrays.

    NaN-safe: the two arrays are assumed pairwise-aligned (same probe
    indexing); we drop indices where either side is NaN before ranking.
    Without this, R's NA-coded underdetermined probes (d = 0 cohorts
    where eBayes drops 12,534 probes) propagate through argsort and
    collapse the overlap to zero.
    """
    n = min(x.size, y.size)
    mask = np.isfinite(x[:n]) & np.isfinite(y[:n])
    if mask.sum() < k:
        return float("nan")
    idx = np.flatnonzero(mask)
    xs = np.abs(x[idx])
    ys = np.abs(y[idx])
    tx = idx[np.argsort(xs)[-k:]]
    ty = idx[np.argsort(ys)[-k:]]
    return float(np.intersect1d(tx, ty).size / float(k))


def _read_python_limma(cohort: str) -> dict[str, dict[str, float]]:
    """{probe_id: {delta_beta, s_sq, s_tilde_sq, t_mod, df_tot, p_value}}."""
    path = DERIVED / f"limma_{cohort}_probes.tsv"
    out: dict[str, dict[str, float]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out[row["probe_id"]] = {
                "delta_beta": float(row["delta_beta"]),
                "s_sq": float(row["s_sq"]),
                "s_tilde_sq": float(row["s_tilde_sq"]),
                "t_mod": float(row["t_mod"]),
                "df_tot": float(row["df_tot"]),
                "p_value": float(row["p_value"]),
            }
    return out


def _read_r_limma(path: Path) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                out[row["probe_id"]] = {
                    "delta_beta": float(row["delta_beta"]),
                    "s_sq": float(row["s_sq"]),
                    "s_tilde_sq": float(row["s_tilde_sq"]),
                    "t_mod": float(row["t_mod"]),
                    "df_tot": float(row["df_tot"]),
                    "p_value": float(row["p_value"]),
                }
            except ValueError:
                continue
    return out


def _write_tsv_for_r(
    probes: list[str],
    samples: list[str],
    betas: np.ndarray,
    group: np.ndarray,
    out_dir: Path,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = out_dir / "betas.tsv"
    print(f"  writing {matrix_path} for R input ...")
    with matrix_path.open("w") as fh:
        fh.write("probe_id\t" + "\t".join(samples) + "\n")
        for i, pid in enumerate(probes):
            row = "\t".join(
                "NA" if not np.isfinite(v) else f"{v:.6g}" for v in betas[i]
            )
            fh.write(f"{pid}\t{row}\n")
    group_path = matrix_path.with_suffix(matrix_path.suffix + ".group")
    group_path.write_text("\n".join(str(int(g)) for g in group) + "\n")
    return matrix_path


def _summary_table(
    py: dict[str, dict[str, float]],
    r: dict[str, dict[str, float]],
) -> dict[str, float | int]:
    common = sorted(set(py) & set(r))
    py_t = np.array([py[k]["t_mod"] for k in common])
    r_t = np.array([r[k]["t_mod"] for k in common])
    py_d = np.array([py[k]["delta_beta"] for k in common])
    r_d = np.array([r[k]["delta_beta"] for k in common])
    py_s = np.array([py[k]["s_tilde_sq"] for k in common])
    r_s = np.array([r[k]["s_tilde_sq"] for k in common])

    sp_t, ps_t = _spearman_pearson(py_t, r_t)
    sp_d, ps_d = _spearman_pearson(py_d, r_d)
    sp_s, ps_s = _spearman_pearson(py_s, r_s)

    # Mask NaN pairs for the absolute-difference diagnostics — R emits NA
    # for probes the eBayes fit drops (n = 2/2 cell-line cohorts have
    # d = 0 residual df and ~12,500 underdetermined probes). Without
    # masking, max/median come back as NaN and obscure the real
    # convergence story.
    diff_mask = np.isfinite(py_t) & np.isfinite(r_t)
    diff_t = np.abs(py_t[diff_mask] - r_t[diff_mask])
    diff_d_mask = np.isfinite(py_d) & np.isfinite(r_d)
    diff_d = np.abs(py_d[diff_d_mask] - r_d[diff_d_mask])

    return {
        "n_python_probes": len(py),
        "n_r_probes": len(r),
        "n_common": len(common),
        "n_common_nonnan_t": int(diff_mask.sum()),
        "n_r_nan_t": int(np.isnan(r_t).sum()),
        "spearman_t_mod": sp_t,
        "pearson_t_mod": ps_t,
        "spearman_delta_beta": sp_d,
        "pearson_delta_beta": ps_d,
        "spearman_s_tilde_sq": sp_s,
        "pearson_s_tilde_sq": ps_s,
        "topk100_overlap_t": _topk_overlap(py_t, r_t, 100),
        "topk1000_overlap_t": _topk_overlap(py_t, r_t, 1000),
        "max_abs_diff_t_mod": float(np.max(diff_t)) if diff_t.size else float("nan"),
        "median_abs_diff_t_mod": float(np.median(diff_t)) if diff_t.size else float("nan"),
        "max_abs_diff_delta_beta": float(np.max(diff_d)) if diff_d.size else float("nan"),
    }


def _format_md(cohort: str, summary: dict[str, float | int]) -> str:
    return "\n".join([
        f"# canonical R `limma::lmFit + eBayes` parity — {cohort}",
        "",
        "Generated by `scripts/run_r_limma_parity.py`. Compares the",
        "pure-Python Smyth (2004) empirical-Bayes moderated-`t`",
        "implementation in `scripts/limma_ebayes.py` against the",
        "canonical R `limma::lmFit + limma::eBayes` pipeline on the",
        "same sample-level β matrix + group vector (extracted via the",
        "shared adapter in `scripts/run_limma_per_cohort.py`).",
        "",
        "## Probe-coverage parity",
        "",
        "| metric | value |",
        "|---|---:|",
        f"| Python probe count | {summary['n_python_probes']:,} |",
        f"| R probe count | {summary['n_r_probes']:,} |",
        f"| Probes in both | {summary['n_common']:,} |",
        f"| Common probes with finite Python and R `t_mod` | {summary['n_common_nonnan_t']:,} |",
        f"| R `t_mod` = NA (probes the eBayes fit dropped) | {summary['n_r_nan_t']:,} |",
        "",
        "## Per-probe statistic agreement",
        "",
        "| statistic | Spearman | Pearson |",
        "|---|---:|---:|",
        f"| `t_mod` (moderated t) | {summary['spearman_t_mod']:.6f} | {summary['pearson_t_mod']:.6f} |",
        f"| `delta_beta` | {summary['spearman_delta_beta']:.6f} | {summary['pearson_delta_beta']:.6f} |",
        f"| `s_tilde_sq` (eBayes posterior σ²) | {summary['spearman_s_tilde_sq']:.6f} | {summary['pearson_s_tilde_sq']:.6f} |",
        "",
        "| top-K |t_mod| Jaccard overlap |",
        "|---|---:|",
        f"| K = 100 | {summary['topk100_overlap_t']:.4f} |",
        f"| K = 1000 | {summary['topk1000_overlap_t']:.4f} |",
        "",
        "## Absolute-difference diagnostics",
        "",
        "| diagnostic | value |",
        "|---|---:|",
        f"| max |Python − R| `t_mod` | {summary['max_abs_diff_t_mod']:.6g} |",
        f"| median |Python − R| `t_mod` | {summary['median_abs_diff_t_mod']:.6g} |",
        f"| max |Python − R| `delta_beta` | {summary['max_abs_diff_delta_beta']:.6g} |",
        "",
        "## Interpretation",
        "",
        "- A Spearman `t_mod` correlation ≥ 0.999 confirms that the",
        "  Python implementation produces the same probe ordering as",
        "  canonical R limma to floating-point precision; downstream",
        "  candidate-mapped AUC and tie-band numbers in PAPER.md §4.3 /",
        "  §5.2.2 are therefore reproducible from the canonical pipeline.",
        "- Sign convention: `delta_beta = β_tumor − β_normal` for both",
        "  pipelines; downstream `limma_candidate_ranking.py` and",
        "  `evidence_class_stratified_benchmark.py` rank by `−t_mod` so",
        "  that targetable sites (β_tumor < β_normal → t < 0 → −t > 0)",
        "  rank high.",
        "",
        "## Reproduce",
        "",
        "```bash",
        f"uv run python scripts/run_r_limma_parity.py --cohort {cohort}",
        "```",
        "",
        "Requires R + Bioconductor `limma` (install via `BiocManager`).",
        "",
    ])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True, choices=sorted(_ADAPTERS.keys()))
    ap.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep the extracted β-matrix TSV in a stable location (~10 GB on full HM450 cohorts).",
    )
    args = ap.parse_args()

    cohort = args.cohort
    print(f"== R limma parity: {cohort} ==", flush=True)

    print(f"extracting (probes, samples, β, group) for {cohort} ...", flush=True)
    probes, samples, betas, group = _ADAPTERS[cohort]()
    print(f"  shape: {betas.shape}; group=1: {int((group==1).sum())}, group=0: {int((group==0).sum())}", flush=True)

    if args.keep_tmp:
        work_dir = REPO / "data" / "tmp" / f"r_limma_parity_{cohort}"
    else:
        tmp = tempfile.mkdtemp(prefix=f"r_limma_parity_{cohort}_")
        work_dir = Path(tmp)

    matrix_path = _write_tsv_for_r(probes, samples, betas, group, work_dir)
    r_out = work_dir / "r_limma_out.tsv"

    rscript = REPO / "scripts" / "r_limma_parity.R"
    print(f"invoking Rscript {rscript.name} ...", flush=True)
    cmd = ["Rscript", str(rscript), str(matrix_path), str(r_out)]
    proc = subprocess.run(cmd, check=True)
    print("R script exited cleanly")

    py = _read_python_limma(cohort)
    r = _read_r_limma(r_out)
    summary = _summary_table(py, r)

    out_md = EXAMPLES / f"r_limma_parity_{cohort}.md"
    out_tsv = EXAMPLES / f"r_limma_parity_{cohort}.tsv"
    out_md.write_text(_format_md(cohort, summary))
    with out_tsv.open("w") as fh:
        fh.write("metric\tvalue\n")
        for k, v in summary.items():
            if isinstance(v, float):
                fh.write(f"{k}\t{v:.6g}\n")
            else:
                fh.write(f"{k}\t{v}\n")
    print(f"wrote {out_md}\nwrote {out_tsv}")

    if not args.keep_tmp:
        shutil.rmtree(work_dir, ignore_errors=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
