#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy"]
# ///
"""Whole-genome V2.5-sigmoid (delta, sigma_fixed) sweep.

This is the journal-hardening sweep for the recommended V2.5-sigmoid axis.
It recomputes

    p_targ * sigmoid((beta_normal - beta_tumor - delta) / sigma_fixed) * p_trust

from the frozen whole-genome scored JSONLs, keeping p_targ and p_trust as
emitted by the shipped V2.5-diff runs. The sweep is intentionally independent
of the older V2.5-diff sigma_floor sweep: delta and sigma_fixed are varied
directly for the sigmoid mode on the final whole-genome denominator.

Outputs:
  - examples/sigmoid_delta_sigma_wg_sweep.tsv
  - examples/sigmoid_delta_sigma_wg_sweep.md
  - examples/tissue_per_positive_wg_ranks.tsv
  - examples/tissue_per_positive_wg_ranks.md
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[1]
DERIVED = REPO / "data" / "derived"
EXAMPLES = REPO / "examples"

OUT_TSV = EXAMPLES / "sigmoid_delta_sigma_wg_sweep.tsv"
OUT_MD = EXAMPLES / "sigmoid_delta_sigma_wg_sweep.md"
TISSUE_TSV = EXAMPLES / "tissue_per_positive_wg_ranks.tsv"
TISSUE_MD = EXAMPLES / "tissue_per_positive_wg_ranks.md"

DELTAS = [0.10, 0.15, 0.20, 0.25, 0.30]
SIGMAS = [0.03, 0.05, math.sqrt(2) * 0.05, 0.10, 0.15, 0.20]
DEFAULT_DELTA = 0.20
DEFAULT_SIGMA = math.sqrt(2) * 0.05

POS_GENES: dict[bytes, str] = {
    b"chr5:38258943+:NNNNCGA": "EGFLAM",
    b"chr6:152011177+:NNNNCGA": "ESR1",
    b"chr10:8087387+:NNNNCGA": "GATA3",
}


@dataclass(frozen=True)
class Cohort:
    label: str
    path: Path
    n_total: int


COHORTS = [
    Cohort("GSE322563 HM450", DERIVED / "scored_gse322563_wg_differential.jsonl", 19_787_820),
    Cohort("GSE322563 native v2", DERIVED / "scored_gse322563_native_wg_differential.jsonl", 35_380_431),
    Cohort("GSE77348", DERIVED / "scored_surrogate_wg_differential.jsonl", 19_787_820),
    Cohort("GSE69914 (tissue)", DERIVED / "scored_gse69914_wg_differential.jsonl", 19_787_820),
]

CID_KEY = b'"candidate_id":"'


def _float_after(line: bytes, key: bytes) -> float:
    i = line.find(key)
    if i < 0:
        return math.nan
    start = i + len(key)
    end = start
    while end < len(line) and line[end] not in b",}":
        end += 1
    raw = line[start:end]
    if raw == b"null":
        return math.nan
    return float(raw)


def _candidate_id(line: bytes) -> bytes:
    start = line.find(CID_KEY)
    if start < 0:
        raise ValueError("candidate_id not found")
    start += len(CID_KEY)
    end = line.find(b'"', start)
    return line[start:end]


def load_numeric_vectors(cohort: Cohort) -> tuple[np.ndarray, np.ndarray, dict[str, int]]:
    """Return (gap, base, positive_indices) for one WG scored JSONL."""
    gap = np.empty(cohort.n_total, dtype=np.float64)
    base = np.empty(cohort.n_total, dtype=np.float64)
    positive_indices: dict[str, int] = {}

    with cohort.path.open("rb") as fh:
        for i, line in enumerate(fh):
            if i >= cohort.n_total:
                raise RuntimeError(f"{cohort.path} has more rows than expected {cohort.n_total:,}")
            bt = _float_after(line, b'"beta_tumor_mean":')
            bn = _float_after(line, b'"beta_normal_mean":')
            p_targ = _float_after(line, b'"p_targetable_tumor":')
            p_trust = _float_after(line, b'"p_observation_trustworthy":')
            if math.isnan(bt) or math.isnan(bn) or math.isnan(p_targ) or math.isnan(p_trust):
                gap[i] = 0.0
                base[i] = 0.0
            else:
                gap[i] = bn - bt
                base[i] = p_targ * p_trust

            cid = _candidate_id(line)
            gene = POS_GENES.get(cid)
            if gene is not None:
                positive_indices[gene] = i

            if (i + 1) % 2_000_000 == 0:
                print(f"  parsed {i + 1:,}/{cohort.n_total:,}", file=sys.stderr)

    if i + 1 != cohort.n_total:
        raise RuntimeError(f"{cohort.path} rows {i + 1:,} != expected {cohort.n_total:,}")
    missing = sorted(set(POS_GENES.values()) - set(positive_indices))
    if missing:
        raise RuntimeError(f"{cohort.label}: missing positives {missing}")
    return gap, base, positive_indices


def score_sigmoid_inplace(gap: np.ndarray, base: np.ndarray, delta: float, sigma: float, out: np.ndarray) -> np.ndarray:
    """Compute base * sigmoid((gap - delta) / sigma) into out."""
    np.subtract(gap, delta, out=out)
    out /= sigma
    np.negative(out, out=out)
    np.exp(out, out=out)
    out += 1.0
    np.reciprocal(out, out=out)
    out *= base
    return out


def auc_from_scores(scores: np.ndarray, pos_indices: dict[str, int]) -> float:
    pos_scores = np.array([scores[pos_indices[g]] for g in ("ESR1", "EGFLAM", "GATA3")], dtype=np.float64)
    n_pos = len(pos_scores)
    n_neg = scores.size - n_pos
    wins = 0.0
    for ps in pos_scores:
        neg_less = int(np.count_nonzero(scores < ps)) - int(np.count_nonzero(pos_scores < ps))
        neg_equal = int(np.count_nonzero(scores == ps)) - int(np.count_nonzero(pos_scores == ps))
        wins += neg_less + 0.5 * neg_equal
    return wins / (n_pos * n_neg)


def tie_band_at_100(scores: np.ndarray) -> int:
    cutoff = np.partition(scores, scores.size - 100)[scores.size - 100]
    return int(np.count_nonzero(scores == cutoff))


def exact_delta_beta_ranks(scored_path: Path, positive_gaps: dict[bytes, float]) -> dict[str, int]:
    """Exact (-delta_beta, candidate_id) ranks for GSE69914 tissue positives."""
    ranks = {gene: 1 for gene in POS_GENES.values()}
    ordered_positive_ids = sorted(positive_gaps)
    with scored_path.open("rb") as fh:
        for i, line in enumerate(fh, start=1):
            cid = _candidate_id(line)
            bt = _float_after(line, b'"beta_tumor_mean":')
            bn = _float_after(line, b'"beta_normal_mean":')
            gap = -math.inf if (math.isnan(bt) or math.isnan(bn)) else bn - bt
            for pcid in ordered_positive_ids:
                pgap = positive_gaps[pcid]
                if gap > pgap or (gap == pgap and cid < pcid):
                    ranks[POS_GENES[pcid]] += 1
            if i % 2_000_000 == 0:
                print(f"  delta-beta rank pass {i:,}", file=sys.stderr)
    return ranks


def percentile(rank: int, n_total: int) -> float:
    return 100.0 * (1.0 - rank / n_total)


def load_rank_artifact(path: Path, axis: str) -> dict[str, int]:
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["axis"] == axis and row["subset"] == "wg":
                return {
                    "ESR1": int(row["ESR1_rank"]),
                    "EGFLAM": int(row["EGFLAM_rank"]),
                    "GATA3": int(row["GATA3_rank"]),
                }
    raise RuntimeError(f"axis {axis!r} not found in {path}")


def load_feature_matched_p() -> dict[str, float]:
    out: dict[str, float] = {}
    with (EXAMPLES / "feature_matched_negative_controls.tsv").open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if row["cohort"] == "GSE69914 (tissue)" and row["scope"] == "wg" and row["axis"] == "v25_sigmoid":
                out[row["gene"]] = float(row["empirical_p"])
    return out


def write_tissue_table(delta_ranks: dict[str, int], n_total: int) -> None:
    gating = EXAMPLES / "gse69914_wg_gating.tsv"
    limma = EXAMPLES / "limma_gse69914_wg.tsv"
    v25_diff = load_rank_artifact(gating, "shipped_v25")
    v25_sigmoid = load_rank_artifact(gating, "gap_sigmoid_0.0707")
    limma_ranks = load_rank_artifact(limma, "limma_t")
    fm_p = load_feature_matched_p()

    rows = []
    for gene in ("ESR1", "EGFLAM", "GATA3"):
        rows.append({
            "gene": gene,
            "delta_rank": delta_ranks[gene],
            "delta_percentile": percentile(delta_ranks[gene], n_total),
            "limma_rank": limma_ranks[gene],
            "limma_percentile": percentile(limma_ranks[gene], n_total),
            "v25_diff_rank": v25_diff[gene],
            "v25_diff_percentile": percentile(v25_diff[gene], n_total),
            "v25_sigmoid_rank": v25_sigmoid[gene],
            "v25_sigmoid_percentile": percentile(v25_sigmoid[gene], n_total),
            "feature_matched_p": fm_p[gene],
        })

    with TISSUE_TSV.open("w") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    md = [
        "# GSE69914 tissue per-positive whole-genome ranks",
        "",
        "Generated by `scripts/sigmoid_delta_sigma_wg_sweep.py`. Ranks are 1-based",
        "against the frozen HM450 whole-genome denominator (N = 19,787,820).",
        "Percentile is `100 * (1 - rank / N)`, higher is better. The",
        "feature-matched p-value is the within-chromosome matched-negative audit",
        "for V2.5-sigmoid from `feature_matched_negative_controls.tsv`.",
        "",
        "| positive | Δβ-only WG %ile | limma WG %ile | V2.5-diff WG %ile | V2.5-sigmoid WG %ile | feature-matched p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        md.append(
            f"| *{row['gene']}* | "
            f"{row['delta_percentile']:.2f}% | {row['limma_percentile']:.2f}% | "
            f"{row['v25_diff_percentile']:.2f}% | {row['v25_sigmoid_percentile']:.2f}% | "
            f"{row['feature_matched_p']:.4f} |"
        )
    TISSUE_MD.write_text("\n".join(md) + "\n")


def format_sigma(sigma: float) -> str:
    return "0.0707" if abs(sigma - DEFAULT_SIGMA) < 1e-8 else f"{sigma:.2f}"


def write_sweep_markdown(rows: list[dict[str, object]]) -> None:
    md = [
        "# Whole-genome V2.5-sigmoid delta/sigma sweep",
        "",
        "Generated by `scripts/sigmoid_delta_sigma_wg_sweep.py` from the frozen",
        "whole-genome scored JSONLs. The recomputed axis is",
        "`p_targ * sigmoid((beta_normal - beta_tumor - delta) / sigma_fixed) * p_trust`.",
        "No hyperparameters are retuned by this sweep; it is a robustness check around",
        "the shipped default `delta = 0.20`, `sigma_fixed = sqrt(2) * 0.05 ~= 0.0707`.",
        "",
    ]
    by_cohort: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_cohort.setdefault(str(row["cohort"]), []).append(row)

    for cohort in [c.label for c in COHORTS]:
        md.extend([f"## {cohort}", "", "AUC:", ""])
        md.append("| delta | " + " | ".join(f"sigma={format_sigma(s)}" for s in SIGMAS) + " |")
        md.append("|---|" + "---:|" * len(SIGMAS))
        for delta in DELTAS:
            cells = []
            for sigma in SIGMAS:
                row = next(r for r in by_cohort[cohort] if r["delta"] == delta and r["sigma_fixed"] == sigma)
                val = f"{float(row['auc']):.3f}"
                if delta == DEFAULT_DELTA and abs(sigma - DEFAULT_SIGMA) < 1e-8:
                    val = f"**{val}**"
                cells.append(val)
            md.append(f"| {delta:.2f} | " + " | ".join(cells) + " |")
        md.extend(["", "`tie_band@100`:", ""])
        md.append("| delta | " + " | ".join(f"sigma={format_sigma(s)}" for s in SIGMAS) + " |")
        md.append("|---|" + "---:|" * len(SIGMAS))
        for delta in DELTAS:
            cells = []
            for sigma in SIGMAS:
                row = next(r for r in by_cohort[cohort] if r["delta"] == delta and r["sigma_fixed"] == sigma)
                val = f"{int(row['tie_band@100']):,}"
                if delta == DEFAULT_DELTA and abs(sigma - DEFAULT_SIGMA) < 1e-8:
                    val = f"**{val}**"
                cells.append(val)
            md.append(f"| {delta:.2f} | " + " | ".join(cells) + " |")
        md.append("")

    md.extend([
        "## Reading the sweep",
        "",
        "The desired robustness pattern is a broad AUC plateau around delta = 0.20",
        "and sigma_fixed ~= 0.05-0.10, with small `tie_band@100` values. The",
        "default need not be the global optimum; the reviewer-facing point is that",
        "it is not a pathological one-point setting.",
    ])
    OUT_MD.write_text("\n".join(md) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", action="append", choices=[c.label for c in COHORTS],
                    help="Restrict to one or more cohort labels.")
    args = ap.parse_args()

    selected = [c for c in COHORTS if not args.cohort or c.label in args.cohort]
    rows: list[dict[str, object]] = []
    tissue_delta_ranks: dict[str, int] | None = None

    for cohort in selected:
        print(f"Loading {cohort.label} from {cohort.path}", file=sys.stderr)
        gap, base, pos_indices = load_numeric_vectors(cohort)
        scores = np.empty_like(gap)
        for delta in DELTAS:
            for sigma in SIGMAS:
                print(f"  {cohort.label}: delta={delta:.2f}, sigma={sigma:.4f}", file=sys.stderr)
                score_sigmoid_inplace(gap, base, delta, sigma, scores)
                rows.append({
                    "cohort": cohort.label,
                    "delta": delta,
                    "sigma_fixed": sigma,
                    "auc": auc_from_scores(scores, pos_indices),
                    "tie_band@100": tie_band_at_100(scores),
                })

        if cohort.label == "GSE69914 (tissue)":
            positive_gaps = {
                pcid: float(gap[pos_indices[gene]])
                for pcid, gene in POS_GENES.items()
            }
            print("  computing exact GSE69914 delta-beta ranks", file=sys.stderr)
            tissue_delta_ranks = exact_delta_beta_ranks(cohort.path, positive_gaps)
            write_tissue_table(tissue_delta_ranks, cohort.n_total)

        del gap, base, scores

    if OUT_TSV.exists() and args.cohort:
        # Merge partial runs with existing rows, replacing selected cohorts.
        existing: list[dict[str, object]] = []
        with OUT_TSV.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            selected_labels = {c.label for c in selected}
            for row in reader:
                if row["cohort"] not in selected_labels:
                    existing.append({
                        "cohort": row["cohort"],
                        "delta": float(row["delta"]),
                        "sigma_fixed": float(row["sigma_fixed"]),
                        "auc": float(row["auc"]),
                        "tie_band@100": int(row["tie_band@100"]),
                    })
        rows = existing + rows

    rows.sort(key=lambda r: (str(r["cohort"]), float(r["delta"]), float(r["sigma_fixed"])))
    with OUT_TSV.open("w") as fh:
        writer = csv.DictWriter(fh, delimiter="\t",
                                fieldnames=["cohort", "delta", "sigma_fixed", "auc", "tie_band@100"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "cohort": row["cohort"],
                "delta": f"{float(row['delta']):.2f}",
                "sigma_fixed": f"{float(row['sigma_fixed']):.6f}",
                "auc": f"{float(row['auc']):.6f}",
                "tie_band@100": int(row["tie_band@100"]),
            })
    write_sweep_markdown(rows)

    if "GSE69914 (tissue)" in [c.label for c in selected] and tissue_delta_ranks is None:
        raise RuntimeError("selected tissue cohort but did not write tissue table")


if __name__ == "__main__":
    main()
