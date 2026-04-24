#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Genome-wide tissue gating experiment: shipped V2.5 vs gap_sigmoid.

Reads a scored JSONL from `thermocas score-cohort` with
`probabilistic_mode: tumor_plus_differential_protection` on the
genome-wide probe-window catalog, and on the chr5/6/10 subset of the
same records. For each axis:

  - shipped V2.5     (the as-scored `p_therapeutic_selectivity` field)
  - gap_sigmoid @ σ_fixed in {0.05, 0.0707, 0.10}
    → `p_targ × sigmoid((Δβ − δ) / σ_fixed) × p_trust`

computes:

  1. AUC at the n = 3 Roth-validated positives (mid-rank Mann-Whitney).
  2. tie_band_size_at_k at K ∈ {20, 100, 1000}.
  3. precision_at_k interval [min, max] at K = 100 (the reviewer-tier test).
  4. per-positive rank under the benchmark's `candidate_id` ascending
     tie-break (1-based). Reports both the genome-wide rank and the
     chr5/6/10-subset rank so the chr-subset / WG robustness is visible.

The scored JSONL is the shipped V2.5 run — we recompute gap_sigmoid
factors directly from the `observation` block (β means + quartiles) on
the fly, keeping `p_targ` and `p_trust` as shipped. Same streaming
pattern as `sigma_floor_sweep.py`.

Inputs:
  --scored PATH         scored JSONL (genome-wide V2.5)
  --positives PATH      positives_roth_validated.txt
  --output-prefix PATH  writes PREFIX.md + PREFIX.tsv

Output: tissue-gating table for §5.2.1 update + tissue decision.
"""

from __future__ import annotations

import argparse
import bisect
import json
import math
from pathlib import Path
from typing import Iterator

_IQR_TO_STDEV = 1.349
DEFAULT_DELTA = 0.2
SIGMA_FLOOR = 0.05
BANDWIDTHS = [0.05, math.sqrt(2) * SIGMA_FLOOR, 0.10]

# δ for the gap-factor recompute. Default 0.2 matches the shipped V2.5 default
# and the cohort YAMLs the V2.5 scoring used; pass --delta 0.1 to test the
# tissue-recommended setting (PAPER.md §5.3.2 / §5.2.2).

CHR_SUBSET = {"chr5", "chr6", "chr10"}

POS_GENE = {
    "chr5:38258943+:NNNNCGA":  "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA":  "GATA3",
}


def sigmoid(x: float) -> float:
    if x >= 0:
        e = math.exp(-x); return 1.0 / (1.0 + e)
    e = math.exp(x); return e / (1.0 + e)


def stream_records(scored_path: Path) -> Iterator[tuple[str, float, float | None, float | None, float | None, float | None, float | None, float | None, float, float, float]]:
    """Yield (cid, chrom_key, β_t, q25_t, q75_t, β_n, q25_n, q75_n, p_t, p_r, shipped_v25)."""
    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            cid = d["candidate"]["candidate_id"]
            chrom = d["candidate"]["chrom"]
            obs = d["observation"]; prob = d["probabilistic"]
            yield (
                cid, chrom,
                obs.get("beta_tumor_mean"),
                obs.get("beta_tumor_q25"), obs.get("beta_tumor_q75"),
                obs.get("beta_normal_mean"),
                obs.get("beta_normal_q25"), obs.get("beta_normal_q75"),
                prob["p_targetable_tumor"],
                prob["p_observation_trustworthy"],
                prob["p_therapeutic_selectivity"],
            )


def gap_sigmoid_score(mu_t, mu_n, sigma_fixed, delta=DEFAULT_DELTA) -> float:
    if mu_t is None or mu_n is None:
        return 0.0
    return sigmoid((mu_n - mu_t - delta) / sigma_fixed)


def compute_all(scored_path: Path, positives: set[str], delta: float = DEFAULT_DELTA) -> dict:
    """One pass; returns per-axis-and-subset score vectors + positive score lookup.

    `delta` is the δ used to recompute the gap_sigmoid axes only. Shipped V2.5
    score in the JSONL is at the YAML's δ at scoring time and is unaffected.
    """
    axes = ["shipped_v25"] + [f"gap_sigmoid_{b:.4f}" for b in BANDWIDTHS]
    # For each (axis, subset ∈ {"wg", "chr5_6_10"}): list of (-score, cid) for sort
    vectors: dict[tuple[str, str], list[tuple[float, str]]] = {
        (a, s): [] for a in axes for s in ("wg", "chr5_6_10")
    }
    pos_scores: dict[tuple[str, str], dict[str, float]] = {
        (a, s): {} for a in axes for s in ("wg", "chr5_6_10")
    }

    n_total_wg = 0; n_total_chr = 0
    for cid, chrom, mu_t, q25_t, q75_t, mu_n, q25_n, q75_n, p_t, p_r, v25_s in stream_records(scored_path):
        n_total_wg += 1
        in_chr = chrom in CHR_SUBSET
        if in_chr:
            n_total_chr += 1
        # shipped V2.5: use stored p_therapeutic_selectivity directly.
        for axis_name, score in (
            ("shipped_v25", v25_s),
            *(
                (f"gap_sigmoid_{b:.4f}",
                 p_t * gap_sigmoid_score(mu_t, mu_n, b, delta=delta) * p_r)
                for b in BANDWIDTHS
            ),
        ):
            vectors[(axis_name, "wg")].append((-score, cid))
            if cid in positives:
                pos_scores[(axis_name, "wg")][cid] = score
            if in_chr:
                vectors[(axis_name, "chr5_6_10")].append((-score, cid))
                if cid in positives:
                    pos_scores[(axis_name, "chr5_6_10")][cid] = score

    return {"axes": axes, "vectors": vectors, "pos_scores": pos_scores,
            "n_total_wg": n_total_wg, "n_total_chr": n_total_chr}


def midrank_auc(sorted_desc: list[tuple[float, str]], positives: set[str]) -> tuple[float, int]:
    """AUC under mid-rank tie handling + n_neg."""
    n_total = len(sorted_desc)
    n_pos = sum(1 for _, c in sorted_desc if c in positives)
    n_neg = n_total - n_pos
    asc_scores = [-k for k, _ in sorted_desc]; asc_scores.reverse()
    asc_cids = [c for _, c in sorted_desc]; asc_cids.reverse()
    midrank: dict[str, float] = {}
    pos_rem = set(positives)
    i = 0
    while i < n_total and pos_rem:
        s = asc_scores[i]; j = i
        while j < n_total and asc_scores[j] == s:
            j += 1
        mid = (i + 1 + j) / 2.0
        for k in range(i, j):
            cid = asc_cids[k]
            if cid in pos_rem:
                midrank[cid] = mid; pos_rem.remove(cid)
                if not pos_rem:
                    break
        i = j
    sum_ranks = sum(midrank.values())
    if n_pos == 0 or n_neg == 0:
        return float("nan"), n_neg
    return (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg), n_neg


def rank_of_positive(sorted_desc: list[tuple[float, str]], positive_cid: str) -> int | None:
    """1-based rank under (-score, candidate_id) ascending tie-break."""
    for idx, (_, cid) in enumerate(sorted_desc, start=1):
        if cid == positive_cid:
            return idx
    return None


def tie_band_at_k(sorted_desc: list[tuple[float, str]], k: int) -> int:
    """Records tied at the K-th position's score."""
    if k > len(sorted_desc):
        return 0
    cutoff = sorted_desc[k - 1][0]
    return sum(1 for s, _ in sorted_desc if s == cutoff)


def precision_at_k_interval(sorted_desc: list[tuple[float, str]], positives: set[str], k: int) -> tuple[float, float, float, int]:
    """(P@K_observed, P@K_min, P@K_max, tie_band)."""
    if k > len(sorted_desc):
        return (float("nan"), float("nan"), float("nan"), 0)
    cutoff = sorted_desc[k - 1][0]
    tied = [(s, c) for s, c in sorted_desc if s == cutoff]
    tie_band = len(tied)
    committed_top_k = sorted_desc[:k]
    committed_non_tied = [c for s, c in committed_top_k if s != cutoff]
    committed_pos = sum(1 for c in committed_non_tied if c in positives)
    k_band = sum(1 for s, c in committed_top_k if s == cutoff)
    band_pos = sum(1 for _, c in tied if c in positives)
    drop_slots = tie_band - k_band
    observed_top_k_pos = sum(1 for _, c in committed_top_k if c in positives)
    p_obs = observed_top_k_pos / k
    p_min = (committed_pos + max(0, band_pos - drop_slots)) / k
    p_max = (committed_pos + min(band_pos, k_band)) / k
    return p_obs, p_min, p_max, tie_band


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored", required=True, type=Path)
    ap.add_argument("--positives", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    ap.add_argument("--delta", type=float, default=DEFAULT_DELTA,
                    help=f"δ for gap-factor recompute (default {DEFAULT_DELTA}); "
                          "shipped V2.5 score in the JSONL is at the YAML's δ "
                          "and is unaffected.")
    args = ap.parse_args()

    positives = set(line.strip() for line in args.positives.read_text().splitlines() if line.strip())
    print(f"positives ({len(positives)}): {sorted(positives)}")
    print(f"scoring: {args.scored}")

    print(f"  δ for gap_sigmoid recompute: {args.delta}")
    r = compute_all(args.scored, positives, delta=args.delta)
    print(f"  N (WG)        : {r['n_total_wg']:,}")
    print(f"  N (chr5/6/10) : {r['n_total_chr']:,}")

    # Sort once per (axis, subset).
    for key in r["vectors"]:
        r["vectors"][key].sort()

    # Build rows.
    rows = []
    for axis in r["axes"]:
        for subset in ("wg", "chr5_6_10"):
            sd = r["vectors"][(axis, subset)]
            n_total = len(sd)
            auc, n_neg = midrank_auc(sd, positives)
            tb20 = tie_band_at_k(sd, 20)
            tb100 = tie_band_at_k(sd, 100)
            tb1000 = tie_band_at_k(sd, 1000)
            p100_obs, p100_min, p100_max, _ = precision_at_k_interval(sd, positives, 100)
            pos_ranks = {POS_GENE[cid]: rank_of_positive(sd, cid) for cid in positives}
            rows.append({
                "axis": axis, "subset": subset, "n_total": n_total, "n_neg": n_neg,
                "auc": auc,
                "tb20": tb20, "tb100": tb100, "tb1000": tb1000,
                "p100_obs": p100_obs, "p100_min": p100_min, "p100_max": p100_max,
                "pos_ranks": pos_ranks,
            })

    # TSV.
    tsv = args.output_prefix.with_suffix(".tsv")
    tsv.parent.mkdir(parents=True, exist_ok=True)
    with tsv.open("w") as fh:
        fh.write("axis\tsubset\tn_total\tn_neg\tauc\ttie_band@20\ttie_band@100\ttie_band@1000"
                 "\tP@100_observed\tP@100_min\tP@100_max"
                 "\tESR1_rank\tEGFLAM_rank\tGATA3_rank\n")
        for r_ in rows:
            pr = r_["pos_ranks"]
            fh.write(
                f"{r_['axis']}\t{r_['subset']}\t{r_['n_total']}\t{r_['n_neg']}\t{r_['auc']:.6f}"
                f"\t{r_['tb20']}\t{r_['tb100']}\t{r_['tb1000']}"
                f"\t{r_['p100_obs']:.4f}\t{r_['p100_min']:.4f}\t{r_['p100_max']:.4f}"
                f"\t{pr.get('ESR1', '')}\t{pr.get('EGFLAM', '')}\t{pr.get('GATA3', '')}\n"
            )
    print(f"wrote {tsv}")

    # Markdown companion.
    md = [f"# Genome-wide tissue gating: shipped V2.5 vs gap_sigmoid on {args.scored.name}",
          "",
          f"Generated by `scripts/genome_wide_tissue_gating.py`. δ = {args.delta} for gap_sigmoid recompute "
          f"(shipped V2.5 score in the JSONL is at the YAML's δ at scoring time).",
          f"Genome-wide N = {r['n_total_wg']:,} candidates; chr5/6/10 subset N = {r['n_total_chr']:,}.",
          "",
          "## AUC and tie-band structure at validated labels (n_pos = 3)",
          "",
          "| axis | subset | AUC | tie@20 | tie@100 | tie@1000 | P@100 [min, max] |",
          "|---|---|---:|---:|---:|---:|---|"]
    for r_ in rows:
        pr = r_["pos_ranks"]
        md.append(
            f"| {r_['axis']} | {r_['subset']} | {r_['auc']:.4f} "
            f"| {r_['tb20']:,} | {r_['tb100']:,} | {r_['tb1000']:,} "
            f"| [{r_['p100_min']:.3f}, {r_['p100_max']:.3f}] |"
        )

    md.extend([
        "",
        "## Per-positive ranks (1-based under `candidate_id` ascending tie-break)",
        "",
        "| axis | subset | ESR1 | EGFLAM | GATA3 |",
        "|---|---|---:|---:|---:|",
    ])
    for r_ in rows:
        pr = r_["pos_ranks"]
        md.append(
            f"| {r_['axis']} | {r_['subset']} | "
            f"{pr.get('ESR1', '—')} | {pr.get('EGFLAM', '—')} | {pr.get('GATA3', '—')} |"
        )

    md_path = args.output_prefix.with_suffix(".md")
    md_path.write_text("\n".join(md))
    print(f"wrote {md_path}")


if __name__ == "__main__":
    main()
