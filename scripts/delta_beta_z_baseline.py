#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Δβ_z = (β_n − β_t) / sqrt(σ_t² + σ_n²) — uncertainty-aware Δβ baseline.

A reviewer concern with the §4.3 axis grid is that the four shipped
axes (Δβ-only, V1, V2 tumor_only, V2.5) do not include any standard
effect-size-with-uncertainty ranker. Δβ_z is the simplest one
constructible from the same inputs as V2.5: it uses the same
IQR-derived σ on each side that V2.5 feeds into `p_diff`, but
combines them as a z-score of the cohort-level mean gap rather than
as a tail probability of the gap distribution.

The point of running this is to isolate what `p_targ` and `p_trust`
contribute to V2.5 over and above an effect-size-with-uncertainty
ranker built from the same per-side dispersion signal. If V2.5 only
matched Δβ_z, we would be carrying two extra factors for no
discriminative gain.

Implementation: streams each scored JSONL, computes Δβ_z per record,
re-ranks, reports AUC + per-positive ranks at the n=3 validated
positives. δ and σ_floor do not enter Δβ_z.

Output: `examples/delta_beta_z_baseline.tsv` + Markdown companion
`examples/delta_beta_z_baseline.md`.

Note on σ floor: we do **not** apply V2.5's σ_floor to the Δβ_z
denominator. The σ_floor is a V2.5-specific implementation choice
that prevents `p_diff` from saturating; Δβ_z would be artificially
inflated near boundary β-values if we floored its denominator the
same way. Records with σ_t = σ_n = 0 are emitted with score = 0
(no signal to rank by, treated as ties at the floor). This is
documented in the script and is the relevant ablation direction.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "delta_beta_z_baseline.tsv"
OUT_MD = REPO / "examples" / "delta_beta_z_baseline.md"

_IQR_TO_STDEV = 1.349

COHORTS = [
    ("GSE322563 HM450",     "scored_gse322563_differential.jsonl"),
    ("GSE322563 native v2", "scored_gse322563_native_differential.jsonl"),
    ("GSE77348",            "scored_surrogate_differential.jsonl"),
    ("GSE69914",            "scored_gse69914_differential.jsonl"),
]

POS_GENE = {
    "chr5:38258943+:NNNNCGA":  "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA":  "GATA3",
}


def per_positive_ranks_and_auc(scored_path: Path, positives: list[str]):
    pos_set = set(positives)
    pos_scores: dict[str, float] = {}
    all_scores: list[tuple[float, str]] = []

    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            obs = d["observation"]
            cid = d["candidate"]["candidate_id"]
            mu_t = obs.get("beta_tumor_mean"); mu_n = obs.get("beta_normal_mean")
            q25_t, q75_t = obs.get("beta_tumor_q25"), obs.get("beta_tumor_q75")
            q25_n, q75_n = obs.get("beta_normal_q25"), obs.get("beta_normal_q75")
            if mu_t is None or mu_n is None:
                z = 0.0
            else:
                sigma_t = (q75_t - q25_t) / _IQR_TO_STDEV if (q25_t is not None and q75_t is not None) else 0.0
                sigma_n = (q75_n - q25_n) / _IQR_TO_STDEV if (q25_n is not None and q75_n is not None) else 0.0
                denom = math.sqrt(sigma_t * sigma_t + sigma_n * sigma_n)
                if denom == 0.0:
                    z = 0.0
                else:
                    z = (mu_n - mu_t) / denom
            all_scores.append((-z, cid))
            if cid in pos_set:
                pos_scores[cid] = z

    n_total = len(all_scores)
    all_scores.sort()
    score_at_pos_100 = -all_scores[99][0]
    tie_band_100 = sum(1 for neg_s, _ in all_scores if -neg_s == score_at_pos_100)

    rank_of: dict[str, int] = {}
    for rank_idx, (_, cid) in enumerate(all_scores, start=1):
        if cid in pos_set:
            rank_of[cid] = rank_idx
            if len(rank_of) == len(pos_set):
                break

    asc_scores = [-neg_s for neg_s, _ in all_scores]; asc_scores.reverse()
    asc_cids = [cid for _, cid in all_scores]; asc_cids.reverse()

    midrank_of: dict[str, float] = {}
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
                midrank_of[cid] = mid; pos_rem.remove(cid)
                if not pos_rem:
                    break
        i = j

    n_pos = len(positives); n_neg = n_total - n_pos
    sum_ranks = sum(midrank_of.get(cid, 0.0) for cid in positives)
    auc = (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    per_positive = [
        (cid, POS_GENE[cid], rank_of[cid], n_total, pos_scores.get(cid, 0.0))
        for cid in positives if cid in rank_of
    ]
    return auc, per_positive, tie_band_100


def main() -> None:
    positives = [line.strip() for line in POSITIVES.read_text().splitlines() if line.strip()]
    print(f"positives: {positives}")

    rows = []
    summary: dict[str, tuple[float, int]] = {}
    for cohort, fname in COHORTS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP missing: {path}")
            continue
        print(f"  Δβ_z scoring {cohort} ...")
        auc, per_pos, tie100 = per_positive_ranks_and_auc(path, positives)
        summary[cohort] = (auc, tie100)
        for cid, gene, rank, n_total, score in per_pos:
            pct = 100.0 * (1.0 - rank / n_total)
            rows.append((cohort, gene, cid, rank, n_total, pct, score))
        print(f"    AUC = {auc:.4f}; tie@100 = {tie100:,}")

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\tgene\tcandidate_id\trank\tn_total\tpercentile\tdelta_beta_z\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {OUT_TSV} ({len(rows)} rows)")

    md = ["# Δβ_z baseline — uncertainty-aware Δβ ranker",
          "",
          "Generated by `scripts/delta_beta_z_baseline.py`.",
          "Δβ_z = (β_n − β_t) / sqrt(σ_t² + σ_n²) using the same",
          "IQR-derived σ as V2.5's `p_diff` denominator. No σ floor",
          "(see script header for rationale). Records with",
          "σ_t = σ_n = 0 are scored 0.",
          "",
          "## Δβ_z AUC and tie-band@100 at the n = 3 validated positives",
          "",
          "| cohort | Δβ_z AUC | Δβ_z tie_band@100 |",
          "|---|---:|---:|"]
    for cohort, _ in COHORTS:
        if cohort in summary:
            auc, tie100 = summary[cohort]
            md.append(f"| **{cohort}** | {auc:.3f} | {tie100:,} |")
    md.extend(["",
               "## Per-positive Δβ_z ranks",
               "",
               "| cohort | gene | rank | N | percentile | Δβ_z |",
               "|---|---|---:|---:|---:|---:|"])
    for cohort, gene, cid, rank, n_total, pct, score in rows:
        md.append(f"| {cohort} | *{gene}* | {rank:,} | {n_total:,} | {pct:.3f}% | {score:.3f} |")
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
