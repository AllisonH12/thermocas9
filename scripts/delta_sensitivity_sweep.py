#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""δ sensitivity sweep on the V2.5 differential-protection composite.

δ is the only V2.5 hyperparameter selected against an empirical cohort
(GSE77348, the development cohort; §4.0). The sweep over
δ ∈ {0.1, 0.2, 0.3, 0.4, 0.5} on **every** cohort tests whether the
matched-cell-line AUC ranking is δ-stable, and whether the
development-cohort selection of δ = 0.2 is also reasonable on the
non-development cohorts.

Implementation parallels `sigma_floor_sweep.py`: streams each
`scored_*_differential.jsonl`, recomputes `p_differential_protection`
at the new δ from the in-record per-side β summaries (σ_floor held at
0.05), re-derives `p_therapeutic_selectivity = p_t × p_diff_new × p_trust`,
re-ranks, and reports AUC + per-positive ranks at the n=3 validated
positives.

Output: `examples/delta_sensitivity_sweep.tsv` + Markdown companion
`examples/delta_sensitivity_sweep.md`.

Reading the table:
  - **Development cohort (GSE77348)**: results at δ ≠ 0.2 are *post hoc*
    sensitivity, since δ = 0.2 was selected against this cohort.
  - **Non-development cohorts**: results across the δ sweep are honest
    sensitivity since δ was not selected against them.
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "delta_sensitivity_sweep.tsv"
OUT_MD = REPO / "examples" / "delta_sensitivity_sweep.md"

_IQR_TO_STDEV = 1.349
DEFAULT_SIGMA_FLOOR = 0.05

DELTAS = [0.1, 0.2, 0.3, 0.4, 0.5]

COHORTS = [
    ("GSE322563 HM450",     "scored_gse322563_differential.jsonl",        False),
    ("GSE322563 native v2", "scored_gse322563_native_differential.jsonl", False),
    ("GSE77348 (dev)",      "scored_surrogate_differential.jsonl",        True),
    ("GSE69914",            "scored_gse69914_differential.jsonl",         False),
]

POS_GENE = {
    "chr5:38258943+:NNNNCGA":  "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA":  "GATA3",
}


def p_diff(mu_t: float, mu_n: float, sigma_t: float, sigma_n: float, delta: float, sigma_floor: float) -> float:
    floor = max(0.0, sigma_floor)
    sigma_sq = max(sigma_t, floor) ** 2 + max(sigma_n, floor) ** 2
    z = (delta - (mu_n - mu_t)) / math.sqrt(sigma_sq)
    return 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))


def per_positive_ranks_and_auc(scored_path: Path, delta: float, positives: list[str]):
    pos_set = set(positives)
    pos_scores: dict[str, float] = {}
    all_new_scores: list[tuple[float, str]] = []

    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            obs = d["observation"]; prob = d["probabilistic"]
            cid = d["candidate"]["candidate_id"]
            mu_t = obs.get("beta_tumor_mean"); mu_n = obs.get("beta_normal_mean")
            q25_t, q75_t = obs.get("beta_tumor_q25"), obs.get("beta_tumor_q75")
            q25_n, q75_n = obs.get("beta_normal_q25"), obs.get("beta_normal_q75")
            p_t = prob["p_targetable_tumor"]; p_r = prob["p_observation_trustworthy"]
            if mu_t is None or mu_n is None:
                new_score = 0.0
            else:
                sigma_t = (q75_t - q25_t) / _IQR_TO_STDEV if (q25_t is not None and q75_t is not None) else 0.0
                sigma_n = (q75_n - q25_n) / _IQR_TO_STDEV if (q25_n is not None and q75_n is not None) else 0.0
                pd = p_diff(mu_t, mu_n, sigma_t, sigma_n, delta, DEFAULT_SIGMA_FLOOR)
                new_score = p_t * pd * p_r
            all_new_scores.append((-new_score, cid))
            if cid in pos_set:
                pos_scores[cid] = new_score

    n_total = len(all_new_scores)
    all_new_scores.sort()
    score_at_pos_100 = -all_new_scores[99][0]
    tie_band_100 = sum(1 for neg_s, _ in all_new_scores if -neg_s == score_at_pos_100)

    rank_of: dict[str, int] = {}
    for rank_idx, (_, cid) in enumerate(all_new_scores, start=1):
        if cid in pos_set:
            rank_of[cid] = rank_idx
            if len(rank_of) == len(pos_set):
                break

    # Mid-rank AUC.
    asc_scores = [-neg_s for neg_s, _ in all_new_scores]
    asc_scores.reverse()
    asc_cids = [cid for _, cid in all_new_scores]; asc_cids.reverse()

    midrank_of: dict[str, float] = {}
    pos_rem = set(positives)
    i = 0
    while i < n_total and pos_rem:
        s = asc_scores[i]
        j = i
        while j < n_total and asc_scores[j] == s:
            j += 1
        mid = (i + 1 + j) / 2.0
        for k in range(i, j):
            cid = asc_cids[k]
            if cid in pos_rem:
                midrank_of[cid] = mid
                pos_rem.remove(cid)
                if not pos_rem:
                    break
        i = j

    n_pos = len(positives)
    n_neg = n_total - n_pos
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
    print(f"δ sweep: {DELTAS}; σ_floor fixed = {DEFAULT_SIGMA_FLOOR}")

    rows = []
    summary = defaultdict(dict)
    for cohort, fname, is_dev in COHORTS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP missing: {path}")
            continue
        for delta in DELTAS:
            print(f"  scoring {cohort} at δ = {delta} ...")
            auc, per_pos, tie100 = per_positive_ranks_and_auc(path, delta, positives)
            summary[cohort][delta] = (auc, tie100)
            for cid, gene, rank, n_total, score in per_pos:
                pct = 100.0 * (1.0 - rank / n_total)
                rows.append((cohort, delta, gene, cid, rank, n_total, pct, score, is_dev))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\tdelta\tgene\tcandidate_id\trank\tn_total\tpercentile\tscore\tis_dev_cohort\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {OUT_TSV} ({len(rows)} rows)")

    md = ["# δ sensitivity sweep — V2.5 differential composite",
          "",
          "Generated by `scripts/delta_sensitivity_sweep.py`.",
          "σ_floor = 0.05 (default) fixed. δ swept over",
          "{0.1, 0.2, 0.3, 0.4, 0.5} on the V2.5 composite.",
          "",
          "**(dev)** marks GSE77348, the cohort against which δ = 0.2",
          "was originally selected; results at δ ≠ 0.2 on (dev) are",
          "*post hoc* sensitivity. Results on the other three cohorts",
          "are honest sensitivity (δ was not selected against them).",
          "",
          "## AUC at the n = 3 validated positives, per δ",
          "",
          "| cohort | δ=0.1 | δ=0.2 *(default)* | δ=0.3 | δ=0.4 | δ=0.5 |",
          "|---|---|---|---|---|---|"]
    for cohort, _, _ in COHORTS:
        cells = []
        for d in DELTAS:
            if d in summary[cohort]:
                auc, tie100 = summary[cohort][d]
                cells.append(f"{auc:.3f} (tie@100={tie100:,})")
            else:
                cells.append("—")
        md.append(f"| **{cohort}** | " + " | ".join(cells) + " |")

    md.extend(["",
               "## Per-positive ranks at each δ",
               ""])
    for cohort, _, _ in COHORTS:
        md.append(f"### {cohort}")
        md.append("")
        md.append("| gene | δ=0.1 rank (%ile) | δ=0.2 rank (%ile) | δ=0.3 rank (%ile) | δ=0.4 rank (%ile) | δ=0.5 rank (%ile) |")
        md.append("|---|---|---|---|---|---|")
        cohort_rows = [r for r in rows if r[0] == cohort]
        for gene in ["ESR1", "EGFLAM", "GATA3"]:
            cells = []
            for d in DELTAS:
                hit = next((r for r in cohort_rows if r[1] == d and r[2] == gene), None)
                if hit:
                    _, _, _, _, rank, n_total, pct, _, _ = hit
                    cells.append(f"{rank:,} ({pct:.2f}%)")
                else:
                    cells.append("—")
            md.append(f"| *{gene}* | " + " | ".join(cells) + " |")
        md.append("")

    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
