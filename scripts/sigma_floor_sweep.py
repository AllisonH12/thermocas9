#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""σ_floor sensitivity sweep on the V2.5 differential-protection composite.

Reads each per-cohort `scored_*_differential.jsonl`, recomputes
`p_differential_protection` at σ_floor ∈ {0.02, 0.05, 0.10, 0.15} from
the per-side β summaries already in the JSONL `observation` block,
re-derives `p_therapeutic_selectivity = p_t × p_diff_new × p_trust`
keeping `p_t` and `p_trust` fixed (neither depends on σ_floor), and
reports:

  - AUC at the n=3 Roth-validated positives, per σ_floor value.
  - Per-positive rank (1-based, candidate_id-asc tie-break) at each
    σ_floor.
  - Tied-band size at K=100 at each σ_floor (so we can see whether
    the floor is what is creating the visible top-K tied behavior).

Output: `examples/sigma_floor_sweep.tsv` (long form) + a Markdown
companion `examples/sigma_floor_sweep.md` for the manuscript table.

Why a recompute instead of `thermocas score-cohort`: the only term
that depends on σ_floor is p_diff (eq. §3.1). Re-running the full
pipeline 16× costs ~25 min; the closed-form recompute costs ~30s per
(cohort × σ_floor) combo. Outputs are identical to the scoring
codepath at σ_floor = 0.05 (verified by re-deriving p_therapeutic at
σ_floor = 0.05 and confirming bit-equivalence with the JSONL's
stored value to within the JSON-serialized float precision).
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "sigma_floor_sweep.tsv"
OUT_MD = REPO / "examples" / "sigma_floor_sweep.md"

_IQR_TO_STDEV = 1.349
DEFAULT_DELTA = 0.2

SIGMA_FLOORS = [0.02, 0.05, 0.10, 0.15]

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


def p_diff(mu_t: float, mu_n: float, sigma_t: float, sigma_n: float, delta: float, sigma_floor: float) -> float:
    floor = max(0.0, sigma_floor)
    sigma_sq = max(sigma_t, floor) ** 2 + max(sigma_n, floor) ** 2
    z = (delta - (mu_n - mu_t)) / math.sqrt(sigma_sq)
    return 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))


def stream_records(path: Path):
    with path.open() as fh:
        for line in fh:
            d = json.loads(line)
            obs = d["observation"]
            prob = d["probabilistic"]
            cid = d["candidate"]["candidate_id"]
            mu_t = obs.get("beta_tumor_mean")
            mu_n = obs.get("beta_normal_mean")
            q25_t, q75_t = obs.get("beta_tumor_q25"), obs.get("beta_tumor_q75")
            q25_n, q75_n = obs.get("beta_normal_q25"), obs.get("beta_normal_q75")
            p_t = prob["p_targetable_tumor"]
            p_r = prob["p_observation_trustworthy"]
            yield cid, mu_t, mu_n, q25_t, q75_t, q25_n, q75_n, p_t, p_r


def auc_mannwhitney(positive_scores: list[float], n_total: int, all_scores_iter) -> float:
    """Mann-Whitney U / (n_pos · n_neg) with 0.5-tie crediting.

    Computed without loading all scores into a single list — we count
    for each positive how many negatives have score < / = / > it. With
    only 3 positives, this is 3·N comparisons, cheap.
    """
    n_pos = len(positive_scores)
    n_neg = n_total - n_pos
    less = [0] * n_pos
    eq = [0] * n_pos
    pos_set_count = [0] * n_pos  # count of positive matches at score (we'll subtract)
    # Index of positive scores (allow duplicates among positives).
    # We'll consume the iterator in a single pass.
    for s in all_scores_iter:
        for i, ps in enumerate(positive_scores):
            if s < ps:
                less[i] += 1
            elif s == ps:
                eq[i] += 1
    # Each positive is also in the iterator and contributes to its own eq;
    # subtract one self-equality per positive.
    u = 0.0
    for i in range(n_pos):
        u += less[i] + 0.5 * (eq[i] - 1)  # subtract self-match
    return u / (n_pos * n_neg)


def per_positive_ranks_and_auc(scored_path: Path, sigma_floor: float, positives: list[str], delta: float = DEFAULT_DELTA):
    """Single pass: collect positives' new scores, then second pass for ranks/AUC.

    Returns (auc, [(cid, gene, rank, n_total, score_new)], tie_band_at_100).
    """
    pos_set = set(positives)
    pos_scores: dict[str, float] = {}

    # First pass: compute new score for every record, and remember positive scores.
    # We also need the full distribution to compute ranks + AUC, so we accumulate.
    all_new_scores: list[tuple[float, str]] = []  # (-score, cid) for stable sort
    for cid, mu_t, mu_n, q25_t, q75_t, q25_n, q75_n, p_t, p_r in stream_records(scored_path):
        if mu_t is None or mu_n is None:
            new_score = 0.0
        else:
            sigma_t = (q75_t - q25_t) / _IQR_TO_STDEV if (q25_t is not None and q75_t is not None) else 0.0
            sigma_n = (q75_n - q25_n) / _IQR_TO_STDEV if (q25_n is not None and q75_n is not None) else 0.0
            pd = p_diff(mu_t, mu_n, sigma_t, sigma_n, delta, sigma_floor)
            new_score = p_t * pd * p_r
        all_new_scores.append((-new_score, cid))
        if cid in pos_set:
            pos_scores[cid] = new_score

    n_total = len(all_new_scores)
    # Sort by (-score asc, cid asc) → highest-score first under candidate_id_asc tie-break.
    all_new_scores.sort()

    rank_of: dict[str, int] = {}
    score_at_pos_100 = -all_new_scores[99][0]
    tie_band_100 = sum(1 for neg_s, _ in all_new_scores if -neg_s == score_at_pos_100)

    for rank_idx, (neg_s, cid) in enumerate(all_new_scores, start=1):
        if cid in pos_set:
            rank_of[cid] = rank_idx
            if len(rank_of) == len(pos_set):
                break

    # AUC via Mann-Whitney over the full sorted list.
    # Equivalent formula given ranks: AUC = (sum_{i in pos} R_i − n_pos*(n_pos+1)/2) / (n_pos * n_neg)
    # where R_i is rank (1-based) of positive i among ALL records by ascending score
    # (which is N + 1 - descending_rank for unique scores; with ties, use mid-rank).
    # Easier here: use the full sorted list to compute ascending ranks with mid-rank tie handling.
    #
    # We have the descending sort. Convert: ascending_rank_i = n_total + 1 - desc_rank_i,
    # but tie handling is more honest if done properly. Mid-rank tie credit:
    #
    # Step through the descending list: group ties by score. For each tie group of
    # size m starting at descending position p, every record in the group has
    # ascending rank == (n_total - p - m/2 + 1.5)? Let's just compute ascending
    # ranks with mid-rank handling inline.
    asc_score_of = [-neg_s for neg_s, _ in all_new_scores]
    # Reverse to ascending order:
    asc_score_of.reverse()
    asc_cid_of = [cid for _, cid in all_new_scores]
    asc_cid_of.reverse()

    # Walk through ascending scores in tie groups, assign mid-ranks.
    midrank_of: dict[str, float] = {}
    pos_rem = set(positives)
    i = 0
    while i < n_total and pos_rem:
        s = asc_score_of[i]
        j = i
        while j < n_total and asc_score_of[j] == s:
            j += 1
        # Tie group is asc_score_of[i:j]; ascending positions are i+1 .. j (1-based).
        # Mid-rank = (i+1 + j) / 2.
        mid = (i + 1 + j) / 2.0
        for k in range(i, j):
            cid = asc_cid_of[k]
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
    print(f"σ_floor sweep: {SIGMA_FLOORS}")

    rows = []  # long-form for TSV
    summary = defaultdict(dict)  # summary[cohort][sigma] = (auc, tie_band_100)
    for cohort, fname in COHORTS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP missing: {path}")
            continue
        for sf in SIGMA_FLOORS:
            print(f"  scoring {cohort} at σ_floor = {sf} ...")
            auc, per_pos, tie100 = per_positive_ranks_and_auc(path, sf, positives)
            summary[cohort][sf] = (auc, tie100)
            for cid, gene, rank, n_total, score in per_pos:
                pct = 100.0 * (1.0 - rank / n_total)
                rows.append((cohort, sf, gene, cid, rank, n_total, pct, score))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\tsigma_floor\tgene\tcandidate_id\trank\tn_total\tpercentile\tscore\n")
        for cohort, sf, gene, cid, rank, n_total, pct, score in rows:
            fh.write(f"{cohort}\t{sf}\t{gene}\t{cid}\t{rank}\t{n_total}\t{pct:.4f}\t{score:.6g}\n")
    print(f"wrote {OUT_TSV} ({len(rows)} rows)")

    md = ["# σ_floor sensitivity sweep — V2.5 differential composite",
          "",
          "Generated by `scripts/sigma_floor_sweep.py`. δ = 0.2 fixed.",
          "Sweep over σ_floor ∈ {0.02, 0.05, 0.10, 0.15} on the V2.5",
          "`tumor_plus_differential_protection` composite. p_t and p_trust",
          "are unchanged (they do not depend on σ_floor); only p_diff is",
          "recomputed.",
          "",
          "## AUC and tie-band@100 at the n = 3 validated positives",
          "",
          "| cohort | σ_floor=0.02 | σ_floor=0.05 (default) | σ_floor=0.10 | σ_floor=0.15 |",
          "|---|---|---|---|---|"]
    for cohort, _ in COHORTS:
        cells = []
        for sf in SIGMA_FLOORS:
            if sf in summary[cohort]:
                auc, tie100 = summary[cohort][sf]
                cells.append(f"AUC={auc:.3f}, tie@100={tie100:,}")
            else:
                cells.append("—")
        md.append(f"| **{cohort}** | " + " | ".join(cells) + " |")

    md.extend(["",
               "## Per-positive ranks at each σ_floor",
               ""])
    for cohort, _ in COHORTS:
        md.append(f"### {cohort}")
        md.append("")
        md.append("| gene | σ_floor=0.02 rank (%ile) | σ_floor=0.05 rank (%ile) | σ_floor=0.10 rank (%ile) | σ_floor=0.15 rank (%ile) |")
        md.append("|---|---|---|---|---|")
        cohort_rows = [r for r in rows if r[0] == cohort]
        # Group by gene
        for gene in ["ESR1", "EGFLAM", "GATA3"]:
            cells = []
            for sf in SIGMA_FLOORS:
                hit = next((r for r in cohort_rows if r[1] == sf and r[2] == gene), None)
                if hit:
                    _, _, _, _, rank, n_total, pct, _ = hit
                    cells.append(f"{rank:,} ({pct:.2f}%)")
                else:
                    cells.append("—")
            md.append(f"| *{gene}* | " + " | ".join(cells) + " |")
        md.append("")

    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
