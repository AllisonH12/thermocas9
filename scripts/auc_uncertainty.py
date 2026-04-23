#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Descriptive AUC uncertainty intervals for the n=3 primary endpoint.

Given the small positive count, classical AUC confidence intervals are
not well-defined. We instead report two descriptive distributions per
cohort × scoring axis:

1. **Permutation null.** Draw n_pos = 3 records uniformly at random from
   the full candidate set, compute their AUC under the same tie-handling
   the benchmark uses (mid-rank Mann-Whitney). Repeat 10,000 times. The
   empirical *p*-value is `Pr(AUC_permuted ≥ AUC_observed)`. This
   answers: "how often would three randomly chosen candidates achieve
   AUC at least this high by chance?" With n_pos=3 the discrete null
   has limited resolution (≈10⁵–10⁶ unique outcomes given the tie
   structure), but the upper-tail probability is still meaningful.

2. **Negative-set bootstrap.** Resample the negative set with
   replacement to size `n_neg`, recompute AUC against the (fixed) three
   real positives, repeat 1,000 times. Report the 2.5–97.5 percentile
   spread. This is a *descriptive* spread of how much AUC moves under
   alternative draws of the negative pool, not a calibrated CI.

Both intervals are reported in the manuscript as descriptive — not as
inferential statistics — because n_pos = 3 makes inferential AUC CIs
unreliable in any framework.

Output: `examples/auc_uncertainty.tsv` + `examples/auc_uncertainty.md`.
"""

from __future__ import annotations

import bisect
import json
import math
import random
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "auc_uncertainty.tsv"
OUT_MD = REPO / "examples" / "auc_uncertainty.md"

N_PERM = 10_000
N_BOOT = 1_000

# (cohort_label, axis_label, scored_path, score_field)
RUNS = [
    ("GSE322563 HM450",     "V1 final_score", "scored_gse322563_differential.jsonl",        "final_score"),
    ("GSE322563 HM450",     "V2.5 diff",      "scored_gse322563_differential.jsonl",        "p_therapeutic"),
    ("GSE322563 HM450",     "Δβ-only",        "scored_gse322563_differential.jsonl",        "delta_beta"),
    ("GSE322563 native v2", "V1 final_score", "scored_gse322563_native_differential.jsonl", "final_score"),
    ("GSE322563 native v2", "V2.5 diff",      "scored_gse322563_native_differential.jsonl", "p_therapeutic"),
    ("GSE322563 native v2", "Δβ-only",        "scored_gse322563_native_differential.jsonl", "delta_beta"),
    ("GSE77348",            "V1 final_score", "scored_surrogate_differential.jsonl",        "final_score"),
    ("GSE77348",            "V2.5 diff",      "scored_surrogate_differential.jsonl",        "p_therapeutic"),
    ("GSE77348",            "Δβ-only",        "scored_surrogate_differential.jsonl",        "delta_beta"),
]


def score_of(record: dict, field: str) -> float:
    if field == "final_score":
        s = record.get("final_score")
    elif field == "p_therapeutic":
        s = record.get("probabilistic", {}).get("p_therapeutic_selectivity")
    elif field == "delta_beta":
        obs = record.get("observation", {})
        bt, bn = obs.get("beta_tumor_mean"), obs.get("beta_normal_mean")
        if bt is None or bn is None:
            return float("-inf")
        return bn - bt
    else:
        raise ValueError(field)
    return s if s is not None else float("-inf")


def load_scores(path: Path, field: str, positives: set[str]) -> tuple[list[float], list[float]]:
    """Return (all_scores_sorted_ascending, positive_scores)."""
    all_scores: list[float] = []
    pos_scores: dict[str, float] = {}
    with path.open() as fh:
        for line in fh:
            d = json.loads(line)
            cid = d["candidate"]["candidate_id"]
            s = score_of(d, field)
            all_scores.append(s)
            if cid in positives:
                pos_scores[cid] = s
    all_scores.sort()
    return all_scores, list(pos_scores.values())


def auc_from_indices(asc_scores: list[float], pos_indices: list[int]) -> float:
    """AUC under mid-rank tie handling, given pre-sorted ascending scores
    and 0-based positive *indices into the sorted array*.

    Mid-rank: a tie group at positions [i, j) gets mid-rank (i+1+j)/2.
    Sum of mid-ranks of positives, U = sum_ranks - n_pos*(n_pos+1)/2,
    AUC = U / (n_pos * n_neg).
    """
    n_total = len(asc_scores)
    n_pos = len(pos_indices)
    n_neg = n_total - n_pos
    sum_ranks = 0.0
    for idx in pos_indices:
        s = asc_scores[idx]
        # Find tie group containing idx.
        lo = bisect.bisect_left(asc_scores, s)
        hi = bisect.bisect_right(asc_scores, s)
        mid = (lo + 1 + hi) / 2.0
        sum_ranks += mid
    return (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def auc_from_scores(asc_scores: list[float], pos_scores: list[float]) -> float:
    """AUC where positives are identified by score (look up tie group via
    bisect). Equivalent to auc_from_indices when each positive sits at a
    unique location, with proper mid-rank for tie groups."""
    n_total = len(asc_scores)
    n_pos = len(pos_scores)
    n_neg = n_total - n_pos
    sum_ranks = 0.0
    for s in pos_scores:
        lo = bisect.bisect_left(asc_scores, s)
        hi = bisect.bisect_right(asc_scores, s)
        mid = (lo + 1 + hi) / 2.0
        sum_ranks += mid
    return (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def permutation_null_auc(asc_scores: list[float], n_pos: int, n_iter: int, observed_auc: float, seed: int = 0) -> tuple[float, float, float, float]:
    """Sample n_pos indices uniformly w/o replacement; compute AUC.
    Return (mean, std, p2.5, p97.5, p_value_one_sided)."""
    rng = random.Random(seed)
    n_total = len(asc_scores)
    aucs: list[float] = []
    n_ge = 0
    for _ in range(n_iter):
        idxs = rng.sample(range(n_total), n_pos)
        a = auc_from_indices(asc_scores, idxs)
        aucs.append(a)
        if a >= observed_auc:
            n_ge += 1
    aucs.sort()
    pval = (n_ge + 1) / (n_iter + 1)  # +1/+1 for unbiasedness
    p_lo = aucs[int(0.025 * n_iter)]
    p_hi = aucs[int(0.975 * n_iter)]
    return p_lo, p_hi, pval


def negative_bootstrap_auc(asc_scores: list[float], pos_scores: list[float], n_iter: int, seed: int = 0) -> tuple[float, float]:
    """Resample negatives with replacement to original size; recompute AUC.

    The negatives are all records whose score is not one of the positive
    scores' candidate ids. We approximate by removing positive scores
    from the score list (drop n_pos entries equal to each positive score
    in the sorted array), then bootstrap from that pool. With unique
    positive scores this is exact; with positive-tied-with-negatives it
    over-removes by at most n_pos entries (negligible at N≈3M).
    """
    rng = random.Random(seed)
    # Build negative pool: copy and remove n_pos positive scores.
    neg_pool = list(asc_scores)
    for ps in pos_scores:
        # Remove the leftmost equal entry (mid-rank semantics survive at scale).
        i = bisect.bisect_left(neg_pool, ps)
        if i < len(neg_pool) and neg_pool[i] == ps:
            neg_pool.pop(i)
    n_neg = len(neg_pool)
    n_pos = len(pos_scores)
    aucs: list[float] = []
    for _ in range(n_iter):
        # Resample n_neg negatives with replacement. Sorting after each
        # sample is O(N log N) — too slow at 3M. Instead, count how many
        # sampled negatives are < / = / > each positive via bisect on
        # the *original* sorted negative pool plus the binomial-equivalent
        # of multinomial sampling.
        # Trick: under sampling with replacement, the empirical
        # distribution of resampled negatives is multinomial over the
        # original neg_pool. So for each positive p, the expected
        # resampled count below p is (i_p / n_neg) * n_neg = i_p exactly
        # in the limit; the actual count is Binomial(n_neg, i_p / n_neg).
        # We sample those binomials directly — O(n_pos) per iter.
        sum_ranks = 0.0
        for ps in pos_scores:
            i_lt = bisect.bisect_left(neg_pool, ps)
            i_le = bisect.bisect_right(neg_pool, ps)
            p_lt = i_lt / n_neg
            p_eq = (i_le - i_lt) / n_neg
            below = _binomial(n_neg, p_lt, rng)
            ties = _binomial(n_neg - below, p_eq / (1 - p_lt) if p_lt < 1.0 else 0.0, rng)
            # Position of this positive in the resampled ascending list:
            #   `below` negatives are strictly less; `ties` negatives are
            #   tied (we'll mid-rank-credit half of them).
            mid = below + 1 + ties / 2.0
            sum_ranks += mid
        # AUC counts only the contribution against negatives (not other
        # positives, since positives are kept as-is). With n_pos = 3
        # the positives' relative order among themselves doesn't move
        # (their scores are fixed), so we use the formula:
        #   AUC = (sum_ranks_against_negatives_only) / (n_pos * n_neg)
        # where sum_ranks here is the count of ranks among negatives
        # (positives don't add to the U statistic). We add back the
        # within-positive ordering by using mid-ranks of positives among
        # themselves, but with n_pos=3 and sorted positive scores this
        # is just an additive constant.
        # For a clean comparison, drop the within-positive ordering term:
        u = sum_ranks - n_pos * (n_pos + 1) / 2.0  # standard correction
        aucs.append(u / (n_pos * n_neg))
    aucs.sort()
    p_lo = aucs[int(0.025 * n_iter)]
    p_hi = aucs[int(0.975 * n_iter)]
    return p_lo, p_hi


def _binomial(n: int, p: float, rng: random.Random) -> int:
    """Approximate Binomial(n, p) sample.

    For large n use the normal approximation; for small n iterate. Cheap
    and sufficient for descriptive bootstrap intervals.
    """
    if p <= 0.0:
        return 0
    if p >= 1.0:
        return n
    if n < 30:
        return sum(1 for _ in range(n) if rng.random() < p)
    mu = n * p
    sigma = math.sqrt(n * p * (1 - p))
    s = int(round(rng.gauss(mu, sigma)))
    return max(0, min(n, s))


def main() -> None:
    positives = set(line.strip() for line in POSITIVES.read_text().splitlines() if line.strip())
    print(f"positives: {sorted(positives)}")

    out_rows = []
    for cohort, axis, fname, field in RUNS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP missing: {path}")
            continue
        print(f"loading {path.name} on {field} ...")
        asc_scores, pos_scores = load_scores(path, field, positives)
        if not pos_scores:
            print(f"  no positives found in {fname}")
            continue
        observed = auc_from_scores(asc_scores, pos_scores)
        print(f"  observed AUC = {observed:.4f}")
        print(f"  permutation null (n={N_PERM:,}) ...")
        p_lo, p_hi, pval = permutation_null_auc(asc_scores, len(pos_scores), N_PERM, observed, seed=0)
        print(f"    null 2.5–97.5%: [{p_lo:.4f}, {p_hi:.4f}]  p_one_sided ≈ {pval:.4g}")
        print(f"  negative bootstrap (n={N_BOOT:,}) ...")
        b_lo, b_hi = negative_bootstrap_auc(asc_scores, pos_scores, N_BOOT, seed=0)
        print(f"    bootstrap 2.5–97.5%: [{b_lo:.4f}, {b_hi:.4f}]")
        out_rows.append((cohort, axis, observed, p_lo, p_hi, pval, b_lo, b_hi))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\taxis\tobserved_auc\tperm_null_p2.5\tperm_null_p97.5\tperm_p_value\tboot_neg_p2.5\tboot_neg_p97.5\n")
        for r in out_rows:
            fh.write("\t".join(str(x) for x in r) + "\n")
    print(f"wrote {OUT_TSV} ({len(out_rows)} rows)")

    md = ["# AUC uncertainty — n_pos = 3 primary endpoint",
          "",
          "Generated by `scripts/auc_uncertainty.py`.",
          f"Permutation null: {N_PERM:,} draws of n_pos = 3 random",
          "candidates (uniform without replacement) over the full",
          "candidate set; AUC computed under benchmark mid-rank tie",
          "handling. Reported as 2.5/97.5-percentile envelope and",
          "one-sided empirical *p*-value `Pr(AUC_perm ≥ AUC_observed)`.",
          f"Negative bootstrap: {N_BOOT:,} resamples of the negative",
          "set with replacement (binomial-multinomial closed form;",
          "positive scores held fixed). Both intervals are descriptive,",
          "not inferential CIs.",
          "",
          "| cohort | axis | observed AUC | null 2.5–97.5% | one-sided p | bootstrap-neg 2.5–97.5% |",
          "|---|---|---:|---|---:|---|"]
    for cohort, axis, observed, p_lo, p_hi, pval, b_lo, b_hi in out_rows:
        md.append(f"| {cohort} | {axis} | {observed:.3f} | [{p_lo:.3f}, {p_hi:.3f}] | {pval:.4g} | [{b_lo:.3f}, {b_hi:.3f}] |")
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
