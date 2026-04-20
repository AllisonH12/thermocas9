"""V3 — cross-validation benchmark harness.

Turns "we found a few good sites" into "the prioritization method
generalizes". The reviewer-facing question this module answers is:

    Given a list of *known* positive sites (literature-validated targets,
    pilot-assay hits, or held-back sites), does our scoring put them at the
    top of the ranked list?

Key operations:
  * `split_by_chrom(candidates, holdout)` — partition into train / test by
    chromosome. The framework's ranking is supposed to generalize, so we
    measure it on chromosomes never seen during scoring.
  * `evaluate_ranking(scored, positives, top_k)` — top-K precision and recall
    plus ROC-AUC via the Mann-Whitney U statistic. No scipy / sklearn dep.
  * `BenchmarkResult` — Pydantic record suitable for JSONL.

The harness is *score-agnostic*: it consumes whatever `final_score` came out
of `score_cohort`, so callers can benchmark deterministic, probabilistic, or
spacer-aware scoring against the same positives list and compare.
"""

from __future__ import annotations

from collections.abc import Iterable

from thermocas.models import (
    BenchmarkResult,
    CandidateSite,
    ScoredCandidate,
)


# ---------- chromosome splits ----------


def split_by_chrom(
    candidates: Iterable[CandidateSite],
    holdout_chroms: set[str],
) -> tuple[list[CandidateSite], list[CandidateSite]]:
    """Partition candidates into (train, test) by chromosome membership.

    Test = candidates whose `chrom` is in `holdout_chroms`. Train = the rest.
    Ordering within each set is preserved.
    """

    train: list[CandidateSite] = []
    test: list[CandidateSite] = []
    for c in candidates:
        if c.chrom in holdout_chroms:
            test.append(c)
        else:
            train.append(c)
    return train, test


# ---------- ranking metrics ----------


def evaluate_ranking(
    scored: Iterable[ScoredCandidate],
    positives: set[str],
    *,
    cohort_name: str,
    top_k: int = 10,
    held_out_chromosomes: list[str] | None = None,
    enforce_holdout: bool = True,
    score_field: str = "final_score",
    missing_score_policy: str = "rank_last",
) -> BenchmarkResult:
    """Measure ranking quality of `scored` against a known-positives set.

    Args:
        scored: stream of ScoredCandidate.
        positives: set of `candidate_id`s considered ground-truth positives.
            Anything else in `scored` counts as a negative.
        cohort_name: passed through to the result for downstream JSONL labeling.
        top_k: K for precision@K and recall@K.
        held_out_chromosomes: list of chromosomes the scoring did NOT see during
            training. When `enforce_holdout=True` (default), only candidates on
            these chromosomes are evaluated — this is the cross-validation
            control the docstring advertises.
        enforce_holdout: if True (default), filter `scored` so only candidates
            whose `candidate.chrom` ∈ `held_out_chromosomes` are evaluated.
            Set False to evaluate the full scored set with the held-out list
            recorded only as metadata.
        score_field: which scalar to rank by. Default is the deterministic
            `final_score`. Pass `"p_therapeutic_selectivity"` to rank by the
            V2 probabilistic composite, or `"spacer_final_score"` for V3.
        missing_score_policy: how to handle candidates that lack the requested
            sub-score (e.g. ranking by `spacer_final_score` when `spacer=None`):
              * `"rank_last"` (default) — assign −∞ so the candidate counts
                toward `n_total` / `n_positives` / `n_negatives` but never
                ranks above any candidate that *does* have the sub-score
              * `"drop"` — silently exclude (the V3.0 behavior; preserves the
                old numerics, hides "I don't know" candidates)
              * `"error"` — raise on the first missing sub-score
    """

    if top_k < 1:
        raise ValueError("top_k must be >= 1")
    if missing_score_policy not in {"rank_last", "drop", "error"}:
        raise ValueError(
            f"missing_score_policy must be rank_last|drop|error, got {missing_score_policy!r}"
        )

    holdout_set: set[str] = set(held_out_chromosomes or [])
    if enforce_holdout and not holdout_set:
        # No holdout requested → trivially "all candidates are in the holdout".
        # Treating this as enforce-False keeps the contract intuitive.
        enforce_holdout = False

    pairs: list[tuple[str, float, bool]] = []
    missing_sentinel = float("-inf")
    for sc in scored:
        if enforce_holdout and sc.candidate.chrom not in holdout_set:
            continue
        score = _extract_score(sc, score_field)
        if score is None:
            if missing_score_policy == "drop":
                continue
            if missing_score_policy == "error":
                raise ValueError(
                    f"candidate {sc.candidate.candidate_id} has no {score_field!r}"
                )
            score = missing_sentinel  # rank_last
        pairs.append((sc.candidate.candidate_id, score, sc.candidate.candidate_id in positives))

    n_total = len(pairs)
    n_pos = sum(1 for _, _, is_pos in pairs if is_pos)
    n_neg = n_total - n_pos

    if n_total == 0 or n_pos == 0:
        return BenchmarkResult(
            cohort_name=cohort_name,
            n_total=n_total, n_positives=n_pos, n_negatives=n_neg,
            held_out_chromosomes=held_out_chromosomes or [],
            top_k=top_k,
            precision_at_k=None, recall_at_k=None, roc_auc=None,
        )

    pairs.sort(key=lambda t: t[1], reverse=True)
    k = min(top_k, n_total)
    top_k_positives = sum(1 for _, _, is_pos in pairs[:k] if is_pos)
    precision = top_k_positives / k
    recall = top_k_positives / n_pos
    auc = _roc_auc_mann_whitney(pairs)

    return BenchmarkResult(
        cohort_name=cohort_name,
        n_total=n_total, n_positives=n_pos, n_negatives=n_neg,
        held_out_chromosomes=held_out_chromosomes or [],
        top_k=top_k,
        precision_at_k=precision,
        recall_at_k=recall,
        roc_auc=auc,
    )


# ---------- internals ----------


def _extract_score(sc: ScoredCandidate, field: str) -> float | None:
    if field == "final_score":
        return sc.final_score
    if field == "p_therapeutic_selectivity":
        return sc.probabilistic.p_therapeutic_selectivity if sc.probabilistic else None
    if field == "spacer_final_score":
        return sc.spacer.final_score if sc.spacer else None
    if field == "naive_selectivity":
        # Baseline ablation: pure (β_normal − β_tumor) signal, with no quantile
        # term, no confidence weighting, no penalties. If this beats
        # `final_score` on a benchmark, the framework's machinery isn't adding
        # value — every additional layer should be auditable against this.
        obs = sc.observation
        if obs.beta_tumor_mean is None or obs.beta_normal_mean is None:
            return None
        return obs.beta_normal_mean - obs.beta_tumor_mean
    raise ValueError(f"unknown score_field: {field!r}")


def _roc_auc_mann_whitney(
    pairs: list[tuple[str, float, bool]],
) -> float | None:
    """Compute ROC-AUC as P(score(positive) > score(negative)) via Mann-Whitney U.

    Tied scores contribute 0.5 each (standard tie handling). Returns None when
    either positives or negatives are absent.
    """

    pos_scores = [s for _, s, p in pairs if p]
    neg_scores = [s for _, s, p in pairs if not p]
    n_pos = len(pos_scores)
    n_neg = len(neg_scores)
    if n_pos == 0 or n_neg == 0:
        return None

    # O(n*m) — fine for the sizes a benchmark sees (tens of thousands max).
    wins = 0.0
    for ps in pos_scores:
        for ns in neg_scores:
            if ps > ns:
                wins += 1.0
            elif ps == ns:
                wins += 0.5
    return wins / (n_pos * n_neg)


__all__ = [
    "evaluate_ranking",
    "split_by_chrom",
]
