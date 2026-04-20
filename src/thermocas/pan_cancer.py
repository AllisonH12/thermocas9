"""Cross-cohort aggregator.

Consumes one stream of `ScoredCandidate`s per cohort and produces
`PanCancerAggregate` records: one per candidate, summarizing how it ranks
across cohorts.

Metric definitions (see `models.PanCancerAggregate`):

  * pan_cancer_score = mean(cohort_scores) over observed cohorts
  * recurrence       = n_cohorts_high_score / max(1, n_cohorts_observed)
  * exclusivity      = stdev(cohort_scores) over observed cohorts
                       (high → cohort-specific, low → pan-cancer recurrent)
  * normal_risk_max  = max(beta_normal_mean) across cohorts
                       (higher → better protected in normal tissue)
"""

from __future__ import annotations

import statistics
from collections.abc import Iterable, Iterator

from thermocas.models import (
    EvidenceClass,
    PanCancerAggregate,
    ScoredCandidate,
)

#: Default threshold for "high score" — candidates above this are considered
#: addressable in a given cohort. Override per-call as needed.
DEFAULT_HIGH_SCORE = 0.30


def aggregate(
    cohorts: dict[str, Iterable[ScoredCandidate]],
    high_score_threshold: float = DEFAULT_HIGH_SCORE,
) -> Iterator[PanCancerAggregate]:
    """Aggregate per-candidate metrics across cohorts.

    Args:
        cohorts: cohort_name → iterable of ScoredCandidate from that cohort.
                 Iterables are consumed exactly once.
        high_score_threshold: cohort-level final_score above which a candidate
                              counts toward `n_cohorts_high_score` and `recurrence`.

    Yields:
        One `PanCancerAggregate` per candidate that appears in at least one cohort.
        Candidates are emitted in deterministic (chrom, critical_c_pos, family) order.
    """

    # collect: candidate_id → (CandidateSite-derived facts, dict[cohort → ScoredCandidate])
    by_candidate: dict[str, dict[str, ScoredCandidate]] = {}
    candidate_meta: dict[str, tuple[str, int, str]] = {}  # cid → (chrom, pos, family)

    for cohort_name, scored_iter in cohorts.items():
        for sc in scored_iter:
            cid = sc.candidate.candidate_id
            by_candidate.setdefault(cid, {})[cohort_name] = sc
            meta = (sc.candidate.chrom, sc.candidate.critical_c_pos, sc.candidate.pam_family)
            existing = candidate_meta.setdefault(cid, meta)
            if existing != meta:
                # Same candidate_id, different (chrom, pos, family) across input
                # JSONL files. This indicates the cohort runs were built from
                # incompatible catalogs; silently merging would emit a
                # mislabeled pan-cancer record. Fail fast.
                raise ValueError(
                    f"candidate_id {cid!r} has conflicting metadata across cohorts: "
                    f"{existing} vs {meta} (cohorts likely scored against different catalogs)"
                )

    # deterministic emission
    for cid in sorted(by_candidate, key=lambda c: candidate_meta[c]):
        per_cohort = by_candidate[cid]
        chrom, pos, fam = candidate_meta[cid]

        observed_scores: dict[str, float] = {}
        normal_means: list[float] = []
        for cohort_name, sc in per_cohort.items():
            obs = sc.observation
            if obs.evidence_class == EvidenceClass.UNOBSERVED:
                continue
            observed_scores[cohort_name] = sc.final_score
            if obs.beta_normal_mean is not None:
                normal_means.append(obs.beta_normal_mean)

        n_observed = len(observed_scores)
        n_high = sum(1 for s in observed_scores.values() if s >= high_score_threshold)
        pan_score = statistics.fmean(observed_scores.values()) if n_observed else 0.0
        recurrence = (n_high / n_observed) if n_observed else 0.0
        exclusivity = statistics.pstdev(observed_scores.values()) if n_observed >= 2 else 0.0
        normal_protection_max = max(normal_means) if normal_means else None
        normal_protection_min = min(normal_means) if normal_means else None

        yield PanCancerAggregate(
            candidate_id=cid,
            chrom=chrom,
            critical_c_pos=pos,
            pam_family=fam,
            n_cohorts_observed=n_observed,
            n_cohorts_high_score=n_high,
            cohort_scores=observed_scores,
            pan_cancer_score=pan_score,
            recurrence=recurrence,
            exclusivity=exclusivity,
            normal_protection_max=normal_protection_max,
            normal_protection_min=normal_protection_min,
        )


def top_recurrent(
    aggregates: Iterable[PanCancerAggregate],
    min_cohorts: int = 2,
    limit: int = 50,
) -> list[PanCancerAggregate]:
    """Return the top-`limit` candidates by `pan_cancer_score` that are
    high-confidence in at least `min_cohorts` cohorts.

    Useful for the headline 'pan-cancer addressable target' table in a paper.
    """

    qualifying = [a for a in aggregates if a.n_cohorts_high_score >= min_cohorts]
    qualifying.sort(key=lambda a: (a.pan_cancer_score, a.recurrence), reverse=True)
    return qualifying[:limit]


def top_exclusive(
    aggregates: Iterable[PanCancerAggregate],
    cohort_name: str,
    high_score_threshold: float = DEFAULT_HIGH_SCORE,
    limit: int = 50,
) -> list[PanCancerAggregate]:
    """Return the top-`limit` candidates that are high-scoring in `cohort_name`
    but low-scoring or unobserved elsewhere — the cohort-specific shortlist."""

    out: list[tuple[float, PanCancerAggregate]] = []
    for a in aggregates:
        s = a.cohort_scores.get(cohort_name)
        if s is None or s < high_score_threshold:
            continue
        # margin = cohort score minus the second-best cohort score (or 0 if alone)
        others = [v for k, v in a.cohort_scores.items() if k != cohort_name]
        runner_up = max(others) if others else 0.0
        margin = s - runner_up
        out.append((margin, a))
    out.sort(key=lambda t: t[0], reverse=True)
    return [a for _, a in out[:limit]]
