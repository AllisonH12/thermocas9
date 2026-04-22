"""Cross-cohort aggregator.

Consumes one stream of `ScoredCandidate`s per cohort and produces
`PanCancerAggregate` records: one per candidate, summarizing how it ranks
across cohorts.

Two paths share the same per-candidate emit contract:

  * `aggregate(...)`           — in-memory; loads every cohort into a
                                 candidate_id → cohort_name → ScoredCandidate
                                 map. Peak RAM is
                                 O(N_unique_candidate_ids × N_cohorts ×
                                 sizeof(ScoredCandidate)). Convenient for
                                 tests and small cohorts; do NOT use at
                                 genome scale.
  * `aggregate_streaming(...)` — k-way merge across cohort iterables that
                                 are pre-sorted by (chrom, pos, family,
                                 candidate_id). Candidate-side memory
                                 grows in N_unique_candidate_ids (a small
                                 `seen_cid → metadata` map used to enforce
                                 cross-cohort metadata parity; ~100 B per
                                 entry) rather than multiplying by
                                 N_cohorts × sizeof(ScoredCandidate). The
                                 only record-level state held is the head
                                 of each input stream plus the records for
                                 the current candidate group.

Metric definitions (see `models.PanCancerAggregate`):

  * pan_cancer_score      = mean(cohort_scores) over observed cohorts
  * recurrence            = n_cohorts_high_score / max(1, n_cohorts_observed)
  * exclusivity           = stdev(cohort_scores) over observed cohorts
                            (high → cohort-specific, low → pan-cancer recurrent)
  * normal_protection_max = max(β_normal_mean) across observed cohorts
                            (best-case protection — at least one cohort methylates)
  * normal_protection_min = min(β_normal_mean) across observed cohorts
                            (worst-case / pan-cancer risk view — least-protected cohort)

Both protection metrics are emitted because the single "max" view overstates
safety in cross-cohort summaries (a candidate with β_normal=0.9 in BRCA and
0.1 in LUAD shouldn't report only the 0.9).
"""

from __future__ import annotations

import heapq
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


def _candidate_sort_key(sc: ScoredCandidate) -> tuple[str, int, str, str]:
    c = sc.candidate
    return (c.chrom, c.critical_c_pos, c.pam_family, c.candidate_id)


def _build_aggregate(
    cid: str,
    chrom: str,
    pos: int,
    family: str,
    per_cohort: dict[str, ScoredCandidate],
    high_score_threshold: float,
) -> PanCancerAggregate:
    """Compute a single PanCancerAggregate from per-cohort observations."""

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

    return PanCancerAggregate(
        candidate_id=cid,
        chrom=chrom,
        critical_c_pos=pos,
        pam_family=family,
        n_cohorts_observed=n_observed,
        n_cohorts_high_score=n_high,
        cohort_scores=observed_scores,
        pan_cancer_score=pan_score,
        recurrence=recurrence,
        exclusivity=exclusivity,
        normal_protection_max=normal_protection_max,
        normal_protection_min=normal_protection_min,
    )


def _check_metadata_consistent(
    cid: str,
    existing: tuple[str, int, str],
    incoming: tuple[str, int, str],
) -> None:
    """Same candidate_id must always carry the same (chrom, pos, family)."""
    if existing != incoming:
        # Same candidate_id, different (chrom, pos, family) across input
        # JSONL files — the cohort runs were built from incompatible
        # catalogs; silently merging would emit a mislabeled pan-cancer
        # record. Fail fast.
        raise ValueError(
            f"candidate_id {cid!r} has conflicting metadata across cohorts: "
            f"{existing} vs {incoming} (cohorts likely scored against different catalogs)"
        )


def aggregate(
    cohorts: dict[str, Iterable[ScoredCandidate]],
    high_score_threshold: float = DEFAULT_HIGH_SCORE,
) -> Iterator[PanCancerAggregate]:
    """Aggregate per-candidate metrics across cohorts (in-memory).

    Args:
        cohorts: cohort_name → iterable of ScoredCandidate from that cohort.
                 Iterables are consumed exactly once.
        high_score_threshold: cohort-level final_score above which a candidate
                              counts toward `n_cohorts_high_score` and `recurrence`.

    Yields:
        One `PanCancerAggregate` per candidate that appears in at least one
        cohort, in ascending (chrom, critical_c_pos, pam_family, candidate_id)
        order.

    Memory characteristic:
        This function is **not single-pass streaming**. Despite returning an
        iterator, it first consumes every input record into a
        candidate_id → cohort_name → ScoredCandidate map, then sorts and
        emits. Peak memory is
        O(N_unique_candidate_ids × N_cohorts × sizeof(ScoredCandidate)).

        For genome-scale runs (tens of millions of candidates × tens of
        cohorts), use `aggregate_streaming(...)` instead — it requires each
        input iterable to be pre-sorted by the natural candidate order and
        its candidate-side memory grows in N_unique_candidate_ids only
        (the seen-cid → metadata parity map, ~100 B per entry), not
        multiplied by N_cohorts × sizeof(ScoredCandidate).
    """

    by_candidate: dict[str, dict[str, ScoredCandidate]] = {}
    candidate_meta: dict[str, tuple[str, int, str]] = {}  # cid → (chrom, pos, family)
    candidate_sort_keys: dict[str, tuple[str, int, str, str]] = {}

    for cohort_name, scored_iter in cohorts.items():
        for sc in scored_iter:
            cid = sc.candidate.candidate_id
            # Intra-cohort duplicate: a malformed cohort JSONL emitting the
            # same candidate_id twice for the same cohort used to silently
            # overwrite with the later record via `.setdefault(cid, {})[name] = sc`.
            # The streaming path rejects this; keep the two paths' malformed-
            # input contracts in sync.
            if cohort_name in by_candidate.get(cid, {}):
                raise ValueError(
                    f"cohort {cohort_name!r} contains duplicate candidate_id {cid!r}"
                )
            by_candidate.setdefault(cid, {})[cohort_name] = sc
            meta = (sc.candidate.chrom, sc.candidate.critical_c_pos, sc.candidate.pam_family)
            existing = candidate_meta.setdefault(cid, meta)
            _check_metadata_consistent(cid, existing, meta)
            candidate_sort_keys[cid] = _candidate_sort_key(sc)

    for cid in sorted(by_candidate, key=lambda c: candidate_sort_keys[c]):
        chrom, pos, fam = candidate_meta[cid]
        yield _build_aggregate(
            cid, chrom, pos, fam, by_candidate[cid], high_score_threshold,
        )


def aggregate_streaming(
    sorted_cohorts: dict[str, Iterable[ScoredCandidate]],
    high_score_threshold: float = DEFAULT_HIGH_SCORE,
    *,
    validate_sort: bool = True,
) -> Iterator[PanCancerAggregate]:
    """Streaming pan-cancer aggregator.

    Performs a k-way merge across cohort iterables, grouping consecutive
    records that share the same (chrom, pos, family, candidate_id) and
    yielding one PanCancerAggregate per group.

    PRECONDITION:
        Each input iterable MUST yield ScoredCandidate records in ascending
        (candidate.chrom, candidate.critical_c_pos, candidate.pam_family,
         candidate.candidate_id) order. The output of `score-cohort` is
        chrom-blocked but not lexicographically sorted by chrom (e.g.
        chr5 → chr6 → chr10), so callers typically need to presort each
        JSONL before invoking this function.

    Memory characteristic:
        Candidate-side memory grows in `N_unique_candidate_ids` rather
        than multiplying by `N_cohorts × sizeof(ScoredCandidate)`. The
        streaming path holds:
          - the head of each input stream (≤ N_cohorts records);
          - the records for the current candidate group (≤ N_cohorts);
          - a `seen_cid → (chrom, pos, family)` map used for cross-group
            metadata-mismatch detection (~100 B per unique cid).
        The in-memory `aggregate()` instead buffers every ScoredCandidate
        into a `cid → cohort_name → ScoredCandidate` dict, so its peak
        RAM scales as `N_unique_candidate_ids × N_cohorts × sizeof(
        ScoredCandidate)` — fine for small cohorts, impractical at
        genome scale.

    Args:
        sorted_cohorts: cohort_name → pre-sorted iterable of ScoredCandidate.
        high_score_threshold: see `aggregate(...)`.
        validate_sort: if True (default), monitor each input stream for sort
                       violations and raise ValueError on the first one. Set
                       False if you trust your inputs and want a marginal
                       speedup at large N.

    Yields:
        PanCancerAggregate records in (chrom, pos, family, candidate_id)
        order — identical to `aggregate(...)` output for the same input set.

    Raises:
        ValueError: cross-cohort `candidate_id` metadata mismatch (a cid
                    that appears with different (chrom, pos, family) in two
                    cohorts, possibly in different merge groups — caught
                    via the `seen_cid → metadata` map);
                    duplicate candidate_id within a single cohort;
                    or sort violation in an input stream (when
                    `validate_sort=True`).
    """

    def _keyed(name: str, it: Iterable[ScoredCandidate]):
        prev_key: tuple | None = None
        for sc in it:
            key = _candidate_sort_key(sc)
            if validate_sort and prev_key is not None and key < prev_key:
                raise ValueError(
                    f"cohort {name!r} input is not sorted: {prev_key} > {key}. "
                    "Presort the JSONL by (chrom, critical_c_pos, pam_family, "
                    "candidate_id) or call aggregate(...) instead."
                )
            prev_key = key
            yield (key, name, sc)

    merged = heapq.merge(
        *(_keyed(name, it) for name, it in sorted_cohorts.items()),
        key=lambda t: t[0],
    )

    # `current_key` is the 4-tuple (chrom, pos, family, candidate_id) that
    # defines the current candidate group. key[0:3] gives (chrom, pos,
    # family), so there's no need to track metadata in a separate variable.
    current_key: tuple[str, int, str, str] | None = None
    current_records: dict[str, ScoredCandidate] = {}

    # cid → first-seen (chrom, pos, family). Required to catch the case
    # where the same candidate_id appears in two cohorts with different
    # metadata: those records land in DIFFERENT merge groups (the merge
    # key includes metadata), so per-group grouping alone would miss it.
    # Without this map the streaming path would silently emit two
    # PanCancerAggregate records both labeled with the same cid — a
    # strictly weaker contract than the in-memory aggregator.
    seen_cid_meta: dict[str, tuple[str, int, str]] = {}

    def _check_seen(cid: str, meta: tuple[str, int, str]) -> None:
        existing = seen_cid_meta.setdefault(cid, meta)
        _check_metadata_consistent(cid, existing, meta)

    for key, name, sc in merged:
        if current_key is None:
            current_key = key
        if key != current_key:
            yield _build_aggregate(
                current_key[3], current_key[0], current_key[1], current_key[2],
                current_records, high_score_threshold,
            )
            current_key = key
            current_records = {}

        cid = sc.candidate.candidate_id
        meta = (sc.candidate.chrom, sc.candidate.critical_c_pos, sc.candidate.pam_family)
        # Cross-group consistency: if cid was seen earlier with different
        # metadata, this raises before we corrupt the aggregate.
        _check_seen(cid, meta)
        # Same cohort emitting the same candidate twice would be a duplicate
        # in the input stream — surface it explicitly rather than silently
        # overwriting.
        if name in current_records:
            raise ValueError(
                f"cohort {name!r} contains duplicate candidate_id {cid!r}"
            )
        current_records[name] = sc

    if current_records and current_key is not None:
        yield _build_aggregate(
            current_key[3], current_key[0], current_key[1], current_key[2],
            current_records, high_score_threshold,
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
