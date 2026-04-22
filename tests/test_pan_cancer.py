"""Pan-cancer aggregator tests."""

from __future__ import annotations

import pytest

from thermocas.models import (
    CandidateSite,
    EvidenceClass,
    MethylationObservation,
    ScoreComponents,
    ScoredCandidate,
    Strand,
)
from thermocas.pan_cancer import (
    aggregate,
    aggregate_streaming,
    top_exclusive,
    top_recurrent,
)


def _candidate(cid: str, pos: int) -> CandidateSite:
    return CandidateSite(
        candidate_id=cid,
        chrom="chr1",
        critical_c_pos=pos,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def _scored(cid: str, cohort: str, score: float, *, observed: bool = True) -> ScoredCandidate:
    cand = _candidate(cid, int(cid.split("_")[-1]) if "_" in cid else 100)
    if observed:
        obs = MethylationObservation(
            candidate_id=cid,
            cohort_name=cohort,
            evidence_class=EvidenceClass.EXACT,
            evidence_distance_bp=0,
            probe_id=f"cg-{cid}",
            beta_tumor_mean=0.10, beta_tumor_q25=0.05, beta_tumor_q75=0.15,
            n_samples_tumor=400,
            beta_normal_mean=0.80, beta_normal_q25=0.75, beta_normal_q75=0.85,
            n_samples_normal=80,
        )
    else:
        obs = MethylationObservation(
            candidate_id=cid,
            cohort_name=cohort,
            evidence_class=EvidenceClass.UNOBSERVED,
        )
    return ScoredCandidate(
        candidate=cand,
        observation=obs,
        components=ScoreComponents(
            sequence_score=1.0,
            selectivity_score=score if observed else 0.0,
            confidence_score=1.0 if observed else 0.0,
        ),
        final_score=score if observed else 0.0,
    )


def test_aggregate_basic_metrics():
    cohorts = {
        "BRCA": [_scored("c_100", "BRCA", 1.5)],
        "LUAD": [_scored("c_100", "LUAD", 1.0)],
        "COAD": [_scored("c_100", "COAD", 0.0, observed=False)],
    }
    out = list(aggregate(cohorts, high_score_threshold=0.5))
    assert len(out) == 1
    a = out[0]
    assert a.candidate_id == "c_100"
    assert a.n_cohorts_observed == 2
    assert a.n_cohorts_high_score == 2
    assert a.recurrence == pytest.approx(1.0)
    assert a.pan_cancer_score == pytest.approx((1.5 + 1.0) / 2)
    assert a.exclusivity > 0  # two cohorts → nonzero stdev
    # Best-case + worst-case protection across cohorts. Both BRCA and LUAD here
    # carry β_normal_mean=0.80, so max == min == 0.80.
    assert a.normal_protection_max is not None and a.normal_protection_max == pytest.approx(0.80)
    assert a.normal_protection_min is not None and a.normal_protection_min == pytest.approx(0.80)


def _scored_with_normal(cid: str, cohort: str, score: float, beta_normal: float) -> ScoredCandidate:
    """Variant that lets the test set per-cohort β_normal explicitly."""
    sc = _scored(cid, cohort, score)
    obs = sc.observation.model_copy(update={
        "beta_normal_mean": beta_normal,
        "beta_normal_q25": max(0.0, beta_normal - 0.05),
        "beta_normal_q75": min(1.0, beta_normal + 0.05),
    })
    return sc.model_copy(update={"observation": obs})


def test_aggregate_reports_both_best_and_worst_normal_protection():
    """Regression: the V3.1 rename adds normal_protection_min so users get the
    worst-case (= pan-cancer risk) view, not just the best-case."""

    cohorts = {
        "BRCA": [_scored_with_normal("c_500", "BRCA", 1.0, 0.90)],
        "LUAD": [_scored_with_normal("c_500", "LUAD", 1.0, 0.10)],
    }
    a = next(iter(aggregate(cohorts, high_score_threshold=0.5)))
    assert a.normal_protection_max == pytest.approx(0.90)  # best-case
    assert a.normal_protection_min == pytest.approx(0.10)  # worst-case (the risk view)


def test_aggregate_rejects_metadata_mismatch():
    """Regression: same candidate_id with different (chrom, pos, family) used
    to silently merge with first-wins metadata, mislabeling the output."""

    sc1 = _scored("same", "BRCA", 1.0)
    sc2 = _scored("same", "LUAD", 1.0)
    sc2 = sc2.model_copy(update={
        "candidate": sc2.candidate.model_copy(update={"chrom": "chr2", "critical_c_pos": 999}),
    })
    cohorts = {"BRCA": [sc1], "LUAD": [sc2]}
    with pytest.raises(ValueError, match="conflicting metadata"):
        list(aggregate(cohorts))


def test_aggregate_rejects_intra_cohort_duplicate_candidate():
    """Regression: a malformed cohort JSONL emitting the same candidate_id
    twice for the same cohort used to silently overwrite with the later
    record via `by_candidate.setdefault(cid, {})[name] = sc`. The streaming
    aggregator rejects this; the in-memory path must reject it too so
    both paths have the same contract under malformed input.
    """

    dup = [
        _scored("c_100", "BRCA", 0.1),
        _scored("c_100", "BRCA", 0.9),  # duplicate cid, same cohort
    ]
    with pytest.raises(ValueError, match="duplicate candidate_id"):
        list(aggregate({"BRCA": dup}))


def test_aggregate_tie_breaks_by_candidate_id_ascending():
    """Regression: in-memory aggregate emission order now breaks
    (chrom, pos, family) ties by candidate_id ascending. Previously the
    order was dict-insertion-order (cohort-iteration-order dependent).
    Construct two distinct candidate_ids at identical (chrom, pos, family)
    and feed them in reverse-cid order — output must still be cid-ascending.
    """

    sc_z = _scored_at("z_cand", "chr1", 100, "BRCA", 0.6)
    sc_a = _scored_at("a_cand", "chr1", 100, "BRCA", 0.6)
    # Insertion order: z first, then a. If the old dict-insertion tie-break
    # were still in effect, output would be ["z_cand", "a_cand"].
    out = list(aggregate({"BRCA": [sc_z, sc_a]}))
    cids = [a.candidate_id for a in out]
    assert cids == ["a_cand", "z_cand"], (
        f"expected cid-ascending tie-break, got {cids}"
    )


def test_aggregate_handles_only_unobserved():
    cohorts = {
        "BRCA": [_scored("c_200", "BRCA", 0.0, observed=False)],
        "LUAD": [_scored("c_200", "LUAD", 0.0, observed=False)],
    }
    out = list(aggregate(cohorts))
    assert len(out) == 1
    a = out[0]
    assert a.n_cohorts_observed == 0
    assert a.pan_cancer_score == 0.0
    assert a.recurrence == 0.0
    assert a.normal_protection_max is None
    assert a.normal_protection_min is None


def test_aggregate_emits_in_deterministic_order():
    cohorts = {
        "BRCA": [
            _scored("c_300", "BRCA", 0.5),
            _scored("c_100", "BRCA", 0.5),
            _scored("c_200", "BRCA", 0.5),
        ],
    }
    out = list(aggregate(cohorts))
    positions = [a.critical_c_pos for a in out]
    assert positions == sorted(positions)


def test_top_recurrent_filters_by_min_cohorts():
    cohorts = {
        "BRCA": [_scored("c_100", "BRCA", 1.0), _scored("c_200", "BRCA", 1.0)],
        "LUAD": [_scored("c_100", "LUAD", 1.0), _scored("c_200", "LUAD", 0.0, observed=False)],
    }
    aggs = list(aggregate(cohorts, high_score_threshold=0.5))
    top = top_recurrent(aggs, min_cohorts=2)
    assert {a.candidate_id for a in top} == {"c_100"}


# ---------- streaming aggregator ----------


def _candidate_at(cid: str, chrom: str, pos: int) -> CandidateSite:
    return CandidateSite(
        candidate_id=cid,
        chrom=chrom,
        critical_c_pos=pos,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def _scored_at(cid: str, chrom: str, pos: int, cohort: str, score: float,
               *, beta_normal: float = 0.80, observed: bool = True) -> ScoredCandidate:
    cand = _candidate_at(cid, chrom, pos)
    if observed:
        obs = MethylationObservation(
            candidate_id=cid, cohort_name=cohort,
            evidence_class=EvidenceClass.EXACT, evidence_distance_bp=0,
            probe_id=f"cg-{cid}",
            beta_tumor_mean=0.10, beta_tumor_q25=0.05, beta_tumor_q75=0.15,
            n_samples_tumor=400,
            beta_normal_mean=beta_normal,
            beta_normal_q25=max(0.0, beta_normal - 0.05),
            beta_normal_q75=min(1.0, beta_normal + 0.05),
            n_samples_normal=80,
        )
    else:
        obs = MethylationObservation(
            candidate_id=cid, cohort_name=cohort,
            evidence_class=EvidenceClass.UNOBSERVED,
        )
    return ScoredCandidate(
        candidate=cand, observation=obs,
        components=ScoreComponents(
            sequence_score=1.0,
            selectivity_score=score if observed else 0.0,
            confidence_score=1.0 if observed else 0.0,
        ),
        final_score=score if observed else 0.0,
    )


def _sort_for_streaming(records: list[ScoredCandidate]) -> list[ScoredCandidate]:
    return sorted(records, key=lambda sc: (
        sc.candidate.chrom, sc.candidate.critical_c_pos,
        sc.candidate.pam_family, sc.candidate.candidate_id,
    ))


def test_aggregate_streaming_matches_inmemory_for_random_input():
    """The streaming aggregator must produce identical PanCancerAggregate
    records (in identical order) as the in-memory aggregator when fed the
    same records in pre-sorted form."""

    brca = [
        _scored_at("c_a", "chr5", 200, "BRCA", 0.7),
        _scored_at("c_b", "chr5", 500, "BRCA", 0.9, beta_normal=0.95),
        _scored_at("c_c", "chr10", 100, "BRCA", 0.4),
        _scored_at("c_d", "chr6", 300, "BRCA", 0.0, observed=False),
    ]
    luad = [
        _scored_at("c_a", "chr5", 200, "LUAD", 0.6, beta_normal=0.30),
        _scored_at("c_c", "chr10", 100, "LUAD", 0.85),
        _scored_at("c_d", "chr6", 300, "LUAD", 0.5),
        _scored_at("c_e", "chr10", 50, "LUAD", 0.95),
    ]
    coad = [
        _scored_at("c_b", "chr5", 500, "COAD", 0.0, observed=False),
        _scored_at("c_e", "chr10", 50, "COAD", 0.6, beta_normal=0.10),
    ]

    inmem = list(aggregate(
        {"BRCA": list(brca), "LUAD": list(luad), "COAD": list(coad)},
        high_score_threshold=0.5,
    ))
    streamed = list(aggregate_streaming(
        {
            "BRCA": iter(_sort_for_streaming(brca)),
            "LUAD": iter(_sort_for_streaming(luad)),
            "COAD": iter(_sort_for_streaming(coad)),
        },
        high_score_threshold=0.5,
    ))
    assert [a.candidate_id for a in streamed] == [a.candidate_id for a in inmem]
    for s, m in zip(streamed, inmem, strict=True):
        assert s.model_dump() == m.model_dump(), (
            f"mismatch for {s.candidate_id}: streamed={s.model_dump()} "
            f"inmem={m.model_dump()}"
        )


def test_aggregate_streaming_handles_empty_cohort():
    """An empty cohort iterable should not break the merge or affect
    output for the other cohorts."""

    brca = [_scored_at("c_a", "chr5", 200, "BRCA", 0.7)]
    streamed = list(aggregate_streaming({
        "BRCA": iter(_sort_for_streaming(brca)),
        "EMPTY": iter([]),
    }))
    assert len(streamed) == 1
    a = streamed[0]
    assert a.candidate_id == "c_a"
    assert a.n_cohorts_observed == 1
    assert "BRCA" in a.cohort_scores and "EMPTY" not in a.cohort_scores


def test_aggregate_streaming_detects_unsorted_input():
    """Sort violation in an input stream must raise ValueError before any
    aggregate is emitted, so the user fixes the input rather than getting
    silently wrong groupings."""

    unsorted = [
        _scored_at("c_b", "chr5", 500, "BRCA", 0.9),  # comes first
        _scored_at("c_a", "chr5", 200, "BRCA", 0.7),  # but smaller pos
    ]
    with pytest.raises(ValueError, match="not sorted"):
        list(aggregate_streaming({"BRCA": iter(unsorted)}))


def test_aggregate_streaming_detects_intra_cohort_duplicate_candidate():
    """A cohort emitting the same candidate_id twice would silently overwrite
    in the in-memory path; in the streaming path we surface it as a hard
    error so corrupted JSONLs don't produce wrong aggregates."""

    dup = [
        _scored_at("c_a", "chr5", 200, "BRCA", 0.7),
        _scored_at("c_a", "chr5", 200, "BRCA", 0.9),  # duplicate
    ]
    with pytest.raises(ValueError, match="duplicate candidate_id"):
        list(aggregate_streaming({"BRCA": iter(dup)}))


def test_aggregate_streaming_aggregates_matched_cross_cohort_records():
    """Sanity: same candidate_id with matching (chrom, pos, family) across
    cohorts groups into a single PanCancerAggregate. This is the happy
    path — distinct from the metadata-mismatch case below."""

    brca = [_scored_at("same", "chr5", 200, "BRCA", 0.7)]
    luad = [_scored_at("same", "chr5", 200, "LUAD", 0.5, beta_normal=0.20)]
    out = list(aggregate_streaming({
        "BRCA": iter(brca), "LUAD": iter(luad),
    }))
    assert len(out) == 1
    assert out[0].candidate_id == "same"
    assert out[0].n_cohorts_observed == 2
    assert set(out[0].cohort_scores.keys()) == {"BRCA", "LUAD"}


def test_aggregate_streaming_rejects_cross_cohort_metadata_mismatch():
    """Regression: same candidate_id at conflicting (chrom, pos, family)
    across cohorts lands in different merge groups (the merge key includes
    metadata), so the per-group consistency check alone misses it. The
    `seen_cid → metadata` map catches it on the second appearance and
    raises — restoring parity with the in-memory aggregator's contract.
    """

    brca = [_scored_at("same", "chr5", 200, "BRCA", 0.7)]
    luad = [_scored_at("same", "chr6", 999, "LUAD", 0.7)]  # divergent meta

    # chr5/200 sorts before chr6/999 so BRCA's record is consumed first
    # and seeds seen_cid_meta. LUAD's record then trips the check.
    with pytest.raises(ValueError, match="conflicting metadata"):
        list(aggregate_streaming({
            "BRCA": iter(brca), "LUAD": iter(luad),
        }))


def test_top_exclusive_picks_cohort_specific_winners():
    cohorts = {
        "BRCA": [_scored("c_100", "BRCA", 1.5), _scored("c_200", "BRCA", 1.0)],
        "LUAD": [_scored("c_100", "LUAD", 0.1), _scored("c_200", "LUAD", 0.95)],
    }
    aggs = list(aggregate(cohorts, high_score_threshold=0.5))
    brca_only = top_exclusive(aggs, "BRCA", high_score_threshold=0.5)
    assert brca_only[0].candidate_id == "c_100"  # margin 1.5 - 0.1 = 1.4 > 1.0 - 0.95 = 0.05
