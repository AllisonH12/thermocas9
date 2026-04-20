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
from thermocas.pan_cancer import aggregate, top_exclusive, top_recurrent


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


def test_top_exclusive_picks_cohort_specific_winners():
    cohorts = {
        "BRCA": [_scored("c_100", "BRCA", 1.5), _scored("c_200", "BRCA", 1.0)],
        "LUAD": [_scored("c_100", "LUAD", 0.1), _scored("c_200", "LUAD", 0.95)],
    }
    aggs = list(aggregate(cohorts, high_score_threshold=0.5))
    brca_only = top_exclusive(aggs, "BRCA", high_score_threshold=0.5)
    assert brca_only[0].candidate_id == "c_100"  # margin 1.5 - 0.1 = 1.4 > 1.0 - 0.95 = 0.05
