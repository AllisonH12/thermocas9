"""Scoring schema tests."""

from __future__ import annotations

import pytest

from thermocas.models import (
    CandidateSite,
    CohortConfig,
    EvidenceClass,
    EvidenceThresholds,
    MethylationObservation,
    PamFamily,
    Penalties,
    Strand,
)
from thermocas.scoring import (
    confidence_score_for_evidence,
    heterogeneity_penalty,
    low_coverage_penalty,
    score_candidate,
    selectivity_score,
)


def _candidate() -> CandidateSite:
    return CandidateSite(
        candidate_id="chr10:8045463+:NNNNCGA",
        chrom="chr10",
        critical_c_pos=8045463,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def _family() -> PamFamily:
    return PamFamily(
        name="NNNNCGA",
        regex="[ACGT][ACGT][ACGT][ACGT]CG[AG]",
        critical_c_offset=4,
        is_cpg=True,
        weight=1.0,
    )


def _cohort() -> CohortConfig:
    return CohortConfig(
        name="TCGA-BRCA",
        tumor_dataset="TCGA-BRCA",
        normal_dataset="TCGA-BRCA-normal",
        platform="HM450",
        evidence_thresholds=EvidenceThresholds(),
        penalties=Penalties(),
    )


def _high_selectivity_observation() -> MethylationObservation:
    return MethylationObservation(
        candidate_id="chr10:8045463+:NNNNCGA",
        cohort_name="TCGA-BRCA",
        evidence_class=EvidenceClass.EXACT,
        evidence_distance_bp=0,
        beta_tumor_mean=0.05,
        beta_tumor_q25=0.02,
        beta_tumor_q75=0.10,
        n_samples_tumor=400,
        beta_normal_mean=0.85,
        beta_normal_q25=0.78,
        beta_normal_q75=0.92,
        n_samples_normal=80,
    )


def test_confidence_weights_strict_ordering():
    """Confidence must strictly decrease with looser evidence."""

    weights = [
        confidence_score_for_evidence(ec)
        for ec in [
            EvidenceClass.EXACT,
            EvidenceClass.PROXIMAL_CLOSE,
            EvidenceClass.PROXIMAL,
            EvidenceClass.REGIONAL,
            EvidenceClass.UNOBSERVED,
        ]
    ]
    assert weights == sorted(weights, reverse=True)
    assert weights[-1] == 0.0


def test_selectivity_combines_mean_and_quantile_terms():
    obs = _high_selectivity_observation()
    s = selectivity_score(obs)
    expected = (0.85 - 0.05) + (0.78 - 0.10)
    assert s == pytest.approx(expected)


def test_selectivity_zero_when_tumor_more_methylated_than_normal():
    obs = _high_selectivity_observation().model_copy(update={
        "beta_tumor_mean": 0.95, "beta_tumor_q25": 0.92, "beta_tumor_q75": 0.98,
        "beta_normal_mean": 0.10, "beta_normal_q25": 0.05, "beta_normal_q75": 0.15,
    })
    assert selectivity_score(obs) == 0.0


def test_selectivity_zero_when_unobserved():
    obs = MethylationObservation(
        candidate_id="x",
        cohort_name="c",
        evidence_class=EvidenceClass.UNOBSERVED,
    )
    assert selectivity_score(obs) == 0.0


def test_heterogeneity_penalty_zero_under_threshold():
    cohort = _cohort()
    obs = _high_selectivity_observation()  # IQR 0.08, below 0.30
    assert heterogeneity_penalty(obs, cohort.penalties) == 0.0


def test_heterogeneity_penalty_proportional_to_excess():
    cohort = _cohort()
    obs = _high_selectivity_observation().model_copy(update={
        "beta_tumor_q25": 0.10, "beta_tumor_q75": 0.60,  # IQR 0.5, excess 0.2
    })
    expected = 0.2 * cohort.penalties.heterogeneity_weight
    assert heterogeneity_penalty(obs, cohort.penalties) == pytest.approx(expected)


def test_low_coverage_penalty_zero_at_threshold():
    cohort = _cohort()
    obs = _high_selectivity_observation()
    assert low_coverage_penalty(obs, cohort.penalties) == 0.0


def test_low_coverage_penalty_scales_with_deficit():
    cohort = _cohort()
    obs = _high_selectivity_observation().model_copy(update={"n_samples_tumor": 6})
    threshold = cohort.penalties.low_coverage_n_threshold
    deficit = (threshold - 6) / threshold
    expected = deficit * cohort.penalties.low_coverage_weight
    assert low_coverage_penalty(obs, cohort.penalties) == pytest.approx(expected)


def test_score_candidate_full_path():
    cand = _candidate()
    obs = _high_selectivity_observation()
    fam = _family()
    cohort = _cohort()
    scored = score_candidate(cand, obs, fam, cohort)
    sel = (0.85 - 0.05) + (0.78 - 0.10)
    expected = 1.0 * sel * 1.0 - 0.0 - 0.0
    assert scored.final_score == pytest.approx(expected)
    assert scored.components.selectivity_score == pytest.approx(sel)
    assert scored.components.confidence_score == 1.0


def test_score_candidate_rejects_id_mismatch():
    cand = _candidate()
    obs = _high_selectivity_observation().model_copy(update={"candidate_id": "different"})
    with pytest.raises(ValueError):
        score_candidate(cand, obs, _family(), _cohort())


def test_score_candidate_propagates_cohort_probabilistic_mode():
    """V2.4 — cohort.probabilistic_mode flows through to the emitted record."""
    cand = _candidate()
    obs  = _high_selectivity_observation()
    fam  = _family()

    cohort_to = _cohort().model_copy(update={"probabilistic_mode": "tumor_only"})
    sc_to = score_candidate(cand, obs, fam, cohort_to, compute_probabilistic=True)
    assert sc_to.probabilistic is not None
    assert sc_to.probabilistic.mode == "tumor_only"
    # tumor_only: p_sel = p_targ × p_trust
    prob = sc_to.probabilistic
    expected_to = prob.p_targetable_tumor * prob.p_observation_trustworthy
    assert prob.p_therapeutic_selectivity == pytest.approx(expected_to)

    cohort_tp = _cohort().model_copy(update={"probabilistic_mode": "tumor_plus_normal_protection"})
    sc_tp = score_candidate(cand, obs, fam, cohort_tp, compute_probabilistic=True)
    assert sc_tp.probabilistic is not None
    assert sc_tp.probabilistic.mode == "tumor_plus_normal_protection"
    prob_tp = sc_tp.probabilistic
    expected_tp = (
        prob_tp.p_targetable_tumor
        * prob_tp.p_protected_normal
        * prob_tp.p_observation_trustworthy
    )
    assert prob_tp.p_therapeutic_selectivity == pytest.approx(expected_tp)


def test_score_candidate_rejects_family_mismatch():
    cand = _candidate()
    obs = _high_selectivity_observation()
    fam = _family().model_copy(update={"name": "NNNNCCA"})
    with pytest.raises(ValueError):
        score_candidate(cand, obs, fam, _cohort())


def test_unobserved_observation_yields_zero_penalties():
    """Regression: an UNOBSERVED candidate must not be penalized for low coverage
    or heterogeneity. Default low_coverage_penalty would otherwise push final_score
    to -0.3, ranking unobserved sites *below* genuinely-zero observed sites."""

    cohort = _cohort()
    obs = MethylationObservation(
        candidate_id="x",
        cohort_name="c",
        evidence_class=EvidenceClass.UNOBSERVED,
    )
    assert heterogeneity_penalty(obs, cohort.penalties) == 0.0
    assert low_coverage_penalty(obs, cohort.penalties) == 0.0


def test_score_candidate_unobserved_final_score_is_zero():
    """End-to-end: UNOBSERVED scored candidate has final_score == 0, not -0.3."""

    cand = _candidate()
    obs = MethylationObservation(
        candidate_id=cand.candidate_id,
        cohort_name="TCGA-BRCA",
        evidence_class=EvidenceClass.UNOBSERVED,
    )
    scored = score_candidate(cand, obs, _family(), _cohort())
    assert scored.final_score == 0.0
    assert scored.components.confidence_score == 0.0
    assert scored.components.selectivity_score == 0.0
