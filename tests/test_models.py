"""Validation tests for the Pydantic data model."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from thermocas.models import (
    CandidateSite,
    EvidenceClass,
    EvidenceThresholds,
    MethylationObservation,
    Penalties,
    ProbabilisticScore,
    ScoreComponents,
    ScoredCandidate,
    Strand,
)


def _make_candidate(cid: str = "chr1:100+:NNNNCGA") -> CandidateSite:
    return CandidateSite(
        candidate_id=cid,
        chrom="chr1",
        critical_c_pos=104,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def _make_observation(cid: str = "chr1:100+:NNNNCGA") -> MethylationObservation:
    return MethylationObservation(
        candidate_id=cid,
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


def test_candidate_site_uppercases_local_seq():
    site = CandidateSite(
        candidate_id="x",
        chrom="chr1",
        critical_c_pos=10,
        strand=Strand.PLUS,
        pam="NNNNCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
        local_seq_100bp="acgtnACGTN",
    )
    assert site.local_seq_100bp == "ACGTNACGTN"


def test_candidate_site_rejects_non_dna():
    with pytest.raises(ValidationError):
        CandidateSite(
            candidate_id="x",
            chrom="chr1",
            critical_c_pos=10,
            strand=Strand.PLUS,
            pam="NNNNCGA",
            pam_family="NNNNCGA",
            is_cpg_pam=True,
            local_seq_100bp="ACXGT",
        )


def test_observation_quantile_ordering():
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.EXACT,
            beta_tumor_mean=0.5,
            beta_tumor_q25=0.7,
            beta_tumor_q75=0.6,  # invalid: q25 > q75
        )


def test_observation_unobserved_must_have_no_beta():
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.UNOBSERVED,
            beta_tumor_mean=0.1,
        )


def test_observation_mean_must_lie_within_quantiles():
    """Regression: the model previously accepted a mean outside [q25, q75]."""

    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.EXACT,
            beta_tumor_mean=0.50,    # outside [0.10, 0.20]
            beta_tumor_q25=0.10,
            beta_tumor_q75=0.20,
        )
    # normal-side too
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.EXACT,
            beta_normal_mean=0.05,   # outside [0.20, 0.40]
            beta_normal_q25=0.20,
            beta_normal_q75=0.40,
        )


def test_observation_unobserved_rejects_normal_side_betas():
    """Regression: UNOBSERVED previously only blocked beta_tumor_mean; it now
    rejects any non-default evidence/sample/beta field."""

    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.UNOBSERVED,
            beta_normal_mean=0.5,   # was previously accepted
        )


def test_observation_unobserved_rejects_probe_metadata():
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.UNOBSERVED,
            probe_id="cg000001",
        )
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.UNOBSERVED,
            evidence_distance_bp=42,
        )


def test_observation_unobserved_rejects_nonzero_sample_counts():
    with pytest.raises(ValidationError):
        MethylationObservation(
            candidate_id="x",
            cohort_name="c",
            evidence_class=EvidenceClass.UNOBSERVED,
            n_samples_tumor=10,
        )


def test_pam_family_rejects_offset_outside_motif():
    """Regression: PamFamily must catch malformed configs like
    `regex='A', critical_c_offset=4` at load time."""

    from thermocas.models import PamFamily

    with pytest.raises(ValidationError, match="outside motif width"):
        PamFamily(
            name="bad",
            regex="A",
            critical_c_offset=4,
            is_cpg=True,
            weight=1.0,
        )


def test_pam_family_rejects_offset_pointing_at_non_C():
    """The position pointed to by critical_c_offset must require a C."""

    from thermocas.models import PamFamily

    with pytest.raises(ValidationError, match="not just C"):
        PamFamily(
            name="bad",
            regex="[ACGT][ACGT][ACGT][ACGT]CG[AG]",
            critical_c_offset=0,  # points at N, not C
            is_cpg=True,
            weight=1.0,
        )


def test_pam_family_rejects_invalid_regex():
    from thermocas.models import PamFamily

    with pytest.raises(ValidationError, match="does not compile"):
        PamFamily(
            name="bad",
            regex="[unterminated",
            critical_c_offset=0,
            is_cpg=True,
            weight=1.0,
        )


def test_pam_family_rejects_mixed_width_regex():
    """Regression: a regex with alternatives of different widths (e.g. `A|...CGA`)
    would let `find_pam_matches` produce out-of-bounds critical_c_pos values for
    the shorter matches. Reject up front."""

    from thermocas.models import PamFamily

    with pytest.raises(ValidationError, match="multiple widths"):
        PamFamily(
            name="bad",
            regex="A|[ACGT][ACGT][ACGT][ACGT]CG[AG]",
            critical_c_offset=4,
            is_cpg=True,
            weight=1.0,
        )


def test_pam_family_accepts_concrete_custom_pams():
    """Regression: validation must use regex semantics, not a hard-coded probe
    set. A concrete PAM like `CCCCCGA` (single matching string) is a legitimate
    member of the NNNNCGA family and must validate."""

    from thermocas.models import PamFamily

    PamFamily(
        name="custom-cga",
        regex="CCCCCGA",
        critical_c_offset=4,
        is_cpg=True,
        weight=1.0,
    )
    PamFamily(
        name="custom-cca",
        regex="CCCCCCA",
        critical_c_offset=4,
        is_cpg=False,
        weight=0.9,
    )


def test_pam_family_accepts_canonical_thermocas9_pams():
    """Sanity: the production PAM definitions must validate."""

    from thermocas.models import PamFamily

    PamFamily(
        name="NNNNCGA",
        regex="[ACGT][ACGT][ACGT][ACGT]CG[AG]",
        critical_c_offset=4,
        is_cpg=True,
        weight=1.0,
    )
    PamFamily(
        name="NNNNCCA",
        regex="[ACGT][ACGT][ACGT][ACGT]CC[AG]",
        critical_c_offset=4,
        is_cpg=False,
        weight=0.9,
    )


def test_scored_candidate_rejects_cross_cohort_probabilistic():
    """Regression: probabilistic.cohort_name must match observation.cohort_name.
    Without this, JSONL records can carry inconsistent cohort attribution."""

    from thermocas.models import (
        CandidateSite,
        CohortConfig,
        EvidenceClass,
        EvidenceThresholds,
        MethylationObservation,
        Penalties,
        ProbabilisticScore,
        ScoreComponents,
        ScoredCandidate,
    )

    cand = CandidateSite(
        candidate_id="x",
        chrom="chr1",
        critical_c_pos=10,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )
    obs = MethylationObservation(
        candidate_id="x",
        cohort_name="BRCA",
        evidence_class=EvidenceClass.EXACT,
        evidence_distance_bp=0,
        probe_id="cg001",
        beta_tumor_mean=0.05, beta_tumor_q25=0.02, beta_tumor_q75=0.10,
        beta_normal_mean=0.85, beta_normal_q25=0.78, beta_normal_q75=0.92,
        n_samples_tumor=400, n_samples_normal=80,
    )
    bad_prob = ProbabilisticScore(
        candidate_id="x",
        cohort_name="LUAD",   # mismatched cohort
        mode="tumor_only",
        p_targetable_tumor=0.9,
        p_protected_normal=0.9,
        p_observation_trustworthy=0.9,
        p_therapeutic_selectivity=0.9 * 0.9,
    )
    components = ScoreComponents(sequence_score=1.0, selectivity_score=0.5, confidence_score=1.0)
    with pytest.raises(ValidationError, match="cohort_name must match"):
        ScoredCandidate(
            candidate=cand,
            observation=obs,
            components=components,
            final_score=0.5,
            probabilistic=bad_prob,
        )


def test_observation_unobserved_minimal_record_is_allowed():
    """The empty UNOBSERVED record is still valid — that's how the synthetic
    pipeline and the evidence classifier signal 'no data here'."""

    obs = MethylationObservation(
        candidate_id="x",
        cohort_name="c",
        evidence_class=EvidenceClass.UNOBSERVED,
    )
    assert obs.evidence_class == EvidenceClass.UNOBSERVED
    assert obs.beta_tumor_mean is None
    assert obs.beta_normal_mean is None
    assert obs.n_samples_tumor == 0
    assert obs.n_samples_normal == 0


def test_scored_candidate_id_consistency():
    c = _make_candidate("a")
    o = _make_observation("b")  # mismatched id
    components = ScoreComponents(
        sequence_score=1.0, selectivity_score=0.5, confidence_score=1.0
    )
    with pytest.raises(ValidationError):
        ScoredCandidate(candidate=c, observation=o, components=components, final_score=0.0)


def test_evidence_thresholds_must_be_monotonic():
    with pytest.raises(ValidationError):
        EvidenceThresholds(exact_bp=10, proximal_close_bp=5, proximal_bp=50, regional_bp=500)


def test_evidence_thresholds_default_exact_is_zero():
    """Regression: V3.1 default was exact_bp=1, which over-classified
    1-bp-away probes as EXACT (1.0 confidence) when they should have been
    PROXIMAL_CLOSE (0.7). EXACT means the probe assays the same CpG."""
    thresholds = EvidenceThresholds()
    assert thresholds.exact_bp == 0


def test_evidence_thresholds_rejects_unknown_keys():
    """V3.1: extra='forbid' so YAML typos error rather than silently default."""
    with pytest.raises(ValidationError, match="extra"):
        EvidenceThresholds(exact_bp=0, exatc_bp=5)  # typo: exatc


def test_cohort_config_rejects_unknown_keys():
    """Regression: `min_samples_tumour` (UK spelling) used to silently fall
    back to the default min_samples_tumor=30 — now errors at load."""
    from thermocas.models import CohortConfig

    with pytest.raises(ValidationError, match="extra"):
        CohortConfig(
            name="T", tumor_dataset="a", normal_dataset="b", platform="HM450",
            min_samples_tumour=5,  # typo: tumour
        )


def test_penalties_rejects_unknown_keys():
    from thermocas.models import Penalties

    with pytest.raises(ValidationError, match="extra"):
        Penalties(heterogeneity_iqr_threshhold=0.5)  # typo: threshhold


def test_penalties_have_sensible_bounds():
    p = Penalties()
    assert 0.0 < p.heterogeneity_iqr_threshold <= 1.0
    assert p.low_coverage_n_threshold > 0


def test_probabilistic_score_stores_composite_and_mode():
    """V2.4 — p_therapeutic_selectivity is now a stored field whose value
    depends on `mode`. Callers MUST pass the precomputed composite; it's no
    longer a property derived at read time."""

    ps = ProbabilisticScore(
        candidate_id="x",
        cohort_name="c",
        mode="tumor_plus_normal_protection",
        p_targetable_tumor=0.9,
        p_protected_normal=0.8,
        p_observation_trustworthy=0.5,
        p_therapeutic_selectivity=0.9 * 0.8 * 0.5,
    )
    assert ps.mode == "tumor_plus_normal_protection"
    assert ps.p_therapeutic_selectivity == pytest.approx(0.9 * 0.8 * 0.5)

    # tumor_only variant
    ps2 = ProbabilisticScore(
        candidate_id="x", cohort_name="c",
        mode="tumor_only",
        p_targetable_tumor=0.9,
        p_protected_normal=0.8,
        p_observation_trustworthy=0.5,
        p_therapeutic_selectivity=0.9 * 0.5,
    )
    assert ps2.mode == "tumor_only"
    assert ps2.p_therapeutic_selectivity == pytest.approx(0.9 * 0.5)
