"""Candidate scoring schema.

The score formula is reviewer-auditable: every component is a separate function
returning a non-negative scalar, and `score_candidate` composes them according to

    final_score = sequence × selectivity × confidence
                − heterogeneity_penalty − low_coverage_penalty

Selectivity rewards class separation, not just shifts in mean (per Roth et al.
ranking guidance for methylation-array data).
"""

from __future__ import annotations

from thermocas.models import (
    CandidateSite,
    CohortConfig,
    EvidenceClass,
    MethylationObservation,
    PamFamily,
    Penalties,
    ScoreComponents,
    ScoredCandidate,
)

# Default confidence weights per evidence class. Override via cohort config later.
_DEFAULT_CONFIDENCE: dict[EvidenceClass, float] = {
    EvidenceClass.EXACT: 1.0,
    EvidenceClass.PROXIMAL_CLOSE: 0.7,
    EvidenceClass.PROXIMAL: 0.4,
    EvidenceClass.REGIONAL: 0.1,
    EvidenceClass.UNOBSERVED: 0.0,
}


def confidence_score_for_evidence(evidence: EvidenceClass) -> float:
    """Map an EvidenceClass to a confidence weight in [0, 1]."""

    return _DEFAULT_CONFIDENCE[evidence]


def selectivity_score(observation: MethylationObservation) -> float:
    """Reward tumor-low / normal-high methylation contrast.

        selectivity = max(0, normal_mean - tumor_mean)
                    + max(0, normal_q25  - tumor_q75)

    Returns 0.0 when either side lacks summary statistics — i.e. unobserved
    candidates do not get credit for "selectivity by absence."
    """

    if (
        observation.beta_tumor_mean is None
        or observation.beta_normal_mean is None
    ):
        return 0.0

    mean_term = max(0.0, observation.beta_normal_mean - observation.beta_tumor_mean)

    if (
        observation.beta_tumor_q75 is None
        or observation.beta_normal_q25 is None
    ):
        quantile_term = 0.0
    else:
        quantile_term = max(
            0.0, observation.beta_normal_q25 - observation.beta_tumor_q75
        )

    return mean_term + quantile_term


def heterogeneity_penalty(observation: MethylationObservation, penalties: Penalties) -> float:
    """Penalize candidates whose tumor methylation is highly variable.

    Estimated via tumor IQR (q75 - q25). When IQR exceeds the configured
    threshold the penalty is `(iqr - threshold) * weight`, otherwise zero.

    UNOBSERVED candidates carry no probe-mapped evidence; this penalty is for
    *observed* sites whose tumor distribution is wide. Don't conflate "no data"
    with "wide data."
    """

    if observation.evidence_class == EvidenceClass.UNOBSERVED:
        return 0.0
    if observation.beta_tumor_q25 is None or observation.beta_tumor_q75 is None:
        return 0.0
    iqr = observation.beta_tumor_q75 - observation.beta_tumor_q25
    excess = iqr - penalties.heterogeneity_iqr_threshold
    if excess <= 0.0:
        return 0.0
    return excess * penalties.heterogeneity_weight


def low_coverage_penalty(observation: MethylationObservation, penalties: Penalties) -> float:
    """Penalize candidates summarized from too few tumor samples.

    Only applies to *observed* candidates. UNOBSERVED candidates already have
    selectivity_score=0 and confidence_score=0, so their final_score is 0;
    applying a coverage penalty on top would push them below genuine zero-score
    observed sites and conflate "no probe coverage" with "low sample count."
    """

    if observation.evidence_class == EvidenceClass.UNOBSERVED:
        return 0.0
    if observation.n_samples_tumor >= penalties.low_coverage_n_threshold:
        return 0.0
    deficit = (penalties.low_coverage_n_threshold - observation.n_samples_tumor) / max(
        1, penalties.low_coverage_n_threshold
    )
    return deficit * penalties.low_coverage_weight


def sequence_score_for(family: PamFamily) -> float:
    """Sequence-compatibility score is the configured PAM family weight."""

    return family.weight


def score_candidate(
    candidate: CandidateSite,
    observation: MethylationObservation,
    family: PamFamily,
    cohort: CohortConfig,
    *,
    compute_probabilistic: bool = False,
    compute_spacer: bool = False,
) -> ScoredCandidate:
    """Compose a `ScoredCandidate` from a candidate, its evidence, and config.

    Pure function — no IO, fully deterministic given inputs. Suitable for direct
    unit testing and for vectorized application across genome-wide candidates.

    Args:
        compute_probabilistic: when True, also computes the V2 ProbabilisticScore
            (P(targetable_tumor) × P(protected_normal) × P(trustworthy)) and
            attaches it to `ScoredCandidate.probabilistic`. Default False keeps
            V1 behavior bit-for-bit identical.
        compute_spacer: when True, also computes the V3 SpacerScore (gRNA
            design-quality heuristics on the 20-nt protospacer) and attaches
            it to `ScoredCandidate.spacer`. Returns None on the spacer field
            for candidates whose `local_seq_100bp` doesn't contain a usable
            20-nt window upstream of the PAM (e.g. chromosome edge).
    """

    if candidate.candidate_id != observation.candidate_id:
        raise ValueError(
            f"candidate_id mismatch: {candidate.candidate_id} vs {observation.candidate_id}"
        )
    if candidate.pam_family != family.name:
        raise ValueError(
            f"PAM family mismatch: candidate is {candidate.pam_family}, model is {family.name}"
        )

    components = ScoreComponents(
        sequence_score=sequence_score_for(family),
        selectivity_score=selectivity_score(observation),
        confidence_score=confidence_score_for_evidence(observation.evidence_class),
        heterogeneity_penalty=heterogeneity_penalty(observation, cohort.penalties),
        low_coverage_penalty=low_coverage_penalty(observation, cohort.penalties),
    )

    final = (
        components.sequence_score
        * components.selectivity_score
        * components.confidence_score
        - components.heterogeneity_penalty
        - components.low_coverage_penalty
    )

    probabilistic = None
    if compute_probabilistic:
        # local import to avoid a top-level circular dependency
        from thermocas.probabilistic import (
            DEFAULT_GAP_SIGMOID_SIGMA_FIXED,
            probabilistic_score,
        )
        sigma_fixed = (
            cohort.sigma_fixed
            if cohort.sigma_fixed is not None
            else DEFAULT_GAP_SIGMOID_SIGMA_FIXED
        )
        probabilistic = probabilistic_score(
            observation,
            mode=cohort.probabilistic_mode,
            differential_delta=cohort.differential_delta,
            sigma_fixed=sigma_fixed,
        )

    spacer = None
    if compute_spacer:
        from thermocas.grna import score_spacer
        spacer = score_spacer(candidate, family)

    return ScoredCandidate(
        candidate=candidate,
        observation=observation,
        components=components,
        final_score=final,
        probabilistic=probabilistic,
        spacer=spacer,
    )
