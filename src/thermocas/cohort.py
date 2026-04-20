"""Cohort adapter — wires a CandidateSite catalog and a MethylationBackend
into a stream of `ScoredCandidate`s for one cohort.

Responsibilities (per the four-layer design):

  * pull every probe from the backend → build an `EvidenceClassifier`
  * for each candidate site, find the nearest probe → assign EvidenceClass
  * for observed sites, pull tumor + normal `BetaSummary` and apply
    `min_samples_*` filters from the cohort config
  * build `MethylationObservation` (preserving the model's strict invariants)
  * score via `scoring.score_candidate`
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator

from thermocas.evidence import EvidenceClassifier
from thermocas.methylation_backend import BetaSummary, MethylationBackend
from thermocas.models import (
    CandidateSite,
    CohortConfig,
    EvidenceClass,
    MethylationObservation,
    ScoredCandidate,
)
from thermocas.pam_model import PamModel
from thermocas.scoring import score_candidate


def score_cohort(
    candidates: Iterable[CandidateSite],
    backend: MethylationBackend,
    cohort: CohortConfig,
    pam_model: PamModel,
    *,
    compute_probabilistic: bool = False,
    compute_spacer: bool = False,
) -> Iterator[ScoredCandidate]:
    """Score every candidate against one cohort, yielding `ScoredCandidate`.

    Streaming so genome-scale catalogs don't have to fit in memory.

    Args:
        compute_probabilistic: forwarded to `score_candidate`. When True, each
            yielded ScoredCandidate carries a ProbabilisticScore alongside the
            deterministic ScoreComponents.
        compute_spacer: forwarded to `score_candidate`. When True, each yielded
            ScoredCandidate carries a SpacerScore for the upstream 20-nt
            protospacer (gRNA design-quality heuristics).
    """

    classifier = EvidenceClassifier(backend.probes(), cohort.evidence_thresholds)

    for candidate in candidates:
        family = pam_model.get(candidate.pam_family)

        evidence_class, probe, distance = classifier.classify(
            candidate.chrom, candidate.critical_c_pos
        )

        observation = _make_observation(
            candidate, evidence_class, probe.probe_id if probe else None,
            distance, backend, cohort,
        )
        yield score_candidate(
            candidate, observation, family, cohort,
            compute_probabilistic=compute_probabilistic,
            compute_spacer=compute_spacer,
        )


def _make_observation(
    candidate: CandidateSite,
    evidence_class: EvidenceClass,
    probe_id: str | None,
    distance: int | None,
    backend: MethylationBackend,
    cohort: CohortConfig,
) -> MethylationObservation:
    """Build a validated `MethylationObservation`. Falls back to UNOBSERVED when
    the cohort fails its `min_samples_*` filter — a methylation summary based
    on too few samples is treated as no usable evidence."""

    if evidence_class == EvidenceClass.UNOBSERVED or probe_id is None:
        return MethylationObservation(
            candidate_id=candidate.candidate_id,
            cohort_name=cohort.name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )

    tumor = backend.tumor_summary(probe_id)
    normal = backend.normal_summary(probe_id)

    if not _has_min_samples(tumor, cohort.min_samples_tumor) or not _has_min_samples(
        normal, cohort.min_samples_normal
    ):
        # downgrade to UNOBSERVED so scoring doesn't reward thin evidence
        return MethylationObservation(
            candidate_id=candidate.candidate_id,
            cohort_name=cohort.name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )

    assert tumor is not None and normal is not None  # for type checkers

    return MethylationObservation(
        candidate_id=candidate.candidate_id,
        cohort_name=cohort.name,
        evidence_class=evidence_class,
        evidence_distance_bp=distance,
        probe_id=probe_id,
        beta_tumor_mean=tumor.mean,
        beta_tumor_q25=tumor.q25,
        beta_tumor_q75=tumor.q75,
        n_samples_tumor=tumor.n_samples,
        beta_normal_mean=normal.mean,
        beta_normal_q25=normal.q25,
        beta_normal_q75=normal.q75,
        n_samples_normal=normal.n_samples,
    )


def _has_min_samples(summary: BetaSummary | None, threshold: int) -> bool:
    return summary is not None and summary.n_samples >= threshold
