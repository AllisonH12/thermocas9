"""Cohort adapter tests."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pytest

from thermocas.cohort import score_cohort
from thermocas.evidence import ProbeRecord
from thermocas.methylation_backend import BetaSummary, MethylationBackend
from thermocas.models import (
    CandidateSite,
    CohortConfig,
    EvidenceClass,
    EvidenceThresholds,
    Penalties,
    Strand,
)
from thermocas.pam_model import PamModel

REPO_ROOT = Path(__file__).resolve().parent.parent
PAM_YAML = REPO_ROOT / "config" / "pam_model.yaml"


class FakeBackend(MethylationBackend):
    """Hand-wired backend for tests — no IO."""

    def __init__(
        self,
        probes: list[ProbeRecord],
        tumor: dict[str, BetaSummary],
        normal: dict[str, BetaSummary],
    ) -> None:
        self._probes = probes
        self._tumor = tumor
        self._normal = normal

    def probes(self) -> Iterable[ProbeRecord]:
        return self._probes

    def tumor_summary(self, probe_id: str) -> BetaSummary | None:
        return self._tumor.get(probe_id)

    def normal_summary(self, probe_id: str) -> BetaSummary | None:
        return self._normal.get(probe_id)


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


def _cohort(min_t: int = 30, min_n: int = 10) -> CohortConfig:
    return CohortConfig(
        name="TEST",
        tumor_dataset="TEST-T",
        normal_dataset="TEST-N",
        platform="HM450",
        min_samples_tumor=min_t,
        min_samples_normal=min_n,
        evidence_thresholds=EvidenceThresholds(),
        penalties=Penalties(),
    )


@pytest.fixture
def pam_model() -> PamModel:
    return PamModel.from_yaml(PAM_YAML)


def test_score_cohort_with_exact_evidence(pam_model: PamModel):
    cands = [_candidate("c1", 100)]
    backend = FakeBackend(
        probes=[ProbeRecord(probe_id="cg001", chrom="chr1", pos=100)],
        tumor={"cg001": BetaSummary("cg001", 400, 0.05, 0.02, 0.10)},
        normal={"cg001": BetaSummary("cg001", 80, 0.85, 0.78, 0.92)},
    )
    scored = list(score_cohort(cands, backend, _cohort(), pam_model))
    assert len(scored) == 1
    s = scored[0]
    assert s.observation.evidence_class == EvidenceClass.EXACT
    assert s.observation.probe_id == "cg001"
    assert s.observation.beta_tumor_mean == pytest.approx(0.05)
    assert s.final_score > 0


def test_score_cohort_unobserved_when_no_probe_within_regional(pam_model: PamModel):
    cands = [_candidate("c1", 100)]
    backend = FakeBackend(
        probes=[ProbeRecord(probe_id="cg001", chrom="chr1", pos=100_000)],
        tumor={"cg001": BetaSummary("cg001", 400, 0.05, 0.02, 0.10)},
        normal={"cg001": BetaSummary("cg001", 80, 0.85, 0.78, 0.92)},
    )
    scored = list(score_cohort(cands, backend, _cohort(), pam_model))
    s = scored[0]
    assert s.observation.evidence_class == EvidenceClass.UNOBSERVED
    assert s.observation.beta_tumor_mean is None  # validator enforces this
    assert s.final_score == 0.0


def test_score_cohort_downgrades_when_below_min_samples_tumor(pam_model: PamModel):
    cands = [_candidate("c1", 100)]
    backend = FakeBackend(
        probes=[ProbeRecord(probe_id="cg001", chrom="chr1", pos=100)],
        tumor={"cg001": BetaSummary("cg001", 5, 0.05, 0.02, 0.10)},  # only 5 < 30
        normal={"cg001": BetaSummary("cg001", 80, 0.85, 0.78, 0.92)},
    )
    s = next(iter(score_cohort(cands, backend, _cohort(min_t=30), pam_model)))
    assert s.observation.evidence_class == EvidenceClass.UNOBSERVED


def test_score_cohort_downgrades_when_below_min_samples_normal(pam_model: PamModel):
    cands = [_candidate("c1", 100)]
    backend = FakeBackend(
        probes=[ProbeRecord(probe_id="cg001", chrom="chr1", pos=100)],
        tumor={"cg001": BetaSummary("cg001", 400, 0.05, 0.02, 0.10)},
        normal={"cg001": BetaSummary("cg001", 2, 0.85, 0.78, 0.92)},  # 2 < 10
    )
    s = next(iter(score_cohort(cands, backend, _cohort(min_n=10), pam_model)))
    assert s.observation.evidence_class == EvidenceClass.UNOBSERVED
