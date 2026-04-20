"""Cohort-config YAML loader."""

from __future__ import annotations

from pathlib import Path

import yaml

from thermocas.models import CohortConfig, EvidenceThresholds, Penalties


def load_cohort_config(path: str | Path) -> CohortConfig:
    """Parse a cohort YAML into a validated `CohortConfig`."""

    with Path(path).open() as f:
        raw = yaml.safe_load(f)

    cohort_block = raw["cohort"]
    if "evidence_thresholds" in raw:
        cohort_block["evidence_thresholds"] = EvidenceThresholds(**raw["evidence_thresholds"])
    if "penalties" in raw:
        cohort_block["penalties"] = Penalties(**raw["penalties"])

    return CohortConfig.model_validate(cohort_block)
