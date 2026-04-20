"""Cohort-config loader test."""

from __future__ import annotations

from pathlib import Path

from thermocas.config import load_cohort_config

REPO_ROOT = Path(__file__).resolve().parent.parent
COHORT_YAML = REPO_ROOT / "config" / "cohorts" / "brca_example.yaml"


def test_load_brca_example_cohort():
    cfg = load_cohort_config(COHORT_YAML)
    assert cfg.name == "TCGA-BRCA"
    assert cfg.platform == "HM450"
    assert cfg.min_samples_tumor == 30
    assert cfg.evidence_thresholds.exact_bp == 1
    assert cfg.evidence_thresholds.regional_bp == 500
    assert cfg.penalties.heterogeneity_iqr_threshold == 0.30
