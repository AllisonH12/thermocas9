"""Cohort-config loader test."""

from __future__ import annotations

from pathlib import Path

import pytest

from thermocas.config import load_cohort_config

REPO_ROOT = Path(__file__).resolve().parent.parent
COHORT_YAML = REPO_ROOT / "config" / "cohorts" / "brca_example.yaml"


def test_load_brca_example_cohort():
    cfg = load_cohort_config(COHORT_YAML)
    assert cfg.name == "TCGA-BRCA"
    assert cfg.platform == "HM450"
    assert cfg.min_samples_tumor == 30
    assert cfg.evidence_thresholds.exact_bp == 0  # V3.1: tightened from 1 → 0
    assert cfg.evidence_thresholds.regional_bp == 500
    assert cfg.penalties.heterogeneity_iqr_threshold == 0.30


def test_load_brca_example_exposes_probabilistic_mode_and_delta():
    """V2.5 — YAML must parse `probabilistic_mode` and `differential_delta`
    into the CohortConfig. The example file pins mode=tumor_only, δ=0.2."""
    cfg = load_cohort_config(COHORT_YAML)
    assert cfg.probabilistic_mode == "tumor_only"
    assert cfg.differential_delta == pytest.approx(0.2)


def test_load_differential_mode_yaml(tmp_path: Path):
    """A cohort YAML that opts into the V2.5 mode with an explicit δ must
    round-trip cleanly through the loader."""
    yaml_text = """
cohort:
  name: demo
  tumor_dataset: demo-tumor
  normal_dataset: demo-normal
  platform: HM450
  min_samples_tumor: 30
  min_samples_normal: 10
  probabilistic_mode: tumor_plus_differential_protection
  differential_delta: 0.35
"""
    p = tmp_path / "demo.yaml"
    p.write_text(yaml_text)
    cfg = load_cohort_config(p)
    assert cfg.probabilistic_mode == "tumor_plus_differential_protection"
    assert cfg.differential_delta == pytest.approx(0.35)


def test_load_rejects_unknown_probabilistic_mode(tmp_path: Path):
    """Typos must fail at load, not silently fall back."""
    from pydantic import ValidationError

    yaml_text = """
cohort:
  name: demo
  tumor_dataset: x
  normal_dataset: y
  platform: HM450
  min_samples_tumor: 30
  min_samples_normal: 10
  probabilistic_mode: tumor_plus_bogus
"""
    p = tmp_path / "bad.yaml"
    p.write_text(yaml_text)
    with pytest.raises(ValidationError):
        load_cohort_config(p)
