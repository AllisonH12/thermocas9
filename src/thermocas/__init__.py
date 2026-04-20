"""Methylome-guided ThermoCas9 target-site discovery framework.

Roth et al., Nature 2026 — DOI 10.1038/s41586-026-10384-z.
"""

__version__ = "0.4.0"

from thermocas.benchmark import evaluate_ranking, split_by_chrom
from thermocas.catalog import CatalogStats, build_catalog, stream_catalog
from thermocas.cohort import score_cohort
from thermocas.config import load_cohort_config
from thermocas.evidence import EvidenceClassifier, ProbeRecord, classify_evidence
from thermocas.grna import SPACER_LEN, extract_spacer, score_spacer
from thermocas.io import (
    iter_fasta,
    read_beta_matrix,
    read_jsonl,
    read_sample_subtypes,
    read_tsv,
    write_jsonl,
    write_jsonl_atomic,
)
from thermocas.methylation_backend import (
    BetaSummary,
    GDCBackend,
    LocalArrayBackend,
    LocalSummaryBackend,
    MethylationBackend,
)
from thermocas.models import (
    BenchmarkResult,
    CandidateSite,
    CohortConfig,
    EvidenceClass,
    EvidenceThresholds,
    MethylationObservation,
    PamFamily,
    PamMatch,
    PanCancerAggregate,
    Penalties,
    ProbabilisticMode,
    ProbabilisticScore,
    ScoreComponents,
    ScoredCandidate,
    SpacerScore,
    Strand,
)
from thermocas.pam_model import PamModel, find_pam_matches
from thermocas.pan_cancer import (
    DEFAULT_HIGH_SCORE,
    aggregate,
    top_exclusive,
    top_recurrent,
)
from thermocas.probabilistic import (
    DEFAULT_METHYLATED_THRESHOLD,
    DEFAULT_TRUST_RAMP_N,
    DEFAULT_UNMETHYLATED_THRESHOLD,
    p_observation_trustworthy,
    p_protected_normal,
    p_targetable_tumor,
    probabilistic_score,
)
from thermocas.scoring import (
    confidence_score_for_evidence,
    score_candidate,
    selectivity_score,
)

__all__ = [
    # core
    "BetaSummary",
    "CandidateSite",
    "CatalogStats",
    "CohortConfig",
    "DEFAULT_HIGH_SCORE",
    "DEFAULT_METHYLATED_THRESHOLD",
    "DEFAULT_TRUST_RAMP_N",
    "DEFAULT_UNMETHYLATED_THRESHOLD",
    "EvidenceClass",
    "EvidenceClassifier",
    "EvidenceThresholds",
    "GDCBackend",
    "LocalArrayBackend",
    "LocalSummaryBackend",
    "MethylationBackend",
    "MethylationObservation",
    "PamFamily",
    "PamMatch",
    "PamModel",
    "PanCancerAggregate",
    "Penalties",
    "ProbabilisticMode",
    "ProbabilisticScore",
    "ProbeRecord",
    "ScoreComponents",
    "ScoredCandidate",
    "Strand",
    # functions
    "aggregate",
    "build_catalog",
    "classify_evidence",
    "confidence_score_for_evidence",
    "find_pam_matches",
    "iter_fasta",
    "load_cohort_config",
    "p_observation_trustworthy",
    "p_protected_normal",
    "p_targetable_tumor",
    "probabilistic_score",
    "read_beta_matrix",
    "read_jsonl",
    "read_sample_subtypes",
    "read_tsv",
    "score_candidate",
    "score_cohort",
    "selectivity_score",
    "stream_catalog",
    "top_exclusive",
    "top_recurrent",
    "write_jsonl",
    "write_jsonl_atomic",
]
