"""Pydantic data model for the methylome-guided ThermoCas9 framework.

Every record is JSON / YAML serializable and validated. The data model is the
contract — analysis code in later layers operates only on these types.
"""

from __future__ import annotations

import re
from enum import Enum
from typing import Annotated, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


# ---------- enums ----------


class Strand(str, Enum):
    """Reference-genome strand."""

    PLUS = "+"
    MINUS = "-"


class EvidenceClass(str, Enum):
    """How directly an assay observation maps onto a candidate's critical PAM cytosine.

    Defined in Roth et al. context: methylation arrays measure probe-associated
    CpGs, not arbitrary cytosines. We promote that limitation to a first-class
    record field so scoring can reflect it.
    """

    EXACT = "exact"  # probe assays the exact CpG containing the critical PAM cytosine
    PROXIMAL_CLOSE = "proximal_close"  # nearest probe within ~25 bp
    PROXIMAL = "proximal"  # nearest probe within ~50 bp
    REGIONAL = "regional"  # broader local context only
    UNOBSERVED = "unobserved"  # no usable evidence


# ---------- PAM model ----------


class PamFamily(BaseModel):
    """One ThermoCas9 PAM family in the configured PAM model.

    Validates at load time that:
      * `regex` compiles
      * `critical_c_offset` lies inside the motif width the regex accepts
      * the position pointed to by `critical_c_offset` actually requires a `C`
        (anything else is almost certainly a misconfigured PAM).

    These checks fail fast at YAML load instead of silently producing
    out-of-bounds `critical_c_pos` values during catalog scanning.
    """

    model_config = ConfigDict(frozen=True)

    name: str
    description: str = ""
    regex: str = Field(description="IUPAC/character-class regex for the PAM motif")
    critical_c_offset: int = Field(
        ge=0,
        description="0-indexed position of the methylation-sensing cytosine within the PAM",
    )
    is_cpg: bool = Field(description="True if the critical C is in CpG dinucleotide context")
    weight: float = Field(default=1.0, ge=0.0, le=1.0)

    @field_validator("regex")
    @classmethod
    def _regex_compiles(cls, v: str) -> str:
        try:
            re.compile(v)
        except re.error as e:
            raise ValueError(f"PAM regex {v!r} does not compile: {e}") from e
        return v

    @model_validator(mode="after")
    def _critical_offset_aligned_to_motif(self) -> PamFamily:
        """Validate the PAM regex by enumerating every ACGT-only sequence it
        accepts up to a maximum motif width.

        The framework requires **fixed-width** PAMs because `find_pam_matches`
        computes `critical_c_pos = start + critical_c_offset`. A regex with
        alternatives of different widths (e.g. `A|[ACGT]{4}CG[AG]`) would
        produce out-of-bounds coordinates for the shorter matches, so we
        reject any regex that matches more than one width.

        Width detection is by exhaustive enumeration of ACGT⁴…ACGT⁸, not by
        a hard-coded probe corpus, so legitimate custom PAMs like `CCCCCGA`
        are accepted.
        """

        from itertools import product

        pat = re.compile(self.regex)
        max_width = 8  # 4^8 = 65,536 — covers SpCas9 (3) through ThermoCas9 (7)

        widths_seen: set[int] = set()
        matches_by_width: dict[int, list[str]] = {}
        for L in range(1, max_width + 1):
            examples: list[str] = []
            for tup in product("ACGT", repeat=L):
                s = "".join(tup)
                if pat.fullmatch(s):
                    examples.append(s)
            if examples:
                widths_seen.add(L)
                matches_by_width[L] = examples

        if not widths_seen:
            raise ValueError(
                f"PAM regex {self.regex!r} matches no ACGT-only sequence of "
                f"length 1..{max_width}; framework requires a fixed-width DNA PAM"
            )

        if len(widths_seen) > 1:
            raise ValueError(
                f"PAM regex {self.regex!r} matches multiple widths {sorted(widths_seen)}; "
                f"framework requires a fixed-width PAM (use one regex per width)"
            )

        width = widths_seen.pop()

        if not (0 <= self.critical_c_offset < width):
            raise ValueError(
                f"critical_c_offset {self.critical_c_offset} is outside motif width "
                f"{width} for regex {self.regex!r}"
            )

        # Every concrete match must have C at critical_c_offset; otherwise the
        # "methylation-sensing cytosine" claim is wrong for some matches.
        chars_at_offset = {s[self.critical_c_offset] for s in matches_by_width[width]}
        if chars_at_offset != {"C"}:
            raise ValueError(
                f"critical_c_offset {self.critical_c_offset} in regex {self.regex!r} "
                f"accepts {sorted(chars_at_offset)}, not just C — the methylation "
                f"sensing position must require a cytosine in every match"
            )
        return self


class PamMatch(BaseModel):
    """Concrete PAM occurrence in a sequence."""

    model_config = ConfigDict(frozen=True)

    family: str
    sequence: str
    start: int = Field(description="0-indexed start coordinate of the PAM in sequence space")
    strand: Strand
    critical_c_pos: int = Field(
        description="0-indexed coordinate of the critical PAM cytosine in sequence space"
    )
    is_cpg: bool
    weight: float = Field(ge=0.0, le=1.0)


# ---------- candidate site ----------


class CandidateSite(BaseModel):
    """One ThermoCas9-compatible PAM site in a reference genome.

    Cohort-independent. Methylation evidence and scores are attached separately.
    """

    candidate_id: str = Field(description="Stable identifier, e.g. 'chr10:8045463+:NNNNCGA'")
    chrom: str
    critical_c_pos: int = Field(
        ge=0, description="0-indexed reference coordinate of the methylation-sensing cytosine"
    )
    strand: Strand
    pam: str = Field(description="The matched PAM sequence as it appears on the indicated strand")
    pam_family: str = Field(description="Name of the PAM family in the configured PAM model")
    is_cpg_pam: bool
    local_seq_100bp: str = Field(
        default="",
        description="±50 bp of reference sequence around the critical C, on the indicated strand",
    )
    nearest_gene: str | None = None
    regulatory_context: str | None = Field(
        default=None,
        description="Free-form annotation, e.g. 'promoter', 'enhancer', 'gene_body'",
    )

    @field_validator("local_seq_100bp")
    @classmethod
    def _seq_uppercase_only(cls, v: str) -> str:
        if not all(c in "ACGTN" for c in v.upper()):
            raise ValueError("local_seq_100bp must contain only ACGTN characters")
        return v.upper()


# ---------- methylation evidence ----------


class MethylationObservation(BaseModel):
    """Tumor-vs-normal methylation summary for one candidate site in one cohort.

    Quantiles are stored separately from the mean because the scoring formula
    uses class separation (q25/q75), not just mean shift.
    """

    candidate_id: str
    cohort_name: str

    evidence_class: EvidenceClass
    evidence_distance_bp: int | None = Field(
        default=None,
        ge=0,
        description="Distance from critical C to nearest informative probe (bp). None if unobserved.",
    )
    probe_id: str | None = None

    # Beta value summaries — None when no samples are available.
    beta_tumor_mean: float | None = Field(default=None, ge=0.0, le=1.0)
    beta_tumor_q25: float | None = Field(default=None, ge=0.0, le=1.0)
    beta_tumor_q75: float | None = Field(default=None, ge=0.0, le=1.0)
    n_samples_tumor: int = Field(default=0, ge=0)

    beta_normal_mean: float | None = Field(default=None, ge=0.0, le=1.0)
    beta_normal_q25: float | None = Field(default=None, ge=0.0, le=1.0)
    beta_normal_q75: float | None = Field(default=None, ge=0.0, le=1.0)
    n_samples_normal: int = Field(default=0, ge=0)

    @model_validator(mode="after")
    def _quantile_ordering(self) -> MethylationObservation:
        for q25, q75, label in (
            (self.beta_tumor_q25, self.beta_tumor_q75, "tumor"),
            (self.beta_normal_q25, self.beta_normal_q75, "normal"),
        ):
            if q25 is not None and q75 is not None and q25 > q75:
                raise ValueError(f"{label} q25 ({q25}) must be <= q75 ({q75})")
        return self

    @model_validator(mode="after")
    def _mean_within_quantiles(self) -> MethylationObservation:
        for mean, q25, q75, label in (
            (self.beta_tumor_mean, self.beta_tumor_q25, self.beta_tumor_q75, "tumor"),
            (self.beta_normal_mean, self.beta_normal_q25, self.beta_normal_q75, "normal"),
        ):
            if mean is None or q25 is None or q75 is None:
                continue
            if not (q25 <= mean <= q75):
                raise ValueError(
                    f"{label} mean ({mean}) must lie within [q25={q25}, q75={q75}]"
                )
        return self

    @model_validator(mode="after")
    def _unobserved_carries_no_evidence(self) -> MethylationObservation:
        """UNOBSERVED records must be empty — no probe, no distance, no betas, no samples.

        Otherwise scoring and serialization can leak impossible state (e.g. an
        UNOBSERVED record reporting a probe_id and a non-zero distance).
        """

        if self.evidence_class != EvidenceClass.UNOBSERVED:
            return self

        forbidden_fields = {
            "probe_id": self.probe_id,
            "evidence_distance_bp": self.evidence_distance_bp,
            "beta_tumor_mean": self.beta_tumor_mean,
            "beta_tumor_q25": self.beta_tumor_q25,
            "beta_tumor_q75": self.beta_tumor_q75,
            "beta_normal_mean": self.beta_normal_mean,
            "beta_normal_q25": self.beta_normal_q25,
            "beta_normal_q75": self.beta_normal_q75,
        }
        set_fields = sorted(k for k, v in forbidden_fields.items() if v is not None)
        if set_fields:
            raise ValueError(
                f"UNOBSERVED records must not set: {set_fields}"
            )
        if self.n_samples_tumor != 0 or self.n_samples_normal != 0:
            raise ValueError(
                "UNOBSERVED records must have n_samples_tumor=n_samples_normal=0; "
                f"got tumor={self.n_samples_tumor}, normal={self.n_samples_normal}"
            )
        return self


# ---------- scoring ----------


Score = Annotated[float, Field(ge=0.0)]


class ScoreComponents(BaseModel):
    """Decomposed components of a candidate's targetability score.

    Stored alongside the final score so reviewers can audit ranking decisions.
    """

    sequence_score: Score = Field(description="From the PAM family weight; reflects motif compatibility")
    selectivity_score: Score = Field(
        description="max(0, normal_mean - tumor_mean) + max(0, normal_q25 - tumor_q75)"
    )
    confidence_score: Score = Field(
        description="In [0, 1]; from the EvidenceClass of the methylation observation"
    )
    heterogeneity_penalty: Score = Field(default=0.0)
    low_coverage_penalty: Score = Field(default=0.0)


class ScoredCandidate(BaseModel):
    """Final per-cohort scored candidate."""

    candidate: CandidateSite
    observation: MethylationObservation
    components: ScoreComponents
    final_score: float = Field(
        description=(
            "sequence × selectivity × confidence − heterogeneity_penalty − low_coverage_penalty. "
            "Can go negative when penalties dominate."
        )
    )
    probabilistic: "ProbabilisticScore | None" = Field(
        default=None,
        description=(
            "Optional V2 probabilistic decomposition. None when score_candidate "
            "is invoked without compute_probabilistic=True."
        ),
    )
    spacer: "SpacerScore | None" = Field(
        default=None,
        description=(
            "Optional V3 protospacer (gRNA) scoring. None when score_candidate "
            "is invoked without compute_spacer=True."
        ),
    )

    @model_validator(mode="after")
    def _candidate_id_consistent(self) -> ScoredCandidate:
        if self.candidate.candidate_id != self.observation.candidate_id:
            raise ValueError(
                "ScoredCandidate.candidate.candidate_id must match observation.candidate_id"
            )
        if self.probabilistic is not None:
            if self.probabilistic.candidate_id != self.candidate.candidate_id:
                raise ValueError(
                    "ScoredCandidate.probabilistic.candidate_id must match candidate.candidate_id"
                )
            if self.probabilistic.cohort_name != self.observation.cohort_name:
                raise ValueError(
                    "ScoredCandidate.probabilistic.cohort_name must match "
                    "observation.cohort_name; got "
                    f"{self.probabilistic.cohort_name!r} vs {self.observation.cohort_name!r}"
                )
        if self.spacer is not None and self.spacer.candidate_id != self.candidate.candidate_id:
            raise ValueError(
                "ScoredCandidate.spacer.candidate_id must match candidate.candidate_id"
            )
        return self


# ---------- gRNA spacer (V3) ----------


class SpacerScore(BaseModel):
    """V3 — protospacer (gRNA) design quality for one candidate.

    Computed from the 20-nt sequence immediately upstream of the PAM on the
    candidate's strand. Each component sits in [0, 1] (higher = better);
    `final_score` is their geometric mean so a single bad component pulls the
    whole spacer down.

    Components capture standard sgRNA design heuristics:
      * gc_content_score    : favors 40–60% GC
      * tm_score            : favors melting temp 60–75 °C
      * runs_score          : penalizes long mononucleotide runs (≥5)
      * hairpin_score       : penalizes 4-bp self-complementary palindromes
    """

    candidate_id: str
    spacer_seq: str = Field(description="20-nt protospacer (5'→3' on strand)")
    gc_fraction: float = Field(ge=0.0, le=1.0)
    tm_C: float = Field(description="Estimated melting temperature in °C")
    longest_run: int = Field(ge=0, description="Longest mononucleotide run length")
    gc_content_score: float = Field(ge=0.0, le=1.0)
    tm_score: float = Field(ge=0.0, le=1.0)
    runs_score: float = Field(ge=0.0, le=1.0)
    hairpin_score: float = Field(ge=0.0, le=1.0)
    final_score: float = Field(
        ge=0.0, le=1.0,
        description="Geometric mean of the four component scores",
    )

    @field_validator("spacer_seq")
    @classmethod
    def _spacer_dna_only(cls, v: str) -> str:
        v = v.upper()
        if not all(c in "ACGTN" for c in v):
            raise ValueError("spacer_seq must contain only ACGTN")
        return v


# ---------- probabilistic (V2) ----------


class ProbabilisticScore(BaseModel):
    """V2 probabilistic decomposition: each factor is a probability in [0, 1].

    P(therapeutic_selectivity) = P(targetable_in_tumor)
                               × P(protected_in_normal)
                               × P(observation_trustworthy)
    """

    candidate_id: str
    cohort_name: str

    p_targetable_tumor: float = Field(ge=0.0, le=1.0)
    p_protected_normal: float = Field(ge=0.0, le=1.0)
    p_observation_trustworthy: float = Field(ge=0.0, le=1.0)

    @property
    def p_therapeutic_selectivity(self) -> float:
        return (
            self.p_targetable_tumor * self.p_protected_normal * self.p_observation_trustworthy
        )


# ---------- cohort config ----------


class EvidenceThresholds(BaseModel):
    """Distance bins (bp) used by the evidence classifier.

    EXACT means the probe assays the *same CpG* as the critical PAM cytosine —
    that's distance 0. The previous default of `exact_bp=1` upgraded
    one-base-away probes from PROXIMAL_CLOSE (0.7 confidence) to EXACT (1.0),
    which doesn't match the documented contract. Tighten to 0.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    exact_bp: int = Field(default=0, ge=0)
    proximal_close_bp: int = Field(default=25, ge=0)
    proximal_bp: int = Field(default=50, ge=0)
    regional_bp: int = Field(default=500, ge=0)

    @model_validator(mode="after")
    def _monotonic(self) -> EvidenceThresholds:
        # exact < proximal_close < proximal < regional
        if not (self.exact_bp <= self.proximal_close_bp <= self.proximal_bp <= self.regional_bp):
            raise ValueError(
                "Evidence thresholds must satisfy exact <= proximal_close <= proximal <= regional"
            )
        return self


class Penalties(BaseModel):
    """Per-cohort penalty weights."""

    model_config = ConfigDict(frozen=True, extra="forbid")

    heterogeneity_iqr_threshold: float = Field(default=0.30, ge=0.0, le=1.0)
    heterogeneity_weight: float = Field(default=0.4, ge=0.0)
    low_coverage_n_threshold: int = Field(default=30, ge=0)
    low_coverage_weight: float = Field(default=0.3, ge=0.0)


Platform = Literal["HM450", "EPIC", "RRBS", "WGBS", "CUSTOM"]


# ---------- pan-cancer aggregation ----------


class PanCancerAggregate(BaseModel):
    """Cross-cohort summary for a single candidate site.

    Produced by `pan_cancer.aggregate(...)` after each cohort has been scored.
    """

    candidate_id: str
    chrom: str
    critical_c_pos: int
    pam_family: str

    n_cohorts_observed: int = Field(ge=0, description="Cohorts where evidence_class != UNOBSERVED")
    n_cohorts_high_score: int = Field(
        ge=0, description="Cohorts where final_score >= high_score_threshold"
    )
    cohort_scores: dict[str, float] = Field(
        default_factory=dict, description="cohort_name → final_score"
    )

    pan_cancer_score: float = Field(
        description="Mean of cohort_scores across observed cohorts; 0 when none observed."
    )
    recurrence: float = Field(
        ge=0.0, le=1.0,
        description="n_cohorts_high_score / max(1, n_cohorts_observed)",
    )
    exclusivity: float = Field(
        ge=0.0,
        description=(
            "Standard deviation of cohort_scores across observed cohorts. "
            "Higher = more cohort-specific (one tumor type addressable, others not). "
            "Lower = pan-cancer recurrent."
        ),
    )
    normal_protection_max: float | None = Field(
        default=None,
        description=(
            "Max β_normal_mean across observed cohorts. **Best-case** protection: "
            "the cohort where normal tissue is most methylated (least exposed). "
            "Use `normal_protection_min` for the worst-case / pan-cancer risk view."
        ),
    )
    normal_protection_min: float | None = Field(
        default=None,
        description=(
            "Min β_normal_mean across observed cohorts. **Worst-case** protection: "
            "the cohort where normal tissue is *least* methylated (most exposed to "
            "off-tumor editing). The honest pan-cancer risk metric."
        ),
    )


# ---------- benchmark (V3) ----------


class BenchmarkResult(BaseModel):
    """V3 — ranking-quality metrics for one cohort run against a positives list.

    Produced by `benchmark.evaluate_ranking(...)`. Captures the standard
    "did the prioritization work?" question that turns the framework from a
    descriptive tool into a methods paper.

    All metrics are well-defined when at least one positive is in the scored
    set; otherwise they're None.
    """

    cohort_name: str
    n_total: int = Field(ge=0, description="Total scored candidates evaluated")
    n_positives: int = Field(ge=0, description="Positives present in the scored set")
    n_negatives: int = Field(ge=0, description="Negatives present in the scored set")
    held_out_chromosomes: list[str] = Field(
        default_factory=list,
        description="Chromosomes excluded from training; empty for full-set eval",
    )
    top_k: int = Field(ge=1, description="K used for precision/recall@K")
    precision_at_k: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description="Fraction of top-K that are positives",
    )
    recall_at_k: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description="Fraction of all positives captured in the top-K",
    )
    roc_auc: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description="Area under the ROC curve (Mann-Whitney U formulation)",
    )


class CohortConfig(BaseModel):
    """One cohort spec — the unit a Cohort Adapter consumes.

    Unknown keys are rejected (`extra='forbid'`). YAML typos like
    `min_samples_tumour` (UK spelling) used to silently fall back to defaults
    while the user thought their override was applied. Now they error at load.
    """

    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    tumor_dataset: str
    normal_dataset: str
    platform: Platform
    subtype_field: str | None = None
    min_samples_tumor: int = Field(default=30, ge=1)
    min_samples_normal: int = Field(default=10, ge=1)
    evidence_thresholds: EvidenceThresholds = EvidenceThresholds()
    penalties: Penalties = Penalties()
