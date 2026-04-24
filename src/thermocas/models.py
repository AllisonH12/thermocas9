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


ProbabilisticMode = Literal[
    "tumor_only",
    "tumor_plus_normal_protection",
    "tumor_plus_differential_protection",
    "tumor_plus_gap_sigmoid",
]


class ProbabilisticScore(BaseModel):
    """V2 probabilistic decomposition: each factor is a probability in [0, 1].

    Per-candidate output. The three core factors (p_targ, p_prot, p_trust) are
    always emitted for auditability. Which combination drives
    `p_therapeutic_selectivity` depends on `mode`:

      * `tumor_only`                          : p_sel = p_targ × p_trust
      * `tumor_plus_normal_protection`        : p_sel = p_targ × p_prot × p_trust
      * `tumor_plus_differential_protection`  : p_sel = p_targ × p_diff × p_trust
        (experimental V2.5 — differential protection with a configurable margin)
      * `tumor_plus_gap_sigmoid`              : p_sel = p_targ × p_gap_sigmoid × p_trust
        where p_gap_sigmoid = sigmoid((Δβ − δ) / σ_fixed) with fixed-bandwidth
        σ_fixed (tissue-recommended mode per PAPER.md §5.2.2 gating experiment)

    `tumor_only` is the framework default because `p_protected_normal` encodes
    `P(β_normal > 0.5)`, which is anti-predictive on cohorts where the normal
    comparator doesn't systematically methylate target promoters (TCGA
    adjacent-normal bulk, and even MCF-10A for gene-body probes — Phase 5b
    ablation AUC 0.384 for p_prot alone). Opt into
    `tumor_plus_normal_protection` per cohort only when the normal biology
    actually supports the assumption.

    `tumor_plus_differential_protection` is the V2.5 experimental mode: replaces
    the threshold-based p_prot with `P(β_normal − β_tumor > δ)` via a normal
    approximation on tumor/normal summary statistics. On the MCF-7/MCF-10A
    surrogate at δ=0.2 it raises AUC from V1 0.657/0.628 to 0.705/0.721
    (loose/tight); P@100 on that cohort is tie-band sensitive due to
    p_trust saturation at n=3 and cannot be robustly claimed. Validated on
    one surrogate cohort only; stays opt-in pending higher-n matched
    validation.
    """

    candidate_id: str
    cohort_name: str
    mode: ProbabilisticMode = Field(
        description="Scoring policy — which factors participate in p_therapeutic_selectivity"
    )

    p_targetable_tumor: float = Field(ge=0.0, le=1.0)
    p_protected_normal: float = Field(ge=0.0, le=1.0)
    p_observation_trustworthy: float = Field(ge=0.0, le=1.0)

    p_differential_protection: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description=(
            "V2.5 experimental factor: P(β_normal − β_tumor > differential_delta) "
            "under a normal approximation. Populated only when mode is "
            "`tumor_plus_differential_protection`; None otherwise."
        ),
    )
    differential_delta: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description=(
            "The δ used to compute p_differential_protection (or p_gap_sigmoid). "
            "Populated when mode is `tumor_plus_differential_protection` or "
            "`tumor_plus_gap_sigmoid`; None otherwise."
        ),
    )
    p_gap_sigmoid: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description=(
            "Fixed-bandwidth sigmoid gap factor: sigmoid((Δβ − δ) / σ_fixed). "
            "Populated only when mode is `tumor_plus_gap_sigmoid`; None otherwise. "
            "Paired with `sigma_fixed` on the same record for auditability."
        ),
    )
    sigma_fixed: float | None = Field(
        default=None, gt=0.0,
        description=(
            "The σ_fixed bandwidth used to compute p_gap_sigmoid. Populated only "
            "when mode is `tumor_plus_gap_sigmoid`; None otherwise. Strictly "
            "positive (σ_fixed = 0 would be a hard step and divides by zero "
            "in the sigmoid)."
        ),
    )

    p_therapeutic_selectivity: float = Field(
        ge=0.0, le=1.0,
        description="Composite score computed at scoring time according to `mode`. "
                    "Stored (not a property) so downstream tools see the actual value "
                    "without needing to know the mode.",
    )

    @model_validator(mode="after")
    def _mode_specific_fields_match_mode(self) -> ProbabilisticScore:
        """Per-mode populated fields, iff semantics:

          * `tumor_plus_differential_protection` → requires `p_differential_protection`
            and `differential_delta`. Forbids `p_gap_sigmoid` / `sigma_fixed`.
          * `tumor_plus_gap_sigmoid` → requires `p_gap_sigmoid`, `sigma_fixed`,
            and `differential_delta` (the δ threshold reused from V2.5). Forbids
            `p_differential_protection`.
          * Other modes → forbid all four mode-specific fields.

        This keeps records self-describing: a reviewer can tell from one row
        which mode's math was in effect without loading the cohort YAML."""

        has_p_diff = self.p_differential_protection is not None
        has_delta = self.differential_delta is not None
        has_p_gap = self.p_gap_sigmoid is not None
        has_sigma = self.sigma_fixed is not None

        if self.mode == "tumor_plus_differential_protection":
            if not (has_p_diff and has_delta):
                raise ValueError(
                    "mode=tumor_plus_differential_protection requires both "
                    "p_differential_protection and differential_delta to be set"
                )
            if has_p_gap or has_sigma:
                raise ValueError(
                    "p_gap_sigmoid/sigma_fixed may not be set when "
                    "mode=tumor_plus_differential_protection"
                )
        elif self.mode == "tumor_plus_gap_sigmoid":
            if not (has_p_gap and has_sigma and has_delta):
                raise ValueError(
                    "mode=tumor_plus_gap_sigmoid requires p_gap_sigmoid, "
                    "sigma_fixed, and differential_delta to all be set"
                )
            if has_p_diff:
                raise ValueError(
                    "p_differential_protection may not be set when "
                    "mode=tumor_plus_gap_sigmoid"
                )
        else:
            if has_p_diff or has_delta or has_p_gap or has_sigma:
                raise ValueError(
                    "p_differential_protection/differential_delta/p_gap_sigmoid/"
                    f"sigma_fixed may not be set when mode={self.mode!r}"
                )
        return self

    @model_validator(mode="after")
    def _composite_matches_factors_and_mode(self) -> ProbabilisticScore:
        """`p_therapeutic_selectivity` must equal the factor product implied by
        `mode`. Closes a real integrity gap: downstream ranking reads the
        stored scalar, not the factors. Without this check a malformed JSONL
        row or a future producer bug could silently change rankings while
        still passing validation.

        Comparison tolerance: `math.isclose` with `rel_tol=1e-9, abs_tol=1e-12`
        — tight enough to catch real divergence, loose enough to tolerate the
        float rounding path through JSON round-trip.
        """

        import math

        if self.mode == "tumor_only":
            expected = self.p_targetable_tumor * self.p_observation_trustworthy
        elif self.mode == "tumor_plus_normal_protection":
            expected = (
                self.p_targetable_tumor
                * self.p_protected_normal
                * self.p_observation_trustworthy
            )
        elif self.mode == "tumor_plus_differential_protection":
            # `_mode_specific_fields_match_mode` has already run; p_diff is not None.
            assert self.p_differential_protection is not None
            expected = (
                self.p_targetable_tumor
                * self.p_differential_protection
                * self.p_observation_trustworthy
            )
        elif self.mode == "tumor_plus_gap_sigmoid":
            # `_mode_specific_fields_match_mode` has already run; p_gap_sigmoid is not None.
            assert self.p_gap_sigmoid is not None
            expected = (
                self.p_targetable_tumor
                * self.p_gap_sigmoid
                * self.p_observation_trustworthy
            )
        else:  # pragma: no cover — exhaustive over ProbabilisticMode
            raise ValueError(f"unknown mode for composite check: {self.mode!r}")

        if not math.isclose(
            self.p_therapeutic_selectivity, expected,
            rel_tol=1e-9, abs_tol=1e-12,
        ):
            raise ValueError(
                f"p_therapeutic_selectivity={self.p_therapeutic_selectivity!r} does not "
                f"match the composite implied by mode={self.mode!r}: expected {expected!r}"
            )
        return self


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
    tie_band_size_at_k: int = Field(
        default=0, ge=0,
        description=(
            "Number of records tied at the K-th position's score. When > 1, "
            "top-K membership is not fully determined by score — "
            "`precision_at_k` and `recall_at_k` are sensitive to whichever "
            "secondary tie-break policy the benchmark applied (the shipped "
            "implementation uses candidate_id lexicographic order). Read "
            "P@K with caution when this is large relative to K."
        ),
    )
    tie_break_policy: str = Field(
        default="candidate_id_asc",
        description=(
            "Deterministic secondary sort key used inside tied score bands. "
            "Default: `candidate_id_asc` — lexicographic ascending on "
            "candidate_id. Recorded on the result so reviewers can verify "
            "ranking determinism across runs."
        ),
    )

    precision_at_k_min: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description=(
            "Worst-case `precision_at_k` under adversarial tie-break: if every "
            "tied-band positive that *could* be pushed out of the top-K is "
            "pushed out, how many positives remain in the top-K. Equals "
            "`precision_at_k` when `tie_band_size_at_k == 1`. Read the "
            "[min, max] interval as the tie-band-aware uncertainty on P@K."
        ),
    )
    precision_at_k_max: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description=(
            "Best-case `precision_at_k` under adversarial tie-break: if every "
            "tied-band positive is pulled into the top-K, up to the available "
            "slots. Equals `precision_at_k` when `tie_band_size_at_k == 1`."
        ),
    )
    recall_at_k_min: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description="Worst-case `recall_at_k` counterpart to `precision_at_k_min`.",
    )
    recall_at_k_max: float | None = Field(
        default=None, ge=0.0, le=1.0,
        description="Best-case `recall_at_k` counterpart to `precision_at_k_max`.",
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

    # V2.4 — which factors participate in the probabilistic composite.
    # Default `tumor_only` because Phase 5b ablation showed p_protected_normal
    # inverts on cohorts where the normal comparator doesn't systematically
    # methylate target promoters. Opt into tumor_plus_normal_protection or
    # tumor_plus_differential_protection per cohort only when the biology
    # (or the empirical benchmark) supports the assumption.
    probabilistic_mode: ProbabilisticMode = "tumor_only"

    # V2.5 — δ margin for `tumor_plus_differential_protection` and
    # `tumor_plus_gap_sigmoid`. Default 0.2 matches the offline experiment's
    # sweet spot. Tissue-regime users may prefer δ = 0.1 (§5.3.2).
    differential_delta: float = Field(default=0.2, ge=0.0, le=1.0)

    # `tumor_plus_gap_sigmoid` fixed-bandwidth σ. Default None means: use the
    # package default (√2 × σ_floor ≈ 0.0707, §5.2.1 bandwidth-robust range).
    # Only consulted when probabilistic_mode = "tumor_plus_gap_sigmoid".
    # Strictly positive: σ_fixed = 0 would produce a hard step (and divide by
    # zero in p_gap_sigmoid); use a very small positive value (e.g. 1e-6) for
    # an effectively-hard threshold.
    sigma_fixed: float | None = Field(default=None, gt=0.0)
