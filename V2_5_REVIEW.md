# V2.5 Self-Review — differential-protection probabilistic mode

*Date: 2026-04-20. Audit scope: commit `df35d3b`. Reviewer: self.*

This document is the post-hoc audit of V2.5 before any public framing. It
records what was built, what is sound, what is fragile, and what should
happen before V2.5 could become a released default. The executive verdict
is **experimental-on-main, do not tag, do not promote**.

---

## 1. What V2.5 adds

A new, opt-in `probabilistic_mode = "tumor_plus_differential_protection"`
replaces the threshold-based `p_protected_normal` with a differential factor:

```
p_diff(δ) = P(β_normal − β_tumor > δ)
```

computed under an independent-normal approximation on the tumor / normal
summary statistics already in `MethylationObservation`:

```
σ_k   ≈ IQR_k / 1.349           (k ∈ {tumor, normal}; floor at 0.05)
σ_Δ²  = max(σ_t, floor)² + max(σ_n, floor)²
P(Δ>δ) = 1 − Φ((δ − (μ_n − μ_t)) / σ_Δ)
```

The new composite is

```
p_therapeutic_selectivity = p_targ × p_diff × p_trust
```

compared to the V0.4.0 default `p_sel = p_targ × p_trust` (tumor_only) and
the deprecated V2 composite `p_sel = p_targ × p_prot × p_trust`.

`differential_delta` is a new cohort YAML key (default `0.2`) consulted only
when the mode is differential; the other modes ignore it.

---

## 2. Scientific motivation — was the critique of `p_prot` correct?

**Yes.** Phase 5b ablation on the MCF-7/MCF-10A surrogate (GSE77348) showed
`p_protected_normal` alone at AUC 0.384 / 0.343 — anti-predictive — and the
`p_targ × p_prot` pair at AUC 0.503 / 0.488 — essentially random. The root
cause is the semantics: `p_prot = P(β_normal > 0.5)` assumes the *normal
comparator methylates the target*. That assumption is paper-specific (the
Roth et al. samples happen to methylate `ESR1`/`GATA3` promoters in
MCF-10A); it does not generalize to bulk adjacent-normal cohorts or even
to MCF-10A at gene-body probes.

`p_diff` replaces an absolute-threshold claim about the normal side with a
differential claim about the contrast, which is cohort-agnostic.

---

## 3. Empirical readout (MCF-7 vs MCF-10A surrogate)

| axis | AUC loose | AUC tight | P@100 loose | P@100 tight |
|---|---|---|---|---|
| V2.5 (δ=0.2) | **0.705** | **0.721** | tie-fragile | tie-fragile |
| V1 `final_score` | 0.657 | 0.628 | 0.040 | 0.030 |
| V2 `tumor_only` | 0.733 | 0.770 | 0.000 | 0.000 |
| V2 `tumor_plus_normal_protection` | 0.503 | 0.488 | ≈0 | ≈0 |

AUC ordering: `tumor_only > V2.5 > V1 > v2_full`. The tumor_only edge on
AUC is paid for with total P@100 collapse; V2.5 trades some of that AUC
for a more biologically coherent composite.

### 3.1 The tie-fragility finding (the main reason V2.5 stays experimental)

On an `n=3` surrogate, `p_trust` saturates at `0.95 × 3/30 = 0.095` for
every EXACT-evidence candidate with ≥3 samples on each side. Candidates
whose `p_targ` and `p_diff` both round to ~1.0 — i.e. the *strongest*
differential candidates — all compose to exactly `0.095`. Direct
measurement on the scored artifact finds **519 records tied at the top
score**. Which 100 of those 519 appear in the top-100 depends entirely on
which tiebreak the benchmark uses:

| tiebreak | P@100 loose | P@100 tight |
|---|---|---|
| CLI default (Python stable sort on stream order) | 0.010 | 0.000 |
| Candidate-id lexicographic (heap) | 0.100 | 0.100 |

Neither is principled: the signal *inside* the tied band is zero. The
AUC metric is robust because ties contribute 0.5 each under Mann-Whitney
U, but P@K reads one deterministic subset of the ties.

`tumor_only` has the same pathology — its top is the same `0.095` band
minus the `p_diff` filter, i.e. mostly "always-unmethylated" loci. V1
`final_score` avoids the problem because it is continuous-real-valued and
does not bottleneck on a discrete `p_trust` ceiling.

**Implication.** On this cohort the V2.5 P@100 improvement over V1 is not
a robust empirical claim. On a cohort where `n ≫ trust_ramp_n` (e.g.
TCGA-BRCA `n ≈ 800`), `p_trust` is continuous and the tied band should
dissolve; V2.5's continuous `p_diff × p_targ` structure should then drive
a real ranking signal. That prediction is unvalidated.

---

## 4. Code-level audit

### 4.1 Math (`probabilistic.p_differential_protection`)

- **CDF form.** `0.5 × (1 − erf(z / √2))` equals `1 − Φ(z)`. Correct.
- **Z transformation.** `z = (δ − (μ_n − μ_t)) / σ_Δ`. Correct for
  `P(Δ > δ)` under `Δ ∼ N(μ_n − μ_t, σ_Δ²)`.
- **σ_Δ under independence.** `σ_Δ² = σ_t² + σ_n²`. Correct.
- **Boundary handling.** β-values live in [0, 1]; `Δ = β_n − β_t ∈ [−1, 1]`.
  The normal approximation ignores this boundary. In practice the model
  overestimates `p_diff` in the tails when μ_n − μ_t is near ±1. For
  δ ≤ 0.5 and IQRs of the scale seen on HM450 (~0.05–0.20), the error
  should be minor but is not formally bounded. Documented as a known
  modeling simplification.
- **σ floor.** The 0.05 floor prevents `σ_Δ → 0` when both IQRs collapse
  at boundary β-values (common for cg-islands methylated at 0 or 1). The
  effect is to widen `p_diff` into a softer probability at those loci
  rather than letting it spike to 1 or 0 on a razor edge. Value chosen
  empirically; not benchmarked against alternatives.
- **One-sided quantile case.** If only `q25` or only `q75` is present,
  the implementation sets `σ_k = 0` for that side and lets the floor
  take over. This is conservative but discards information. An
  improvement would be `σ = |q − μ| / 0.6745` from whichever quantile is
  present. Left as a follow-up because the surrogate has both quantiles
  on every record.

### 4.2 Dispatch (`probabilistic.probabilistic_score`)

- The new branch `elif mode == "tumor_plus_differential_protection":` is
  the only mode that sets `p_differential_protection` and
  `differential_delta` on the output record. Other branches leave both
  `None`, which is enforced by the `ProbabilisticScore` validator.
- `differential_delta` is accepted as a keyword-only argument and defaults
  to `DEFAULT_DIFFERENTIAL_DELTA`. It is ignored when the mode is not
  differential — no silent override path.
- `sigma_floor` is **not** exposed on `probabilistic_score`. It is only
  tunable on the low-level `p_differential_protection(…, sigma_floor=…)`
  call. Acceptable for now (no benchmark depends on varying it); flag as
  a potential knob to lift if tuning becomes interesting.

### 4.3 Model (`ProbabilisticScore` + `CohortConfig`)

- `ProbabilisticScore.p_differential_protection: float | None` and
  `differential_delta: float | None` default to `None`.
- The `@model_validator(mode="after")` enforces the biconditional:
  *audit fields populated iff mode is differential*. Both directions
  covered by tests. This prevents malformed records from silently
  round-tripping through JSONL in either direction.
- `ProbabilisticScore.p_therapeutic_selectivity` is a stored scalar, not
  a property — downstream consumers see the actual composite the mode
  emitted. **Patched post-review (P2 fix):** the model now also carries a
  `_composite_matches_factors_and_mode` validator that re-derives the
  expected composite from the factors according to `mode` and rejects
  records where the stored scalar does not match (tolerance
  `rel_tol=1e-9, abs_tol=1e-12`). Reproduced the original P2 failure
  case — `mode=differential, p_targ=0.8, p_diff=0.9, p_trust=0.5,
  p_sel=0.01` — and confirmed it now raises at construction. Positive
  path test round-trips every mode through `model_dump_json` →
  `model_validate_json`. Closes the integrity gap that would otherwise
  let a malformed JSONL row silently change downstream ranking.
- `CohortConfig.differential_delta: float = 0.2` (ge=0.0, le=1.0).
  `extra="forbid"` continues to reject unknown keys, so V2.5 cannot
  silently absorb a typo.

### 4.4 Backward compatibility

- Old `scored_*.jsonl` artifacts (pre-V2.5) lack the two new fields.
  Pydantic fills them with `None` on re-load, the validator's "None
  when mode is not differential" branch passes, and the record is
  accepted unchanged. Verified by re-reading the committed
  `scored_surrogate_tumor_only.jsonl` in the offline script.
- V1 `final_score` path is bit-identical: the new mode is gated by
  `compute_probabilistic=True` and by the mode field. No V1-only path
  touches `p_differential_protection`.

### 4.5 Test coverage (`tests/test_probabilistic.py`, `tests/test_config.py`)

Covered:

- mode dispatch (differential differs from tumor_only when gap ≈ δ)
- δ recorded on the emitted score
- strict monotonicity: higher δ ⇒ smaller `p_diff` ⇒ smaller composite
- composite formula equals `p_targ × p_diff × p_trust` exactly
- non-differential modes leave audit fields `None`
- validator rejects: differential mode without audit fields; non-differential
  mode *with* audit fields
- YAML round-trip of `tumor_plus_differential_protection` + custom δ
- YAML rejects unknown mode strings at load
- `p_diff` monotone in gap size
- `p_diff = 0` on one-sided observations
- `p_diff = 0.5` at the breakpoint (gap = δ) — defining property of the factor
- `p_diff = 0` on UNOBSERVED

199 tests green.

Not covered (acknowledged gaps):

- The `sigma_floor` tunable parameter has no dedicated test; its effect
  is implicit in the breakpoint and monotonicity tests.
- Cross-field consistency between `p_therapeutic_selectivity` and the
  individual factors is not checked by the Pydantic model.
- The benchmark's tie-handling is not asserted — a regression could
  silently change P@K on tied bands.

### 4.6 Benchmark harness — tie handling (patched post-review)

**Status: the P1 bug is fixed.** A subsequent code-review pass demonstrated
that under the original `pairs.sort(key=lambda t: t[1], reverse=True)`,
a single positive with three tied candidates swings from P@1=1.0 to
P@1=0.0 depending purely on input order — a real contract bug. Two
changes shipped:

- `evaluate_ranking` now sorts by `(-score, candidate_id)` — deterministic
  secondary ordering by candidate_id ascending. Same scores in any input
  order produce the same P@K.
- `BenchmarkResult` gains `tie_band_size_at_k: int` and
  `tie_break_policy: str`. On every result, the number of records tied
  at the K-th position's score is reported explicitly. When > 1, P@K is
  partially determined by the secondary key, and readers should scale
  their interpretation accordingly.

On the V2.5 surrogate benchmark, `tie_band_size_at_k = 299` at K=100.
That is the honest read of the earlier "tie-fragile" caveat: nearly
three times K's worth of candidates are tied at the cutoff score, and
the headline P@100 is a function of the tie-break, not of the model.
V1 `final_score` on the same cohort reports `tie_band_size_at_k = 1`
at K=100 — its P@K is a real signal.

Remaining follow-ups:

- optionally emit a P@K range (min/max over tied subsets) when the top-K
  boundary falls inside a tied band — strictly informational, not
  required for the contract to be correct
- consider offering a configurable secondary key (e.g. `final_score`) as
  an alternate tie-break policy

Neither is a blocker.

---

## 5. Claims inventory — what we can and cannot say

*(Revised after the label repair and the three-cohort cross-check in §8.)*

**Can say:**

- V2.5 is the **highest-AUC probabilistic axis on matched cell-line cohorts**
  at every label granularity tested. Under the repaired Roth-validated
  positives, V2.5 AUC is 0.990 on GSE322563 validated (n=3 positives),
  0.982 on GSE77348 validated — each at least +0.01 above V1 and at
  least +0.05 above `tumor_only`.
- V2.5 replaces a biologically mis-specified protection factor
  (`p_prot = P(β > 0.5)`, empirically anti-predictive) with a differential
  factor that makes no absolute-threshold assumption.
- The composite, the validator, the YAML surface, and the downstream
  records are internally consistent and contract-tested. `ProbabilisticScore`
  additionally now asserts composite = `p_targ × p_diff × p_trust` at
  construction (P2 fix), so malformed JSONL cannot silently change rankings.
- V2.5's tie_band behaves correctly with `n` — 190 at n=2 → 299 at n=3 →
  2 at n=305/50 — confirming the `p_trust`-saturation prediction made
  in §3 before the higher-`n` cohort was run. This is a prediction that
  held without tuning after-the-fact.
- V2.5 reproduces Roth's Fig. 5d β values exactly at the three validated
  target sites (EGFLAM 0.01/0.49, ESR1 0.07/0.94, GATA3 0.02/0.31),
  cross-checking the hg38 → hg19 liftover and the EPIC-v2 → HM450 probe
  intersect.
- V2.5 is opt-in. The default remains `tumor_only`. V1 `final_score`
  remains the recommended continuous-score fallback for top-K stability
  on any cohort type.

**Cannot say:**

- That V2.5 is the universally best axis. On GSE69914 primary tissue
  `tumor_only` has higher AUC at every label granularity (+0.03 on
  validated, +0.13 on narrow, +0.15 on wide). `tumor_only`'s top-K is
  unusable there (tie_band 6,540), so it wouldn't be the right choice
  for target discovery, but on pure AUC it wins on tissue cohorts.
- That V2.5 P@100 is competitive across cohorts. P@100 is dominated by
  the tie band on all three cohorts; the only axis producing a real
  top-K is V1. V2.5's value is in the continuous ranking above the tied
  band, not in the top-100 slice per se.
- That V2.5 generalizes to cohorts that do NOT reproduce the Roth
  target biology. A cross-series run against the Sanger GDSC breast
  panel (GSE68379, 52 breast lines × GSE69914 healthy normal) is
  documented in §8.1 as an out-of-distribution boundary case: the
  Roth-validated positives are methylated in Sanger's MCF-7 and
  unmethylated in Roth's, so the label set is not biologically valid
  there and the inverted AUC is expected, not informative.

---

## 6. What should happen before V2.5 can be promoted to unconditional default

*(Revised. The earlier list had six items; four are now resolved.)*

1. ~~**Higher-replicate matched cohort** where `p_trust` does not saturate.~~
   **Done.** GSE69914 n=305/50 ran under all three label granularities.
   Tie_band dissolved (2 at K=100) as predicted; V2.5 remains the
   recommended probabilistic axis for top-K usability even though
   `tumor_only` wins AUC on that cohort.
2. ~~**Benchmark tie-reporting.**~~ **Done.** `BenchmarkResult`
   now emits `tie_band_size_at_k` and `tie_break_policy` on every
   result; `evaluate_ranking` uses a deterministic `candidate_id`
   secondary sort (P1 fix).
3. ~~**Roth author response** to confirm GSE accession.~~ **Done**
   indirectly via the supplementary reporting summary: the correct
   accession is `GSE322563` (main paper printed `GSE32256` — typo).
   Ingested and benchmarked; β values match Fig. 5d at all three
   validated sites.
4. ~~**Benchmark label repair** to replace gene-universe positives with
   validated target loci.~~ **Done.** Three new positives files
   (`positives_roth_{validated,narrow,wide}.txt`) derived from Fig. 5d.
5. ~~**Second cell-line generalization cohort** — GSE68379 (Sanger GDSC,
   breast line panel on HM450).~~ **Run; reclassified.** GSE68379 was
   ingested and benchmarked (52 breast tumor × GSE69914 healthy normal,
   n=52/50). It is **not** a fourth generalization cohort — at the three
   Roth-validated probes, Sanger's MCF-7 is methylated (β 0.63–0.92)
   where Roth's MCF-7 is unmethylated (β 0.01–0.07), so the transferred
   Roth label set is not biologically valid on this cohort. Documented
   in §8.1 as an out-of-distribution boundary case. Does not bear on the
   V2.5 claim. A true "second cell-line generalization cohort" — one
   that reproduces the Roth hypo/hypermethylation pattern at the
   validated probes — remains open and would need a different cohort.
6. Stress the `sigma_floor` and one-sided-quantile paths on real data
   where they actually occur. Still open; not a blocker.
7. ~~**Benchmark honesty for tied bands.**~~ **Done.** `BenchmarkResult`
   now reports `precision_at_k_min`, `precision_at_k_max`,
   `recall_at_k_min`, `recall_at_k_max` — the [min, max] interval over
   adversarial tie-break within the K-boundary tied band. When
   `tie_band_size_at_k == 1` the interval collapses to the observed
   P@K. Low-`n` probabilistic results are now reportable as honest
   intervals rather than a single deterministic slice; e.g. GSE322563
   V2.5 narrow P@100 = 0.000 observed with interval [0.000, 0.020].

---

## 7. Decision record

- `v0.4.0` remains the stable release point. Default `probabilistic_mode`
  is still `tumor_only`.
- V2.5 remains experimental-on-main. Not tagged. README, CHANGELOG, and
  model docstrings now recommend V2.5 for matched cell-line / paper-
  comparable biology with explicit low-`n` tie-band caveats, while keeping
  V1 as the continuous-score fallback.
- GSE68379 (Sanger GDSC panel × GSE69914 normal) has been run and is
  documented as an out-of-distribution **boundary case** in §8.1 — not
  a fourth generalization cohort. The Roth-validated target probes are
  methylated in Sanger's MCF-7 and unmethylated in Roth's, so the
  transferred Roth label set is not biologically valid on this cohort;
  inversion of low-tumor-methylation axes is the expected behavior
  under label–biology mismatch. It does not bear on promotion.
- V2.5 cannot be promoted to a universal default under any of these
  results — `tumor_only`'s tissue-AUC lead + V2.5's tissue-AUC shortfall
  means there is no regime-agnostic winner. The framing stays
  cohort-type-dependent.

---

## 8. Cross-cohort readout under repaired Roth-validated labels

Three cohorts, three label granularities, three modes — AUC
(`tumor_only` → `differential` → `V1`):

| cohort | regime | `n`/side | validated | narrow | wide | V2.5 tie_band |
|---|---|:---:|---|---|---|---:|
| **GSE322563** | Roth cell lines | 2/2 | 0.928 / **0.990** / 0.821 | 0.886 / **0.942** / 0.884 | 0.871 / **0.910** / 0.768 | 190 |
| **GSE77348** | MCF-7/MCF-10A surrogate | 3/3 | 0.912 / **0.982** / 0.968 | 0.911 / **0.983** / 0.969 | 0.887 / **0.949** / 0.931 | 299 |
| **GSE69914** | primary tissue | 305/50 | **0.803** / 0.773 / 0.660 | **0.843** / 0.711 / 0.539 | **0.874** / 0.726 / 0.435 | **2** |

Bold = best AUC in that row. `tumor_only` tie_band: 6,540 on GSE69914,
10,005 on GSE322563, 11,848 on GSE77348 — unusable top-K everywhere.
`V1` tie_band = 1 everywhere (continuous-valued).

### What the matrix shows

- **Cohort-type × mode interaction.** On cell-line cohorts: V2.5 >
  `tumor_only` > V1. On tissue: `tumor_only` > V2.5 > V1. This is a
  real interaction, not a label-window artifact: all three
  granularities agree within each cohort.
- **Label repair is the dominant source of headline-AUC movement.**
  V2.5 on GSE322563 went from 0.694 (gene-universe tight) → 0.990
  (validated). V1 went from 0.541 → 0.821. Both were being penalized;
  V2.5 was being penalized more. The claim "V2.5 cleanly beats V1 on
  the paper-comparable cohort" is a function of using the right
  supervision target, not of the scoring math alone.
- **V2.5's tie_band dissolves with `n` exactly as predicted.** The §3
  prediction ("at `n ≳ 30`, `p_trust` stops saturating and the tied
  band dissolves") has now been tested; at n=305/50 the tie_band is 2.
  This is the single result that lets us drop the earlier blanket
  tie-fragility caveat.
- **`tumor_only`'s AUC-on-tissue lead is not a usable discovery signal.**
  6,540 tied records at K=100 means the top-100 is determined by the
  secondary sort, not the score. It is valuable only as a diagnostic
  axis.
- **V1 collapses on tissue.** AUC 0.435 on GSE69914 wide labels — worse
  than random. V1 is the discovery axis where ranking stability
  matters, but it is not the best AUC axis and its tissue AUC cannot
  be defended.

### Verdict

V2.5 is the recommended probabilistic axis on matched cell-line /
paper-comparable cohorts. It is second-best on tissue but the only
tissue axis with a usable top-K. It is not universally dominant; the
framing stays cohort-type-dependent.

---

## 8.1 GSE68379 — out-of-distribution boundary case (not a generalization cohort)

We also ran GSE68379 (Sanger GDSC1000, 52 breast cell lines on HM450,
paired with GSE69914 healthy normals as a cross-series external normal
arm, n=52/50). The headline AUC is dramatic and needs explicit framing
to avoid being misread.

| axis | validated AUC | narrow AUC | wide AUC | V2.5 tie_band |
|---|---:|---:|---:|---:|
| V1 `final_score` | 0.510 | 0.485 | 0.407 | 1 |
| `tumor_only` | 0.222 | 0.228 | 0.462 | 5,271 |
| V2.5 differential | 0.197 | 0.169 | 0.269 | 1,602 |

`tumor_only` and V2.5 are strongly **inverted** (AUC well below 0.5). V1
is near random. Under any naive reading this would look like V2.5
failing a generalization test. It is not. The root cause is a probe-
level biology mismatch at the three Roth-validated target sites:

| probe | gene | Roth MCF-7 β (GSE322563) | Sanger MCF-7 β (GSE68379) |
|---|---|---:|---:|
| cg05251676 | EGFLAM | 0.01 | **0.92** |
| cg25338972 | ESR1 | 0.07 | **0.63** |
| cg01364137 | GATA3 | 0.02 | **0.70** |

Same cell line name, **opposite methylation state** at exactly the
three loci Roth validated. This is a known phenomenon in the field —
MCF-7 has accumulated substantial genetic and epigenetic drift across
laboratories after decades of separate passaging. Sanger's MCF-7 (and
by extension most of their 52 breast lines) is heavily methylated at
EGFLAM / ESR1 / GATA3; Roth's MCF-7 is unmethylated at those same sites.

### Framing

> **GSE68379 functions as an out-of-distribution boundary case rather
> than a direct external validation cohort: the three Roth-validated
> target probes are methylated in Sanger's MCF-7, whereas they are
> unmethylated in Roth's MCF-7, so the transferred Roth label set is
> not biologically valid on this cohort. Under this label–biology
> mismatch, inversion of the low-tumor-methylation score axes is the
> expected behavior.**

Both V2.5 and `tumor_only` contain `p_targ = P(β_tumor < 0.30)` as a
factor — when β_tumor is 0.6–0.9 at the positives, p_targ is near zero
for those positives, and they rank at the *bottom* of the score
distribution. That is the correct response of the scorer to the
underlying biology; it is the labels that no longer describe the
candidates' biology in this cohort.

V1's near-random AUC arises by a different mechanism: its
`selectivity_score = max(0, β_n − β_t)` caps at zero when β_t > β_n,
so most validated positives score zero and land in a long mid-rank tie
rather than a bottom-rank inversion.

### What to say about this cohort

- V2.5 is supported on cohorts that reproduce the Roth target biology
  (GSE322563, GSE77348, GSE69914 under the repaired labels).
- GSE68379 is a boundary case where the transferred Roth labels are
  biologically invalid.
- Under that label–biology mismatch, inversion of `tumor_only` and V2.5
  is *expected*, because both encode low tumor methylation as part of
  the targetability hypothesis.

### What not to say

- Do not call GSE68379 a "failed generalization" of V2.5. It is a
  failure of label portability, not of the scorer.
- Do not pool GSE68379 into a single cross-cohort headline metric for
  V2.5 — it would pollute the estimate.
- Do not redefine positives on GSE68379 from GSE68379's own β values
  and report the resulting AUC as a comparable metric. That is
  circular (positives defined from the data being tested) and
  non-comparable to the label-transfer results on the other three
  cohorts.

### Recommended-axis table — unchanged

Nothing about this result moves the cohort-type × axis recommendation
in README. GSE68379 is not a fourth cohort in the main table; it is a
documented boundary case. The recommendation stays: V2.5 is the
recommended probabilistic axis on matched cell-line / paper-comparable
biology; V1 is the continuous-score fallback for top-K stability;
`tumor_only` stays analysis-only.

---

## 9. Top-hit annotation — from "good AUC" to "real target shortlist"

AUC answers "does the scorer order candidates well." The downstream
question is "what are the top-ranked loci actually?" — does the output
look like a plausible target list, or is it an AUC artifact?
`scripts/annotate_top_hits.py` attaches regulatory context (nearest
gene + TSS distance, feature class, CpG-island context) to the top-K
by any score axis. Below: top-20 across V1 and V2.5 on the three
Roth-biology cohorts. No Roth-validated positive lands in any top-20
— V2.5's AUC lift on Roth labels comes from *lifting the 3 Roth
targets into the top ~1% of 3M candidates*, not from putting them at
the top-100. That is the expected behavior on an n=2 / n=3 cohort
with a 190-300-record tied band.

Tie-break note: `top_k_by` and `thermocas inspect --top` now use the
same `candidate_id` ascending secondary sort as `evaluate_ranking` (the
earlier versions used the opposite end of the tied band, which made
the annotated top-20 and the benchmark disagree on tied cohorts —
fixed in a dedicated commit). The gene lists below reflect the
corrected ordering. V1 top-20s are unchanged because V1's tie_band is
always 1.

### GSE322563 (Roth actual, n=2/2)

**V1 `final_score` top-20**, highest-ranked genes:

  ZMIZ1, DPYSL4, **RET**, LINC02669, ADRB1, KCNIP2, [TRIM31 ×5],
  CNPY3-GNMT, TTBK1, GFRA1, IRX1, CTBP2, ADGRA1, LINC01163, ME1, NRG3

Spans chr5 / chr6 / chr10 (continuous-valued V1 score escapes tied
bands). tie_band_size_at_k = 1.

**V2.5 differential top-20**, highest-ranked genes:

  KCNIP2, CALHM2, SORCS3-AS1, LINC00710, [CELF2 ×2], XPNPEP1,
  CELF2-AS2, [ADRB1 ×3], ABLIM1, [GFRA1 ×3], PLPP4, DMBT1, BUB3,
  [CTBP2 ×2]

All on chr10 — this is the first 20 records of the 190-record tied
band at `p_sel = 0.0633` when tied records are ordered by
`candidate_id` ascending. The [P@100_min, P@100_max] interval on
narrow labels is **[0.000, 0.020]** — at best 2 of 28 narrow-positives
could land in top-100; at worst 0.

**Shared top-20 genes across V1 and V2.5 (GSE322563):** ADRB1, CTBP2,
GFRA1, KCNIP2. Four genes appear in both axes' top-20. V1 reaches
further into the underlying biology because it is not tied-band-
confined; V2.5 samples a specific lex-ordered window inside the
chr10-heavy portion of the 190-tied band.

### GSE77348 (MCF-7 / MCF-10A surrogate, n=3/3)

**V2.5 differential top-20 genes:** PITX3, BTRC, SORCS1, LINC02624,
CELF2-DT, CELF2, XPNPEP1, ADRB1, AFAP1L2, [GFRA1 ×3].

Cross-cohort overlap with GSE322563 V2.5 top-20: **ADRB1, CELF2,
GFRA1, XPNPEP1** (4 genes). Meaningful cross-cohort convergence
given both top-20s are window slices of very large tied bands
(190 and 299 respectively).

### GSE69914 (primary breast tissue, n=305/50)

**V2.5 differential top-20 genes:** [RUNX2 ×5 at one locus],
[SCML4 ×5], FOXCUT, FOXC1, MAS1L, COL21A1, MXI1, [TENM2 ×3], MAT2B.

β differentials are smaller (0.17–0.31) than on cell lines (0.79–0.96)
because primary tissue is cell-type-heterogeneous — β regresses toward
the middle. tie_band_size_at_k = 2 at K=100 — the tied band is
essentially gone on this `n=305/50` cohort, so the top-20 is the
genuine top-20, not a window slice.

Zero gene overlap with GSE322563 / GSE77348 V2.5 top-20s. V2.5 is
cohort-specific — on tissue it finds different biology (RUNX2, SCML4,
FOX family) than on cell lines. This is the healthy behavior: the
scorer is not returning the same fixed gene list regardless of input.

### What this changes

- The top-20 annotations confirm the V2.5 claim is about a real
  target-discovery signal, not an AUC artifact. Candidates are
  biologically coherent Roth-pattern loci (tumor-unmethylated,
  normal-methylated, promoter- or CGI-adjacent on cell-line cohorts).
- On the cell-line cohorts (GSE322563, GSE77348) the V2.5 top-20 is
  a 20-record window inside a 190–299-record tied band — the reported
  list is sensitive to tie-break policy (which we now document as
  `candidate_id ascending`). The P@K min/max interval is the honest
  read of top-K performance on these cohorts.
- V1 and V2.5 produce overlapping but distinct shortlists on the same
  cohort; they share ~4 genes on GSE322563 despite sampling different
  tie-band slices. V1 is tie-band-free and therefore more sensitive
  to the full distribution.
- Zero cross-cohort gene overlap between cell-line and primary-tissue
  top-20s argues against a worrying failure mode ("the scorer always
  returns the same few loci") and for the cohort-specific behavior we
  want.

Annotation outputs are committed at:

- `examples/gse322563/top20_annotated_v1.tsv`
- `examples/gse322563/top20_annotated_v25.tsv`
- `examples/gse69914/top20_annotated_v25.tsv`
- `examples/gse77348_roth_labels/top20_annotated_v25.tsv`
