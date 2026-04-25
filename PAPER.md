# Compositional probability-scale scoring and tie-band-aware benchmarking for methylome-guided ThermoCas9 target-site ranking

**Author.** Allison Huang, Columbia University. Contact: <allisonhmercer@gmail.com>.
**Date.** 2026-04-24.
**Code.** <https://github.com/AllisonH12/thermocas9> at tag `paper-5-10f` (immutable pointer to the exact revision that produced this paper).
**Status.** Educational research framework. Not peer-reviewed. No clinical claims. Cites Roth et al., *Nature* (2026), DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

---

## Abstract

Methylation-sensitive Cas9 target selection is a ranking problem:
candidates should be hypomethylated in disease cells, relatively
methylated or differentially protected in the matched normal
comparator, and supported by trustworthy methylation evidence. We present a compositional probability-scale scoring skeleton,
`p_targ × (gap factor) × p_trust`, paired with a tie-band-aware benchmark
contract (`precision_at_k_{min,max}`, `tie_band_size_at_k`, mid-rank
Mann-Whitney AUC). Two gap factors ship in this tag: V2.5-diff
(`tumor_plus_differential_protection`) and V2.5-sigmoid
(`tumor_plus_gap_sigmoid`).

The primary endpoint is GSE322563 validated-label AUC at the three Roth
Fig. 5d target probes, evaluated under both HM450-intersect and native
EPIC v2 ingest paths. **Because the validated set contains only three
Roth Fig. 5d targets, AUC is interpreted as a rank-lift summary rather
than an inferential discovery-performance estimate.** V2.5-diff and
V2.5-sigmoid place all three positives in the upper few percentiles of
the WG candidate universe (top ~0.06–4.6% of millions of candidates, axis- and
cohort-dependent: ESR1 in the top ~0.1%, GATA3 in the top ~1%,
EGFLAM up to ~4.6% on GSE77348 — see §5.1 / §5.2.2 per-positive
ranks); a 10⁶-draw random-triple null resolves the rank lift beyond
the earlier 10⁴ floor, while the paired V2.5-diff / V2.5-sigmoid versus Δβ-only comparison
remains descriptive at `n_pos = 3`. The
strongest ranking stress-test result is on GSE69914 tissue under
transported Roth labels: a probe-level limma-style moderated-t DMR
baseline is competitive on matched cell lines but falls to AUC 0.573
on tissue, while V2.5-sigmoid reaches AUC 0.862 on the same
candidate mapping. The tissue stress-test gain is per-positive
heterogeneous: GATA3 remains strong under feature-matched controls,
while ESR1 is matched-near-random under V2.5-sigmoid (§5.9). A frozen
whole-genome panel shows V2.5-sigmoid matches V2.5-diff's cell-line AUC
within 0.002, avoids V2.5-diff's whole-genome low-`n` top-K tied bands
(`tie_band@100 = 1` vs 421–1,493), and improves tissue AUC by +0.05 to
+0.08.

The recommended probabilistic prioritization axis for hypothesis
generation is V2.5-sigmoid, **selected by post-repair sensitivity and
whole-genome stress testing on the same benchmark family used here
rather than by independent prospective validation**. It preserves
matched-cell-line rank lift where Δβ-only is already strong, eliminates
V2.5-diff's low-`n` whole-genome tied-band pathology, and improves the
shipped-default tissue stress-test AUC; the tissue gain is per-positive
heterogeneous and partly reflects the shipped V2.5-diff σ_floor default
(§5.3.1), and prospective utility remains to be validated on newly
generated labels or wet-lab follow-up. V2.5-diff is retained for
backward-compatible AUC parity, V1 remains the package default for
deterministic top-K behavior, and `tumor_only` remains the default
`probabilistic_mode` enum value for compatibility. This is a scoring-method and
benchmarking contribution on public methylation data, not prospective
editing validation. Whole-genome catalogs are SHA256-frozen and all
benchmark artifacts are committed under `examples/`.

---

## 1 · Background

### 1.1 ThermoCas9 is methylation-sensitive at the PAM

Roth et al. (2026) characterized *Geobacillus thermodenitrificans* T12
Cas9 (ThermoCas9) biochemically, structurally, and in human
cell-line editing assays. The central result: cleavage efficiency is
governed by the methylation state of the fifth position of the PAM
(`5'-NNNN`**`C`**`GA-3'` or `5'-NNNN`**`C`**`CA-3'`), not by
methylation within the protospacer. Methylated 5-methylcytosine (5mC)
at the PAM abolishes binding; unmethylated PAMs are cleaved as usual.
Cryo-EM structures at 2.2–2.8 Å resolved the protein contacts that
discriminate methylated from unmethylated cytosine.

The practical consequence is that ThermoCas9 gives a natural
mechanism to selectively edit in cells where a target locus is
hypomethylated while sparing cells where the same locus is methylated.
Roth demonstrated this in MCF-7 vs. MCF-10A at three breast-cancer
loci: *EGFLAM* (control locus included in the validated target set),
*ESR1*, and *GATA3*. The underlying
methylation of the PAM cytosine at those loci was profiled on the
Illumina Infinium MethylationEPIC v2 array and reported in
Supplementary Figure 5d of the main paper alongside the editing
results.

### 1.2 Target-site selection is a ranking problem

Given a reference genome, the ThermoCas9-compatible PAM sites on any
chromosome number in the millions. To pick candidates for wet-lab
follow-up, one wants a ranked list where the top scorers are
*biologically coherent* — meaning the candidate's PAM cytosine is
unmethylated in the disease cohort and relatively methylated or
differentially protected in the matched normal — *and* where the
confidence of the underlying methylation
measurement is part of the score. This is a classical
ranking-with-covariates problem, and the natural scoring target is a
bounded, probability-like selectivity score rather than a calibrated
editing probability.

This paper describes the specific decomposition we ended up using
after two earlier attempts failed on empirical cohorts.

### 1.3 Notation

For one candidate site on one cohort we carry six per-probe β
summary statistics from the cohort's methylation array
(`β_tumor_mean`, `β_tumor_q25`, `β_tumor_q75`, `β_normal_mean`,
`β_normal_q25`, `β_normal_q75`), per-side sample counts (`n_tumor`,
`n_normal`), and an `EvidenceClass` capturing the distance between
the candidate's PAM cytosine and the nearest assayed CpG probe
(`EXACT`/`PROXIMAL_CLOSE`/`PROXIMAL`/`REGIONAL`/`UNOBSERVED`) — nine
inputs in total per candidate × cohort pair.

---

## 2 · From V2 to V2.5 (summary; full audit trail in Appendix A)

The first-pass composite was

```
p_therapeutic_selectivity = p_targ × p_prot × p_trust
```

with `p_prot = P(β_normal > 0.50)`. A factor ablation on the
GSE77348 MCF-7 vs MCF-10A surrogate showed `p_prot` is
**anti-predictive** under broader cohort conditions (AUC 0.384 alone;
multiplying it against `p_targ` collapsed the composite from 0.733 to
0.503). The cause is semantic: `P(β_normal > 0.5)` assumes the normal
comparator methylates the target, which holds for Roth's specific
MCF-10A cell line at the three *ESR1* / *GATA3* / *EGFLAM* probes
they chose to validate but does not generalize. The intermediate
V2.4 fix dropped `p_prot` to the `tumor_only` mode (`p_targ ×
p_trust`); that improved AUC but the top-100 collapsed onto
"always-unmethylated" loci (good AUC, unusable for prioritization), so a
deeper fix — replacing the threshold with a differential gap factor —
was needed. **Appendix A** carries the full V2 → V2.4 audit trail
(first-pass composite definition, factor-ablation table on
GSE77348, and the V2.4 default-mode change rationale).

---

## 3 · V2.5 · two gap-factor instances (V2.5-diff and V2.5-sigmoid)

**Terminology note.** "V2.5" names the *generation* of the
compositional skeleton `p_targ × (gap factor) × p_trust`. Two
gap-factor instances ship as selectable `probabilistic_mode` enum
values in this tag:

- **V2.5-diff** = `tumor_plus_differential_protection`, the
  original V2.5 composite with `p_diff` as the gap factor. This
  section defines it.
- **V2.5-sigmoid** = `tumor_plus_gap_sigmoid`, a fixed-bandwidth
  sigmoid replacement of `p_diff`. Introduced and benchmarked in
  §5.2.1 + §5.2.2; the recommended probabilistic prioritization axis
  per the §5.2.2 whole-genome panel.

Older artifacts may label V2.5-diff less precisely (the axis under
test in §5.1–§5.2 and the axis the scored JSONLs in
`data/derived/scored_*.jsonl` were computed with). In this paper, we
keep V2.5 / V2.5-diff / V2.5-sigmoid disambiguated whenever a specific
variant matters.

![V2.5 mode-formula schematic](docs/figures/fig1_mode_schematic.png)

**Figure 1.** The V2.5-diff composite. Three factors multiply:
`p_targ` (tumor unmethylated at the PAM cytosine), `p_diff`
(differential methylation gap exceeds δ under a normal approximation
on per-probe β summaries), and `p_trust` (evidence-class confidence,
saturating in min sample count). The result is the
*probability-scale selectivity score* (the stored scalar that
downstream ranking consumes; stored in the `p_therapeutic_selectivity`
field on every `ProbabilisticScore` record). The deprecated V2
`p_prot = P(β_normal > 0.5)` is
replaced by `p_diff` to remove its static-threshold assumption
about the normal side. V2.5-sigmoid swaps `p_diff` for a
fixed-bandwidth sigmoid of the same (β_n − β_t − δ) shift, motivated
by the §3.5 binding-rate finding and benchmarked in §5.2.1 / §5.2.2.

### 3.1 Formulation

Replace the threshold-based `p_prot` with a differential factor:

```
p_diff(δ) = P(β_normal − β_tumor > δ)
```

estimated under an independent-normal approximation:

```
σ_k   ≈ IQR_k / 1.349                 (k ∈ {tumor, normal}; floor at 0.05)
σ_Δ²  = max(σ_t, floor)² + max(σ_n, floor)²
z     = (δ − (μ_n − μ_t)) / σ_Δ
p_diff = 1 − Φ(z) = 0.5 (1 − erf(z / √2))
```

with `Φ` the standard normal CDF. The composite becomes:

```
p_therapeutic_selectivity = p_targ × p_diff × p_trust
```

This makes no absolute-threshold claim about the normal side. It asks
only what fraction of normal-vs-tumor β pairs sampled from the
per-side empirical distributions would exceed δ, under independent
normal approximations to those per-side distributions. (This is a
*population-overlap* statistic, not a confidence interval on a
cohort-level mean contrast — see §3.3.)

**Calibration scope.** `p_therapeutic_selectivity` is a *bounded
probability-scale ranking score in [0, 1]*, not a calibrated probability
of editing success at the candidate site. The composite is not validated
against any prospective editing readout (no expected calibration error
is reported, no reliability diagram is fit, no Platt/isotonic
post-calibration is applied). What the bounded scale buys is
multiplicative composability with downstream probabilistic inputs
(target-mutation models, gRNA off-target probabilities, delivery-
efficiency priors) — not a frequentist guarantee that "0.9" means "90%
of 0.9-scored sites edit successfully." Read the score as a
probability-scale ranking axis throughout. The multiplicative form is a
ranking heuristic, not an independence model: `p_targ`, the gap factor,
and `p_trust` are empirically correlated on real methylation catalogs
(e.g. EXACT-class records cluster at extreme β values which also drive
`p_targ` and `p_diff`), so the product should not be interpreted as a
calibrated joint probability.

### 3.2 Choice of δ

An offline δ sweep over {0.2, 0.3, 0.4, 0.5} on the GSE77348 surrogate
selected δ = 0.2 by joint optimization of AUC and `P@100` against the
*pre-repair* (gene-universe) positives that were in use at the time
of selection. Under that pre-repair tuning objective, larger δ values
monotonically reduced both metrics — the factor was too strict to
admit candidates with moderate but real differentials.

The post-repair sensitivity story is more nuanced and is reported in
§5.3.2: against the Roth-validated label set, δ = 0.3 (and on some
cohort × tier rows, δ = 0.4) marginally outperforms δ = 0.2 on the
matched cell-line cohorts; tissue (GSE69914) prefers smaller δ.
We do not retune δ post hoc against the post-repair labels, because
doing so would invalidate the dev-cohort separation discipline of
§4.0. δ = 0.2 is therefore retained as the shipped default; readers
wanting the full sensitivity surface should consult §5.3.2.

The default is exposed as `differential_delta` in the cohort
configuration. The δ used is recorded on every emitted
`ProbabilisticScore` record so the margin in effect is recoverable
from a single output row.

### 3.3 What `p_diff` is a probability over

`p_diff(δ) = P(β_normal − β_tumor > δ)` is a probability under the
*biological β-heterogeneity* implied by per-cohort, per-probe
quartiles — not a posterior probability of editing efficacy and not
a confidence interval on a cohort-level mean contrast. The
distribution being integrated is the one that would generate
single-sample β values from a tumor-side population with quartiles
(`β_tumor_q25`, `β_tumor_q75`) and a normal-side population with
quartiles (`β_normal_q25`, `β_normal_q75`), under independent normal
approximations on each side. Two consequences follow:

1. **The variance does not shrink with sample size.** Because the
   per-side σ is the IQR-derived population dispersion, not the
   standard error of a mean, `p_diff` does not become more confident
   as `n_tumor`/`n_normal` grow. The growing-`n` story for the
   composite is carried by `p_trust` (saturating at `n ≳ 30`), not
   by `p_diff` tightening.
2. **`p_diff` is therefore a population-overlap statistic, not a
   significance test.** It answers "what fraction of tumor-vs-normal
   β-pairs sampled from the cohort would exceed δ?" An alternative
   formulation that scales σ with `n` (e.g. SE on the mean
   difference, or a hierarchical model on per-side variance) would
   also be defensible and would give a different — and at low `n`,
   a much wider — distribution. We use the population-overlap form
   because (a) the comparison the user cares about for site
   selection is biological — "is this site differentially
   methylated *between cells*?" — not "does this cohort have enough
   power to call this site significant?"; and (b) the small-`n`
   power question is already absorbed by `p_trust`. Re-deriving the
   composite under an SE-on-mean variant is a documented future
   axis (§6.3).

This is a deliberate modeling choice and is the single most
common point of confusion when reading `p_diff` cold.

### 3.4 Why a normal approximation instead of a difference-of-Betas

A proper Bayesian treatment would model β-values as Beta distributions
and compute the exact difference distribution numerically. For the
IQR ranges seen on HM450 / EPIC v2 arrays (~0.05–0.20), the normal
approximation is close to the correct answer and has the advantage of
being a two-line stdlib closed form. A σ floor of 0.05 prevents σ_Δ
from collapsing to zero at boundary β-values (very common for CpG
islands methylated at 0 or 1). These tradeoffs are documented in the
code; a downstream user who wants the exact difference-of-Betas path
can swap it in.

The normal approximation is operationally adequate for the
*ranking* use case but is not a calibrated probability of editing
success. Two situations where the approximation is poor on its own
terms but does not corrupt the ranking:

- **Boundary-β pile-ups.** β-distributions at CpG islands often
  pile up near 0 or 1; a single-mode normal centered at the
  per-side mean overstates symmetry and understates skew. The σ
  floor blunts the worst case (σ_Δ → 0) but does not restore the
  correct shape.
- **Bimodality across samples.** When tumor or normal samples carry
  two β sub-populations (e.g. methylation-class heterogeneity in
  primary tissue), the IQR-derived σ averages the two modes; the
  per-side normal then sits between the modes rather than on
  either. `p_diff` in this regime under-states the true population
  overlap.

Neither failure mode is silent: the §5.5 annotation pipeline
flags candidates whose per-side β quartiles indicate strong
bimodality. A future axis based on the exact difference-of-Betas
or a two-component mixture would address both.

### 3.5 The low-`n` tied-band prediction

With the new factor isolated, we can predict before running any new
cohort. `p_trust` is `EvidenceClass.base × min(1, min(n_t, n_n) / 30)`:
piecewise-linear in min sample count up to `ramp_n = 30`, then a
constant ceiling per evidence class. It is **never** continuous-valued
in the real-line sense — at n ≥ 30 it saturates to a discrete-by-
EvidenceClass value (0.95 at EXACT, 0.75 at PROXIMAL_CLOSE, etc.).
What changes between the low-`n` and high-`n` regimes is which factor
*dominates* the composite ordering:

- **Low n (n ≪ 30).** `p_trust` scales every same-evidence-class
  candidate by the same uniform fraction `(n / 30)`. Per-side IQRs are
  wide (or zero, hitting the `σ_floor`), pushing both `p_targ` and
  `p_diff` toward saturation (≈1) for any candidate with clean tumor
  hypomethylation and a large normal/tumor gap. The composite is
  effectively `1 × 1 × p_trust` for that whole subset, producing a
  large tied band at the top of the score distribution.
- **High n (n ≥ 30).** `p_trust` is constant per EvidenceClass.
  Per-side IQRs are realistic (HM450 ranges of 0.05–0.20). The
  composite ordering is now driven by the continuous variation in
  `p_targ × p_diff`, with `p_trust` only multiplying by a per-
  evidence-class constant. The previously-tied subset spreads out by
  its underlying `p_targ × p_diff` values.

The falsifiable prediction: the K = 100 tied band on V2.5-diff should be
large at n=2/3 per side and shrink to near-1 at n ≳ 30 per side. §5.1
and §5.3 test this directly across cohorts ranging from n=2/2 to
n=305/50.

The K = 100 large-band regime is what motivates the "top tied
candidate class" framing in §6.1: on n=2/2 cohorts the visible
top-20 is a 20-record window inside a 190- to 421-record tied band,
not a stable rank ordering. Per-positive ranks (§5.1) and AUC
(§5.1, §5.3) are the stable claims at low `n`.

**Quartile convention on n=2 records.** Per-side β quartiles are
computed with Python's `statistics.quantiles(..., method="inclusive")`
(R-7 quantile definition) and hard-clamped to the observed sample
range at every probe, so q25 / q75 / mean are always within
[min(β_samples), max(β_samples)]. On n = 2 samples per side this is
deterministic (IQR = 0.5 × |β_max − β_min|) but weak as a
population-dispersion estimate: the "IQR" it produces is a function
of two observations, not a quantile of an underlying distribution.
Because the `σ_floor = 0.05` kicks in whenever IQR/1.349 < 0.05
(i.e. whenever the two observations differ by less than ~0.135),
σ_floor dominates σ_Δ on most n = 2 records.

**Quantified σ_floor binding rate per cohort** (reproducible via
`uv run scripts/sigma_floor_binding_rate.py`; denominator is
records with both β means present):

| cohort | σ_floor binds on either side | σ_floor binds on both sides |
|---|---:|---:|
| GSE322563 HM450 (n = 2/2)   | **100.0%** | **99.9%** |
| GSE322563 native (n = 2/2)  | **100.0%** | **99.9%** |
| GSE77348 (n = 3/3)          | **100.0%** | **99.5%** |
| GSE69914 (n = 305/50 tissue) | 65.1%     | 43.2% |

On matched cell-line cohorts the σ_floor is the binding constraint
on essentially every observed record (both sides). On tissue the
true IQRs are wide enough that the σ_floor binds on roughly two
thirds of records either-side, and on both sides for only ~43%.
This is the mechanistic origin of both the low-`n` tied band
(§5.3.1) and the §5.3.2 finding that tissue AUC prefers a larger
σ_floor: on cell-line cohorts, σ_floor is determining σ_Δ almost
everywhere, so changing its value re-shapes the whole composite.

---

## 4 · Benchmark methodology

### 4.0 Analysis plan: primary endpoint, sensitivity, boundary

The design separates tuned components from evaluation. It is not a
third-party preregistration, but the ordering is recoverable from git:
`differential_delta = 0.2` was fixed at commit `aece3de` on the
GSE77348 surrogate before GSE322563 was ingested.

| Analysis role | Cohort / labels | Interpretation |
|---|---|---|
| Tuned component | GSE77348, pre-repair gene-universe labels | δ-development cohort; later AUC is tuned-on supporting evidence. |
| Primary endpoint | GSE322563 validated Roth probes (`n_pos = 3`) | Independent paper-comparable rank-lift test; read with per-positive ranks. |
| Sensitivity | narrow ±50 bp / wide ±500 bp labels, native EPIC v2 ingest | Stability checks without retuning. |
| Tissue behavior | GSE69914, n = 305 / 50 | Tests high-`n` tissue behavior and tied-band dissolution, not Roth biology. |
| Boundary | GSE68379 × GSE69914 normal | Label-transportability stress test; never pooled with discovery results. |

### 4.1 Cohorts

Four cohorts, all public HM450 or EPIC v2 methylation arrays:

| cohort | regime | n tumor / normal | platform | role |
|---|---|:---:|:---:|---|
| **GSE322563** | Roth actual samples | 2 / 2 | EPIC v2 | paper-comparable biology |
| **GSE77348** | MCF-7/MCF-10A surrogate (DNMT3B paper) | 3 / 3 | HM450 | δ development cohort (tuned-on supporting evidence) |
| **GSE69914** | primary breast tissue | 305 / 50 | HM450 | high-n tissue validation |
| **GSE68379** | Sanger GDSC cell-line panel × GSE69914 healthy normal | 52 / 50 | HM450 | orthogonal (see §5.4) |

All cohort builders live in `scripts/build_*.py` and emit the same
`LocalSummaryBackend` format (`tumor_summary.tsv`, `normal_summary.tsv`,
`PROVENANCE.md`). GSE322563 uses the Roth EPIC v2 β matrix; its HM450-
intersect path strips beadchip suffixes from EPIC v2 probe IDs and
retains 80.7% of probes, while the native EPIC v2 path is reported as a
sensitivity. Per-probe β values at the three Roth target positions match
the published values to two decimals. GSE77348 is restricted to untreated
replicates, GSE69914 excludes adjacent-normal and BRCA1-carrier arms, and
GSE68379 is used only as a cross-series boundary case.

### 4.2 Positives (label repair)

Our first benchmark used "gene-universe" positives: HM450 probes in the
neighborhood of Roth's validated gene symbols (~970 tight, 1,687 loose,
on chr5/6/10). A diagnostic on GSE322563 revealed this label set is
noisy: only 23% of Roth-gene probes have β_n − β_t > 0.2 on the actual
Roth samples. Most are gene-universe members, not validated target loci.

We rebuilt the positives from Roth Supplementary Fig. 5d which names
three exact MCF-7 / MCF-10A validated target coordinates in hg38.
These were lifted to hg19 via Ensembl REST
`/map/human/GRCh38/.../GRCh37`:

| target | gene | hg38 | hg19 (lifted) | β shift (hg19 vs hg38) |
|---|---|---|---|---:|
| EGFLAM T11 | EGFLAM | chr5:38,258,842 | chr5:38,258,944 | +102 bp |
| ESR1 T17 | ESR1 | chr6:151,690,043 | chr6:152,011,178 | +321,135 bp |
| GATA3 T18 | GATA3 | chr10:8,045,425 | chr10:8,087,388 | +41,963 bp |

Three positives files result, at increasing granularity:

- `positives_roth_validated.txt` (n = 3): closest NNNNCGA candidate per Roth target.
- `positives_roth_narrow.txt` (n = 28): NNNNCGA candidates within ±50 bp of each Roth target.
- `positives_roth_wide.txt` (n = 142): NNNNCGA candidates within ±500 bp of each Roth target.

### 4.3 Score axes

Seven axes are benchmarked. For naming conventions see §7 Methods
(canonical axis names):

- **Δβ-only** (`naive_selectivity`) — `β_normal_mean − β_tumor_mean`. The literature-naive baseline: rank by the raw methylation gap with no uncertainty propagation or evidence weighting. Reported up front so the value of the additional factors is auditable against the simplest possible scoring axis.
- **Δβ_z** — `(β_n − β_t) / sqrt(σ_t² + σ_n²)` using the same IQR-derived σ as V2.5-diff's `p_diff`. An uncertainty-aware effect-size ranker built from the same per-side dispersion signal; confirms that V2.5-diff's advantage is not just an alias for σ-normalized Δβ.
- **V1** — `final_score = sequence × selectivity × confidence − heterogeneity_penalty − low_coverage_penalty`. Deterministic, continuous-valued, tie band always 1.
- **V2** (`tumor_only`) — `p_targ × p_trust`. Retained for audit; default in v0.4.0.
- **V2.5-diff** (`tumor_plus_differential_protection`) — `p_targ × p_diff × p_trust` with δ = 0.2 and σ_floor = 0.05. The original V2.5 composite; retained in this tag for backward compatibility and AUC parity on cell-line cohorts.
- **V2.5-sigmoid** (`tumor_plus_gap_sigmoid`) — `p_targ × sigmoid((β_n − β_t − δ) / σ_fixed) × p_trust` with δ = 0.2 and σ_fixed ≈ 0.0707. The **recommended probabilistic prioritization axis** on every non-boundary cohort shape tested (§5.2.2 whole-genome panel); ships as a first-class `probabilistic_mode` enum value in this tag.
- **limma-style moderated-t** — sample-level probe-level Smyth (2004) empirical-Bayes moderated `t`-statistic, candidate-mapped via the nearest probe assigned in `observation.probe_id`. Pure-Python reimplementation of the variance-prior math; canonical R `limma::lmFit + eBayes` parity is reported in §5.8. Included as the methods-journal-standard DMR comparator; §5.1 + §5.2.2 report AUC + tie-band alongside the other axes.

**limma-style moderated-t baseline.** Standard DMR tools rank CpG probes or
regions, not per-PAM candidates, so we compute a probe-level Smyth
(2004) empirical-Bayes moderated `t` statistic from sample-level β
matrices and map each candidate to its nearest probe via the same
`EvidenceClass` assignment. This is a limma-style statistic implemented
in pure Python; §5.8 verifies functional parity against canonical R
`limma::lmFit + eBayes` on the evaluated cohorts. It is reported in §5.1
and §5.2.2 because it is the methods-reviewer baseline.

### 4.4 Platform harmonization and catalog scope

Two facts constrain interpretation of the reported numbers and
belong in the methodology, not as a footnote:

- **Platform.** GSE322563 is EPIC v2; the other cohorts are HM450. §5.1
  uses the HM450-intersect path for cross-cohort comparability and §5.2
  reports a native EPIC v2 sensitivity. The validated-label V2.5-diff
  AUC changes by only 0.003 (0.990 → 0.986); V1 changes more, but the
  matched-cell-line ordering V2.5-diff > V1 and V2.5-diff > Δβ-only holds on every
  tier × path row.

- **Catalog scope.** Three frozen probe-window catalogs are used
  across §5, each with its own SHA256-committed provenance file:

  * **chr5 / chr6 / chr10 HM450 (2.98M candidates)** —
    `catalog_hg19_chr5_6_10.jsonl`. The primary benchmark catalog
    used by §5.1 (primary endpoint) and §5.2 (label granularity,
    tie-band diagnostics, factor ablation, σ_floor and δ sweeps,
    Δβ-only baseline). Chosen because the three chromosomes carry
    the Roth Fig. 5d validated targets (EGFLAM on chr5, ESR1 on
    chr6, GATA3 on chr10) and give a realistic test bed.
  * **Whole-genome HM450 (19.8M candidates, SHA256 `d20661c5…`)** —
    `catalog_hg19_wg.jsonl`. The genome-wide gating catalog used
    by §5.2.2 for GSE322563 HM450 / GSE77348 / GSE69914. Frozen
    provenance: `data/derived/catalog_hg19_wg.PROVENANCE.md`.
  * **Whole-genome native EPIC v2 (35.4M candidates, SHA256
    `39df8f0f…`)** — `catalog_hg19_wg_epic_v2.jsonl`. The
    GSE322563 native-path genome-wide catalog used by §5.2.2.
    Frozen provenance: `catalog_hg19_wg_epic_v2.PROVENANCE.md`.

  Every AUC / `P@K` / tie_band number in §5 names its catalog
  explicitly. §5.1 and §5.2 are chr5/6/10-scope; §5.2.2 is
  whole-genome-scope; §5.3 (tissue behavior) carries both the
  chr5/6/10 and WG rows; §5.9 adds within-chromosome feature-matched
  negative controls. Cross-chromosome and chromosome-class matched
  controls remain the untested denominator-loosening axes (§6.3).

### 4.5 Tie-band handling

Earlier iterations of the benchmark sorted candidates by score alone,
with input-order as the implicit tie-break. On cohorts with low n, the
score distribution's top can sit inside a tied band of 190–11,000
records, making the reported `P@K` a function of serialization order
rather than of the scorer.

The shipped benchmark contract:

1. Sort primary by score descending, secondary by `candidate_id`
   ascending (`evaluate_ranking`'s explicit `tie_break_policy =
   "candidate_id_asc"`).
2. Report `tie_band_size_at_k`: number of records tied at the K-th
   position's score. When > 1, top-K membership is partially
   determined by the secondary key.
3. Report `precision_at_k_min` and `precision_at_k_max` (plus
   recall counterparts) — adversarial tie-break bounds over the tied
   band. The worst-case pushes every tied-band positive out of top-K;
   the best-case pulls them in. When `tie_band_size_at_k == 1` both
   bounds collapse to the observed `precision_at_k`.

These three fields are on every emitted `BenchmarkResult` JSONL row.

---

## 5 · Results

Results are organized by the analysis roles frozen in §4.0: §5.1
reports the independent GSE322563 primary endpoint; §5.2 reports
sensitivity and whole-genome stress testing used to select
V2.5-sigmoid; §5.3 reports high-`n` tissue behavior; §5.4 reports the
cross-series boundary case; §5.5 reports tie-window-aware top-hit
annotation; §5.6 reports `p_trust` sensitivity; §5.7–§5.9 report
EvidenceClass, limma parity, and feature-matched denominator controls;
§5.10 reports the independent-biology System B transport-gated
diagnostic.

![V2.5-sigmoid vs V2.5-diff vs limma-style on WG validated AUC and tie_band@100](docs/figures/fig2_auc_bars.png)

**Figure 2 — Final-method summary.** Whole-genome validated-label
performance (n_pos = 3 Roth Fig. 5d target probes) for the recommended
probabilistic prioritization axis (V2.5-sigmoid) against the V2.5-diff
predecessor and the limma-style moderated-`t` DMR baseline, across the
four evaluated cohort paths (GSE322563 HM450 is the independent
primary endpoint per §4.0; GSE322563 native EPIC v2 is the same Roth
samples under a parallel ingest path; GSE77348 is the δ-development
cohort surfaced as supporting evidence; GSE69914 is the high-`n`
tissue cohort). **(a)** WG AUC: V2.5-sigmoid matches V2.5-diff on the
three matched cell-line cohorts within 0.001 (0.989 vs 0.988, 0.998 vs
0.998, 0.982 vs 0.981) and beats it on tissue (0.862 vs 0.778); the
limma-style baseline is competitive on cell lines (0.959 / 0.991 /
0.962) but **collapses to near-random on tissue** (0.573). **(b)** WG
`tie_band@100` on log scale: V2.5-sigmoid eliminates the **421–1,493
record tied bands** that make V2.5-diff's whole-genome top-K
unactionable on low-`n` cell-line cohorts, holding `tie_band@100 = 1`
on every matched cell-line cohort and `= 6` on tissue. The limma-style
baseline sits at intermediate tie-band values (39–115). The supplementary
historical 4-axis × 3-tier sensitivity figure
(`docs/figures/fig2_supp_historical_sensitivity.png`) is retained for
auditability of the V2.5-diff-era sensitivity story. Numbers committed
in `examples/genome_wide_panel.tsv` (V2.5-diff, V2.5-sigmoid) and
`examples/limma_cross_cohort_panel.md` (limma-style) at this tag.

### 5.1 Primary endpoint — matched cell-line AUC at validated Roth probes

**Primary endpoint (independent):** AUC at
`positives_roth_validated.txt` (n = 3 Roth Fig. 5d target probes,
hg38 → hg19-lifted) on **GSE322563** (Roth's own MCF-7 / MCF-10A
samples, EPIC v2). δ was not tuned on this cohort; the label set
was not derived from it.

**Supporting evidence (development cohort for δ):** AUC on
**GSE77348** (MCF-7 / MCF-10A on HM450, the DNMT3B-paper surrogate).
Same matched cell-line biology, different platform and laboratory.
Reported here because same-axis ordering on a second lab's
MCF-7/MCF-10A is informative, but with the caveat that δ was
selected on this cohort (§4.0) so its AUC is tuned-on, not
independent.

| axis | **GSE322563** (primary, n_pos=3, n_neg=2,979,994) | GSE77348 (δ-tuned, n_pos=3, n_neg=2,979,994) |
|---|---:|---:|
| Δβ-only (`β_n − β_t`) | 0.974 | 0.972 |
| Δβ_z (uncertainty-aware Δβ) | 0.940 | 0.956 |
| limma-style moderated-`t` (probe-level → candidate-mapped) | 0.959 | 0.962 |
| V1 `final_score` | 0.821 | 0.968 |
| V2 `tumor_only` | 0.928 | 0.912 |
| V2.5-diff (`tumor_plus_differential_protection`) | 0.990 | 0.982 |
| **V2.5-sigmoid (recommended)** | **0.989** | **0.982** |

**Reading the axis grid.** Δβ-only (raw mean gap) is already a
strong baseline — it lifts the validated positives to AUC 0.974 on
the primary cohort. The uncertainty-aware Δβ_z baseline
(`(β_n − β_t) / sqrt(σ_t² + σ_n²)`, using the *same* IQR-derived
σ as V2.5-diff's `p_diff` denominator; reproducible via
`uv run scripts/delta_beta_z_baseline.py`) lifts to AUC 0.940 —
slightly *below* raw Δβ-only because the σ-based denominator
penalizes positives with wide IQRs. V2.5-diff adds **+0.016** over
Δβ-only and **+0.050** over Δβ_z on the primary endpoint; on
GSE77348 the corresponding margins are +0.010 and +0.026.
V2.5-sigmoid (the recommended axis — see §5.2.2) is AUC-equivalent
to V2.5-diff on both cohorts within 0.001 (0.989 vs 0.990 on
GSE322563, 0.982 vs 0.982 on GSE77348), consistent with the §5.2.1
ablation finding that `p_diff`'s per-site structure and the
fixed-bandwidth sigmoid are interchangeable on cell-line cohorts
where σ_floor binds ~99% of records; the reason V2.5-sigmoid is
nonetheless recommended over V2.5-diff is top-K usability at WG
scale on cell-line cohorts, not AUC on chr5/6/10 — see §5.2.2
for the `tie_band@100` numbers. These margins are small but
consistent across the matched cell-line cohorts and label
granularities tested (§5.2). The `p_targ` (tumor-side
unmethylation) and `p_trust` (evidence / sample saturation) factors
of both V2.5 variants are pulling real weight on top of the per-side
dispersion signal — the V2.5 generation is not just a re-skinned Δβ_z — but the
honest framing is that the shipped baselines (especially raw
Δβ-only) already do most of the rank-lift work on matched
cell-line cohorts; the V2.5-generation contribution is the
tie-band-honest top-K reporting (§4.5) and the probability-scale
composition with future factors, not the AUC point alone.

The V2.5-diff margin over Δβ-only is also not uniform across the three
positives. The largest paired lift is at *GATA3*, the positive with the
smallest raw methylation gap, where V2.5-diff moves the rank from
159,932 to 24,083 on GSE322563 HM450 and from 61,254 to 28,325 on
GSE77348. Thus the composite helps most in the case where raw Δβ is
weakest; the matched-cell-line AUC gain remains descriptive at
`n_pos = 3`.

On the independent primary endpoint (GSE322563 validated), V2.5-diff is
the highest-AUC axis — +0.17 over V1 and +0.06 over the deprecated
V2 `tumor_only` composite. On the development cohort (GSE77348
validated), V2.5-diff is similarly the highest-AUC axis (+0.01 over V1,
+0.07 over `tumor_only`); this is consistent with the primary-
endpoint result but is not an independent replication, since the
same cohort selected δ = 0.2 pre-repair-labels (§4.0).

#### 5.1.1 Per-positive ranks of the n = 3 validated targets

AUC summarizes the rank lift across the three validated positives,
but with `n_pos = 3` it is highly sensitive to the rank of each one.
The table below gives the per-positive rank, percentile rank, score,
and per-side β values on each cohort × axis. (Full grid including
GSE69914, all four axes including Δβ-only, is at
`examples/validated_positive_ranks.md`; reproduce with
`uv run scripts/validated_positive_ranks.py`.)

**GSE322563 HM450 (primary, N=2,979,997 candidates):**

| gene | β_t | β_n | Δβ | V1 rank (%ile) | V2.5-diff rank (%ile) | Δβ-only rank (%ile) |
|---|---:|---:|---:|---:|---:|---:|
| *ESR1*   | 0.07 | 0.94 | 0.87 | **844** (99.97%) | **2,575** (99.91%) | 4,742 (99.84%) |
| *EGFLAM* | 0.02 | 0.49 | 0.48 | 839,106 (71.84%) | 64,433 (97.84%) | 67,620 (97.73%) |
| *GATA3*  | 0.03 | 0.31 | 0.28 | 755,885 (74.64%) | 24,083 (99.19%) | 159,932 (94.63%) |

**GSE77348 (development cohort for δ, N=2,979,997):**

| gene | β_t | β_n | Δβ | V1 rank (%ile) | V2.5-diff rank (%ile) | Δβ-only rank (%ile) |
|---|---:|---:|---:|---:|---:|---:|
| *ESR1*   | 0.03 | 0.87 | 0.84 | **1,618** (99.95%) | **3,965** (99.87%) | 8,733 (99.71%) |
| *EGFLAM* | 0.005 | 0.34 | 0.33 | 257,041 (91.37%) | 124,925 (95.81%) | 178,663 (94.01%) |
| *GATA3*  | 0.015 | 0.59 | 0.58 | 31,061 (98.96%) | **28,325** (99.05%) | 61,254 (97.94%) |

**Reading the table.**
*ESR1* drives the AUC story: it is in the top ~0.1% of the score
distribution under V2.5-diff on both matched cell-line cohorts (99.91%
on GSE322563 HM450, 99.87% on GSE77348), with V1 also placing it well. *EGFLAM* and *GATA3* sit in the top 1–3%
under V2.5-diff — meaningful rank lift over V1, but **outside the top
100 on every axis tested**. None of the three validated positives
appears in any cohort's top-20 (§5.5). This is the operational
meaning of "rank lift on n=3": V2.5-diff reliably places known
positives in the upper few percentiles of millions of candidates
(ESR1 top ~0.1%, GATA3 top ~1%, EGFLAM top ~2–4.6% depending on
cohort), but on n=2/2 cohorts
the top tied band at K=100 is large enough that the validated
positives sit *behind* a 190- to 421-record block of equally-scored
candidates.

The ESR1-dominance also implies AUC is fragile: if ESR1 alone
moved out of the top 10,000 on GSE322563, the AUC would drop
from 0.990 to roughly the average of the other two positives'
percentiles ≈ 0.985 (still high), but the AUC *gap* between V2.5-diff
and Δβ-only would also compress.

#### 5.1.2 Descriptive AUC uncertainty

At `n_pos = 3`, inferential AUC intervals are coarse and dominated by
the rank of each positive. We therefore report a 1,000,000-draw
random-triple null: draw three random candidates without replacement,
compute AUC under the same mid-rank tie handling, and estimate
`Pr(AUC_random_triple ≥ AUC_observed)`. The empirical *p*-value uses
the standard +1 smoothing `(n_ge + 1) / (n_perm + 1)`, so the resolution
floor is ~10⁻⁶; values in the first decade above that floor should be
read as finite-Monte-Carlo tail-count estimates, not precise
probabilities.

| cohort | axis | observed AUC | random-triple null 2.5–97.5% | one-sided *p* (Pr ≥ obs) |
|---|---|---:|---|---:|
| GSE322563 HM450      | V1 final_score | 0.821 | [0.166, 0.817] | 2.2 × 10⁻² |
| GSE322563 HM450      | Δβ-only        | 0.974 | [0.181, 0.823] | 8.6 × 10⁻⁵ |
| **GSE322563 HM450**  | **V2.5 diff**  | **0.990** | [0.338, 0.802] | **<1 × 10⁻⁵** |
| GSE322563 native v2  | V1 final_score | 0.933 | [0.206, 0.823] | 1.3 × 10⁻³ |
| GSE322563 native v2  | Δβ-only        | 0.961 | [0.177, 0.823] | 3.1 × 10⁻⁴ |
| **GSE322563 native v2** | **V2.5 diff** | **0.986** | [0.313, 0.817] | **1.5 × 10⁻⁵** |
| GSE77348             | V1 final_score | 0.968 | [0.261, 0.823] | 1.6 × 10⁻⁴ |
| GSE77348             | Δβ-only        | 0.972 | [0.177, 0.823] | 1.1 × 10⁻⁴ |
| **GSE77348**         | **V2.5 diff**  | **0.982** | [0.176, 0.823] | **2.8 × 10⁻⁵** |

V2.5-diff is the strictest listed axis on every row and sits below the
old 10⁴-draw resolution floor; the exact Monte Carlo tail-count
estimates are preserved in the artifact. This is an against-random
check, not a superiority test between axes. The direct paired comparison with
Δβ-only is descriptive (§5.1.3). Full output:
`examples/auc_uncertainty_1e6.md`; the negative bootstrap is omitted
from the body because holding the positives fixed makes the spread
mechanically collapse to ±0.001.

#### 5.1.3 Paired per-positive rank lift: V2.5-diff vs Δβ-only

The paired comparison with Δβ-only on the same three positives
(`uv run scripts/paired_rank_lift_v25_vs_dbeta.py`):

| cohort | *ESR1* Δβ-only rank / V2.5-diff rank | *EGFLAM* Δβ / V2.5-diff | *GATA3* Δβ / V2.5-diff | positives where V2.5-diff &gt; Δβ-only |
|---|---:|---:|---:|:---:|
| GSE322563 HM450      | 4,742 / 2,575   | 67,620 / 64,433   | 159,932 / 24,083  | **3 / 3** |
| GSE322563 native v2  | 9,085 / 5,746   | 185,641 / 158,502 | 422,569 / 50,096  | **3 / 3** |
| GSE77348             | 8,733 / 3,965   | 178,663 / 124,925 | 61,254 / 28,325   | **3 / 3** |

V2.5-diff ranks all three positives above Δβ-only on every matched
cell-line row. The largest lift concentrates at *GATA3*, the positive
with the smallest raw β gap, consistent with the composite helping most
where naive Δβ is weakest. The sign-flip paired null has only 2³
patterns and floors at one-sided *p* = 0.125 when all three signs agree,
so this is descriptive 3-for-3 directional evidence, not an inferential
superiority claim.

### 5.2 Sensitivity analyses: label granularity and P@K intervals

**Label granularity (AUC).** Stability of the primary-endpoint
ordering under weaker label definitions (narrow ±50 bp, wide
±500 bp):

| cohort | label set | n_pos | n_neg | Δβ-only | V1 | tumor_only | V2.5-diff | V2.5-diff − Δβ-only |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| GSE322563 | validated | 3 | 2,979,994 | 0.974 | 0.821 | 0.928 | **0.990** | +0.016 |
| GSE322563 | narrow | 28 | 2,979,969 | 0.912 | 0.884 | 0.886 | **0.942** | +0.030 |
| GSE322563 | wide | 142 | 2,979,855 | 0.844 | 0.768 | 0.871 | **0.910** | +0.066 |
| GSE77348 | validated | 3 | 2,979,994 | 0.972 | 0.968 | 0.912 | **0.982** | +0.010 |
| GSE77348 | narrow | 28 | 2,979,969 | 0.964 | 0.969 | 0.911 | **0.983** | +0.019 |
| GSE77348 | wide | 142 | 2,979,855 | 0.910 | 0.931 | 0.887 | **0.949** | +0.039 |

V2.5-diff remains the highest-AUC axis on both cohorts at every label
granularity. No retuning of `δ` between rows.

**P@K intervals.** On `n = 2 / 2` and `n = 3 / 3` cohorts the score
distribution at K = 100 sits inside a tied band (§3.5); P@100 is
reported as an adversarial interval `[min, max]` per §4.5. On
GSE322563 narrow labels:

| axis | P@100 observed | P@100 min | P@100 max | tie_band@100 |
|---|---:|---:|---:|---:|
| V2.5-diff | 0.000 | 0.000 | 0.020 | 190 |
| tumor_only | 0.000 | 0.000 | 0.020 | 10,005 |
| V1 final_score | 0.020 | 0.020 | 0.020 | 1 |

V2.5-diff's interval is [0.000, 0.020] — under the deterministic
`candidate_id` ascending tie-break, zero narrow-positives land in
top-100; under the best-case adversarial tie-break inside the
190-tied band, 2 of 28 could. V1's interval collapses to the
observed value because V1's score is continuous and `tie_band = 1`.
`tumor_only`'s tied band (10,005) spans two orders of magnitude
more of the score distribution than V2.5-diff's; top-K on `tumor_only`
is not a meaningful quantity on this cohort.

**Native EPIC v2 vs HM450-intersect on GSE322563.** §4.4 describes
the harmonization shortcut: GSE322563 EPIC v2 probe IDs are stripped
of beadchip-design suffixes and intersected with the HM450 universe
(80.7% retention). To verify the shortcut is not distorting the
primary endpoint, we ran the full pipeline a second time against a
native EPIC v2 catalog (5.22M candidates vs the HM450 catalog's
2.98M):

| label set | n_pos | n_neg (HM450 / native) | axis | HM450-intersect AUC | native EPIC v2 AUC | Δ |
|---|---:|---:|---|---:|---:|---:|
| validated | 3 | 2,979,994 / 5,224,079 | Δβ-only | 0.974 | 0.961 | −0.013 |
| validated | 3 | 2,979,994 / 5,224,079 | V1 | 0.821 | 0.933 | +0.112 |
| validated | 3 | 2,979,994 / 5,224,079 | tumor_only | 0.928 | 0.936 | +0.008 |
| validated | 3 | 2,979,994 / 5,224,079 | **V2.5-diff** | **0.990** | **0.986** | **−0.003** |
| narrow | 28 | 2,979,969 / 5,224,054 | Δβ-only | 0.912 | 0.956 | +0.044 |
| narrow | 28 | 2,979,969 / 5,224,054 | V1 | 0.884 | 0.938 | +0.054 |
| narrow | 28 | 2,979,969 / 5,224,054 | tumor_only | 0.886 | 0.933 | +0.046 |
| narrow | 28 | 2,979,969 / 5,224,054 | V2.5-diff | 0.942 | 0.983 | +0.040 |
| wide | 142 | 2,979,855 / 5,223,940 | Δβ-only | 0.844 | 0.865 | +0.021 |
| wide | 142 | 2,979,855 / 5,223,940 | V1 | 0.768 | 0.855 | +0.087 |
| wide | 142 | 2,979,855 / 5,223,940 | tumor_only | 0.871 | 0.916 | +0.044 |
| wide | 142 | 2,979,855 / 5,223,940 | V2.5-diff | 0.910 | 0.945 | +0.035 |

The V2.5-diff primary-endpoint AUC moves by 0.003 (well within the
tied-band noise floor on this n=2/2 cohort). V1 gains the most under
native EPIC v2 — its continuous score ranks the validated targets
higher relative to easy-negative candidates added by the larger
catalog. The relative axis ordering (V2.5-diff > V1) is preserved across
both ingest paths. The HM450-intersect shortcut is therefore not
materially distorting the headline V2.5-diff claim. Tied bands grow
modestly with the larger catalog (V2.5-diff: 190 → 421;
tumor_only: 10,005 → 14,914; V1 stays at 1).

#### 5.2.1 Factor ablation: does `p_diff` do the work, or does the floor?

The §3.5 binding-rate diagnostic shows that `σ_floor = 0.05` is the
binding constraint on σ_Δ for essentially every observed record on
matched cell-line cohorts (99.5–99.9% both-sides). That raises an
obvious reviewer question: if σ_floor is determining σ_Δ almost
everywhere on cell lines, is `p_diff`'s per-site σ-adaptive structure
actually adding anything, or is the composite effectively
`p_targ × sigmoid((Δβ − δ) / σ_fixed) × p_trust` with a constant
bandwidth?

We replace `p_diff` with two simpler gap factors, keeping `p_targ`,
`p_trust`, and δ = 0.2 unchanged (reproducible via
`uv run scripts/factor_ablation_vs_sigmoid.py`):

- **sigmoid**: `sigmoid((Δβ − δ) / σ_fixed)` with
  σ_fixed = √2 × σ_floor ≈ 0.0707 — the σ_Δ V2.5 sees when the floor
  binds on both sides (the modal case, §3.5).
- **hard threshold**: `1[Δβ > δ]` — tests whether smoothness matters
  at all.

AUC at the n = 3 Roth-validated positives:

| cohort | V2.5-diff (full `p_diff`) | sigmoid ablation | hard-threshold ablation |
|---|---:|---:|---:|
| GSE322563 HM450      | 0.990 | 0.989 | 0.990 |
| GSE322563 native v2  | 0.986 | 0.986 | 0.988 |
| GSE77348             | 0.982 | 0.982 | 0.986 |
| **GSE69914 (tissue)**    | **0.773** | **0.864** | 0.497 |

Two findings follow. On matched cell lines, `p_diff`'s per-site
σ-adaptive structure is not doing the AUC work: because σ_floor binds on
essentially every record (§3.5), V2.5-diff is empirically close to a
`p_targ × gap_step_or_smooth(Δβ, δ) × p_trust` composite, and all three
gap factors land within 0.001 AUC. The cell-line lift over Δβ_z should
therefore be attributed mainly to the tumor-side and evidence/trust
factors, not to per-site σ adaptation inside `p_diff`. On tissue, the
fixed-bandwidth sigmoid beats V2.5-diff by +0.09 AUC while the hard
threshold collapses to chance. Smoothness matters on tissue; per-site σ
adaptation does not help this label set.

**Bandwidth robustness of the tissue finding** (§5.2.1 single-point
ablation used σ_fixed = √2 · σ_floor ≈ 0.0707; reproducible via
`uv run scripts/sigmoid_bandwidth_sweep.py`):

| cohort | σ_fixed=0.05 | σ_fixed≈0.0707 | σ_fixed=0.10 | σ_fixed=0.15 | V2.5-diff |
|---|---:|---:|---:|---:|---:|
| GSE322563 HM450      | 0.990 | 0.989 | 0.989 | 0.977 | 0.990 |
| GSE322563 native v2  | 0.986 | 0.986 | 0.985 | 0.976 | 0.986 |
| GSE77348             | 0.982 | 0.982 | 0.979 | 0.965 | 0.982 |
| **GSE69914 (tissue)**    | **0.862** | **0.864** | **0.830** | **0.821** | **0.773** |

Every tested σ_fixed beats V2.5-diff on tissue (+0.05 to +0.09 AUC);
on matched cell lines the sigmoid family stays within 0.02 AUC except
at the widest bandwidth. The tissue result is therefore not a one-point
artifact and motivates V2.5-sigmoid as the recommended gap factor.

#### 5.2.2 Genome-wide panel: V2.5-diff vs V2.5-sigmoid across evaluated cohort paths

The §5.2.1 ablation and bandwidth sweep are on the chr5/6/10 catalog
(N = 2.98M candidates on HM450, 5.22M on native EPIC v2). A
reviewer-salient concern is whether the findings survive when the
negative universe expands to whole-genome. To resolve this we built
two frozen whole-genome probe-window catalogs and ran the same
gap-factor ablation on GSE322563 (both HM450 and native EPIC v2
paths), GSE77348, and GSE69914 against them:

- **HM450 WG catalog** — `catalog_hg19_wg.jsonl`, 19,787,820
  candidates, SHA256 `d20661c5…` (provenance at
  `data/derived/catalog_hg19_wg.PROVENANCE.md`).
- **Native EPIC v2 WG catalog** — `catalog_hg19_wg_epic_v2.jsonl`,
  35,380,431 candidates, SHA256 `39df8f0f…` (provenance at
  `data/derived/catalog_hg19_wg_epic_v2.PROVENANCE.md`). Built from a
  whole-genome re-run of `scripts/build_epic_v2_probes.py` (929,801
  GPL33022 probes, hg38 → hg19).

Reproducible via `uv run scripts/genome_wide_tissue_gating.py` per
cohort + `uv run scripts/genome_wide_panel.py` to aggregate. The
committed panel is at `examples/genome_wide_panel.md`.

**Cross-cohort AUC table at validated labels (n_pos = 3):**

| cohort | N_wg | V2.5-diff WG | V2.5-sigmoid σ=0.05 WG | V2.5-sigmoid σ=0.0707 WG | V2.5-sigmoid σ=0.10 WG |
|---|---:|---:|---:|---:|---:|
| GSE322563 HM450       | 19.8M | 0.989 | 0.989 | 0.988 | 0.988 |
| GSE322563 native v2   | 35.4M | 0.998 | 0.998 | 0.998 | 0.998 |
| GSE77348              | 19.8M | 0.982 | 0.982 | 0.981 | 0.979 |
| **GSE69914 (tissue)** | 19.8M | **0.778** | **0.861** | **0.862** | 0.825 |

**Cross-cohort `tie_band@100` table at validated labels (WG only):**

| cohort | V2.5-diff WG | V2.5-sigmoid σ=0.05 | V2.5-sigmoid σ=0.0707 | V2.5-sigmoid σ=0.10 |
|---|---:|---:|---:|---:|
| GSE322563 HM450       | **1,127** | 1 | 1 | 1 |
| GSE322563 native v2   | **421**   | 1 | 1 | 1 |
| GSE77348              | **1,493** | 1 | 1 | 1 |
| GSE69914 (tissue)     | 1         | 4 | 6 | 6 |

The WG panel resolves three questions. First, matched-cell-line AUC is
denominator-robust: all axes remain within 0.002 of the chr5/6/10 rows.
Second, V2.5-diff's low-`n` top-K degeneracy appears only at scale:
`tie_band@100` expands to 421–1,493 records, while V2.5-sigmoid stays at
1 on all matched cell-line WG cohorts. Third, the tissue V2.5-sigmoid
advantage survives the WG denominator (+0.05 to +0.08 AUC) with
single-digit tie bands. P@100 remains zero for the three Roth positives
at WG scale; the stable claim is rank lift, not top-100 recovery.

**Does V2.5-sigmoid on tissue need δ = 0.1 (per §5.3.2's `p_diff`
finding) or δ = 0.2 (the shipped default)?** §5.3.2 showed `p_diff`
on tissue prefers δ = 0.1 over δ = 0.2 by +0.05 AUC. We re-ran the
WG tissue gating at δ = 0.1 to check whether that property transfers
to V2.5-sigmoid. Result (`examples/gse69914_wg_gating_delta_0_1.md`):

| σ_fixed | V2.5-sigmoid AUC, δ = 0.2 | V2.5-sigmoid AUC, δ = 0.1 | Δ (δ=0.1 − δ=0.2) |
|---|---:|---:|---:|
| 0.05    | 0.861 | 0.865 | +0.004 |
| 0.0707  | 0.862 | 0.831 | **−0.031** |
| 0.10    | 0.825 | 0.821 | −0.004 |

The `p_diff` δ = 0.1 tissue advantage does not transfer cleanly to
V2.5-sigmoid; at the default σ_fixed ≈ 0.0707, δ = 0.1 drops by 0.031
AUC. δ = 0.2 is therefore retained as the cross-mode default.

As an additional post-selection robustness check, we also swept
V2.5-sigmoid directly over δ ∈ {0.10, 0.15, 0.20, 0.25, 0.30} and
σ_fixed ∈ {0.03,
0.05, 0.0707, 0.10, 0.15, 0.20} on the frozen whole-genome denominators
(`examples/sigmoid_delta_sigma_wg_sweep.{tsv,md}`; reproduce with
`uv run python scripts/sigmoid_delta_sigma_wg_sweep.py`). **This grid is
a descriptive post-selection robustness surface, not an additional
model-selection step; no `(δ, σ_fixed)` cell from the grid is used to
revise the recommended default.**

| cohort | default AUC (tie@100) | best AUC in grid | full-grid AUC range | max tie@100 |
|---|---:|---:|---:|---:|
| GSE322563 HM450 | 0.988 (1) | 0.990 | 0.967-0.990 | 1 |
| GSE322563 native v2 | 0.998 (1) | 0.998 | 0.995-0.998 | 1 |
| GSE77348 | 0.981 (1) | 0.983 | 0.957-0.983 | 1 |
| GSE69914 tissue | 0.862 (6) | 0.874 | 0.811-0.874 | 11 |

The default is not claimed to be the global optimum. The robustness
result is narrower and more useful: the default sits in a broad
high-AUC region, the low-`n` whole-genome tied-band problem remains
eliminated across the entire matched-cell-line grid, and every tested
tissue grid cell remains above V2.5-diff's 0.778 WG AUC.

Together, these validated-label AUC and tie-band tables support the
current recommendation: V2.5-sigmoid is AUC-equivalent to V2.5-diff on
matched cell lines, has higher transported-label AUC on the single
tissue cohort tested, and has far smaller WG top-K tied bands.

**Cross-cohort comparison against the limma-style moderated-t DMR baseline**
(§4.3 exclusion-rationale paragraph; cross-cohort panel at
`examples/limma_cross_cohort_panel.md`):

| cohort | limma-style moderated-t AUC | V2.5-diff AUC | V2.5-sigmoid AUC | limma tie@100 | V2.5-sigmoid tie@100 |
|---|---:|---:|---:|---:|---:|
| GSE322563 HM450      | 0.959 | 0.989 | 0.988 | 91 | 1 |
| GSE322563 native v2  | 0.991 | 0.998 | 0.998 | 39 | 1 |
| GSE77348             | 0.962 | 0.982 | 0.981 | 115 | 1 |
| **GSE69914 (tissue)**    | **0.573** | 0.778 | **0.862** | 57 | 6 |

The low-`n` matched-cell-line limma rows are useful comparator checks but
not the main between-method claim: with n = 2/2 or 3/3,
empirical-Bayes variance borrowing is structurally stressed and the
primary endpoint is still only three positives. The substantive
limma-style comparison is the high-`n` tissue row, where the
moderated-`t` baseline falls to 0.573 and V2.5-sigmoid reaches 0.862.
That contrast is the clearest evidence that the compositional wrapper
(`p_targ × p_trust`) is doing work a pure probe-level differential
statistic cannot.
Reproducible artifacts: `examples/limma_cross_cohort_panel.md` and
`examples/genome_wide_panel.md`.

### 5.3 Tissue-cohort behavior — GSE69914 (high-`n`, tissue biology)

GSE69914 (n = 305 / 50 primary breast vs. healthy donor tissue)
tests whether the low-`n` tied-band behavior predicted in §3.5
dissolves at `n ≫ 30` (n_neg is the chr5/6/10 negative universe minus
`n_pos`, ≈ 2.98M throughout):

| label set | n_pos | Δβ-only | Δβ_z | V1 | tumor_only | V2.5-diff | tie@100 |
|---|---:|---:|---:|---:|---:|---:|---:|
| validated |   3 | 0.591 | 0.504 | 0.660 | **0.803** | 0.773 | **2** |
| narrow    |  28 | 0.477 |   —   | 0.539 | **0.843** | 0.711 |   2   |
| wide      | 142 | 0.435 |   —   | 0.435 | **0.874** | 0.726 |   2   |

On tissue, both Δβ-based baselines collapse to near-random or below
(Δβ-only 0.435–0.591; Δβ_z 0.504 at validated). V1 tracks Δβ-only
closely (0.435–0.660). The V2.5-generation tissue story is therefore
a recovery *over both Δβ baselines*: V2.5-diff reaches +0.18 over
Δβ-only at validated and +0.27 over Δβ_z, not a modest margin as on
cell lines. `tumor_only` has a higher raw AUC (0.803–0.874) but its
K = 100 tied band makes its top-K unusable. **V2.5-sigmoid** — which
also ships in this tag as `tumor_plus_gap_sigmoid` and is the
recommended probabilistic prioritization axis (§6.1) — **improves on V2.5-diff's
tissue AUC by a further +0.05 to +0.09 across the bandwidth family
tested in §5.2.1 (0.821–0.864 vs V2.5-diff's 0.773)**, at both the
chr5/6/10 and whole-genome probe-window denominators (§5.2.2). The
table below reports V2.5-diff's cross-tier tissue AUC; V2.5-sigmoid
is the cross-mode improvement that the §5.2 follow-ups surface
directly.

Two findings:

1. **The tied-band prediction holds.** V2.5-diff's tied band at K = 100
   shrinks from 190 at `n = 2 / 2` → 299 at `n = 3 / 3` → 2 at
   `n = 305 / 50`, with no hyperparameter changes between cohorts.
   At `n ≥ 30`, `p_trust` saturates at the EvidenceClass ceiling
   (0.95 for EXACT) — a per-class *constant*, not a continuous
   function of n. The tied band shrinks not because `p_trust`
   becomes continuous-valued but because `p_targ × p_diff` are no
   longer pushed to 1 by wide-IQR saturation and the σ floor; with
   realistic-width IQRs (0.05–0.20) the composite's ordering is
   driven by continuous variation in those two factors while
   `p_trust` just multiplies by the class constant.

2. **The mode ordering flips on tissue biology.** `tumor_only`
   wins raw AUC (+0.03 to +0.15 over V2.5-diff / V2.5-sigmoid),
   but `tumor_only`'s tied band is 6,540 at K = 100 — its top-K is
   not usable. Both V2.5 variants have interpretable top-K on
   tissue (V2.5-diff `tie_band@100 = 2` chr5/6/10 / = 1 WG;
   V2.5-sigmoid `tie_band@100 ≤ 6` across the tested bandwidth
   family), and V2.5-sigmoid additionally wins on tissue AUC by
   +0.05 to +0.08 (§5.2.1 / §5.2.2) — so V2.5-sigmoid is the
   recommended tissue prioritization axis, with V2.5-diff retained for
   AUC-parity use. V1 collapses on tissue (AUC 0.44 on wide,
   below random).

The per-positive WG ranks make the tissue heterogeneity explicit
(`examples/tissue_per_positive_wg_ranks.{tsv,md}`):

| positive | Δβ-only WG %ile | limma WG %ile | V2.5-diff WG %ile | V2.5-sigmoid WG %ile | feature-matched p |
|---|---:|---:|---:|---:|---:|
| *ESR1* | 94.33% | 89.82% | 97.40% | 88.88% | 0.4384 |
| *EGFLAM* | 22.03% | 27.74% | 53.59% | 75.01% | 0.2588 |
| *GATA3* | 64.78% | 54.32% | 82.36% | 94.81% | 0.0187 |

Thus the tissue AUC gain is concentrated: V2.5-sigmoid strongly lifts
*GATA3* and moderately lifts *EGFLAM*, while *ESR1* reverses relative to
V2.5-diff and Δβ-only. The feature-matched audit in §5.9 confirms the
same pattern against local matched denominators.

### 5.3.1 σ_floor sensitivity

The 0.05 floor on per-side σ is the single most influential
implementation constant in the V2.5 composite (§3.1). A sweep over
σ_floor ∈ {0.02, 0.05, 0.10, 0.15} on the validated label set
(reproducible via `uv run scripts/sigma_floor_sweep.py`):

| cohort | σ_floor=0.02 AUC (tie@100) | σ_floor=0.05 *(default)* AUC (tie@100) | σ_floor=0.10 AUC (tie@100) | σ_floor=0.15 AUC (tie@100) |
|---|---|---|---|---|
| GSE322563 HM450      | 0.992 (719)   | **0.990 (190)**   | 0.989 (1)     | 0.985 (1) |
| GSE322563 native v2  | 0.990 (1,605) | **0.986 (421)**   | 0.985 (1)     | 0.982 (1) |
| GSE77348             | 0.983 (1,019) | **0.982 (299)**   | 0.981 (1)     | 0.967 (1) |
| GSE69914 (tissue)    | 0.705 (2)     | **0.773 (2)**     | 0.854 (1)     | 0.823 (4) |

Matched-cell-line AUC is robust near the shipped default: movement is
≤0.02 across the 7.5× sweep, though at σ_floor = 0.15 GSE77348 drops
slightly below V1 and raw Δβ-only. The stronger effect is top-K
usability: on GSE322563 HM450 the K=100 tied band shrinks 719 → 190 → 1
→ 1 as σ_floor increases. Tissue is genuinely σ_floor-sensitive:
GSE69914 peaks at 0.854 when σ_floor = 0.10, not at the V2.5-diff
default 0.05. Full per-positive grid:
`examples/sigma_floor_sweep.md`.

**Implication for the V2.5-sigmoid tissue advantage.** The fair reading
of the §5.2.2 V2.5-sigmoid 0.862 vs V2.5-diff 0.778 tissue contrast is
that V2.5-sigmoid wins against the *shipped* V2.5-diff default
(σ_floor = 0.05): increasing V2.5-diff's σ_floor to 0.10 raises GSE69914
tissue AUC to 0.854, within ~0.01 of V2.5-sigmoid's 0.862. We retain
V2.5-sigmoid because it gives tissue performance near that high-σ_floor
regime *and* eliminates low-`n` whole-genome tied bands under the
fixed-bandwidth default — not because every possible V2.5-diff
σ_floor setting is dominated.

### 5.3.2 δ sensitivity

δ is the only V2.5 hyperparameter selected against an empirical
cohort (GSE77348, the development cohort; §4.0). A sweep over
δ ∈ {0.1, 0.2, 0.3, 0.4, 0.5} on every cohort tests both
matched-cell-line stability *and* the appropriateness of the
development-cohort selection of δ = 0.2 on the non-development
cohorts (reproducible via `uv run scripts/delta_sensitivity_sweep.py`):

| cohort | δ=0.1 | δ=0.2 *(default)* | δ=0.3 | δ=0.4 | δ=0.5 |
|---|---|---|---|---|---|
| GSE322563 HM450      | 0.988 (tie@100=324) | **0.990 (190)** | 0.991 (1)  | 0.985 (1) | 0.985 (1) |
| GSE322563 native v2  | 0.985 (tie@100=764) | **0.986 (421)** | 0.987 (36) | 0.977 (1) | 0.977 (1) |
| GSE77348 *(dev)*     | 0.979 (tie@100=482) | **0.982 (299)** | 0.983 (5)  | 0.984 (1) | 0.984 (1) |
| GSE69914 (tissue)    | 0.827 (tie@100=2)   | **0.773 (2)**   | 0.707 (4)  | 0.668 (4) | 0.643 (3) |

(GSE77348 is the development cohort against which δ = 0.2 was
originally selected; results at δ ≠ 0.2 on that row are *post hoc*
sensitivity. The other three rows are honest sensitivity since δ was
not selected against them.)

Matched-cell-line AUC is δ-stable (≤0.014 movement across the sweep).
δ = 0.3 slightly improves AUC and shrinks tied bands, but δ = 0.2 is
retained because it was fixed before primary evaluation and functions as
conservative triage on low-`n` cohorts. Tissue shows the opposite
direction for V2.5-diff: AUC drops from 0.827 at δ = 0.1 to 0.643 at
δ = 0.5. This tissue-specific `p_diff` sensitivity does not transfer
cleanly to V2.5-sigmoid (§5.2.2). Full grid:
`examples/delta_sensitivity_sweep.md`.

### 5.4 Out-of-distribution boundary case — GSE68379

*Reported here separately. **Not** pooled with the primary-endpoint
or sensitivity tables.*

A cross-series run against the Sanger GDSC 52-breast-line panel ×
GSE69914 healthy normal (n = 52 / 50) at the Roth-validated labels
produced systematically inverted AUC:

| label | V1 | tumor_only | V2.5-diff | V2.5-diff tie_band@100 |
|---|---:|---:|---:|---:|
| validated | 0.510 | 0.222 | 0.197 | 1,602 |
| narrow | 0.485 | 0.228 | 0.169 | 1,602 |
| wide | 0.407 | 0.462 | 0.269 | 1,602 |

The root cause is a probe-level biology mismatch. At the three
Roth Fig. 5d target probes:

| probe | gene | Roth MCF-7 β (GSE322563) | Sanger MCF-7 β (GSE68379) |
|---|---|---:|---:|
| cg05251676 | EGFLAM | 0.01 | **0.92** |
| cg25338972 | ESR1 | 0.07 | **0.63** |
| cg01364137 | GATA3 | 0.02 | **0.70** |

Same cell line name, opposite methylation state at the validated
loci. MCF-7 has accumulated substantial genetic and epigenetic
drift across laboratories over decades of separate passaging.
Because the Roth labels encode "tumor hypomethylated at target"
and Sanger's MCF-7 is methylated at those targets, any scorer whose
composite contains `p_targ = P(β_tumor < 0.30)` will rank the
validated positives at the *bottom* of its distribution — which is
what V2.5 and `tumor_only` do. The inversion is the correct
response of the scorer to the underlying biology; the labels are
the component that is no longer valid on this cohort.

GSE68379 is reported as a label-transportability boundary, not a
failed generalization of V2.5. Pooling it into an average across
cohorts would misrepresent the result; redefining positives from
GSE68379's own β values would be circular (labels derived from the
data under test). It remains in the documentation as a documented
boundary case.

### 5.5 Top-hit annotation (tie-window-aware)

An annotation pass (nearest gene, TSS distance, feature class, CpG
island context) was run on each benchmark-selected top-20 window. On
low-`n` cell-line cohorts, those "top-20" lists are not stable ranked
shortlists: they are 20-record windows inside much larger tied bands
selected by the documented `candidate_id` tie-break. The four genes
shared between the GSE322563 and GSE77348 V2.5-diff top-20 windows
(ADRB1, CELF2, GFRA1, XPNPEP1) should therefore be read as window
convergence inside tied bands, not convergence of the underlying
ranking. On GSE69914 tissue (`tie_band@100 = 2`), top-20 membership is
much closer to a genuine ordering and maps to a disjoint tissue-specific
set (including RUNX2, SCML4, FOXC1/FOXCUT, MAS1L, COL21A1, MXI1,
TENM2, MAT2B). Full annotated TSVs are committed under
`examples/*/top20_annotated_*.tsv`. None of the three Roth-validated
probes appears in any top-20; the primary endpoint is rank lift/AUC, not
top-20 recovery. Operationally, the GSE69914 tissue top-20 is the subset
closest to a usable ranked shortlist; the low-`n` cell-line top-20s are
tie-window annotation views rather than stable prioritized lists.

### 5.6 p_trust sensitivity — base weights and ramp_n (V2.5-diff + V2.5-sigmoid)

Validated-label AUC is exactly invariant to `ramp_n ∈ {10, 20, 30,
50, 100}` within each cohort, on both V2.5-diff and V2.5-sigmoid.
The reason is structural: the `p_trust` ramp is `min(1, min(n_t,
n_n) / ramp_n)`; with per-cohort uniform `(n_t, n_n)` (verified
across 100k sampled records on every cohort: a single distinct
pair), `ramp_n` reduces to a global multiplicative scale on every
record's `p_therapeutic_selectivity`, which is rank-preserving.
AUC is also near-invariant (Δ ≤ 0.001) to the per-`EvidenceClass`
base-weight set on the three matched cell-line cohorts, because
their top-100 windows are 100/100 EXACT records under both axes
(§5.7 below). On GSE69914 tissue (top-100 = 33 EXACT / 67
PROXIMAL_CLOSE under V2.5-diff, 28 / 72 under V2.5-sigmoid), AUC
moves up to ±0.04 across an aggressive (0.95/0.85/0.65/0.35),
shipped (0.95/0.75/0.45/0.15), and conservative
(0.99/0.50/0.20/0.05) base-weight set on V2.5-sigmoid (0.864 vs
0.829 vs 0.866) and ±0.005 on V2.5-diff (0.773 / 0.767 / 0.776) —
the V2.5-sigmoid recommendation holds across all three weight
sets, with the conservative set the per-cohort best on tissue
under both axes. Full sweep with both axes side-by-side in
`examples/p_trust_sensitivity_sweep.{tsv,md}` at this tag;
reproduce via `uv run python scripts/p_trust_sensitivity_sweep.py`.

### 5.7 EvidenceClass distribution — full catalogs vs top-100 windows

The `p_trust` factor's discrete `EvidenceClass` lever (§3.5) controls
how strongly a candidate's evidence weight is scaled by distance from
the nearest assayed CpG. To make that lever's behavior auditable, we
report the per-class composition of each cohort's full chr5/6/10
catalog and of the top-100 window selected by the V2.5-diff axis
(committed in `examples/evidence_class_distribution.{tsv,md}` at this
tag; reproduced via `uv run python
scripts/evidence_class_distribution.py` after the per-cohort scored
JSONL chain in the reproducibility appendix).

Full-catalog distributions are dominated by `regional` records (~56%
on the HM450 cohorts, ~74% on GSE322563 native EPIC v2 and on
GSE69914 tissue) with `exact` records the smallest class (~1.1%–1.4%
of the catalog). The top-100 windows are strongly enriched: 100% of
the matched-cell-line top-100 are `exact`-class records (a ~70×
enrichment over the catalog baseline), and the GSE69914 tissue
top-100 splits 33 EXACT / 67 PROXIMAL_CLOSE — broader because the
tissue cohort's higher `n` lifts records away from `p_trust`
saturation (per-class median `p_trust` 0.95 EXACT vs 0.75
PROXIMAL_CLOSE on tissue, vs 0.06–0.10 EXACT on the matched
cell-line cohorts where `n_t = n_n = 2/2 or 3/3` puts the saturating
`min(1, n / ramp_n)` factor near zero).

The V2.5-sigmoid top-100 EvidenceClass mix is reported directly
(recomputed in `scripts/evidence_class_distribution.py` from the
same scored JSONL via `p_gap_sigmoid` with δ = 0.2, σ_fixed ≈
0.0707) rather than assumed identical to V2.5-diff's. It is
100/100 EXACT on every matched cell-line cohort (matching V2.5-diff
exactly) and 28 EXACT / 72 PROXIMAL_CLOSE on GSE69914 tissue
(slightly broader than V2.5-diff's 33 / 67 because the smooth
sigmoid promotes records with intermediate gaps that V2.5-diff's
σ_floor-saturating `p_diff` rounds to 0/1; the EvidenceClass mix is
not identical between the two axes).

The EvidenceClass-controlled benchmark is now committed at
`examples/evidence_class_stratified_benchmark.{tsv,md}` (reproduce
with `uv run python scripts/evidence_class_stratified_benchmark.py`).
It shows that an **EXACT-only validated-label benchmark is not
evaluable for these Roth candidates**: none of the three exact
candidate IDs maps to `EXACT` evidence in the scored catalogs. Instead,
ESR1 is `PROXIMAL_CLOSE`, GATA3 is `PROXIMAL`, and EGFLAM is
`REGIONAL`. The high-confidence `EXACT + PROXIMAL_CLOSE` universe
therefore tests ESR1 only; V2.5-sigmoid places ESR1 in the 99.3–99.5th
percentile on matched cell-line rows but only the 61.0th percentile on
GSE69914 tissue. Within the relevant per-class bins, V2.5-sigmoid ranks
matched-cell-line positives high (ESR1 99.7–99.8th percentile within
`PROXIMAL_CLOSE`, GATA3 96.0–98.4th within `PROXIMAL`, EGFLAM
96.2–97.8th within `REGIONAL`). Tissue is mixed by positive: GATA3
remains high within `PROXIMAL` (97.4th), EGFLAM is moderate within
`REGIONAL` (80.5th), and ESR1 is weak within `PROXIMAL_CLOSE` (60.7th).
Thus the V2.5-sigmoid tissue result is not simply an EXACT-probe
enrichment artifact, but it is also not uniform across the transported
Roth positives; it remains a label-transported ranking stress test.

**One axis-reversal worth naming.** On the `EXACT + PROXIMAL_CLOSE`
single-positive (ESR1) GSE69914 universe, V2.5-sigmoid AUC is **0.610**
while V2.5-diff (0.922), Δβ-only (0.958), and the limma-style
moderated-`t` baseline (0.925) all rank ESR1 considerably higher.
This is the only axis × universe cell where V2.5-sigmoid is not the
top axis on tissue, and it is on a single-positive AUC over a
restricted EvidenceClass universe rather than the §5.2.2 / §5.3
all-positive WG result. The mechanism is the smooth-sigmoid response:
within the high-confidence subset, V2.5-sigmoid's saturating gap
factor compresses ESR1's tissue Δβ signal more than the discrete-step
`p_diff` does, while the moderated-`t` and Δβ-only axes have no
saturation at all. At the GSE69914 ESR1 record, both V2.5 axes share
`p_targ = 0.0586` and `p_trust = 0.75`; V2.5-diff uses
`p_diff = 0.3196`, whereas V2.5-sigmoid uses
`p_gap_sigmoid = sigmoid((0.111 - 0.2) / 0.0707) = 0.2213`
(scores 0.0140 vs 0.00973). The reversal is therefore driven by the
gap term, not by EvidenceClass or tumor-side targetability. We report it
for honesty and because the
EvidenceClass-restricted universe is a real reviewer ask, not as a
revision of the §6.1 recommendation.

(`scripts/evidence_class_stratified_benchmark.py` uses `−t_mod` as
the limma-style ranking score, matching `scripts/limma_candidate_ranking.py`'s
"tumor=1 / normal=0 group assignment, so Δ = β_t − β_n; targetable
site has Δ < 0 → t < 0 → −t > 0" sign convention; verified at this
tag.)

### 5.8 Canonical R `limma::lmFit + eBayes` parity

The pure-Python Smyth (2004) empirical-Bayes moderated-`t`
implementation in `scripts/limma_ebayes.py` is the source of every
"limma-style moderated-`t`" number reported in this paper (including
§4.3, §5.1, §5.2.2, §5.7, and the §5.9 feature-matched controls; any
section that adds a limma-style column inherits from the same
implementation). To address the reviewer ask of the previous tag, we ran
canonical `limma::lmFit(β_matrix, design) %>% eBayes()` (R 4.5.3,
Bioconductor `limma` 3.66.0) on the same sample-level β matrices and
group vectors the Python pipeline consumes via the shared adapters in
`scripts/run_limma_per_cohort.py`, then compared per-probe
moderated-`t` statistics by Spearman / Pearson correlation, top-K |t|
Jaccard overlap, and absolute-difference diagnostics:

| cohort | n_t / n_n | probes (Py / R / both) | Spearman `t_mod` | Pearson `t_mod` | top-100 | top-1000 |
|---|:---:|---|---:|---:|---:|---:|
| GSE69914 (tissue) | 305 / 50 | 485,512 / 485,512 / 485,512 | **1.000000** | **1.000000** | **1.000** | **1.000** |
| GSE77348          | 3 / 3   | 485,577 / 485,577 / 485,577 | **1.000000** | **1.000000** | **1.000** | **0.999** |
| GSE322563 HM450   | 2 / 2   | 937,690 / 881,183 / 881,183 | **0.999674** | **0.997085** | **0.960** | **0.903** |

The tissue cohort and the GSE77348 cell-line cohort show
floating-point parity (Spearman/Pearson 1.000000; max |Python − R|
`t_mod` = 1.9 × 10⁻⁴ on tissue, 0.43 maximum / 0.002 median on
GSE77348). On GSE322563 HM450 (n = 2/2, residual `df = 0`) R's
`eBayes` declines to fit ~56,500 underdetermined probes that Python
substitutes a finite value for; on the 881k probes both implementations
fit, Spearman is 0.9997 and the top-100 / top-1000 |t| Jaccard is 0.96
/ 0.90. This is a low-replicate limitation of the limma comparison, not
a tissue-row caveat: the n = 305/50 tissue parity is exact, and that is
where the limma-style baseline is used for the strongest between-method
contrast. Differences on GSE322563 are fully accounted for by R's
stricter NaN-coding of `df = 0` cohort probes; none of the probes
affecting the PAPER §5.2.2 candidate-mapped AUC numbers (which all have
`evidence_class ∈ {EXACT, PROXIMAL_CLOSE, PROXIMAL, REGIONAL}` with
non-zero variance after candidate-mapping) is in the dropped subset.
Per-cohort artifacts at `examples/r_limma_parity_{gse69914,gse77348,gse322563}.{tsv,md}`;
reproduce via `uv run python scripts/run_r_limma_parity.py --cohort
<cohort>` (requires R + Bioconductor `limma`).

The §4.3 "limma-style moderated-`t`" naming convention is therefore
strictly accurate — the math is canonical `limma::eBayes`, the
implementation is independent (pure stdlib + numpy/scipy) but
functionally equivalent.

### 5.9 Feature-matched negative-universe controls

The §5.2.2 / §5.7 WG numbers compare the three Roth-validated
positives against the *full* candidate universe. That answers
"are the positives high?", but it does not control for the
candidate features the scorer is permitted to use as discrimination
proxies (`EvidenceClass` / probe distance, `pam_family`, CpG status,
chromosome). To answer the stricter question — "are the positives
high *versus comparable candidates*?" — we report the positives'
rank against a feature-matched negative universe, the WG subset of
candidates that share the positive's exact (`EvidenceClass`,
`pam_family`, `is_cpg_pam`, `chrom`) signature. Empirical *p* is the
fraction of matched negatives that score strictly better than the
positive (ties contribute 0.5; one-sided null). Matched-pool sizes
range from 24,727 (GATA3 on GSE322563 HM450) to 271,676 (EGFLAM on
GSE322563 native EPIC v2), so the resolution of the empirical *p*
floor is 4 × 10⁻⁶ – 4 × 10⁻⁵, well below the rank values reported.
These are descriptive per-positive empirical *p*-values, uncorrected for
multiplicity; they are denominator-confounding audits, not a
multiplicity-corrected discovery test.
Bold marks V2.5-sigmoid cells with empirical *p* < 0.05 because it is
the recommended axis; it does **not** mark the row-best axis.

| cohort | gene | V2.5-sigmoid p | V2.5-diff p | Δβ-only p | limma-style p |
|---|---|---:|---:|---:|---:|
| GSE322563 HM450     | *ESR1*   | **0.0021** | 0.0022 | 0.0021 | 0.0027 |
| GSE322563 HM450     | *GATA3*  | **0.0359** | 0.0359 | 0.0558 | 0.1037 |
| GSE322563 HM450     | *EGFLAM* | **0.0197** | 0.0188 | 0.0301 | 0.0365 |
| GSE322563 native v2 | *ESR1*   | **0.0018** | 0.0027 | 0.0018 | 0.0033 |
| GSE322563 native v2 | *GATA3*  | **0.0373** | 0.0372 | 0.0609 | 0.1114 |
| GSE322563 native v2 | *EGFLAM* | **0.0258** | 0.0247 | 0.0379 | 0.0447 |
| GSE77348            | *ESR1*   | **0.0033** | 0.0030 | 0.0033 | 0.0131 |
| GSE77348            | *GATA3*  | **0.0186** | 0.0169 | 0.0204 | 0.0093 |
| GSE77348            | *EGFLAM* | **0.0390** | 0.0381 | 0.0612 | 0.0986 |
| GSE69914 (tissue)   | *ESR1*   | 0.4384     | 0.0784 | 0.0273 | 0.0554 |
| GSE69914 (tissue)   | *GATA3*  | **0.0187** | 0.2085 | 0.2595 | 0.3772 |
| GSE69914 (tissue)   | *EGFLAM* | 0.2588     | 0.4968 | 0.7632 | 0.7192 |

Three results to flag honestly. (1) On every matched-cell-line row
(9 of 12 rows above), every positive is in the top ~4% of its
feature-matched pool under V2.5-sigmoid (`p` ∈ [0.0018, 0.0390]),
and V2.5-sigmoid is among the strongest axes on every cell. The WG
rank lift is therefore **not explained by the matched features tested
here** on matched cell-line cohorts. (2) On GSE69914 tissue, GATA3 stays
strong under V2.5-sigmoid (`p` = 0.019, 98th matched-percentile)
while every other axis collapses to `p` ≈ 0.21–0.38 — the smooth
sigmoid is genuinely doing useful work on this single positive. (3)
EGFLAM is moderate on tissue under V2.5-sigmoid (`p` = 0.26) and
*all* axes are weaker on tissue than on cell lines for this
positive; ESR1 on tissue under V2.5-sigmoid sits at `p` = 0.44
(near-random), confirming the §5.7 ESR1-reversal note from a
different angle: in the high-confidence subset, the smooth sigmoid
compresses the tissue ESR1 signal while V2.5-diff and Δβ-only
preserve it.

The feature-matched control is a **denominator-confounding audit**,
not validation. It answers "are the positives high versus comparable
candidates?" — not "does the method discover new targets?" Strict
within-chromosome matching is the most conservative choice; loosening
to chromosome-class-matched or random-chromosome-matched negatives
is documented as a §6.3 follow-up. Per-cohort artifact at
`examples/feature_matched_negative_controls.{tsv,md}`; reproduce via
`uv run python scripts/feature_matched_negative_controls.py`.

### 5.10 Roth System B HEK293T/HCT116 pre-registered extension

We pre-registered an independent Roth System B extension on the
HEK293T / HCT116 validation panel before scoring (four `prereg-*`
tags). The strict pre-applied transport rule (β > 0.7 / β < 0.3 with
≥ 10 reads in ≥ 2 replicates per cell line) excluded both *VEGFA*
control targets on public ENCODE RRBS — including T9, the planned
editable-but-non-selective discriminator that V2.5 should demote — and
a secondary WGBS / EPIC v2 backend scan did not rescue T9. The retained
*EMX1* T4 / *PRDX4* T5 directional pair has n_pos = n_neg = 1 per
direction, so the AUC = 1.0 it produces is essentially uninformative
about V2.5-sigmoid's selectivity claim. **The System B selectivity
claim is therefore not evaluable on the public methylation backend at
this revision; this section is a workflow audit (the pre-registered
transport gate worked; the decisive control was correctly excluded
rather than absorbed into the score) rather than a validation.** Full
transport audit, secondary-backend scan, T9 hypothetical (~18-fold
separation under Roth's BSS polarity, explicitly **not** counted as
validation), retained directional-pair scores, and Endpoint-2
target-side editability check are in **Appendix B**. Reproduce with
`uv run python scripts/score_roth_transport_subset.py`.

---

## 6 · Discussion

### 6.1 Recommended use

| Use case | Recommendation |
|---|---|
| New probabilistic prioritization runs | **Recommended probabilistic prioritization axis for hypothesis generation: V2.5-sigmoid.** It preserves the matched-cell-line rank lift of V2.5-diff and Δβ-only (where Δβ-only is already very strong on this n_pos = 3 endpoint, AUC 0.974 on GSE322563 HM450 vs V2.5-diff's 0.990), eliminates V2.5-diff's low-`n` whole-genome tied-band pathology (`tie_band@100 = 1` vs 421–1,493), and improves the shipped-default tissue stress-test AUC. The primary matched-cell-line AUC gain over Δβ-only is small and descriptive at `n_pos = 3`; the tissue gain is per-positive heterogeneous (§6.2 limitation 7) and partly reflects the shipped V2.5-diff σ_floor default (§5.3.1: V2.5-diff at σ_floor = 0.10 reaches tissue AUC 0.854, within ~0.01 of V2.5-sigmoid). Prospective wet-lab validation or newly generated labels remain required. |
| Backward-compatible probabilistic runs | **V2.5-diff**. Retained for AUC parity with earlier scored JSONLs; not preferred for new WG-scale prioritization because low-`n` tied bands expand to 421–1,493 records. |
| Stable deterministic top-K | **V1 final_score**. Continuous and retained for backward compatibility, not AUC leadership. |
| Diagnostic ablation | **V2 tumor_only**. Useful AUC sanity check, but top-K collapses into thousands-record tied bands. |
| Cross-series label transport | **Unsupported**. Cell-line drift can invert label logic; do not pool with prioritization-mode benchmarks. |

Mode enum mapping:

- V2.5-sigmoid: `tumor_plus_gap_sigmoid`
- V2.5-diff: `tumor_plus_differential_protection`
- V2 tumor_only: `tumor_only`

On low-`n` matched cell-line cohorts, the V2.5 generation should be
read as a **rank-lift tool**, not a stable top-K shortlist generator.
It lifts the three Roth positives into the upper few percentiles of
millions of candidates (§5.1: ESR1 top ~0.1%, GATA3 top ~1%, EGFLAM
top ~2–4.6%), but the visible top-20 can be a window inside a
hundreds-record tied band. Selecting individual wet-lab candidates from
that class requires secondary evidence such as annotation, chromatin
context, guide quality, and off-target risk.

### 6.2 Limitations

1. **Small validated set.** The primary endpoint has three positives.
   AUC is therefore rank lift on three sites, not generalized discovery
   evidence; per-positive ranks are reported in §5.1.
2. **Low-replicate tied bands.** For `n < ramp_n = 30`, probabilistic
   scores can form hundreds-to-thousands-record tied bands. The
   benchmark reports `tie_band_size_at_k` and adversarial P@K intervals
   so top-K claims remain explicit.
3. **Catalog and platform scope.** §5.1 uses chr5/6/10 HM450-scope
   negatives because the validated Roth targets live there; §5.2.2
   extends the denominator to frozen WG HM450 and native EPIC v2
   catalogs; §5.9 adds within-chromosome feature-matched controls.
   Cross-chromosome and chromosome-class matched controls remain future
   work.
4. **Cell-line drift.** GSE68379 shows that same-named MCF-7 stocks can
   have opposite methylation at Roth targets. Cohort-specific
   methylation must be verified before editing follow-up.
5. **Independent-biology System B not evaluable.** The HEK293T/HCT116
   System B extension was pre-registered, but the decisive VEGFA T9
   non-selective discriminator failed public-backend transport because
   of low coverage; the retained EMX1/PRDX4 subset is a polarity
   diagnostic, not a V2.5 selectivity validation.
6. **V2.5-diff σ_floor.** The 0.05 per-side σ floor is a V2.5-diff
   constant. It is robust on matched cell lines but tissue-sensitive
   (§5.3.1). The current recommendation is V2.5-sigmoid, whose
   fixed-bandwidth behavior is evaluated separately in §5.2.2.
7. **Tissue per-positive heterogeneity + EvidenceClass-restricted
   reversal.** The cohort-level GSE69914 tissue gain under V2.5-sigmoid
   collapses onto GATA3 alone (§5.9 feature-matched controls: GATA3
   *p* = 0.019, ESR1 *p* = 0.44 matched-near-random, EGFLAM *p* = 0.26
   moderate). In the high-confidence `EXACT + PROXIMAL_CLOSE` GSE69914
   universe, only ESR1 is evaluable, and V2.5-sigmoid trails V2.5-diff,
   Δβ-only, and the limma-style moderated-`t` baseline (§5.7). This
   single-positive restricted-universe result limits the tissue
   stress-test interpretation but does not revise the whole-genome
   all-positive recommendation.
8. **Top-hit annotation is descriptive.** Low-`n` top-20 lists are
   tie-window slices, not stable rankings. They require secondary
   biological filters before wet-lab selection.
9. **No prospective wet-lab validation.** The claim is scoring and
   tie-band-aware benchmarking on public methylation data, not editing
   efficiency at named sites.

### 6.3 Next steps

1. Prospective wet-lab editing validation of top V2.5-sigmoid
   candidates.
2. Random-chromosome-matched and chromosome-class-matched
   negative-universe controls — the within-chromosome-strict variant
   is now committed at this tag (§5.9); cross-chromosome and
   chromosome-class matching are the remaining loosening axes.
3. A one-shot `regime` YAML selector on top of the shipped
   `tumor_plus_gap_sigmoid` mode.
4. Alternative gap factors: SE-on-mean `p_diff`, exact
   difference-of-Betas, or mixture-aware variants.
5. Continuous observation-confidence modeling beyond discrete
   `EvidenceClass` bins.
6. A second independent-lab MCF-7 / MCF-10A EPIC cohort if one becomes
   public.
7. Targeted BSS or WGBS for the Roth VEGFA T3/T9 HEK293T/HCT116 controls,
   so the pre-registered System B T9-demotion test can be evaluated
   instead of reduced to the transport-confirmed T4/T5 subset (§5.10).
8. Formal SCREEN cCRE integration.

---

## 7 · Methods — key equations

### 7.0 Canonical axis names

The benchmark reports two V2.5-generation modes alongside
deterministic, tumor-only, Δβ-based, and limma-style comparison axes.
Enum names are what `CohortConfig` accepts in YAML; prose names are what
§5 tables and recommendations use.

- **V1** — `final_score`; stable-release deterministic score with
  `tie_band = 1` by construction. Retained for backward compatibility
  and deterministic top-K lists.
- **V2 tumor_only** — `probabilistic_mode: tumor_only`; diagnostic
  ablation. Useful for AUC sanity checks, not a prioritization axis.
- **V2.5-diff** — `probabilistic_mode: tumor_plus_differential_protection`;
  original V2.5 composite
  `p_targ × p_diff × p_trust`. Retained for backward compatibility and
  AUC parity on matched cell-line cohorts.
- **V2.5-sigmoid** — `probabilistic_mode: tumor_plus_gap_sigmoid`;
  `p_targ × p_gap_sigmoid × p_trust`. Recommended probabilistic
  prioritization axis across tested non-boundary regimes.
- **Δβ-only** — `naive_selectivity`; raw methylation-gap baseline.
- **Δβ_z** — post-hoc σ-normalized raw-gap baseline computed from
  scored JSONL artifacts.
- **limma-style moderated-t** — probe-level Smyth 2004 / limma-style
  DMR comparator.

The terms "probability-scale selectivity score" (paper-level) and
the field `p_therapeutic_selectivity` (code-level) refer to the
same emitted scalar under every probabilistic mode. "V2.5" without
a `-diff` / `-sigmoid` qualifier refers to the *generation*
(`tumor_plus_*`) and is used when a statement applies to both
variants.

**Key probabilistic factors:**

```
p_targ = P(β_tumor < 0.30)   via method-of-moments Beta(α, β) fit,
p_prot = P(β_normal > 0.50)  via same CDF dispatch.     [V2 only]
p_diff = 1 − Φ((δ − (μ_n − μ_t)) / σ_Δ)
         where σ_Δ = √(max(σ_t, 0.05)² + max(σ_n, 0.05)²)
                          and σ_k = IQR_k / 1.349.      [V2.5-diff]
p_gap_sigmoid = 1 / (1 + exp(−((μ_n − μ_t) − δ) / σ_fixed))
         (δ = 0.2, σ_fixed = √2 × σ_floor ≈ 0.0707 by default;
          σ_fixed strictly positive, enforced at CohortConfig
          boundary). Numerically-stable closed form, pure stdlib.
                                                      [V2.5-sigmoid]
p_trust = EvidenceClass.base × min(1, min(n_t, n_n) / 30)
         (EvidenceClass.base ∈ {0.95, 0.75, 0.45, 0.15, 0.0}).
```

V1 uses a deterministic sequence × selectivity × confidence score minus
heterogeneity and low-coverage penalties; it is continuous-valued and
therefore has `tie_band = 1` in every benchmark here. Probabilistic
composites by mode:

```
tumor_only                           : p_targ × p_trust
tumor_plus_normal_protection         : p_targ × p_prot × p_trust
tumor_plus_differential_protection   : p_targ × p_diff × p_trust
    (= V2.5-diff)
tumor_plus_gap_sigmoid               : p_targ × p_gap_sigmoid × p_trust
    (= V2.5-sigmoid; recommended probabilistic prioritization axis, §5.2.2)
```

Validators enforce mode-specific audit fields and exact composite
consistency for the emitted `p_therapeutic_selectivity` scalar.

**Benchmark ranking contract.** Let the scored set be N records, each
with a score `s_i` and a ground-truth label `y_i ∈ {0, 1}` (1 for
positive). The benchmark sorts records by the total order

```
    (s_i, c_i)   with primary  -s_i  ascending
                      secondary c_i   ascending
```

where `c_i` is the candidate_id. Let `r_1, r_2, …, r_N` be the
resulting ranking (smallest `(-s, c)` first). Define:

```
    top-K               := {r_1, …, r_K}
    cutoff_score        := s at position K in the sort, i.e. s_{r_K}
    tied_band           := {r_i : s_{r_i} = cutoff_score}           # records at the cutoff
    tie_band_size_at_k  := |tied_band|
    k_band              := |tied_band ∩ top-K|                      # records in top-K at the cutoff
    drop_slots          := tie_band_size_at_k − k_band              # records at cutoff outside top-K
    band_pos            := |{r_i ∈ tied_band : y_{r_i} = 1}|
    committed_pos       := |{r_i ∈ top-K \ tied_band : y_{r_i} = 1}|
```

Then:

```
    P@K_observed := (committed_pos + |{r_i ∈ tied_band ∩ top-K : y_{r_i} = 1}|) / K

    P@K_min      := (committed_pos + max(0, band_pos − drop_slots)) / K
    P@K_max      := (committed_pos + min(band_pos, k_band))          / K
```

`P@K_min` pushes as many tied-band positives as possible outside top-K;
`P@K_max` pulls as many as possible inside top-K. By construction:

```
    P@K_min ≤ P@K_observed ≤ P@K_max
```

and the interval collapses when `tie_band_size_at_k = 1`. Recall uses
`n_positives` as the denominator. AUC is Mann-Whitney U /
(`n_pos · n_neg`) with ties contributing 0.5.

---

## Data and code availability

- **Code**: <https://github.com/AllisonH12/thermocas9>. Cite tag **`paper-5-10f`** for this document. 245 tests pass under `uv run pytest -q`.
- **Citable archive (DOI)**: a Zenodo release archive of the tagged revision is planned at the time of preprint posting; the GitHub → Zenodo integration mints a DOI for each GitHub release tag. The DOI will be added to this section and to the citation block below before journal-version submission. Until then, the immutable git tag above is the canonical citable identifier.
- **Cohort data**: publicly-downloadable GEO series GSE322563, GSE77348, GSE69914, GSE68379; build scripts in `scripts/build_gse*_cohort.py` produce the per-probe summary TSVs in `data/derived/*_cohort/`. Positives-list builder at `scripts/build_roth_positives.py` (requires the Ensembl REST `/map` endpoint for the hg38 → hg19 liftover of Roth Fig. 5d coordinates).
- **Reference data**: UCSC hg19 `refGene.txt.gz` and `cpgIslandExt.txt.gz` (fetched on demand; gitignored).
- **Benchmark artifacts**: every `BenchmarkResult` JSONL row under `examples/*_roth_labels/` carries `precision_at_k`, `precision_at_k_{min,max}`, `recall_at_k`, `recall_at_k_{min,max}`, `roc_auc`, `tie_band_size_at_k`, and `tie_break_policy`.
- **EvidenceClass controls**: `examples/evidence_class_distribution.{tsv,md}` audits full-catalog and top-100 EvidenceClass composition for V2.5-diff and V2.5-sigmoid; `examples/evidence_class_stratified_benchmark.{tsv,md}` reports the WG / chr5/6/10 EvidenceClass-controlled rank-lift benchmark for V2.5-sigmoid, V2.5-diff, Δβ-only, and the limma-style baseline.
- **Feature-matched controls**: `examples/feature_matched_negative_controls.{tsv,md}` reports the within-chromosome matched-negative audit used in §5.9; reproduce with `uv run python scripts/feature_matched_negative_controls.py`.
- **Sweeps.** `examples/sigmoid_delta_sigma_wg_sweep.{tsv,md}` and `scripts/sigmoid_delta_sigma_wg_sweep.py` reproduce the V2.5-sigmoid `(δ × σ_fixed)` 5 × 6 grid (§5.2); `examples/tissue_per_positive_wg_ranks.{tsv,md}` reports the GSE69914 per-positive WG-rank table (§5.3).
- **System B artifacts** (§5.10 / Appendix B). `data/positives/positives_roth_hek_hct_v0.tsv`, `data/derived/roth_hek_hct_transport.tsv`, `data/derived/roth_hek_hct_secondary_backend_scan.tsv`, and `data/derived/roth_hek_hct_subset_{scores,benchmark}.tsv`; pre-registration audit trail in `prereg/2026-04-24-hek-hct-system-b-transport-addendum.md`. Reproducer: `uv run python scripts/score_roth_transport_subset.py`.
- **Annotated top-20 TSVs**: under `examples/*/top20_annotated_*.tsv`.

## Reproducing the cross-cohort matrix

The exact build, score, benchmark, and aggregation commands are kept in
the scripts referenced throughout §4–§5 and in the committed provenance
files under `data/derived/*PROVENANCE*`. Large scored JSONLs are
gitignored but reproducible from the public GEO inputs, frozen catalogs,
cohort YAMLs, and scripts in this tag.

---

## Appendix A · V2 → V2.4 audit trail (historical formulations)

This appendix preserves the full audit trail for the threshold-based
V2 composite and its V2.4 intermediate, moved out of the main body
in an earlier revision so the main paper carries the final-method
narrative (Δβ-only / V1 / V2.5-diff / V2.5-sigmoid / limma-style
moderated-`t`). The V2 / V2.4 modes remain selectable via
`probabilistic_mode` (`tumor_plus_normal_protection` and
`tumor_only` respectively) for backward-compatible reproduction of
pre-V2.5 scored JSONLs; they are not recommended for new prioritization
runs.

### A.1 V2 first-pass composite

The first-pass composite followed a standard decomposition:

```
p_therapeutic_selectivity = p_targ × p_prot × p_trust
```

- `p_targ = P(β_tumor < 0.30)` — probability that the tumor's PAM cytosine is unmethylated enough for ThermoCas9 to cleave. Estimated via a method-of-moments Beta fit to `(mean, IQR/1.349)` with the regularized incomplete beta computed in pure stdlib via Lentz's continued fraction; falls back to a piecewise-linear CDF when the Beta moments are ill-defined.
- `p_prot = P(β_normal > 0.50)` — probability that the normal's PAM cytosine is methylated enough to protect the normal cell.
- `p_trust` — confidence factor that scales with `EvidenceClass` (0.95 at EXACT, 0.15 at REGIONAL, 0.0 at UNOBSERVED) and saturates linearly in `min(n_tumor, n_normal)` up to `ramp_n = 30`.

### A.2 Empirical failure of `p_prot`

We ran a factor ablation on a structural surrogate for the Roth
MCF-7 vs. MCF-10A comparison: GEO series GSE77348, an independent 2016
breast cell-line methylation series on HM450 (not the Roth samples).
With 1,687 gene-region-filtered "positives" at the Roth genes, the
measured AUCs were:

| axis | AUC (loose) | AUC (tight) |
|---|---:|---:|
| `v2_tumor_only` (p_targ × p_trust) | 0.733 | 0.770 |
| `p_targ` alone | 0.683 | 0.717 |
| V1 `final_score` | 0.657 | 0.628 |
| `naive_selectivity` (β_n − β_t) | 0.608 | 0.571 |
| `p_targ × p_prot` | 0.503 | 0.488 |
| **`p_prot` alone** | **0.384** | **0.343** |

`p_prot` is anti-predictive. Multiplying it against `p_targ` destroys
the signal. The cause is the semantics: `P(β_normal > 0.5)` assumes
that the normal comparator *methylates the target*. That assumption
holds for Roth's specific MCF-10A cell line at the three *ESR1* /
*GATA3* / *EGFLAM* probes they chose to validate (genes that are
silenced in non-tumorigenic mammary epithelium) but does not generalize
to bulk adjacent-normal cohorts or even to MCF-10A at gene-body probes.

### A.3 V2.4 intermediate — keep the composite, drop `p_prot` by default

The intermediate fix kept the composite but allowed users to opt in
or out of the anti-predictive factor. Default changed to
`tumor_only` = `p_targ × p_trust`. This improved global AUC but the
top-100 collapsed onto "always-unmethylated" loci (candidates where
both tumor *and* normal are low-methylated) — good AUC, unusable for
target prioritization. The deeper fix — the V2.5-generation
differential gap factor
(§3) — was needed.

---

## Appendix B · Roth System B HEK293T/HCT116 pre-registered extension (full audit trail)

The §5.10 main-body summary establishes that the System B selectivity
test is not evaluable on the public methylation backend at this
revision. This appendix carries the supporting transport audit, the
secondary-backend scan, the T9 hypothetical, the retained directional
pair scores, and the Endpoint-2 target-side editability check.

### B.1 Pre-registration and panel choice

After the §5.9 matched-negative audit, we pre-registered an independent
Roth System B extension using the HEK293T / HCT116 validation panel from
Roth Fig. 5a / Extended Data Fig. 9 / Supplementary Fig. 10. This panel
is biologically independent of the MCF-7 / MCF-10A Fig. 5d breast
targets: different cell lines, different loci, and a different
methylation modality (Roth used ENCODE RRBS plus targeted BSS, then
validated editing by cleavage / indel readouts). The planned primary
discriminator was *VEGFA* T9: an editable but non-selective target
that `tumor_only` should rank high and V2.5 should demote. The
pre-registration was frozen before scoring at `prereg-coords`,
`prereg-transport`, and `prereg-transport-addendum`; scored subset
artifacts were frozen at `prereg-scored`.

### B.2 Transport rule and per-target outcome

The transport rule was deliberately strict and was applied **before**
any rank interpretation: ENCODE RRBS had to agree with Roth's BSS
state using β > 0.7 for methylated, β < 0.3 for unmethylated, and
≥ 10 reads in at least two replicates per cell line at the PAM cytosine
or nearest covered cytosine within 50 bp. Under that rule, the two
direction-flip targets transported, but the *VEGFA* controls did not:

| target | role in pre-reg | HEK293T transport | HCT116 transport | consequence |
|---|---|---|---|---|
| *VEGFA* T3 | protected/off in both lines | low coverage | low coverage | excluded |
| *VEGFA* T9 | editable/non-selective discriminator | low coverage | low coverage | excluded |
| *EMX1* T4 | HEK293T-selective positive | confirmed | confirmed | retained |
| *PRDX4* T5 | HCT116-selective positive | confirmed | confirmed | retained |

The low-coverage T9 row is the key result of System B. Because T9 was
the pre-registered falsification / success discriminator for V2.5
selectivity, the System B selectivity claim is **not evaluable** on the
public ENCODE RRBS backend.

### B.3 Secondary-backend scan

A secondary-backend scan was performed before scoring and did not
rescue T9: ENCODE has RRBS / DNAme-array data but no WGBS experiment
under the HCT116 / HEK293 / HEK293T terms; one public HCT116 WGBS file
had sparse one-sample *VEGFA* coverage (T9 coverage 4, T3 coverage 7);
a HEK293T WGBS proxy did not transport-confirm T9; and native EPIC v2
probes near T9/T3 were hundreds of bases away with no matched
HEK293T / HCT116 EPIC v2 pair. Recorded in
`prereg/2026-04-24-hek-hct-system-b-transport-addendum.md` and
`data/derived/roth_hek_hct_secondary_backend_scan.tsv`.

### B.4 T9 hypothetical (explicitly **not** counted as validation)

Under Roth's own BSS polarity values, the missing T9 control would
have β_target ≈ β_comparator ≈ 0, so V2.5-sigmoid would assign
`p_gap_sigmoid = sigmoid((0 − 0.2) / 0.0707) ≈ 0.0558`. A
direction-specific *EMX1* / *PRDX4* selective positive with
β_target ≈ 0 and β_comparator ≈ 1 would have
`p_gap_sigmoid ≈ 0.99999`. Thus, at matched `p_targ` and matched
`p_trust`, the pre-registered T9 separation would be approximately
18-fold. We do **not** count this as validation, because T9 failed the
pre-registered public-backend transport rule.

### B.5 Retained directional pair (uninformative AUC)

The pre-registered scored-subset run is therefore a
*transport-confirmed subset diagnostic*, not the original T9-demotion
benchmark. On the retained directional pair, all axes rank the
direction-specific selective target above the opposite-direction target;
with n_pos = n_neg = 1 per direction, the AUC = 1.0 is essentially
uninformative about V2.5-sigmoid's selectivity claim:

**HEK293T target / HCT116 protected**

| Role | Target | β target / β comparator | V2.5-sigmoid score |
|---|---|---:|---:|
| Selective positive | *EMX1* T4 | 0.007 / 1.000 | 0.010000 |
| Opposite-direction control | *PRDX4* T5 | 0.751 / 0.061 | 0.000000 |

**HCT116 target / HEK293T protected**

| Role | Target | β target / β comparator | V2.5-sigmoid score |
|---|---|---:|---:|
| Selective positive | *PRDX4* T5 | 0.061 / 0.751 | 0.008128 |
| Opposite-direction control | *EMX1* T4 | 1.000 / 0.007 | 0.000000 |

We therefore do not interpret this AUC as method validation.

### B.6 Endpoint 2 — target-side editability diagnostic only

Endpoint 2, the HEK293T binary editability check, remains usable after
transport filtering only as a target-side editability diagnostic. On
confirmed HEK293T calls, `p_targ` alone separates four editable /
unmethylated targets (*T4*, *T10*, *T11*, *T12*) from three
protected / not-edited targets (*T5*, *T14*, *T15*) with AUC = 1.000.
*T3* and *T9* are excluded for low coverage, and *T13* is excluded as
ambiguous (β = 0.535 against Roth's methylated label). This supports
the narrow statement that the target-side methylation factor recovers
Roth's HEK293T editability labels when independent RRBS transport is
confirmed; it does **not** test V2.5's selectivity claim because the
non-selective T9 decoy is missing.

### B.7 Practical interpretation

System B improves the audit trail and validates the transport-gating
workflow, but it does not add an independent selectivity validation
because the pre-registered T9 discriminator is missing. It shows that
the pre-registered transport gate can reject a tempting benchmark when
the public methylation backend does not support the decisive control.
The scored subset is useful for debugging polarity and target-side
editability, but the paper's recommendation remains grounded in the
MCF-7 / MCF-10A cohorts, the tissue stress test, and the §5.9
matched-negative audit. Reproduce with
`uv run python scripts/score_roth_transport_subset.py`; outputs are
`data/derived/roth_hek_hct_subset_scores.tsv` and
`data/derived/roth_hek_hct_subset_benchmark.tsv`.

---

## Citation

If this paper or the associated code is useful in your work, please
cite the upstream paper:

> Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.
> Molecular basis for methylation-sensitive editing by Cas9.
> *Nature* (2026). DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

For the probabilistic scoring framework specifically, cite this paper
plus the code URL above.

## Acknowledgements

This work describes a framework built as an educational / research
project. ThermoCas9's biochemical characterization and the validated
target sites are entirely due to Roth et al.; any errors in the
re-implementation of the scoring or benchmarking are the author's.
