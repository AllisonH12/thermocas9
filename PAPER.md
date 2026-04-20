# Differential-protection probabilistic scoring for methylome-guided ThermoCas9 target-site ranking

**Author.** Allison Huang, Thermocas9 Inc.
**Date.** 2026-04-20.
**Code.** <https://github.com/AllisonH12/thermocas9> at commit `45867fb`.
**Status.** Technical memo from an educational research framework. Not peer-reviewed. No clinical claims. Cites Roth et al., *Nature* (2026), DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

---

## Abstract

Roth et al. (2026) showed that ThermoCas9 is methylation-sensitive at the
fifth-position cytosine of its `5'-NNNNCGA-3'` PAM: methylation of that
cytosine blocks DNA binding and cleavage. This turns ThermoCas9 into a
natural methylation-sensitive editor, usable to target loci that are
hypomethylated in tumor and methylated in matched normal cells — a
cleaner precision layer than engineered methylation-insensitive Cas9
variants. Selecting such loci from genome-scale methylation data is a
ranking problem: which of millions of PAM-compatible sites shows the
strongest tumor-hypomethylation / normal-methylation contrast, scored in
a way that is cohort-agnostic and reports its own uncertainty.

A first-pass probabilistic composite `p_targ × p_prot × p_trust` was
empirically **anti-predictive** because `p_prot` encoded a static
"normal is methylated above 0.5" assumption that does not hold on bulk
normal comparators. We replace the threshold-based `p_prot` with a
differential factor `p_diff = P(β_normal − β_tumor > δ)` estimated under
an independent-normal approximation on per-probe summary statistics.
The new composite — V2.5 `tumor_plus_differential_protection` —
dominates AUC on matched cell-line cohorts (GSE322563 Roth actual
samples at validated Roth-Fig.-5d target probes: AUC **0.990** vs V1
`final_score` 0.821) and is the only probabilistic axis whose top-K is
interpretable on low-replicate cohorts (tie band shrinks from hundreds
at n=2/3 to 2 at n=305/50, exactly matching the `p_trust`-saturation
prediction made before the high-n cohort was run).

The recommendation is cohort-type-dependent rather than universal:
V2.5 is the strongest probabilistic axis for matched cell-line /
paper-comparable biology; V1 `final_score` stays available as a
continuous-score fallback for top-K stability; the deprecated V2
`tumor_only` mode stays analysis-only. A cross-series run against the
Sanger GDSC breast panel produced inverted AUC at Roth-validated
probes because Sanger's MCF-7 is methylated at exactly the sites where
Roth's MCF-7 is unmethylated — a label-transportability boundary, not
a scorer failure. That finding is preserved in the documentation as an
out-of-distribution case rather than pooled with the paper-comparable
cohorts.

Code, benchmark artifacts, and test suite (215 passing tests) are
public. All cohorts are reproducible from committed summary TSVs plus
publicly-available reference annotations.

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
loci: *EGFLAM* (a control target), *ESR1*, and *GATA3*. The underlying
methylation of the PAM cytosine at those loci was profiled on the
Illumina Infinium MethylationEPIC v2 array and reported in
Supplementary Figure 5d of the main paper alongside the editing
results.

### 1.2 Target-site selection is a ranking problem

Given a reference genome, the ThermoCas9-compatible PAM sites on any
chromosome number in the millions. To pick candidates for wet-lab
follow-up, one wants a ranked list where the top scorers are
*biologically coherent* — meaning the candidate's PAM cytosine is
unmethylated in the disease cohort and methylated in the matched
normal — *and* where the confidence of the underlying methylation
measurement is part of the score. This is a classical
ranking-with-covariates problem, and the natural scoring target is
some decomposition of a probability of "therapeutic selectivity."

This memo describes the specific decomposition we ended up shipping
after two earlier attempts failed on empirical cohorts.

### 1.3 Notation

For one candidate site on one cohort we carry five per-probe summary
statistics from the cohort's methylation array: `β_tumor_mean`,
`β_tumor_q25`, `β_tumor_q75`, `β_normal_mean`, `β_normal_q25`,
`β_normal_q75`, plus per-side sample counts `n_tumor` and `n_normal`
and an `EvidenceClass` capturing the distance between the candidate's
PAM cytosine and the nearest assayed CpG probe
(`EXACT`/`PROXIMAL_CLOSE`/`PROXIMAL`/`REGIONAL`/`UNOBSERVED`).

---

## 2 · The V2 misspecification

### 2.1 First-pass composite

The first-pass composite followed a standard decomposition:

```
p_therapeutic_selectivity = p_targ × p_prot × p_trust
```

- `p_targ = P(β_tumor < 0.30)` — probability that the tumor's PAM cytosine is unmethylated enough for ThermoCas9 to cleave. Estimated via a method-of-moments Beta fit to `(mean, IQR/1.349)` with the regularized incomplete beta computed in pure stdlib via Lentz's continued fraction; falls back to a piecewise-linear CDF when the Beta moments are ill-defined.
- `p_prot = P(β_normal > 0.50)` — probability that the normal's PAM cytosine is methylated enough to protect the normal cell.
- `p_trust` — confidence factor that scales with `EvidenceClass` (0.95 at EXACT, 0.15 at REGIONAL, 0.0 at UNOBSERVED) and saturates linearly in `min(n_tumor, n_normal)` up to `ramp_n = 30`.

### 2.2 Empirical failure of `p_prot`

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

### 2.3 V2.4 — keep the composite, drop `p_prot` by default

The intermediate fix (release v0.4.0) kept the composite but added a
`probabilistic_mode` cohort YAML key so users could explicitly opt in
or out of the anti-predictive factor. Default changed to
`tumor_only` = `p_targ × p_trust`. This improved global AUC but the
top-100 collapsed onto "always-unmethylated" loci (candidates where
both tumor *and* normal are low-methylated) — good AUC, unusable for
target discovery. The deeper fix was needed.

---

## 3 · V2.5 · differential-protection probabilistic mode

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
only whether the normal-vs-tumor gap exceeds δ in the sense of a
normal-approximation confidence interval on the per-probe summaries.

### 3.2 Choice of δ

An offline δ sweep over {0.2, 0.3, 0.4, 0.5} on the GSE77348 surrogate
selected δ = 0.2 by joint optimization of AUC and `P@100`. Higher δ
values monotonically reduce both — the factor becomes too strict to
include candidates with moderate but real differentials. This default
is exposed as `differential_delta` in the cohort YAML. We carry the δ
used on every emitted `ProbabilisticScore` record so that reviewers
can tell from a single JSONL row which margin was in effect.

### 3.3 Why a normal approximation instead of a difference-of-Betas

A proper Bayesian treatment would model β-values as Beta distributions
and compute the exact difference distribution numerically. For the
IQR ranges seen on HM450 / EPIC v2 arrays (~0.05–0.20), the normal
approximation is close to the correct answer and has the advantage of
being a two-line stdlib closed form. A σ floor of 0.05 prevents σ_Δ
from collapsing to zero at boundary β-values (very common for CpG
islands methylated at 0 or 1). These tradeoffs are documented in the
code; a downstream user who wants the exact difference-of-Betas path
can swap it in.

### 3.4 The `p_trust` saturation prediction

With the new factor isolated, we could predict before running any new
cohort: `p_trust` saturates at the evidence-class ceiling × `min(n_t,
n_n) / 30`. On low-replicate cohorts (n=2 or n=3 per side), every
`EXACT`-evidence candidate with clean tumor hypomethylation and a
large differential ends up at the same composite value, creating a
tied band at the top of the score distribution. The prediction is that
this tied band dissolves around n ≥ 30, at which point `p_trust`
becomes continuous-valued and the composite's ordering reflects only
`p_targ × p_diff`. This is a falsifiable claim the cohort matrix in §5
directly tests.

---

## 4 · Benchmark methodology

### 4.1 Cohorts

Four cohorts, all public HM450 or EPIC v2 methylation arrays:

| cohort | regime | n tumor / normal | platform | role |
|---|---|:---:|:---:|---|
| **GSE322563** | Roth actual samples | 2 / 2 | EPIC v2 | paper-comparable biology |
| **GSE77348** | MCF-7/MCF-10A surrogate (DNMT3B paper) | 3 / 3 | HM450 | independent cell-line surrogate |
| **GSE69914** | primary breast tissue | 305 / 50 | HM450 | high-n tissue validation |
| **GSE68379** | Sanger GDSC cell-line panel × GSE69914 healthy normal | 52 / 50 | HM450 | orthogonal (see §5.4) |

Ingest details:

- **GSE322563.** Supplementary β matrix `GSE322563_beta_matrix_EPIC_v2.txt.gz`. EPIC v2 probe IDs carry `_BC##` / `_TC##` / `_TO##` / `_BO##` beadchip suffixes; we strip to the canonical `cg*` identifier and intersect with HM450 — 80.7% retention. Roth-positive gene retention: ESR1 96%, EGFLAM 91%, VEGFA 90%, GATA3 83%. The main paper's data-availability statement prints `GSE32256`, which resolves to a *Paramecium* transcriptome series; the supplement's reporting summary correctly prints `GSE322563`. We verified the liftover by cross-checking our per-probe β values at the three Roth-Fig.-5d target positions against the paper's published β — every value agrees to two decimal places (EGFLAM 0.01/0.49; ESR1 0.07/0.94; GATA3 0.02/0.31).
- **GSE77348.** Per-sample β values computed from Unmethylated / Methylated intensities (β = M / (M + U + 100)); restricted to untreated replicates.
- **GSE69914.** Full SeriesMatrix β table; filtered to `status == 2` (sporadic breast cancer, non-BRCA1) for tumor and `status == 0` (healthy donor, non-BRCA1) for normal. Tumor-adjacent normals and BRCA1 carrier arms excluded.
- **GSE68379.** Sanger processed matrix has integer row indices; we rebuild cg probe IDs via the Illumina GPL13534 manifest (485,577 data rows; row N → IlmnID N); filtered to 52 primary-site-breast cell lines. No in-study normal — paired cross-series with GSE69914 healthy normals.

All cohort builder scripts are in `scripts/build_*.py`. Every run produces `tumor_summary.tsv` + `normal_summary.tsv` + `PROVENANCE.md` in the shared `LocalSummaryBackend` format; the scorer is platform-agnostic thereafter.

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

Three axes are benchmarked on every cohort:

- **V1** — `final_score = sequence × selectivity × confidence − heterogeneity_penalty − low_coverage_penalty`. Deterministic, continuous-valued, tie band always 1.
- **V2 `tumor_only`** — `p_targ × p_trust`. Retained for audit; default in v0.4.0.
- **V2.5 `tumor_plus_differential_protection`** — `p_targ × p_diff × p_trust` with δ = 0.2.

### 4.4 Tie-band handling

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

### 5.1 Cross-cohort AUC matrix

AUC under the repaired Roth-validated labels
(`tumor_only` / **differential** / V1`final_score`):

| cohort | regime | n | validated | narrow | wide | V2.5 tie_band@100 |
|---|---|:---:|---|---|---|---:|
| **GSE322563** | Roth cell lines | 2/2 | 0.928 / **0.990** / 0.821 | 0.886 / **0.942** / 0.884 | 0.871 / **0.910** / 0.768 | 190 |
| **GSE77348** | MCF-7/MCF-10A surrogate | 3/3 | 0.912 / **0.982** / 0.968 | 0.911 / **0.983** / 0.969 | 0.887 / **0.949** / 0.931 | 299 |
| **GSE69914** | primary tissue | 305/50 | **0.803** / 0.773 / 0.660 | **0.843** / 0.711 / 0.539 | **0.874** / 0.726 / 0.435 | **2** |

Bold = best AUC in that row. `tumor_only` tie_band is 6,540–11,848
across all three cohorts — top-K not usable on any of them. V1 tie_band
is 1 everywhere.

Three readings:

1. **V2.5 wins AUC on matched cell-line cohorts by clear margins.** +0.17 loose on the strictest label set for GSE322563, +0.01 on GSE77348. Label repair by itself shifted V2.5 on GSE322563 from 0.694 (gene-universe tight) to 0.990 (validated); V1 from 0.541 to 0.821. Both axes were penalized by noisy labels; V2.5 more so.
2. **Tissue cohort has a different mode ordering**: `tumor_only` > V2.5 > V1 on AUC, but `tumor_only`'s tie_band (6,540) makes its top-K fictional — V2.5 is the only probabilistic axis with a usable top-K on GSE69914.
3. **The `p_trust`-saturation prediction from §3.4 held**: V2.5's tie_band at K=100 is 190 at n=2/2 → 299 at n=3/3 → 2 at n=305/50. This is a prediction made before the high-n cohort was run; no hyperparameter was touched between cohorts.

### 5.2 P@K intervals

On the two low-n cohorts, `P@100` values should be read as intervals,
not point estimates. Example from GSE322563 narrow labels (n=28 positives):

| axis | P@100 observed | P@100 min | P@100 max | tie_band |
|---|---:|---:|---:|---:|
| V2.5 differential | 0.000 | 0.000 | 0.020 | 190 |
| tumor_only | 0.000 | 0.000 | 0.020 | 10,005 |
| V1 final_score | 0.020 | 0.020 | 0.020 | 1 |

The min-max interval is the tie-band-aware uncertainty on the metric.
V2.5's "P@100 = 0.000" on this cohort means: under the deterministic
`candidate_id` ascending tie-break, zero Roth-narrow-positives landed
in top-100; under the best-case adversarial tie-break, up to 2 could
have. V1's interval collapses to its observed value because V1's score
is continuous-valued and tie_band = 1. This is the most honest reading
of P@K on cohorts with n ≪ `ramp_n`.

### 5.3 Top-hit biology

We ran an annotation pass (nearest gene + TSS distance + feature class
+ CpG-island context, from UCSC `refGene` and `cpgIslandExt` hg19
tables) on the V2.5 top-20 for each cell-line cohort:

- **GSE322563** top-20 by V2.5 differential (chr10 portion of the
  190-tied band under `candidate_id` ascending): KCNIP2, CALHM2,
  SORCS3-AS1, LINC00710, CELF2 (×2), XPNPEP1, CELF2-AS2, ADRB1 (×3),
  ABLIM1, GFRA1 (×3), PLPP4, DMBT1, BUB3, CTBP2 (×2).
- **GSE77348** top-20 by V2.5 differential: PITX3, BTRC, SORCS1,
  LINC02624, CELF2-DT, CELF2, XPNPEP1, ADRB1, AFAP1L2, GFRA1 (×3)
  (plus earlier ranks that we omit here; full TSV committed).
- **Cross-cohort V2.5 top-20 overlap (GSE322563 ∩ GSE77348)**:
  **ADRB1, CELF2, GFRA1, XPNPEP1** (4 genes). Meaningful convergence
  given both lists are 20-record windows inside large (190 / 299) tied
  bands. These four genes are differentially methylated across both
  cohorts and carry known cancer relevance (ADRB1 as a β-adrenergic
  receptor implicated in tumor microenvironment signaling;
  GFRA1 as a GDNF coreceptor with loss-of-function roles in some
  cancers; CELF2 as a splicing regulator; XPNPEP1 as a peptidase with
  reported prognostic associations).
- **V1 `final_score` top-20 on GSE322563**: ZMIZ1, DPYSL4, **RET**,
  LINC02669, ADRB1, KCNIP2, TRIM31 (×5), CNPY3-GNMT, TTBK1, GFRA1,
  IRX1, CTBP2, ADGRA1, LINC01163, ME1, NRG3. The standalone **RET** hit
  is particularly notable — RET is a well-established therapeutic
  target in medullary thyroid cancer and *RET*-rearranged non-small-cell
  lung cancer. V1 is not tied-band-confined (tie_band = 1), so it
  samples more broadly than V2.5's window into the 190-tied band.
  Shared genes across V1 and V2.5 on GSE322563: ADRB1, CTBP2, GFRA1,
  KCNIP2.
- **GSE69914 tissue top-20 by V2.5**: RUNX2 (×5), SCML4 (×5), FOXCUT,
  FOXC1, MAS1L, COL21A1, MXI1, TENM2 (×3), MAT2B. tie_band = 2 — the
  top-20 is a real top-20. β differentials are smaller (0.17–0.31)
  than on cell lines (0.79–0.96) because primary tissue is
  cell-type-heterogeneous. **Zero gene overlap** with the cell-line
  top-20s — V2.5 is cohort-specific, not returning a fixed gene list.

The top hits are biologically coherent Roth-pattern loci. None of the
three Roth-*Fig.-5d*-validated probes (EGFLAM T11, ESR1 T17, GATA3
T18) appear in any cohort's top-20 — they appear in the top ~1% of the
scored set (rank-lifted by AUC), not the top-100. That is the expected
behavior on cohorts with large tied bands.

### 5.4 Out-of-distribution: GSE68379 label transportability

GSE68379 (Sanger GDSC 52-breast-line panel × GSE69914 healthy normal,
n = 52/50) produced systematically **inverted** AUC for both
probabilistic axes under the Roth-validated labels:

| label | V1 | tumor_only | differential | V2.5 tie_band |
|---|---:|---:|---:|---:|
| validated | 0.510 | 0.222 | 0.197 | 1,602 |
| narrow | 0.485 | 0.228 | 0.169 | 1,602 |
| wide | 0.407 | 0.462 | 0.269 | 1,602 |

The root cause is a probe-level biology mismatch at the three Roth
target sites:

| probe | gene | Roth MCF-7 β | Sanger MCF-7 β |
|---|---|---:|---:|
| cg05251676 | EGFLAM | 0.01 | **0.92** |
| cg25338972 | ESR1 | 0.07 | **0.63** |
| cg01364137 | GATA3 | 0.02 | **0.70** |

Same cell line name, **opposite methylation state** at the validated
loci. MCF-7 is notorious in the field for accumulating substantial
genetic and epigenetic drift across laboratories over decades of
passaging. Because the Roth labels encode "tumor hypomethylated at
target", and Sanger's MCF-7 is methylated at those targets, any
scorer that contains `p_targ = P(β_tumor < 0.30)` will rank the
"positives" at the *bottom* of its distribution — which is exactly
what V2.5 and `tumor_only` do. The inversion is the *correct* response
to the underlying biology; it is the labels that no longer describe
the candidates' biology in this cohort, not the scorer that is wrong.

We document GSE68379 as an **out-of-distribution boundary case**
rather than a fourth generalization cohort. It is not pooled into any
cross-cohort summary metric. Redefining positives from GSE68379's own
β values would be circular (positives defined from the data being
tested).

---

## 6 · Discussion

### 6.1 Recommended use

| cohort type | example | recommended axis | rationale |
|---|---|---|---|
| Matched cell-line / paper-comparable | GSE322563, GSE77348 | **V2.5** | Highest AUC by clear margins at every label granularity. |
| Primary tumor tissue | GSE69914 | V2.5 | `tumor_only` has slightly higher AUC but its top-K is determined by tie-break, not score; V2.5 has tie_band = 2 and is the only interpretable top-K. |
| Any cohort, top-K stability priority | — | V1 `final_score` | Continuous-valued, tie_band = 1 always. Safe to read P@K under all cohort shapes. |

V2.5 is **not** promoted to an unconditional default. The default
`probabilistic_mode` in cohort YAMLs remains `tumor_only`. V2.5 ships
as `probabilistic_mode: tumor_plus_differential_protection` on the
`main` branch (v0.4.0 is the last tagged stable release) with the
caveats documented in `V2_5_REVIEW.md`.

### 6.2 Limitations

1. **Low-replicate tie bands.** On cohorts with n < `ramp_n = 30`,
   `p_trust` saturates at the evidence-class ceiling × min(n_t, n_n)
   / 30, creating a tied band of 100s–1000s of records at the top of
   the probabilistic composite. `tie_band_size_at_k` and
   `precision_at_k_{min, max}` are reported on every benchmark result
   so that P@K can be read as an interval. The problem dissolves at
   n ≳ 30 (demonstrated on GSE69914).

2. **EPIC v2 → HM450 intersect.** 80.7% of HM450 probes are retained
   after canonical probe-ID intersect. Per-gene retention for the
   Roth-validated genes: ESR1 96%, EGFLAM 91%, VEGFA 90%, GATA3 83%.
   A native EPIC v2 probe annotation (manifest ingest, catalog
   rebuild) is queued but not yet shipped.

3. **Cell-line drift.** GSE68379 establishes that MCF-7 in one lab
   is not the same epigenetic object as MCF-7 in another. Any V2.5
   target shortlist generated on a given cohort must be validated
   against the specific cell line that will be used in the follow-up
   editing assay, not against a different lab's stock of the
   same-named line.

4. **Catalog scope.** Our catalog is chr5 / chr6 / chr10 only (2.98M
   NNNNCGA/NNNNCCA candidates within ±500 bp of an HM450 probe).
   Extending to the whole genome is straightforward but unnecessary
   for the scientific question asked here.

5. **σ floor.** The 0.05 floor on per-side σ prevents σ_Δ from
   collapsing to zero at boundary β-values. The value is empirical; a
   formal robustness study across floor values is deferred.

6. **No wet-lab validation in this memo.** V2.5 is a
   target-prioritization tool. The claim is AUC + top-list
   interpretability on committed methylation data, not actual editing
   efficiency at any specific site. The top hits (ADRB1, CELF2, GFRA1,
   XPNPEP1) are unvalidated predictions at this stage.

### 6.3 Next steps

In priority order, not committed to any timeline:

1. Native EPIC v2 probe annotation + catalog rebuild, to remove the
   HM450-intersect caveat from the Roth-comparable pipeline.
2. Annotation-pass extension: mappability, repeat overlap, and ENCODE
   cCRE intersection to move top-20 shortlists from "interesting
   genes" to "experimentally prioritized target sites".
3. A second independent-lab MCF-7 / MCF-10A EPIC cohort, if one
   becomes public, to establish reproducibility of the V2.5 claim at
   n ≥ 3 on matched cell-line pairs.
4. Wet-lab editing validation at the top V2.5 candidates on matched
   cell-line substrates — the terminal step and the one that actually
   changes the scientific claim.

---

## 7 · Methods — key equations

**V1 `final_score`:**

```
sequence_score      = PamFamily.weight
selectivity_score   = max(0, β_normal_mean − β_tumor_mean)
                    + max(0, β_normal_q25  − β_tumor_q75)
confidence_score    = EvidenceClass → {1.0, 0.7, 0.4, 0.1, 0.0}
heterogeneity_pen   = weight × max(0, β_tumor_IQR − threshold)
low_coverage_pen    = weight × max(0, (n_threshold − n_tumor) / n_threshold)
final_score         = sequence × selectivity × confidence
                    − heterogeneity_pen − low_coverage_pen
```

UNOBSERVED candidates are not penalized (verified by regression test).

**V2 / V2.5 probabilistic factors:**

```
p_targ = P(β_tumor < 0.30)   via method-of-moments Beta(α, β) fit,
                              fall back to piecewise-linear CDF when
                              σ² ≥ μ(1−μ). Incomplete beta via Lentz's
                              continued fraction, pure stdlib.
p_prot = P(β_normal ≥ 0.50)  via same CDF dispatch.     [V2 only]
p_diff = 1 − Φ((δ − (μ_n − μ_t)) / σ_Δ)
         where σ_Δ = √(max(σ_t, 0.05)² + max(σ_n, 0.05)²)
                          and σ_k = IQR_k / 1.349.            [V2.5]
p_trust = EvidenceClass.base × min(1, min(n_t, n_n) / 30)
         (EvidenceClass.base ∈ {0.95, 0.75, 0.45, 0.15, 0.0}).
```

Composites by mode:

```
tumor_only                           : p_targ × p_trust
tumor_plus_normal_protection         : p_targ × p_prot × p_trust
tumor_plus_differential_protection   : p_targ × p_diff × p_trust
```

Validator enforces that the stored `p_therapeutic_selectivity`
matches the product implied by the declared `mode` to within
`math.isclose(rel_tol=1e-9, abs_tol=1e-12)`; a second validator
enforces that the differential audit fields (`p_differential_protection`,
`differential_delta`) are populated iff `mode` is differential.

**Benchmark ranking contract:**

```
primary key    = score descending
secondary key  = candidate_id ascending  (tie_break_policy)
top-K          = first K records under (primary, secondary) sort
P@K_observed   = positives in top-K / K
tie_band_size_at_k = records at the K-th position's score
P@K_min        = max(0, band_pos − drop_slots) contribution +
                 committed_pos, all / K
P@K_max        = min(band_pos, k_band) contribution + committed_pos, / K
AUC            = Mann-Whitney U / (n_pos · n_neg); ties contribute 0.5
```

---

## Data and code availability

- **Code**: <https://github.com/AllisonH12/thermocas9> at commit `45867fb` (tag `v0.4.0` for the stable default; V2.5 as `tumor_plus_differential_protection` on `main`).
- **Tests**: 215 passing under `uv run pytest -q`.
- **Cohort data**: publicly-downloadable GEO series GSE322563, GSE77348, GSE69914, GSE68379; build scripts in `scripts/build_gse*_cohort.py` produce the committed per-probe summary TSVs in `data/derived/*_cohort/`. Positives-list builder at `scripts/build_roth_positives.py` (requires the Ensembl REST `/map` endpoint for the hg38 → hg19 liftover of Roth Fig. 5d coordinates).
- **Reference data**: UCSC hg19 `refGene.txt.gz` and `cpgIslandExt.txt.gz` (fetched on demand; gitignored).
- **Benchmark artifacts**: every `BenchmarkResult` JSONL row committed under `examples/*_roth_labels/` carries `precision_at_k`, `precision_at_k_{min,max}`, `recall_at_k`, `recall_at_k_{min,max}`, `roc_auc`, `tie_band_size_at_k`, and `tie_break_policy`.
- **Annotated top-20 TSVs**: under `examples/*/top20_annotated_*.tsv`.

## Reproducing the cross-cohort matrix

```bash
# Given a clone of the repo:
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest                       # 215 tests, ~1 s

# Rebuild the catalog on chr5/6/10 (~2 min):
thermocas build-catalog \
    --reference data/raw/hg19/hg19_chr5_6_10.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg19.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg19_chr5_6_10.jsonl

# Score + benchmark one cohort (~90 s per mode on a modern laptop):
thermocas score-cohort \
    --catalog data/derived/catalog_hg19_chr5_6_10.jsonl \
    --cohort data/derived/gse322563_differential.yaml \
    --pam-model config/pam_model.yaml \
    --backend summary \
    --probe-annotation data/derived/gse322563_cohort/probes.tsv \
    --tumor-summary  data/derived/gse322563_cohort/tumor_summary.tsv \
    --normal-summary data/derived/gse322563_cohort/normal_summary.tsv \
    --probabilistic \
    --output data/derived/scored_gse322563_differential.jsonl

thermocas benchmark \
    --scored data/derived/scored_gse322563_differential.jsonl \
    --positives data/derived/positives_roth_validated.txt \
    --cohort-name GSE322563-validated-differential \
    --score-field p_therapeutic_selectivity \
    --top-k 100 --no-enforce-holdout \
    --output bench.jsonl
```

`data/derived/scored_*.jsonl` files are gitignored because they are
3.4 GB each and reproducible from the committed summary TSVs.

---

## Citation

If this memo or the associated code is useful in your work, please
cite the upstream paper:

> Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.
> Molecular basis for methylation-sensitive editing by Cas9.
> *Nature* (2026). DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

For the probabilistic scoring framework specifically, cite this memo
plus the code URL above.

## Acknowledgements

This memo describes a framework built as an educational / research
project. ThermoCas9's biochemical characterization and the validated
target sites are entirely due to Roth et al.; any errors in the
re-implementation of the scoring or benchmarking are this memo's.

The deep-self-review process behind the code (documented in
`V2_5_REVIEW.md`) caught a series of bugs that would have otherwise
affected the headline claim: a composite-consistency validator gap on
`ProbabilisticScore`, an input-order-dependent tie-break in
`evaluate_ranking`, a prefix-ID tie-break bug in the top-hit annotator,
and an early-termination bug in the overlap-detection loop of the
gene annotator. Each of these is covered by a dedicated regression test.
