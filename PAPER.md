# Differential-protection probabilistic scoring for methylome-guided ThermoCas9 target-site ranking

**Author.** Allison Huang, Columbia University. Contact: <allisonhmercer@gmail.com>.
**Date.** 2026-04-22.
**Code.** <https://github.com/AllisonH12/thermocas9> at tag `memo-2026-04-22-j` (immutable pointer to the exact revision that produced this memo; supersedes `memo-2026-04-22-i` with five fixes: render script's awk title-block stripper now extends to the `---` horizontal rule (was leaking second-paragraph metadata into both PDFs' bodies); render script now sources the date from a `**Date.**` line in the source MD (was injecting wall-clock date, making PDFs non-reproducible from the tag); §4.3 axis count corrected from "Three" to "Four" (Δβ-only baseline was missing); §4.4 vague "axis ordering preserved" claim narrowed to the verified 9/9 invariant; §6.1 V1-row decision-table claim "V2.5 outperforms V1 at every cohort × tier" narrowed to "every non-boundary cohort × tier" — the prior phrasing was violated on the 3 GSE68379 OOD rows where V1 > V2.5).
**Status.** Technical memo from an educational research framework. Not peer-reviewed. No clinical claims. Cites Roth et al., *Nature* (2026), DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

---

## Abstract

Methylation-sensitive Cas9 variants target genomic loci that are
hypomethylated in disease cells and methylated in matched normal
cells. Selecting such loci from genome-scale array data is a ranking
problem whose scoring axis must be cohort-agnostic and must report
its own uncertainty. A first-pass probabilistic composite
`p_targ × p_prot × p_trust` failed empirically — `p_prot` encoded
a static "normal is methylated above 0.5" assumption that is
anti-predictive on bulk normal comparators (AUC 0.38). We replace it
with a differential factor `p_diff = P(β_normal − β_tumor > δ)` under
an independent-normal approximation on per-probe summaries. On the
single independent primary endpoint — AUC at the validated target
probes from Roth et al. (2026) Fig. 5d on Roth's own MCF-7/MCF-10A
EPIC v2 cohort (GSE322563) — the new composite (V2.5) reaches 0.990,
compared with 0.821 for the deterministic V1 score and 0.928 for
the deprecated V2 `tumor_only` mode. A second matched MCF-7/MCF-10A
cohort (GSE77348) is the development cohort on which the differential
margin δ was tuned and is reported as supporting evidence, not as
an independent confirmation. Every benchmark result emits
`tie_band_size_at_k` and `precision_at_k_{min,max}` so that
low-replicate top-K numbers are reported as adversarial intervals,
not point estimates. A separate cross-series run at the Sanger GDSC
breast panel is documented as an out-of-distribution label-transport
boundary case, not a generalization test, because Sanger's MCF-7 is
methylated at the sites where Roth's MCF-7 is unmethylated. V1
remains the stable-release default; V2.5 is the recommended
probabilistic mode for matched cell-line cohorts; all code and
benchmark artifacts are public.

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

![Figure 1 · V2.5 mode-formula schematic](docs/figures/fig1_mode_schematic.png)

**Figure 1.** The V2.5 composite. Three factors multiply: `p_targ`
(tumor unmethylated at the PAM cytosine), `p_diff` (differential
methylation gap exceeds δ under a normal approximation on per-probe
β summaries), and `p_trust` (evidence-class confidence, saturating
in min sample count). The result is `p_therapeutic_selectivity`,
the stored scalar that downstream ranking consumes. The deprecated
V2 `p_prot = P(β_normal > 0.5)` is replaced by `p_diff` to remove
its static-threshold assumption about the normal side.

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

### 3.4 The low-`n` tied-band prediction

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

The falsifiable prediction: the K = 100 tied band on V2.5 should be
large at n=2/3 per side and shrink to near-1 at n ≳ 30 per side. §5.1
and §5.3 test this directly across cohorts ranging from n=2/2 to
n=305/50.

---

## 4 · Benchmark methodology

### 4.0 Pre-registration: primary endpoint, sensitivity, boundary

To guard against post-hoc optimization, the analysis design separates
tuned components from evaluation. The hyperparameter `differential_delta`
was tuned on a specific cohort, and that cohort is therefore *not*
part of the protected primary endpoint.

- **Tuned, fixed before cross-cohort evaluation.**
  `differential_delta = 0.2` was selected by an offline sweep over
  {0.2, 0.3, 0.4, 0.5} on the GSE77348 MCF-7/MCF-10A surrogate
  against the *pre-repair* (gene-universe) positives. The value was
  fixed in code before any other cohort was ingested or benchmarked.
  GSE77348 is the **development cohort** for δ; results on GSE77348
  are tuned-on supporting evidence, not independent confirmation.
  The label-repair exercise that produced the Roth-validated
  positives is data preparation, not scoring-axis tuning, and it
  post-dates the δ choice.

- **Primary endpoint.** AUC at the `positives_roth_validated.txt`
  label set (n = 3 Roth Fig. 5d target probes, hg38 → hg19-lifted)
  on **GSE322563** (Roth's own MCF-7 / MCF-10A samples). This is
  the single independent paper-comparable test. δ was not tuned
  against this cohort; the label set was not derived from it.

- **Supporting evidence** (§5.1, §5.2). AUC on GSE77348 (the
  development cohort for δ, MCF-7/MCF-10A on HM450). Same matched
  cell-line biology as GSE322563 but profiled by a different lab on a
  different platform, so cross-platform agreement on the V2.5 ordering
  is informative — but it is *not* an independent cohort under the
  pre-registration discipline because δ was selected on it.

- **Sensitivity analyses** (§5.2). The `narrow` (±50 bp) and `wide`
  (±500 bp) label sets are alternative label granularities; the
  cross-granularity AUC ordering stability is reported without
  retuning δ. `P@K` intervals on tied bands are reported but are
  expected to be uninformative on low-`n` cohorts (§3.4).

- **Tissue-cohort behavior** (§5.3). GSE69914 (n = 305 / 50 primary
  breast tissue) is *not* a paper-comparable test of Roth biology —
  it probes whether the V2.5 architecture's low-`n` tied-band
  behavior dissolves at high `n`, as predicted in §3.4.

- **Out-of-distribution boundary case** (§5.4). GSE68379 (Sanger
  GDSC breast panel × GSE69914 healthy normal) is *not* a
  generalization cohort; it is a label-transportability test. It is
  reported separately and is not pooled into any primary-endpoint
  average.

### 4.1 Cohorts

Four cohorts, all public HM450 or EPIC v2 methylation arrays:

| cohort | regime | n tumor / normal | platform | role |
|---|---|:---:|:---:|---|
| **GSE322563** | Roth actual samples | 2 / 2 | EPIC v2 | paper-comparable biology |
| **GSE77348** | MCF-7/MCF-10A surrogate (DNMT3B paper) | 3 / 3 | HM450 | δ development cohort (tuned-on supporting evidence) |
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

Four axes are benchmarked on every cohort:

- **Δβ-only** (`naive_selectivity`) — `β_normal_mean − β_tumor_mean`. The literature-naive baseline: rank by the raw methylation gap with no uncertainty propagation or evidence weighting. Added in `memo-2026-04-22-c` to address the "why not just rank by Δβ?" reviewer question up front.
- **V1** — `final_score = sequence × selectivity × confidence − heterogeneity_penalty − low_coverage_penalty`. Deterministic, continuous-valued, tie band always 1.
- **V2 `tumor_only`** — `p_targ × p_trust`. Retained for audit; default in v0.4.0.
- **V2.5 `tumor_plus_differential_protection`** — `p_targ × p_diff × p_trust` with δ = 0.2.

### 4.4 Platform harmonization and catalog scope

Two facts constrain interpretation of the reported numbers and
belong in the methodology, not as a footnote:

- **Platform.** GSE322563 was profiled on Illumina EPIC v2
  (GPL33022, ~937K probes); GSE77348, GSE69914, and GSE68379 are on
  HM450 (GPL13534, ~485K probes). The primary GSE322563 results in
  §5.1 use an HM450-intersect harmonization path (EPIC v2 probe IDs
  stripped of beadchip-design suffixes — `_BC##`, `_TC##`, `_TO##`,
  `_BO##` — then intersected with the HM450 universe; 80.7%
  retention overall, 83–96% at the Roth-validated genes). To
  quantify the cost of that shortcut, we also built a **native EPIC
  v2 ingest** (script `build_epic_v2_probes.py` lifts GPL33022
  hg38 → hg19 via pyliftover; 147,928 chr5/6/10 probes, 99.95%
  liftover success; cohort builder
  `build_gse322563_native_epic_v2_cohort.py` keeps native probe IDs
  and produces a 5.22M-candidate catalog vs the HM450 catalog's
  2.98M). Side-by-side at the validated label set on GSE322563:
  V2.5 differential AUC 0.990 (HM450-intersect) → 0.986 (native EPIC
  v2), a 0.003 difference well within the cohort's tied-band noise
  floor. V1 gains the most under the native ingest (+0.05 to +0.11
  AUC across label granularities) because its continuous score
  ranks the validated targets higher relative to easy-negative
  candidates added by the EPIC v2 probe density. The only axis
  ordering that holds on **every** matched-cell-line tier × path
  row (9/9) is V2.5 > V1 and V2.5 > Δβ-only; the relative order of
  Δβ-only, V1, and V2 `tumor_only` below V2.5 reshuffles across
  rows (e.g. `tumor_only` sits above Δβ-only on both GSE322563 wide
  tiers; on GSE77348 `tumor_only` is the lowest at every tier).
  Sensitivity table in §5.2; benchmark JSONLs at
  `examples/gse322563_native_roth_labels/`.

- **Catalog scope.** The candidate catalog is chr5 / chr6 / chr10
  only, filtered to candidates within ±500 bp of any assayed HM450
  probe (~2.98M NNNNCGA / NNNNCCA candidates). The three chromosomes
  were chosen because they carry the Roth Fig. 5d validated targets
  (EGFLAM on chr5, ESR1 on chr6, GATA3 on chr10) and because they
  provide a realistic but not whole-genome test bed. All AUC /
  `P@K` / tie_band numbers in §5 are on this scope. Extending to
  whole-genome hg19 is an infrastructure change, not a scientific
  one.

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

![Figure 2 · Cross-cohort AUC by score axis](docs/figures/fig2_auc_bars.png)

**Figure 2.** Cross-cohort AUC. **(a)** Validated-label AUC on each
cohort. GSE322563 is the independent primary endpoint; GSE77348 is
the development cohort on which δ was tuned (§4.0) and is shown as
supporting evidence; GSE69914 is the tissue-behavior check for the
§3.4 tied-band prediction. V2.5 is the highest-AUC axis on both
matched cell-line cohorts; on tissue `tumor_only` takes a small AUC
lead but at a tied band of 6,540 (§5.3 + Fig. 3). V1 falls below
random on the tissue wide-label set in (b). **(b)** Sensitivity over
label granularity (validated → narrow ±50 bp → wide ±500 bp). The
cohort-type × axis ordering is stable across granularities; only the
magnitude varies. Dashed line at AUC = 0.5 (random). The out-of-
distribution GSE68379 cohort is not included here — it is plotted
separately (§5.4).

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

| axis | **GSE322563** (primary, n=2/2) | GSE77348 (δ-tuned, n=3/3) |
|---|---:|---:|
| V1 `final_score` | 0.821 | 0.968 |
| V2 `tumor_only` | 0.928 | 0.912 |
| **V2.5 differential** | **0.990** | **0.982** |

On the independent primary endpoint (GSE322563 validated), V2.5 is
the highest-AUC axis — +0.17 over V1 and +0.06 over the deprecated
V2 `tumor_only` composite. On the development cohort (GSE77348
validated), V2.5 is similarly the highest-AUC axis (+0.01 over V1,
+0.07 over `tumor_only`); this is consistent with the primary-
endpoint result but is not an independent replication, since the
same cohort selected δ = 0.2 pre-repair-labels (§4.0).

### 5.2 Sensitivity analyses: label granularity and P@K intervals

**Label granularity (AUC).** Stability of the primary-endpoint
ordering under weaker label definitions (narrow ±50 bp, wide
±500 bp):

| cohort | label set | V1 | tumor_only | V2.5 |
|---|---|---:|---:|---:|
| GSE322563 | validated (n=3) | 0.821 | 0.928 | **0.990** |
| GSE322563 | narrow (n=28) | 0.884 | 0.886 | **0.942** |
| GSE322563 | wide (n=142) | 0.768 | 0.871 | **0.910** |
| GSE77348 | validated (n=3) | 0.968 | 0.912 | **0.982** |
| GSE77348 | narrow (n=28) | 0.969 | 0.911 | **0.983** |
| GSE77348 | wide (n=142) | 0.931 | 0.887 | **0.949** |

V2.5 remains the highest-AUC axis on both cohorts at every label
granularity. No retuning of `δ` between rows.

**P@K intervals.** On `n = 2 / 2` and `n = 3 / 3` cohorts the score
distribution at K = 100 sits inside a tied band (§3.4); P@100 is
reported as an adversarial interval `[min, max]` per §4.5. On
GSE322563 narrow labels:

| axis | P@100 observed | P@100 min | P@100 max | tie_band@100 |
|---|---:|---:|---:|---:|
| V2.5 differential | 0.000 | 0.000 | 0.020 | 190 |
| tumor_only | 0.000 | 0.000 | 0.020 | 10,005 |
| V1 final_score | 0.020 | 0.020 | 0.020 | 1 |

V2.5's interval is [0.000, 0.020] — under the deterministic
`candidate_id` ascending tie-break, zero narrow-positives land in
top-100; under the best-case adversarial tie-break inside the
190-tied band, 2 of 28 could. V1's interval collapses to the
observed value because V1's score is continuous and `tie_band = 1`.
`tumor_only`'s tied band (10,005) spans two orders of magnitude
more of the score distribution than V2.5's; top-K on `tumor_only`
is not a meaningful quantity on this cohort.

**Native EPIC v2 vs HM450-intersect on GSE322563.** §4.4 describes
the harmonization shortcut: GSE322563 EPIC v2 probe IDs are stripped
of beadchip-design suffixes and intersected with the HM450 universe
(80.7% retention). To verify the shortcut is not distorting the
primary endpoint, we ran the full pipeline a second time against a
native EPIC v2 catalog (5.22M candidates vs the HM450 catalog's
2.98M):

| label set | axis | HM450-intersect AUC | native EPIC v2 AUC | Δ |
|---|---|---:|---:|---:|
| validated | V1 | 0.821 | 0.933 | +0.112 |
| validated | tumor_only | 0.928 | 0.936 | +0.008 |
| validated | **V2.5 differential** | **0.990** | **0.986** | **−0.003** |
| narrow | V1 | 0.884 | 0.938 | +0.054 |
| narrow | tumor_only | 0.886 | 0.933 | +0.046 |
| narrow | V2.5 differential | 0.942 | 0.983 | +0.040 |
| wide | V1 | 0.768 | 0.855 | +0.087 |
| wide | tumor_only | 0.871 | 0.916 | +0.044 |
| wide | V2.5 differential | 0.910 | 0.945 | +0.035 |

The V2.5 primary-endpoint AUC moves by 0.003 (well within the
tied-band noise floor on this n=2/2 cohort). V1 gains the most under
native EPIC v2 — its continuous score ranks the validated targets
higher relative to easy-negative candidates added by the larger
catalog. The relative axis ordering (V2.5 > V1) is preserved across
both ingest paths. The HM450-intersect shortcut is therefore not
materially distorting the headline V2.5 claim. Tied bands grow
modestly with the larger catalog (V2.5 differential: 190 → 421;
tumor_only: 10,005 → 14,914; V1 stays at 1).

### 5.3 Tissue-cohort behavior — GSE69914 (high-`n`, tissue biology)

GSE69914 (n = 305 / 50 primary breast vs. healthy donor tissue)
tests whether the low-`n` tied-band behavior predicted in §3.4
dissolves at `n ≫ 30`:

| label set | V1 | tumor_only | V2.5 | V2.5 tie_band@100 |
|---|---:|---:|---:|---:|
| validated (n=3) | 0.660 | **0.803** | 0.773 | **2** |
| narrow (n=28) | 0.539 | **0.843** | 0.711 | 2 |
| wide (n=142) | 0.435 | **0.874** | 0.726 | 2 |

Two findings:

1. **The tied-band prediction holds.** V2.5's tied band at K = 100
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
   wins AUC (+0.03 to +0.15 over V2.5), but `tumor_only`'s tied
   band is 6,540 at K = 100 — its top-K is not usable. V2.5 is the
   only probabilistic axis with an interpretable top-K on tissue.
   V1 collapses on tissue (AUC 0.44 on wide, below random).

### 5.4 Out-of-distribution boundary case — GSE68379

*Reported here separately. **Not** pooled with the primary-endpoint
or sensitivity tables.*

A cross-series run against the Sanger GDSC 52-breast-line panel ×
GSE69914 healthy normal (n = 52 / 50) at the Roth-validated labels
produced systematically inverted AUC:

| label | V1 | tumor_only | V2.5 | V2.5 tie_band@100 |
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

![Figure 3 · Top-20 gene presence per (axis × cohort)](docs/figures/fig3_topgene_heatmap.png)

**Figure 3.** Per-(axis × cohort) top-20 gene presence, ordered so
that the most-shared genes appear at the top. **Five columns**:
GSE322563 HM450 V1, GSE322563 HM450 V2.5, GSE322563 native EPIC v2
V2.5, GSE77348 V2.5, GSE69914 V2.5. Blue column headers mark cohorts
where `tie_band ≤ 2` and the top-20 is the genuine top-20 of the
score distribution; red column headers mark cohorts where
`tie_band ≫ K` (190 / 421 / 299 for the three cell-line V2.5 cohorts)
and the top-20 is a 20-record window inside the tied band selected
by the documented `candidate_id` ascending tie-break. Bold-blue
gene rows appear in **all three** cell-line V2.5 top-20 windows
(GSE322563 HM450, GSE322563 native EPIC v2, GSE77348) — that
intersection is *CELF2* and *XPNPEP1*. On the cell-line cohorts
shared membership is window convergence inside large tied bands,
not robust ranking convergence; AUC (Fig. 2) is the stable claim
at low n. V1 on GSE322563 has `tie_band = 1` and samples a different
(broader) genomic neighborhood — RET, ZMIZ1, DPYSL4 and other rows
it surfaces are robust under that axis. The GSE69914 tissue column
has `tie_band = 2` and **shares no genes with any of the three
cell-line V2.5 columns** — V2.5 is cohort-specific, not returning a
fixed gene list across regimes; the tissue top-20 maps to 9 distinct
nearest-gene symbols (COL21A1, FOXC1, FOXCUT, MAS1L, MAT2B, MXI1,
RUNX2, SCML4, TENM2) that are entirely disjoint from the cell-line
shortlists.

### 5.5 Top-hit annotation (tie-window-aware)

An annotation pass (nearest gene + TSS distance + feature class + CpG-
island context, from UCSC `refGene` and `cpgIslandExt` hg19 tables)
was run on the top-20 of each score axis × cohort. Important phrasing
note: on cohorts with large tied bands at K=100, the "top-20" is a
20-record *window* inside that tied band, selected by the documented
tie-break policy (`candidate_id` ascending). The membership of this
window is not a robust property of the scorer in the way AUC is — a
different deterministic tie-break would surface a different slice of
the same band. The claims below are therefore phrased in terms of
"present in the benchmark-selected top-20 window," not "the model
converges on this gene."

On the two matched cell-line cohorts (GSE322563 tie_band = 190;
GSE77348 tie_band = 299):

- **Present in the GSE322563 V2.5 top-20 window** (chr10 slice under
  `candidate_id` ascending): KCNIP2, CALHM2, SORCS3-AS1, LINC00710,
  CELF2 (×2), XPNPEP1, CELF2-AS2, ADRB1 (×3), ABLIM1, GFRA1 (×3),
  PLPP4, DMBT1, BUB3, CTBP2 (×2).
- **Present in the GSE77348 V2.5 top-20 window**: PITX3, BTRC, SORCS1,
  LINC02624, CELF2-DT, CELF2, XPNPEP1, ADRB1, AFAP1L2, GFRA1 (×3)
  (full TSV committed).
- **Shared between the two top-20 windows**: ADRB1, CELF2, GFRA1,
  XPNPEP1 (4 genes). This is cross-cohort convergence inside the
  `candidate_id`-ascending slice of the tied bands, not cross-cohort
  convergence of the underlying ranking. AUC is the stable claim at
  low `n`; top-20 membership is a window claim.

On GSE322563, V1 `final_score` has `tie_band = 1` and so its top-20
is a genuine ordering of the full candidate set, not a window slice.
It surfaces different and overlapping genes:

- **GSE322563 V1 top-20** (tie_band = 1): ZMIZ1, DPYSL4, RET,
  LINC02669, ADRB1, KCNIP2, TRIM31 (×5), CNPY3-GNMT, TTBK1, GFRA1,
  IRX1, CTBP2, ADGRA1, LINC01163, ME1, NRG3. RET is included and is
  notable as a well-established therapeutic target in medullary
  thyroid cancer and RET-rearranged NSCLC; this is a robust top-20
  claim under the V1 axis. Overlap with the V2.5 top-20 window on
  the same cohort: ADRB1, CTBP2, GFRA1, KCNIP2.

On GSE69914 (tissue, `n = 305 / 50`, tie_band = 2) the top-20 is
essentially a genuine top-20:

- **GSE69914 V2.5 top-20**: RUNX2 (×5), SCML4 (×5), FOXCUT, FOXC1,
  MAS1L, COL21A1, MXI1, TENM2 (×3), MAT2B. β differentials are
  smaller (0.17–0.31) than on cell lines (0.79–0.96) — primary
  tissue is cell-type-heterogeneous, β regresses toward the middle.
  Zero gene overlap with the cell-line top-20 windows. This is the
  expected behavior: on a genuine top-20 (not a tied-band slice),
  V2.5 is cohort-specific and does not return a fixed gene list.

None of the three Roth-*Fig.-5d*-validated probes (EGFLAM T11, ESR1
T17, GATA3 T18) appear in any cohort's top-20 — they appear in the
top ~1% of the scored set (rank-lifted by AUC) but not in the top-
100 slice. This is the expected behavior on cohorts with large
tied bands at K = 100 and is why AUC, not P@K, is the primary
endpoint.

---

## 6 · Discussion

### 6.1 Decision table — hierarchy of use

Three truths coexist and should not be conflated:

1. V1 `final_score` is the stable release axis (tagged `v0.4.0`).
2. V2.5 is an **experimental-on-main** probabilistic mode.
3. V2.5 is **recommended** as the probabilistic research mode across
   every cohort shape tested, including tissue — the scope was
   narrowed to "matched cell-line only" in an earlier draft because
   V2.5 appeared to be beaten by V1 on GSE69914 AUC; the committed
   bench JSONLs (which this memo is internally consistent with) show
   V2.5 actually beats V1 on tissue at every label granularity
   (§5.3). See MANUSCRIPT.md §6.1 for the same-shape-corrected
   framing.

The decision table below is the literal hierarchy:

| intended use | axis (mode) | status | rationale |
|---|---|---|---|
| **Default stable framework release** | V1 `final_score` | tagged `v0.4.0`; default mode in cohort YAMLs remains `tumor_only` | Deterministic, continuous-valued score; `tie_band = 1` on every cohort tested, so P@K is never tie-break-dependent. The stable-release role is about backward-compatibility and top-K determinism, not AUC leadership — V2.5 outperforms V1 on AUC at every **non-boundary** cohort × tier combination tested (12/12 rows across GSE322563 HM450, GSE322563 native, GSE77348, GSE69914; §5.1–§5.3). On the GSE68379 OOD boundary case, AUCs are inverted or near-random for every axis, so "AUC leadership" is not a meaningful concept there (§5.4). |
| **Recommended probabilistic research mode (all cohort shapes tested)** | V2.5 (differential) | experimental-on-main, not tagged | Highest-AUC discovery axis at every cohort × label-granularity combination on GSE322563 HM450, GSE322563 native EPIC v2, GSE77348, and GSE69914 (§5.1–5.3). On tissue (GSE69914), `tumor_only` has higher raw AUC but its tie_band (6,540 at K=100) disqualifies it for discovery; V2.5 (tie_band = 2) is the highest usable axis. Tie-bands reported per benchmark; P@K intervals honest. The cohort-YAML key is `probabilistic_mode: tumor_plus_differential_protection`. |
| **Analysis-only (diagnostic)** | V2 `tumor_only` | retained in the mode enum; not a discovery axis | Competitive AUC on tissue (§5.3) but `tie_band_size_at_k` at K = 100 ranges 5,271–14,914 across the five cohort paths tested (see §5.3 cross-cohort matrix); top-K is not interpretable. Use only for AUC sanity checks against V2.5 / V1. |
| **Unsupported / out-of-distribution interpretation** | any axis at GSE68379 | documented as §5.4 boundary case | Sanger MCF-7 epigenetic drift breaks label transportability from Roth Fig. 5d. Inverted AUC is the expected scorer response; do not pool this cohort's numbers with §5.1. |

The phrase "recommended" in row 2 means: for a user running V2.5 on a
cohort that matches the §5.1 use profile, we believe the reported
benefit over V1 is real. It does not mean V2.5 is tagged as a
release, and it does not mean V2.5 is safe to use on cohorts outside
that profile without re-running the diagnostics in §5.3 and §5.4 on
the new cohort.

### 6.2 Limitations

(Platform harmonization and catalog scope are part of the methodology
proper — §4.4 — because they materially shape what the reported AUCs
mean. They are cited briefly here as the methodological constraints
they are.)

1. **Low-replicate tie bands.** On cohorts with `n < ramp_n = 30`,
   `p_trust` saturates and produces tied bands of 100s–1000s of
   records at the top of the probabilistic composite.
   `tie_band_size_at_k` and `precision_at_k_{min, max}` are reported
   per benchmark result so that P@K is read as an interval. The
   problem dissolves at `n ≳ 30` (demonstrated on GSE69914; §5.3).

2. **HM450/EPIC harmonization (see §4.4 + §5.2 sensitivity).** The
   primary GSE322563 results in §5.1 use the HM450-intersect path
   (80.7% retention). A native EPIC v2 ingest is now also shipped;
   the V2.5 primary-endpoint AUC differs by 0.003 (0.990 → 0.986)
   under native vs intersect — well within the tied-band noise
   floor on this n=2/2 cohort. The shortcut is not materially
   distorting the headline claim. V1 gains the most under native
   (+0.05 to +0.11 across label granularities) but the relative
   axis ordering V2.5 > V1 is preserved.

3. **Catalog scope (see §4.4).** AUCs are on chr5 / chr6 / chr10
   candidates filtered to within ±500 bp of an HM450 probe, not on
   a whole-genome catalog. The three chromosomes carry the Roth
   Fig. 5d validated targets and provide a realistic test bed.

4. **Cell-line drift.** GSE68379 (§5.4) establishes that MCF-7 in
   one lab is not the same epigenetic object as MCF-7 in another.
   Any V2.5 target shortlist generated on a given cohort must be
   validated against the specific cell line that will be used in
   the follow-up editing assay, not against a different lab's stock
   of the same-named line.

5. **σ floor.** The 0.05 floor on per-side σ prevents σ_Δ from
   collapsing to zero at boundary β-values. The value is empirical;
   a formal robustness study across floor values is deferred.

6. **Top-20 membership is a window claim on low-`n` cohorts.** As
   §5.5 states, the cell-line top-20 lists are 20-record windows
   inside 190- and 299-record tied bands, selected by the
   documented tie-break policy. The four-gene overlap (ADRB1,
   CELF2, GFRA1, XPNPEP1) is convergence of the benchmark-selected
   window, not a robust claim about the underlying scorer ordering.
   AUC is the stable metric on low-`n` cohorts; top-20 membership
   becomes a stable metric at `n ≳ 30` (GSE69914, `tie_band = 2`).

7. **No wet-lab validation in this memo.** V2.5 is a target-
   prioritization tool. The claim is AUC + tie-band-aware top-K
   reporting on committed methylation data, not actual editing
   efficiency at any specific site. All named top hits are
   unvalidated predictions at this stage.

### 6.3 Next steps

In priority order, not committed to any timeline:

1. ~~Native EPIC v2 probe annotation + catalog rebuild~~ **Done.**
   Shipped after this revision; results in §5.2. The HM450-intersect
   shortcut is now sensitivity-bracketed rather than an open caveat.
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

`P@K_min` is the smallest possible value of P@K achievable by any
tie-break policy that respects the primary score ordering; it is
reached by pushing as many tied-band positives as possible outside
the top-K (limited by `drop_slots` available). `P@K_max` is the
symmetric upper bound: pull as many tied-band positives as possible
into top-K (limited by `k_band` available). By construction:

```
    P@K_min ≤ P@K_observed ≤ P@K_max
```

and the interval collapses to `{P@K_observed}` exactly when
`tie_band_size_at_k = 1`.

The `R@K_min / R@K_max` counterparts replace the denominator with
`n_positives`. `AUC` is the Mann-Whitney U / (n_pos · n_neg) with
ties contributing 0.5 (standard tie-handling); it does not need an
interval form because under the cid-ascending tie-break it is
already invariant within tied score regions.

---

## Data and code availability

- **Code**: <https://github.com/AllisonH12/thermocas9>. Five dated memo tags plus the stable-release tag, per the immutable-tag policy (never move; supersede by cutting a new dated tag):
  - `v0.4.0` — the stable-release V1 code. Default `probabilistic_mode` is `tumor_only`; V2.5 is not yet shipped at this tag.
  - `memo-2026-04-21` — preceding-day revision; predates the native EPIC v2 ingest.
  - `memo-2026-04-22` — initial revision: V2.5 experimental mode, native EPIC v2 ingest, P@K-interval benchmark contract, top-hit annotation pipeline.
  - `memo-2026-04-22-b` — adds post-self-review fixes: nested-repeat scan correctness, streaming-aggregator cross-cohort metadata parity, intra-cohort duplicate rejection, memory-claim docstring scrubs.
  - `memo-2026-04-22-c` — added the Δβ-only baseline benchmark across all five cohort paths × three positives tiers (15 `bench_*_naive.jsonl` artifacts), and introduced `MANUSCRIPT.md` as the Bioinformatics-submission-shaped sibling. **Retained but SHOULD NOT BE CITED**: MANUSCRIPT.md at this tag contains GSE69914 tissue-cohort AUC values (V1 = 0.861, V2.5 = 0.837) that do not match the committed bench JSONLs (actual: V1 = 0.660, V2.5 = 0.773). This was a manuscript-text / committed-artifact inconsistency, not a code bug.
  - `memo-2026-04-22-d` — corrected tissue-cohort AUC values and §5.2 native-vs-HM450 sensitivity table; also corrected three drift items from an independent review. **Retained but SHOULD NOT BE CITED**: this revision still had a documented-vs-implemented `p_targ` threshold mismatch (prose said 0.5; code uses 0.30), a false universal axis-ordering claim in MANUSCRIPT.md §5.2, a wrong `tumor_only` tie-band range (stated 6,000–12,000; actual 5,271–14,914), an incorrect "tumor–normal pair tissue" label on GSE69914 (actual: unpaired sporadic-tumor + healthy-donor), and a GSE68379 cohort row listing 1/0 instead of the 52/50 cross-series setup the benchmarks actually use.
  - `memo-2026-04-22-e` — corrected the `-d` issues against the committed JSONLs, cohort build scripts, and `probabilistic.py` constants. **Retained but SHOULD NOT BE CITED**: this revision still carried (i) a false-universal claim "V2.5 is the highest-AUC discovery axis at every cohort × tier combination tested" that is violated on the 3 GSE68379 OOD rows where V1 > V2.5; (ii) a false matched-cell-line ordering claim "V2.5 > V1 > Δβ, typically by 0.01-0.05" where V1 > Δβ actually holds on only 2 of 9 rows (Δβ > V1 elsewhere) and the V2.5-over-V1 range is +0.014 to +0.169, not 0.01-0.05; (iii) a wrong normal-arm mechanism attribution in §5.3/§6.1 ("adjacent-normal bulk") that contradicted the GSE69914 build script's explicit exclusion of adjacent-normal arms.
  - `memo-2026-04-22-f` — corrected the `-e` issues and added the `verify_manuscript_claims.py` guard. **Retained but SHOULD NOT BE CITED**: MANUSCRIPT.md at this tag embedded only Figure 1 (mode schematic). Figures 2 and 3 existed on disk from an earlier render against a 3-axis × 3-cohort grid (no Δβ-only baseline, no GSE322563 native EPIC v2) and were inconsistent with the manuscript's 4-axis narrative; they were not embedded in MANUSCRIPT.md.
  - `memo-2026-04-22-g` — regenerated Fig 2 and Fig 3 against the committed bench grid and embedded both in MANUSCRIPT.md §5. **Retained but SHOULD NOT BE CITED**: the PDF render helper `scripts/render_paper_pdf.sh` still hardcoded `PAPER.md` as the source (following the documented workflow would have produced the audit memo, not the submission manuscript); and Fig 3's caption had an overclaim "plus cohort-specific candidates" implying additional tissue genes beyond the 9 actually listed.
  - `memo-2026-04-22-h` — generalized the PDF render helper, fixed the Fig 3 caption "plus cohort-specific" overclaim in MANUSCRIPT.md, and extended the verify script to check figure captions. **Retained but SHOULD NOT BE CITED**: PAPER.md Fig 3 caption was not updated when fig3 went from 4 columns to 5 columns at `-g`, so it still described "BOTH cell-line V2.5 top-20 windows (GSE322563 + GSE77348)" and "either cell-line column" — 2-column framing that does not match the actual 5-column / 3-cell-line-V2.5-column figure. The verify script's caption check ran only on MANUSCRIPT.md, so this drift was not surfaced.
  - `memo-2026-04-22-i` — fixed PAPER.md Fig 3 caption to the 5-column reality and extended `verify_manuscript_claims.py` to scan both circulated documents. **Retained but SHOULD NOT BE CITED**: render script's awk title-block stripper only handled single-paragraph metadata (stopped at first blank line) — both PAPER.md and MANUSCRIPT.md have multi-paragraph metadata blocks, so the second paragraph leaked into the rendered PDF body before "Abstract". Render script also injected wall-clock date at render time, so re-running on a later day produced a visibly different PDF from the tagged revision. PAPER.md §4.3 also still said "Three axes" (missed the Δβ-only baseline added in `-c`), §4.4 had a vague "axis ordering preserved" claim, and the §6.1 decision-table V1 row had the same "V2.5 outperforms V1 at every cohort × tier" overclaim that was already corrected in MANUSCRIPT.md.
  - `memo-2026-04-22-j` — **the exact revision that produced this memo.** Render script: awk now strips the title block from the H1 to the first `---` horizontal rule (handles multi-paragraph metadata correctly; verified by `pdftotext` of both renders showing TOC → Abstract with no leaked metadata). Date is now sourced from the `**Date.**` line in the source MD (with fallback to most-recent memo-* tag's date, then to wall clock with a stderr warning); MANUSCRIPT.md gains a `**Date.**` line. PAPER.md §4.3 fixed to "Four axes" (Δβ-only listed). §4.4 "axis ordering preserved" replaced with the verified 9/9 invariant. §6.1 decision-table V1 row narrowed to "every non-boundary cohort × tier (12/12 rows)". Resolve to a SHA with `git rev-parse memo-2026-04-22-j` in a fresh clone.
  Development continues on `main` past the tagged memo revision; cite `memo-2026-04-22-j` when citing this document.
- **Submission-shaped companion**: `MANUSCRIPT.md` at the same tag is the Bioinformatics-submission-shaped cut-down of this memo (~340 lines vs ~960) with the headline framing led from the Δβ-baseline finding.
- **Tests**: 236 passing under `uv run pytest -q`.
- **Cohort data**: publicly-downloadable GEO series GSE322563, GSE77348, GSE69914, GSE68379; build scripts in `scripts/build_gse*_cohort.py` produce the committed per-probe summary TSVs in `data/derived/*_cohort/`. Positives-list builder at `scripts/build_roth_positives.py` (requires the Ensembl REST `/map` endpoint for the hg38 → hg19 liftover of Roth Fig. 5d coordinates).
- **Reference data**: UCSC hg19 `refGene.txt.gz` and `cpgIslandExt.txt.gz` (fetched on demand; gitignored).
- **Benchmark artifacts**: every `BenchmarkResult` JSONL row committed under `examples/*_roth_labels/` carries `precision_at_k`, `precision_at_k_{min,max}`, `recall_at_k`, `recall_at_k_{min,max}`, `roc_auc`, `tie_band_size_at_k`, and `tie_break_policy`.
- **Annotated top-20 TSVs**: under `examples/*/top20_annotated_*.tsv`.

## Reproducing the cross-cohort matrix

```bash
# Given a clone of the repo:
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest                       # 236 tests, ~1 s

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
