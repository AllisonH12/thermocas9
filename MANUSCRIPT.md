# thermocas: probabilistic ranking of methylation-protected ThermoCas9 target sites with tie-aware benchmarking

**Allison Huang** · Columbia University · <allisonhmercer@gmail.com>

**Date.** 2026-04-22.

Code at <https://github.com/AllisonH12/thermocas9>, immutable tag `memo-2026-04-22-be`.

---

## Abstract

**Motivation.** Methylation of the fifth PAM cytosine protects genomic loci from ThermoCas9 cleavage (Roth et al., *Nature* 2026), enabling tumor-selective targeting at loci that diverge in methylation between cancer and matched normal tissue. Translating this mechanism into a target-discovery shortlist requires a scorer that is cohort-agnostic, honest about evidence quality, and explicit about ties in its output. No such scorer was publicly available. The obvious baseline — rank candidates by the raw methylation gap Δβ = β_normal − β_tumor — was not benchmarked against any proposed alternative in the open literature.

**Results.** We present `thermocas`, an open probabilistic scoring framework built around a *compositional skeleton* `p_targ × (gap factor) × p_trust` paired with a tie-band-aware benchmarking contract (`precision_at_k_{min, max}`, `tie_band_size_at_k`, mid-rank Mann-Whitney AUC on every `BenchmarkResult`). The gap-factor slot admits several instances; two ship in this tag as selectable `probabilistic_mode` enum values: **V2.5-diff** (`tumor_plus_differential_protection`, `p_diff = P(β_normal − β_tumor > δ)` under an independent-normal approximation on per-probe β summaries with σ_floor = 0.05) and **V2.5-sigmoid** (`tumor_plus_gap_sigmoid`, a fixed-bandwidth logistic `sigmoid((β_n − β_t − δ) / σ_fixed)` with σ_fixed ≈ 0.0707). We benchmark across four public methylation cohorts — GSE322563 (Roth MCF-7/MCF-10A, primary endpoint; both HM450-intersect and native EPIC v2 ingest paths), GSE77348 (δ-development cohort, n=3/3), GSE69914 (n=305 primary tumor + 50 healthy-donor tissue, unpaired), GSE68379 (Sanger panel, OOD boundary) — against five baselines: raw Δβ-only, an uncertainty-aware Δβ_z, the deterministic V1 score, the V2 `tumor_only` diagnostic, and a sample-level limma-style moderated-t DMR baseline (Smyth 2004 empirical-Bayes variance prior, probe-level first then candidate-mapped). **On matched cell-line cohorts, Δβ-only is already strong** (validated AUC 0.96–0.97) and V2.5 variants add small but consistent margins (+0.010 to +0.080 across nine tier × path rows). **A 1,000,000-draw permutation null on the n=3 primary endpoint** places V2.5-diff at one-sided random-triple `p = 4 × 10⁻⁶ / 1.5 × 10⁻⁵ / 2.8 × 10⁻⁵` across the three primary cohorts — lower random-null *p*-values than Δβ-only's `8.6 × 10⁻⁵ / 3.1 × 10⁻⁴ / 1.1 × 10⁻⁴` on the same cohorts (against-random check, not a formal pairwise superiority test; the paired V2.5-vs-Δβ comparison floors at sign-flip `p = 0.125` with `n_pos = 3` and is reported separately in PAPER.md §5.1.3). **Because the validated set contains only three Roth Fig. 5d targets, AUC is interpreted as a rank-lift summary rather than an inferential discovery-performance estimate.** **The headline tissue ranking stress-test result (under transported Roth labels, not tissue-specific ThermoCas9 validation).** A probe-level limma-style moderated-t baseline (Smyth 2004 empirical-Bayes, sample-level β matrices, then candidate-mapped) is competitive on cell lines (0.959 / 0.991 / 0.962) but **collapses to near-random on bulk tissue (GSE69914, AUC 0.573)**, while V2.5-sigmoid **recovers tissue to AUC 0.862** — a +0.29 swing on the same probe-level statistical inputs, confirming that the `p_targ × p_trust` wrapping is the load-bearing contribution, not the gap factor's specific shape. *The tissue stress-test gain is per-positive heterogeneous: GATA3 remains strong under feature-matched controls, while ESR1 is matched-near-random under V2.5-sigmoid (PAPER.md §5.9).* **A frozen whole-genome panel** (HM450 SHA256 `d20661c5…`, 19.8M candidates; EPIC v2 SHA256 `39df8f0f…`, 35.4M candidates; PAPER.md §5.2.2) further shows V2.5-sigmoid matches V2.5-diff's AUC on matched cell-line cohorts within 0.002 but strictly beats V2.5-diff on top-K usability — `tie_band@100 = 1` vs 421 / 1,127 / 1,493 for V2.5-diff under the WG denominator — and improves tissue AUC over V2.5-diff by +0.05 to +0.08 across a bandwidth-robust sweep. The shipped recommendation is: **V1 as overall package default** (backward compatibility + `tie_band = 1`); **V2.5-sigmoid as recommended probabilistic prioritization axis** for new non-boundary runs (hypothesis generation; not a discovery claim absent prospective wet-lab validation), **selected by post-repair sensitivity and whole-genome stress testing — its prospective ranking utility remains to be validated on newly generated labels or wet-lab follow-up**; V2.5-diff retained for audit/backward-compatibility and AUC parity; `tumor_only` as `probabilistic_mode`-enum default (opt-in-required for V2.5-diff/V2.5-sigmoid via explicit cohort-YAML configuration). The package also adds a k-way-merge pan-cancer aggregator for genome-scale atlas builds and a per-candidate annotation pipeline (nearest gene, CpG-island context, RepeatMasker overlap, ENCODE DNase-HS cluster breadth) with a Markdown shortlist aimed at experimental collaborators.

**Availability and implementation.** `thermocas` is open source at <https://github.com/AllisonH12/thermocas9>, tagged `memo-2026-04-22-be` for the version evaluated here. 245 unit tests pass under `uv run pytest -q`. Python 3.11+, BSD-3.

**Contact.** <allisonhmercer@gmail.com>

---

## 1 Introduction

Roth et al. (2026) characterized *Geobacillus thermodenitrificans* T12 Cas9 (ThermoCas9) biochemically and structurally. The central result is that cleavage efficiency is governed by the methylation state of the fifth PAM cytosine — 5-methylcytosine at that position abolishes binding, while an unmethylated PAM cleaves as usual. The practical consequence: ThermoCas9 permits selective editing at loci that are unmethylated in disease cells and relatively methylated or differentially protected in matched normal cells. Roth demonstrated this with three breast-cancer targets in MCF-7 vs MCF-10A on the Illumina MethylationEPIC v2 platform.

Turning that mechanism into a genome-scale target-discovery shortlist is a ranking problem with three awkward properties. First, per-probe methylation summaries are noisy (arrays commonly report n = 2–3 replicates per side for cell-line cohorts) and the scorer must carry its own uncertainty rather than emit a scalar that hides it. Second, the "normal" arm is inconsistent: adjacent-normal tissue in bulk methylation studies is not the same object as a matched cell-line pair, and any assumption that baked in a fixed β_normal threshold would fail across cohort types. Third, on low-replicate cohorts the top-K of any probabilistic composite can sit inside a large tied band; reporting a single Precision@K value without a tie interval misleads a reader about the score's discriminative power.

No public tool addressed these three problems together when we started. The closest prior art — generic CRISPR guide-scoring and Δβ ranking on methylation arrays — handles neither uncertainty propagation nor tie-band reporting. `thermocas` is the first open framework to target all three.

**Contributions.**

1. The V2.5 probabilistic composite, where a differential-protection factor replaces a failed threshold-based factor. The decomposition is `p_targ × p_diff × p_trust`, each factor has a closed-form derivation, and the final score is a bounded probability-scale ranking composite — not a calibrated posterior probability of editing efficacy and not a hypothesis test (see Methods §2 and PAPER.md §3.3).
2. A benchmarking protocol that reports Precision@K and Recall@K as adversarial intervals (`_min`, `_max`) together with the tie-band size at K, so reviewers are not misled by arbitrary tie-break choices on low-replicate cohorts.
3. A benchmark against three baselines — a Δβ-only ranker, V1 deterministic scoring, and a rejected V2 formulation — across four public cohorts with positives tiered from Roth's own validated targets to a wide genomic window.
4. A k-way-merge pan-cancer aggregator with O(N_cohorts + N_unique_candidate_ids) memory and cross-cohort metadata-parity enforcement, for genome-scale atlas builds.
5. An experiment-facing shortlist pipeline: TSV + Markdown companion that attach nearest gene, CpG-island context, RepeatMasker overlap, and ENCODE DNase-HS breadth, with rule-based triage flags (STRONG / CAUTION / NOTE) deterministically derived from the annotated fields.

**Scope.** This paper reports methods development and benchmarking on public methylation cohorts. It does not include prospective wet-lab validation of predicted target sites; the framework is an open educational research tool, not a clinical decision-support system. Per-site p-values are not reported because `p_observation_trustworthy` saturates to a discrete per-`EvidenceClass` value at high replicate counts; the composite provides a ranking axis, not a hypothesis test.

---

## 2 Methods — Formulation

For each candidate site on each cohort, we carry six per-probe β-value summaries (`β_tumor_mean`, `β_tumor_q25`, `β_tumor_q75`, `β_normal_mean`, `β_normal_q25`, `β_normal_q75`), sample counts `n_tumor` and `n_normal`, and an `EvidenceClass` capturing the distance between the PAM cytosine and the nearest assayed CpG (`EXACT` / `PROXIMAL_CLOSE` / `PROXIMAL` / `REGIONAL` / `UNOBSERVED`).

### 2.1 V2.5 generation — two shipped composites

The *compositional skeleton* is
`probability-scale selectivity score = p_targ × (gap factor) × p_trust`
(stored in the `p_therapeutic_selectivity` field on every
`ProbabilisticScore` record). The gap-factor slot admits two shipped
instances:

```
V2.5-diff    (tumor_plus_differential_protection)  — p_targ × p_diff × p_trust
V2.5-sigmoid (tumor_plus_gap_sigmoid)              — p_targ × p_gap_sigmoid × p_trust
             — RECOMMENDED probabilistic prioritization axis (§5.2.2 WG panel)
```

Both ship as selectable `probabilistic_mode` enum values in this
tag. V2.5-diff is retained for backward compatibility and AUC parity
on cell-line cohorts; V2.5-sigmoid is the recommended axis on every
tested non-boundary cohort shape per the PAPER.md §5.2.2 WG panel.

**Calibration scope.** `p_therapeutic_selectivity` is a *bounded
probability-scale ranking score in [0, 1]*, not a calibrated probability
of editing success at the candidate site. The composite is not
validated against any prospective editing readout (no expected
calibration error, no reliability diagram, no Platt/isotonic
post-calibration). Read it as a probability-scale ranking axis, not as
"`0.9` means 90% of `0.9`-scored sites edit successfully". The bounded
scale buys multiplicative composability with downstream probabilistic
inputs (target-mutation models, gRNA off-target probabilities,
delivery-efficiency priors), which the V1 weighted-sum heuristic and
the raw-Δβ baseline both lack — see PAPER.md §3.1, §3.3.

**Shared factors.**

- **`p_targ = P(β_tumor < τ_u)`** with `τ_u = 0.30` (the `DEFAULT_UNMETHYLATED_THRESHOLD` constant in `src/thermocas/probabilistic.py`) — the probability that the tumor arm is unmethylated at the PAM cytosine. Estimated by the same IQR-normal approximation used below. The 0.30 threshold is conservative relative to the methylation-convention "0.5 split between hypo- and hyper-methylated": we require the tumor to be *confidently* hypomethylated, not just below the mid-point, to avoid scoring-inflation at borderline probes. Reported benchmark results assume `τ_u = 0.30`; a downstream user who wants a different operational cutoff can pass `unmethylated_threshold` to `p_targetable_tumor(...)`.
- **`p_trust = EvidenceClass.base × min(1, min(n_t, n_n) / 30)`** — evidence-class confidence with a linear ramp up to n = 30, saturating to a per-class constant (`EXACT`: 0.95, `PROXIMAL_CLOSE`: 0.75, etc.).

**V2.5-diff gap factor.**

- **`p_diff(δ) = P(β_normal − β_tumor > δ)`** — the probability that the normal-minus-tumor methylation gap exceeds a configurable threshold δ (default 0.2). Estimated by:
  ```
  σ_k  ≈ IQR_k / 1.349                (k ∈ {tumor, normal}; floor at 0.05)
  σ_Δ² = max(σ_t, floor)² + max(σ_n, floor)²
  z    = (δ − (μ_n − μ_t)) / σ_Δ
  p_diff = 1 − Φ(z)
  ```
  The σ_floor prevents σ_Δ from collapsing to zero at boundary β-values (CpGs methylated at 0 or 1, which are common at islands).

**V2.5-sigmoid gap factor.**

- **`p_gap_sigmoid(δ, σ_fixed) = sigmoid((β_normal − β_tumor − δ) / σ_fixed)`** — a fixed-bandwidth logistic transform of the same gap, with δ = 0.2 and σ_fixed ≈ 0.0707 (≈ √2 · σ_floor; strictly positive, enforced at the `CohortConfig` boundary). Motivated by the PAPER.md §3.5 binding-rate finding (σ_floor binds on 99.5–99.9% of records in low-`n` matched cell-line cohorts, so `p_diff` collapses to a fixed-bandwidth response in practice there); the V2.5-sigmoid formulation operationalizes this directly and avoids the tie-band blow-up V2.5-diff experiences under whole-genome denominators on low-`n` cell-line cohorts (PAPER.md §5.2.2).

Figure 1 summarizes V2.5-diff (which ships as the historical composite and is used as the axis under test in §5.1; V2.5-sigmoid swaps `p_diff` for `p_gap_sigmoid` and is benchmarked alongside it in §5.2 and PAPER.md §5.2.2).

![V2.5 mode-formula schematic](docs/figures/fig1_mode_schematic.png)

**Figure 1.** The V2.5-diff composite. Three factors multiply: `p_targ` (tumor unmethylated at the PAM cytosine), `p_diff` (differential methylation gap exceeds δ under a normal approximation on per-probe β summaries), and `p_trust` (evidence-class confidence, saturating in min sample count). The result is the *probability-scale selectivity score* (stored on every `ProbabilisticScore` record as the `p_therapeutic_selectivity` field; downstream ranking consumes that scalar directly). V2.5-sigmoid replaces `p_diff` with `p_gap_sigmoid = sigmoid((β_n − β_t − δ) / σ_fixed)` inside the same skeleton.

### 2.2 Baselines

| name | score | rationale |
|---|---|---|
| **Δβ-only** (naive) | `β_normal_mean − β_tumor_mean` | Literature-naive: rank by the raw methylation gap. No uncertainty propagation, no evidence weighting. The "obvious thing." |
| **V1 `final_score`** | weighted sum of sequence, selectivity, confidence, heterogeneity, and coverage components (continuous-valued) | The framework's stable-release axis, retained for top-K interpretability. |
| **V2 `tumor_only`** | `p_targ × p_trust` (no protection factor) | Intermediate formulation after we dropped V2's `p_prot` factor. Retained as a diagnostic; not a prioritization axis. |
| **V2** (deprecated) | `p_targ × p_prot × p_trust` with `p_prot = P(β_normal > 0.5)` | First-pass composite with a static threshold on the normal side. Empirically anti-predictive (AUC 0.38 on tissue cohort); rejected. Not carried forward as a recommended axis but reported in the comparison table for completeness. |

### 2.3 Why a normal approximation

A proper Bayesian treatment would model each β-value as a Beta distribution and compute the exact difference distribution numerically. For the IQR ranges observed on HM450 and EPIC v2 arrays (≈0.05–0.20), the independent-normal approximation is close to the correct answer and has the advantage of being a two-line closed form on summary statistics. A downstream user who wants the exact difference-of-Betas path can swap it in — the relevant code is isolated behind a single function.

### 2.4 Low-`n` tied-band prediction

A consequence of the composite structure: on low-replicate cohorts (`n ≪ 30`), `p_trust` scales every same-`EvidenceClass` candidate by the same uniform factor `(n / 30)`. Per-side IQRs are wide or zero, pushing `p_targ` and `p_diff` toward saturation for any candidate with clean tumor hypomethylation and a large differential. The composite is effectively `1 × 1 × p_trust` for that entire subset, producing a large tied band at the top of the score distribution. At `n ≳ 30`, `p_trust` is a per-class constant and the continuous variation in `p_targ × p_diff` drives the ordering. The prediction is falsifiable: the tied band at K = 100 should be large at n = 2/3 per side and shrink to near-1 at n ≳ 30 per side. §4.1 and §4.3 test this directly.

---

## 3 Methods — Benchmarking protocol

### 3.1 Tie-band-aware Precision@K and Recall@K

For score `s` and candidate_id `c`, the benchmarking harness sorts by `(−s, c)` — score descending, `candidate_id` ascending as tie-break. This defines a deterministic top-K window. For any K:

- **`tie_band_size_at_k`**: the number of candidates sharing the score of the K-th record. If K = 100 and the 100th, 101st, 102nd records all score `0.0633`, the band size is at least 3. We compute it by scanning forward until the score changes.
- **`precision_at_k_min`**: Precision@K under the worst-case tie-break (positives pushed out of the window).
- **`precision_at_k_max`**: Precision@K under the best-case tie-break (positives pulled into the window).
- **`recall_at_k_{min, max}`**: analogous.

These four quantities together tell a reader whether the observed P@K is a meaningful discrimination or an artifact of where the deterministic tie-break happens to fall.

### 3.2 ROC-AUC by Mann–Whitney U

AUC is computed as `P(score(positive) > score(negative))` via Mann–Whitney U with 0.5 contribution for ties (standard tie handling). Returns `None` when either positives or negatives are empty.

### 3.3 Held-out chromosomes

All benchmarks report `held_out_chromosomes`. When enforced (`--enforce-holdout`, default), only candidates on held-out chromosomes contribute to the evaluation, so no candidate that informed δ selection appears in the test set. The primary-endpoint cohort (GSE322563) is independent of δ tuning (done on GSE77348), so held-out enforcement is less load-bearing there, but the field is reported on every BenchmarkResult for auditability.

---

## 4 Methods — Cohorts and positives

### 4.1 Cohorts

| cohort | platform | n_tumor / n_normal | role | source |
|---|---|---:|---|---|
| **GSE322563** | EPIC v2 | 2 / 2 | Independent primary endpoint (Roth's own MCF-7/MCF-10A samples) | GEO |
| GSE77348 | HM450 | 3 / 3 | Development cohort — δ tuned here (tumor_plus_differential_protection, δ-sweep over {0.2, 0.3, 0.4, 0.5}) | GEO |
| GSE69914 | HM450 | 305 / 50 | Sporadic breast tumor vs healthy-donor breast tissue — **unpaired**; adjacent-normal and BRCA1-carrier arms excluded by the build script. Tests predicted `p_trust` desaturation at high n. | GEO |
| GSE68379 | HM450 | 52 / 50 (cross-series) | OOD boundary case: 52 Sanger GDSC1000 breast cell lines paired with 50 external healthy-donor normals from GSE69914 status=0. Cross-series lab/scanner batch effects are a known caveat. Documented as a boundary, not a generalization test. | GEO (tumor) + GEO (normal) |

### 4.2 Positives (label tiers)

Ground-truth positives are derived from Roth et al. Fig. 5d — the three validated target probes. We define three tiers to probe label-granularity sensitivity:

- **validated** (n = 3): exact Roth Fig. 5d probes, hg38 → hg19 lifted via Ensembl REST.
- **narrow** (n = 28): all candidate PAMs within ±50 bp of a validated probe.
- **wide** (n = 142): all candidate PAMs within ±500 bp of a validated probe.

Native EPIC v2 positives lists for GSE322563 use the native EPIC v2 probe annotation (§5.2); HM450-intersect positives use `cg*` ID-stripping.

### 4.3 Platform harmonization

Our first-draft pipeline harmonized GSE322563 EPIC v2 probe IDs to HM450 by stripping the `_BC##` / `_TC##` / `_TO##` / `_BO##` beadchip-design suffix and intersecting with the HM450 universe (~80.7 % retention). We subsequently re-ran the analysis against the native EPIC v2 probe set (147,928 probes on chr5/6/10, hg38 → hg19 lifted from GPL33022). Headline AUCs are reported under both paths as a sensitivity analysis (§5.2).

Catalog scope is chr5, chr6, and chr10 — the three chromosomes carrying Roth's Fig. 5d validated targets. Candidates are restricted to within ±500 bp of a methylation probe to avoid extrapolation beyond array coverage.

---

## 5 Results

![V2.5-sigmoid vs V2.5-diff vs limma-style on WG validated AUC and tie_band@100](docs/figures/fig2_auc_bars.png)

**Figure 2 — Final-method summary.** Whole-genome validated-label performance (n_pos = 3 Roth Fig. 5d targets) for the recommended probabilistic prioritization axis (V2.5-sigmoid) against the V2.5-diff predecessor and the limma-style moderated-`t` DMR baseline, on the four primary cohorts. **(a)** WG AUC: V2.5-sigmoid matches V2.5-diff on the three matched cell-line cohorts within 0.001 (0.989 vs 0.988, 0.998 vs 0.998, 0.982 vs 0.981) and beats it on tissue (0.862 vs 0.778); the limma-style baseline is competitive on cell lines (0.959 / 0.991 / 0.962) but **collapses to near-random on tissue** (0.573). **(b)** WG `tie_band@100` (log scale): V2.5-sigmoid eliminates the **421–1,493-record tied bands** that make V2.5-diff's whole-genome top-K unactionable on low-`n` matched cell-line cohorts, holding `tie_band@100 = 1` on every matched cell-line cohort and `= 6` on tissue. The limma-style baseline sits at intermediate tie-band values (39–115). The historical 4-axis × 3-tier sensitivity figure (Δβ-only, V1, V2 `tumor_only`, V2.5-diff) is retained as supplementary `docs/figures/fig2_supp_historical_sensitivity.png` for the V2.5-diff-era audit trail. Numbers committed in `examples/genome_wide_panel.tsv` (V2.5-diff, V2.5-sigmoid) and `examples/limma_cross_cohort_panel.md` (limma-style) at this tag.

### 5.1 Primary endpoint — matched cell-line AUC at validated Roth probes

Held-out primary endpoint: AUC at `positives_roth_validated.txt` on GSE322563, using HM450-intersected probes. δ was not tuned on this cohort; the label set was not derived from it.

| axis | **GSE322563** (primary, n = 2/2) | GSE77348 (δ-tuned, n = 3/3) |
|---|---:|---:|
| Δβ-only (naive) | 0.974 | 0.972 |
| V1 `final_score` | 0.821 | 0.968 |
| V2 `tumor_only` | 0.928 | 0.912 |
| **V2.5 differential** | **0.990** | **0.982** |

On the independent primary endpoint, **V2.5-diff** is the highest-AUC shipped axis (V2.5-sigmoid matches it within 0.001 — PAPER.md §5.2.2 WG panel — and is the recommended prioritization axis once top-K usability is also considered). The margin over the next-best baseline is narrow: +0.02 over Δβ-only, +0.06 over V2 `tumor_only`, +0.17 over V1. **Δβ-only is an unexpectedly strong baseline on matched cell-line cohorts** (0.974 on GSE322563 validated, 0.972 on GSE77348 validated). This is an honest finding the framework's complexity has to justify: the V2.5 variants' AUC advantage on matched cell-line data is small. The V2.5 generation's structural advantages — a probability-scale interpretation, tie-band-honest top-K reporting, and principled uncertainty propagation via `p_trust` — are what earn the complexity on low-replicate cohorts; the AUC gap over Δβ widens on tissue (§5.3), and V2.5-sigmoid further widens the tissue gap over V2.5-diff by +0.05 to +0.08 AUC (PAPER.md §5.2.2).

### 5.2 Sensitivity — label granularity, P@K intervals, platform path

**Label granularity.** Ordering is stable under weaker label definitions:

| cohort | tier (n_pos) | Δβ-only | V1 | V2 tumor_only | V2.5 |
|---|---|---:|---:|---:|---:|
| GSE322563 | validated (3) | 0.974 | 0.821 | 0.928 | **0.990** |
| GSE322563 | narrow (28) | 0.912 | 0.884 | 0.886 | **0.942** |
| GSE322563 | wide (142) | 0.844 | 0.768 | 0.871 | **0.910** |
| GSE77348 | validated (3) | 0.972 | 0.968 | 0.912 | **0.982** |
| GSE77348 | narrow (28) | 0.964 | 0.969 | 0.911 | **0.983** |
| GSE77348 | wide (142) | 0.910 | 0.931 | 0.887 | **0.949** |

V2.5-diff (the historically-shipped instance used as the axis under test in §5.1) is the highest-AUC shipped axis at every tier on both matched cell-line cohorts (HM450-intersect path); V2.5-sigmoid matches V2.5-diff on AUC within 0.002 (PAPER.md §5.2.2 WG panel) and is the recommended axis on top-K-usability grounds. The margin over Δβ-only ranges from **+0.010 to +0.066**; the margin over V1 ranges from **+0.014 to +0.169**. On GSE77348 narrow, Δβ-only (0.964) sits between V1 (0.969) and V2 `tumor_only` (0.911) — a tight cluster that underscores the "V2.5 generation earns a small but consistent margin" story, not a dominance claim.

**P@K intervals (GSE322563, narrow labels).** On `n = 2/2` the score distribution at K = 100 sits inside a tied band; P@100 is an adversarial interval, not a point estimate:

| axis | P@100 observed | P@100 min | P@100 max | tie_band@100 |
|---|---:|---:|---:|---:|
| V2.5 differential | 0.000 | 0.000 | 0.020 | 190 |
| tumor_only | 0.000 | 0.000 | 0.020 | 10,005 |
| V1 `final_score` | 0.020 | 0.020 | 0.020 | 1 |

V2.5's band at K = 100 is 190 (interval [0.00, 0.02]); V2 `tumor_only`'s band is 10,005, spanning two orders of magnitude more of the distribution — P@K on `tumor_only` is not a meaningful quantity at this cohort size. V1's band is 1 (continuous score), and the interval collapses to the observed value.

**Native EPIC v2 vs HM450-intersect (GSE322563).** Re-running under the native EPIC v2 probe set (no HM450 intersect) gives:

| tier | Δβ-only (HM450 / native) | V1 (HM450 / native) | V2 tumor_only (HM450 / native) | V2.5 (HM450 / native) |
|---|---:|---:|---:|---:|
| validated | 0.974 / 0.961 | 0.821 / 0.933 | 0.928 / 0.936 | 0.990 / 0.986 |
| narrow | 0.912 / 0.956 | 0.884 / 0.938 | 0.886 / 0.933 | 0.942 / 0.983 |
| wide | 0.844 / 0.865 | 0.768 / 0.855 | 0.871 / 0.916 | 0.910 / 0.945 |

The V2.5 primary-endpoint AUC differs by 0.004 under native vs intersect (0.990 vs 0.986) — within the tied-band noise floor. V1 and V2 `tumor_only` gain more under native on average (V1 +0.054 to +0.112; V2 +0.008 to +0.047), consistent with the broader probe coverage exercising more of the continuous-score dynamic range. The only axis-ordering invariant that holds at every matched-cell-line tier × path row (9/9) is **V2.5 > Δβ-only** and **V2.5 > V1**; the relative order of Δβ, V1, and V2 `tumor_only` reshuffles across tiers (e.g. on GSE322563 HM450 wide and native wide, `tumor_only` sits above Δβ; on GSE77348 all three tiers, `tumor_only` is the lowest). The HM450-intersect shortcut does not distort the headline claim: V2.5 is the top-scoring axis on every matched-cell-line row on both paths.

### 5.3 Tissue cohort — GSE69914 (n = 305/50)

At high replicate count, `p_trust` desaturates and the composite's ordering is driven by continuous `p_targ × p_diff`. On GSE69914 (numbers straight from the committed `examples/gse69914_roth_labels/bench_*.jsonl`):

| axis | validated AUC | narrow AUC | wide AUC | tie_band@100 |
|---|---:|---:|---:|---:|
| Δβ-only (naive) | 0.591 | 0.477 | 0.435 | n/a (K=20) |
| V1 `final_score` | 0.660 | 0.539 | 0.435 | 1 |
| V2 `tumor_only` | **0.803** | **0.843** | **0.874** | 6,540 |
| V2.5 differential | 0.773 | 0.711 | 0.726 | 2 |

- **The tied-band prediction from §2.4 holds.** V2.5's band at K = 100 shrinks from 190 (GSE322563 `n = 2/2`) to 299 (GSE77348 `n = 3/3`) to **2** here (`n = 305/50`) with no hyperparameter changes. At `n ≳ 30`, `p_trust` saturates to the per-`EvidenceClass` ceiling (a constant, not a continuous function of n); the tied-band shrinkage comes from `p_targ × p_diff` no longer being pushed to 1 by wide-IQR σ-floor saturation.
- **Δβ-only collapses to near-random** (0.435–0.591). Intuition: both arms of GSE69914 are bulk-tissue β values — 305 sporadic breast tumors and 50 healthy-donor breast tissues — and each sample is a mixture of cell types (epithelial, stromal, immune, adipose) rather than a pure population. A CpG that is cleanly hypomethylated in a pure MCF-7 culture can appear partially methylated in a bulk breast sample because the cell type carrying the hypomethylation is diluted at the same locus by other cell types in the mix. Large matched-cell-line gaps on the order of 0.8 attenuate to smaller, noisier gaps of 0.1–0.3 on bulk tissue, which destroys the ranking signal for a raw-magnitude Δβ ranker. V2.5 recovers the ordering by weighting each candidate's tumor–normal gap by its per-probe uncertainty (§2.1) rather than by magnitude alone.
- **V1 degrades sharply on tissue.** V1 reaches 0.435 on wide (at chance) and 0.539 on narrow — worse than V2.5 on every tier. The V1 continuous-score advantage (`tie_band = 1` by construction) does not translate into AUC leadership on tissue; it just makes V1's failure mode *deterministic*.
- **`tumor_only` wins raw AUC but is excluded from prioritization recommendations.** Its tie_band of 6,540 at K = 100 means its top-100 window sits inside a massively tied region — ordering within the top-K is a coin flip under any tie-break. Retained as a diagnostic axis (AUC sanity check), not a prioritization axis.
- **V2.5-sigmoid is the recommended tissue prioritization axis.** V2.5-diff reaches +0.113 / +0.172 / +0.291 AUC over V1 on validated / narrow / wide with stable top-K (`tie_band = 2`); V2.5-sigmoid additionally beats V2.5-diff by +0.05 to +0.09 AUC across the §5.2.1 bandwidth family and at the WG denominator in PAPER.md §5.2.2 (0.862 vs 0.773 at σ_fixed ≈ 0.0707 on the validated label set). Both ship as selectable modes in this tag; V2.5-sigmoid is the recommended default and V2.5-diff is retained for AUC-parity use. The tissue-specific gain of V2.5-sigmoid over V2.5-diff is a direct consequence of the PAPER.md §3.5 binding-rate picture: on cell lines `σ_floor` binds almost everywhere and V2.5-diff's `p_diff` collapses to a fixed-bandwidth response anyway (so the two modes are equivalent there); on tissue the per-site σ varies and hurts the ranking, so the explicit fixed-bandwidth V2.5-sigmoid wins.

**The cross-cohort picture after whole-genome validation supports `gap_sigmoid` for new prioritization runs.** PAPER.md §5.2.2 cross-cohort panel: AUC within 0.002 on matched cell-line (cell-line AUC parity), +0.05 to +0.08 transported-label AUC on the single tissue cohort tested, and crucially `tie_band@100 = 1` for `gap_sigmoid` vs 421–1,493 for shipped V2.5 across the three WG matched cell-line cohorts — V2.5's top-K is unusable at WG scale on n = 2/2 or 3/3 cohorts, while `gap_sigmoid`'s is clean. `gap_sigmoid` ships as `probabilistic_mode: tumor_plus_gap_sigmoid` in this tag. Shipped V2.5 (`tumor_plus_differential_protection`) is retained as a selectable mode for backward compatibility but is no longer the recommended prioritization axis for new non-boundary runs. V1 remains the stable-release default because of its continuous-score guarantees (deterministic `tie_band = 1`), not because it wins AUC anywhere — §6.2 explains the shipping logic.

### 5.4 Out-of-distribution boundary — GSE68379

GSE68379 is the 52-line Sanger GDSC1000 breast cell-line panel, paired cross-series with 50 external healthy-donor normals (from GSE69914 status=0) to form the OOD test cohort. Sanger's MCF-7 is a member of this panel, and the panel's consensus methylation at the Roth Fig. 5d validated sites is shifted relative to Roth's own MCF-7 — well-documented cell-line epigenetic drift. Applying Roth-derived positives across this cross-series cohort should fail: the validated sites are hypermethylated across most Sanger breast lines where they are hypomethylated in Roth MCF-7. AUC inverts below 0.5 on both Δβ-only and V2.5:

| axis | validated AUC | narrow AUC | wide AUC |
|---|---:|---:|---:|
| Δβ-only (naive) | 0.201 | 0.150 | 0.194 |
| V2.5 differential | 0.197 | 0.169 | 0.269 |

Both the literature-naive baseline and V2.5 invert correctly below 0.5 at every tier: the scorers are consistently identifying that this label set does not transport from Roth MCF-7 to the Sanger panel. This is documentation of a known boundary (cross-series cell-line drift + cross-lab batch), not a negative generalization result, and is not pooled with the §5.1 primary-endpoint numbers.

### 5.5 Top-hit annotation

For each cohort, the top-20 candidates are annotated with:

- nearest gene symbol, TSS distance, and feature class (promoter / gene_body / intergenic),
- CpG-island context (island / shore / shelf / open_sea) and distance to the nearest CGI boundary,
- RepeatMasker overlap (class / family / name) from UCSC rmsk,
- ENCODE DNase-HS cluster membership and breadth (number of cell types with the peak) from `wgEncodeRegDnaseClusteredV3`.

A rule-based flag set (STRONG: island-localized promoter / active promoter with DNase support; CAUTION: overlaps repeat / sparse evidence / small Δβ; NOTE: gene body without DNase support) is emitted in both the TSV and a Markdown companion file. On GSE69914's top-20, 6/20 carry at least one STRONG flag and 1/20 is caution-free — a distribution that is directly useful for experimental triage. Example annotated shortlists are committed under `examples/*/top20_annotated_v25.md`.

The formal SCREEN cCRE Registry is gated behind a JavaScript challenge on `screen.wenglab.org` and could not be fetched directly; DNase-HS clusters serve as the v1 regulatory-activity proxy. This is documented in the script's docstring; a user with SCREEN access can swap in the formal Registry with a one-line change.

![Top-20 gene presence per axis × cohort](docs/figures/fig3_topgene_heatmap.png)

**Figure 3.** Tie-window-aware top-20 gene presence. Columns are per-cohort per-axis top-20 shortlists (ordered by the deterministic tie-break documented in §3.1); rows are the union of nearest-gene symbols sorted by cross-column presence. Column-header colors: **blue** = deterministic top-K under the current score distribution (`tie_band = 1` for V1, `tie_band = 2` for V2.5 on high-`n` tissue); **red** = top-K is a 20-record window inside a large tied band (`tie_band = 190`, `299`, `421` for V2.5 on low-`n` cell-line cohorts), so membership is the documented tie-break policy's choice, not the scorer's discrimination. Genes bolded in blue on the row labels are present in **all three** cell-line V2.5 shortlists (GSE322563 HM450, GSE322563 native, GSE77348) — cross-platform / cross-laboratory convergence gives two such genes: *CELF2* and *XPNPEP1*. **The tissue cohort (GSE69914 V2.5) has zero nearest-gene overlap with any of the three cell-line V2.5 shortlists** — its top-20 candidates map to exactly 9 distinct nearest-gene symbols (*COL21A1, FOXC1, FOXCUT, MAS1L, MAT2B, MXI1, RUNX2, SCML4, TENM2*), all of which are disjoint from every cell-line V2.5 top-20. This is a direct visualization of the cell-line-vs-tissue biological shift that drives the AUC patterns in Figure 2.

---

## 6 Discussion

### 6.1 What V2.5 is and is not

On matched cell-line cohorts — where the PAM cytosine's methylation state cleanly separates the tumor and normal arms — raw Δβ is already a strong baseline and V2.5 improves on it by a consistent but small margin (+0.010 to +0.080 across all nine tier × path rows on GSE322563 and GSE77348). **V2.5 does not dominate Δβ by a large AUC margin on easy cohorts, and we do not claim it does** — the practical case for V2.5 over Δβ on such cohorts is structural, not a matter of AUC supremacy. What V2.5 *does* do, relative to all three baselines, is summarized below:

1. **V2.5-diff matches the highest-AUC prioritization axis at every matched-cell-line and tissue cohort × tier combination tested. V2.5-sigmoid matches V2.5-diff on AUC *on matched cell lines* (within 0.002 across GSE322563 HM450, GSE322563 native EPIC v2, and GSE77348; PAPER.md §5.2.2 WG panel) and strictly beats V2.5-diff on WG top-K usability — `tie_band@100 = 1` vs 421 / 1,127 / 1,493 under the WG denominator. *On tissue (GSE69914), V2.5-sigmoid beats V2.5-diff by +0.05 to +0.08 transported-label AUC* (0.862 vs 0.778 at the validated tier), so the parity statement is cell-line-only and the tissue win is separate.** V2.5-sigmoid ships as `tumor_plus_gap_sigmoid` in this tag and is the recommended probabilistic prioritization axis per PAPER.md §5.2.2 whole-genome panel. (The OOD boundary case GSE68379 correctly inverts for all axes and is reported separately in §5.4; it is not pooled with the prioritization-mode comparisons here.) Row-by-row margins against the shipped baselines:
    - V2.5 over Δβ-only on 9 matched-cell-line rows: **+0.010 to +0.080**, consistently positive.
    - V2.5 over V1 on 9 matched-cell-line rows: **+0.014 to +0.169**, small on GSE77348 (~+0.015) and moderate-to-large on GSE322563 (up to +0.169 on HM450 validated). No single "typical" margin — the magnitude depends on the cohort.
    - V2.5 over V1 on the 3 tissue rows: **+0.113, +0.172, +0.291** (validated, narrow, wide).
    - No general order invariant across Δβ and V1 — both above and below V2.5 reshuffle relative to each other across rows (§5.2). The only invariants that hold 9/9 on matched cell-line are *V2.5 > Δβ* and *V2.5 > V1*.
    - The raw-AUC winner on tissue is V2 `tumor_only` (0.803–0.874), but its `tie_band@100 = 6,540` puts its top-100 inside a massively tied region, disqualifying it as a prioritization axis (retained as a diagnostic for AUC cross-checking).
2. **V2.5 rescues tissue-cohort ranking that V1 and Δβ both lose.** V1 collapses on GSE69914 wide (0.435, at chance) and degrades sharply on narrow (0.539); Δβ-only collapses across tiers (0.435–0.591). The mechanism: GSE69914 is bulk breast tissue on both arms (sporadic tumors + healthy-donor tissue), and bulk-sample β values are cell-type mixtures rather than pure populations, which attenuates clean cell-line differentials to small noisy gaps at the same loci (§5.3). V2.5 recovers the ordering because its IQR-derived per-probe σ (meaningful at n = 305 tumor / n = 50 normal) lets `p_diff` weight each candidate's gap by its confidence — probes with large-but-noisy gaps get down-weighted versus probes with moderate-but-tight gaps. Raw Δβ has no such mechanism, and V1's weighted-sum heuristic does not propagate per-probe uncertainty. The tie-band prediction from §2.4 also holds: V2.5's band at K = 100 shrinks from 190 (`n = 2/2`) → 299 (`n = 3/3`) → 2 (`n = 305/50`) with no hyperparameter changes.
3. **V2.5 emits tie-band-honest top-K reporting.** `tie_band_size_at_k` and `precision_at_k_{min, max}` let a reader see when the top-K window sits inside a saturated region of the composite (the low-replicate case). Δβ and V1 have `tie_band = 1` by construction (continuous scores) so they do not need the explicit reporting; V2.5 needs it because of its discrete `EvidenceClass` saturation, and reporting the intervals is what makes its low-`n` top-K honest.
4. **V2.5 produces a probability-scale interpretation.** Each per-site *probability-scale selectivity score* ∈ [0, 1] is composable with downstream probabilistic inputs (target-mutation models, gRNA off-target probabilities, delivery-efficiency priors). Δβ is a signed similarity and V1 is a weighted-sum heuristic; neither is directly multiplied into such a chain.
5. **On the OOD boundary (GSE68379), both probabilistic axes and Δβ-only invert correctly** below 0.5 — the scorers are consistently identifying that the Roth label set does not transport to Sanger MCF-7. This is a sanity check, not a generalization claim.

This is the methodologically honest framing: Δβ is a strong baseline that should be reported alongside any new methylation-guided Cas9 scorer; `gap_sigmoid` (shipped in this tag as `tumor_plus_gap_sigmoid`) is the right research mode across every non-boundary cohort shape we tested (matched cell-line *and* tissue) after the whole-genome denominator check in §5.2.2 revealed that shipped V2.5's top-K blows up at WG scale on n = 2/2 and 3/3 cohorts while `gap_sigmoid`'s stays clean; shipped V2.5 is retained as a selectable mode for backward compatibility but dominated by `gap_sigmoid` everywhere we measured. V1 remains the stable-release default for reasons of backward compatibility and its deterministic `tie_band = 1` top-K — *not* because V1 wins AUC anywhere. §6.2 below states the shipping logic explicitly.

### 6.2 Decision hierarchy

Three truths coexist:

1. **V1 `final_score` is the stable-release default** (tagged `v0.4.0`). Its continuous deterministic score gives `tie_band = 1` on every cohort tested — top-K is never tie-break-dependent, regardless of replicate count or cohort shape. V1 is the default because of this property and for backward-compatibility with existing users, not because it wins AUC anywhere.
2. **The recommended probabilistic axis is `gap_sigmoid`** (`probabilistic_mode: tumor_plus_gap_sigmoid`, ships in this tag). PAPER.md §5.2.2 cross-cohort WG panel shows `gap_sigmoid` is AUC-equivalent to shipped V2.5 on matched cell lines (GSE322563 HM450 0.989 vs 0.989, native v2 0.998 vs 0.998, GSE77348 0.982 vs 0.982), has higher transported-label AUC on the single tissue cohort tested (GSE69914 0.862 vs 0.778), and critically has `tie_band@100 = 1` at WG scale on every matched cell-line cohort (vs V2.5's 421–1,493 on the same WG catalogs, where V2.5's top-K becomes unusable). Shipped V2.5 (`probabilistic_mode: tumor_plus_differential_protection`) is retained as a selectable mode for backward compatibility and AUC parity on cell-line cohorts but is no longer recommended over `gap_sigmoid` for new prioritization runs. For GSE68379 (OOD cross-series), no axis is supported.
3. **V2 `tumor_only` is retained as a diagnostic**, not a prioritization axis. Its tied bands at K = 100 are **5,271–14,914 across the five cohort paths tested** (GSE68379 5,271; GSE69914 6,540; GSE322563 HM450 10,005; GSE77348 11,848; GSE322563 native 14,914) — in every case large enough that ordering within the top-100 is determined by tie-break, not by score. Its raw AUC can lead (as on GSE69914 tissue, 0.803–0.874), but its top-K is not usable for target-shortlist construction.

The decision table below is the literal recommended hierarchy.

### 6.3 Limitations

1. **Low-replicate tie bands.** On cohorts with `n < 30`, top-K is reported as an interval via `precision_at_k_{min, max}` and `tie_band_size_at_k`. The problem dissolves at `n ≳ 30` (demonstrated on GSE69914).
2. **Cell-line drift.** GSE68379 confirms that "MCF-7" is not a single epigenetic object across labs. Any V2.5 shortlist must be validated against the specific cell line intended for the follow-up editing assay.
3. **V2.5-diff σ_floor = 0.05** is empirical. A sweep over σ_floor ∈ {0.02, 0.05, 0.10, 0.15} (PAPER.md §5.3.1; `examples/sigma_floor_sweep.*`) shows the matched-cell-line AUC ranking is robust (≤0.02 movement across a 7.5×-range sweep) but tissue AUC peaks at σ_floor = 0.10, not the shipped default. This σ_floor sensitivity is a V2.5-diff property only; V2.5-sigmoid uses a fixed bandwidth `σ_fixed` (default ≈0.0707) and is not σ_floor-dependent. Formal regime-specific σ_floor default selection for V2.5-diff is deferred.
4. **Catalog scope.** Primary tables in §5.1 use chr5/6/10 candidates within ±500 bp of a methylation probe (the chromosomes carrying Roth Fig. 5d targets). The frozen whole-genome panel cited in the abstract (PAPER.md §5.2.2) extends the matched cell-line + tissue evidence to 19.8M HM450 (SHA256 `d20661c5…`) and 35.4M EPIC v2 (SHA256 `39df8f0f…`) candidates, and PAPER.md §5.9 adds within-chromosome feature-matched negative controls. Cross-chromosome and chromosome-class matched controls remain future work.
5. **SCREEN cCRE Registry** is not integrated; DNase-HS clusters are a proxy.
6. **No prospective wet-lab validation.** The claim is AUC + tie-band-aware top-K on public methylation data. All named top hits are unvalidated predictions.

### 6.4 Related work

Two bodies of prior art are adjacent but non-overlapping. Generic CRISPR guide-scoring tools (e.g. Azimuth, DeepCRISPR, CRISPRitz) score sequence properties and off-target risk; none of them incorporate methylation state at the PAM cytosine, because none of them target a methylation-sensitive Cas9 variant. Methylation-aware variant analysis tools (e.g. methylKit, minfi) provide per-probe differential analysis on arrays but do not produce ranked target lists; their output feeds a scorer, not replaces one. `thermocas` is the first open framework we are aware of that joins the two — a ranking scorer that consumes methylation array summaries and outputs tumor-selective ThermoCas9 target shortlists with honest uncertainty.

### 6.5 Next steps

1. Prospective wet-lab validation of 3–5 top V2.5-sigmoid candidates on MCF-7 vs MCF-10A (β-pyrosequencing + ThermoCas9 editing readout). This is the terminal step that changes the scientific claim.
2. **Regime preset selector** on top of the shipped `tumor_plus_gap_sigmoid` mode (PAPER.md §6.3): a one-shot `regime` YAML field that picks (gap factor, σ_fixed, δ) defaults per cohort regime. Per §6.1, all non-boundary presets recommend V2.5-sigmoid (the `matched_cell_line` and `primary_tissue` presets differ only in the `EvidenceClass` ramp-n implied by sample count, not in the gap factor). A `backward_compat` preset that opts into V2.5-diff is retained for users needing AUC parity with pre-ag scored JSONLs.
3. **Cross-chromosome denominator controls** — the WG panel (§5.2.2) and within-chromosome feature-matched audit (PAPER.md §5.9) are complete; chromosome-class and random-chromosome matched controls remain as loosening axes.
4. **Canonical R `limma` parity** — completed at `memo-2026-04-22-ba`. Per-probe Spearman/Pearson `t_mod` ≥ 0.9997 across all three primary cohorts versus R `limma::lmFit + eBayes` 3.66.0 on identical sample-level β + group inputs; full results in PAPER.md §5.8 with per-cohort artifacts under `examples/r_limma_parity_*.{tsv,md}`.
5. A second independent-lab matched MCF-7/MCF-10A EPIC v2 cohort if one becomes public, to replicate the V2.5 headline at `n ≥ 3` rather than `n = 2`.
6. Formal SCREEN cCRE Registry integration when a non-challenged download path becomes available.

---

## 7 Conclusion

`thermocas` offers a *compositional probability-scale scoring skeleton* (`p_targ × (gap factor) × p_trust`) paired with a *tie-band-aware benchmarking contract* for methylome-guided ThermoCas9 target prioritization. The compositional skeleton is the durable contribution: it decouples three independently-replaceable questions ("is the tumor side targetable?", "is the tumor-vs-normal gap large enough?", "is the observation trustworthy?"), and the §5.2.1 ablation + §5.2.2 cross-cohort whole-genome panel exercise that decomposition by swapping the gap factor between V2.5's `p_diff` and a simpler fixed-bandwidth `gap_sigmoid`. The current shipped recommendation that emerges from those experiments is `gap_sigmoid` (`probabilistic_mode: tumor_plus_gap_sigmoid`, ships in this tag): AUC-equivalent to the earlier `tumor_plus_differential_protection` (V2.5) on matched cell lines, higher transported-label AUC on the single tissue cohort tested, and substantially better low-`n` WG top-K usability. Shipped V2.5 is retained as a selectable mode for backward compatibility and AUC parity on cell-line cohorts but is no longer the recommended prioritization axis. V1 remains the stable-release default for its deterministic `tie_band = 1` top-K guarantee. Δβ-only is a strong matched-cell-line baseline that any methylome-guided Cas9 scorer should report alongside its own scoring axis, and we report it on every cohort × tier in §5.2. The combination is offered as an open research framework; groups extending Roth et al. to other cell-type pairs are the immediate target audience.

---

## Data and code availability

- **Code**: <https://github.com/AllisonH12/thermocas9>, tagged `memo-2026-04-22-be` for the exact revision evaluated here. BSD-3 licensed. `git rev-parse memo-2026-04-22-be` resolves to a SHA in a fresh clone. Supersedes all prior same-day dated memo tags; see the PAPER.md tag ledger for per-tag provenance (what each tag added and what was corrected in the next one). A `scripts/verify_manuscript_claims.py` guard enforces a set of *selected known-risk* numerical claims — constants (`src/thermocas/probabilistic.py`), universal-quantifier wording, artifact counts, figure captions, and the tag-span claims most recently caught in reviewer cycles — against the committed bench JSONLs and package source. It is not a full table-AUC verifier; new numerical claims should be audited by hand or added to the `check_*` function set. Run it before cutting any subsequent tag.
- **Tests**: 245 passing under `uv run pytest -q`.
- **Cohorts**: public GEO series GSE322563, GSE77348, GSE69914, GSE68379. Build scripts in `scripts/build_gse*_cohort.py`.
- **Reference annotations**: UCSC hg19 `refGene.txt.gz`, `cpgIslandExt.txt.gz`, `rmsk.txt.gz`, and `wgEncodeRegDnaseClusteredV3.txt.gz` (fetched on demand from hgdownload.soe.ucsc.edu).
- **Positives**: `data/derived/positives_roth_{validated, narrow, wide}.txt` (HM450 path) and `data/derived/epic_v2_positives/positives_roth_*.txt` (native EPIC v2 path). Derived from Roth et al. Fig. 5d via Ensembl REST hg38 → hg19 liftover.
- **Benchmark artifacts**: every `BenchmarkResult` JSONL under `examples/*_roth_labels/` carries `precision_at_k`, `precision_at_k_{min, max}`, `recall_at_k`, `recall_at_k_{min, max}`, `roc_auc`, `tie_band_size_at_k`, and `tie_break_policy`.
- **Feature-matched controls**: `examples/feature_matched_negative_controls.{tsv,md}` reports the within-chromosome matched-negative audit used in PAPER.md §5.9.
- **Annotated top-20 TSV + Markdown shortlists** under `examples/*/top20_annotated_v25.{tsv,md}`.

## Competing interests

The author declares no competing interests.

## Acknowledgements

This work would not exist without Roth et al.'s characterization of ThermoCas9 and their transparent release of the GSE322563 matched cell-line methylation data.

---

## Appendix A · Reproducibility

From a fresh clone of the repository at `memo-2026-04-22-be`:

```bash
# One-time env setup
uv sync
uv run pytest -q           # 245 passing

# Rebuild the EPIC v2 probe annotation on chr5/6/10 (~2 min, needs pyliftover)
uv run python scripts/build_epic_v2_probes.py \
  --soft data/raw/epic_v2/GPL33022_family.soft.gz \
  --output data/raw/probes_hg19_epic_v2.tsv

# Rebuild the GSE322563 native EPIC v2 cohort (~30 s)
uv run python scripts/build_gse322563_native_epic_v2_cohort.py \
  --beta-matrix data/raw/gse322563/GSE322563_beta_matrix_EPIC_v2.txt.gz \
  --epic-v2-probes data/raw/probes_hg19_epic_v2.tsv \
  --output-dir data/derived/gse322563_native_epic_v2_cohort/

# Score + benchmark one mode on one cohort (~90 s on a modern laptop)
uv run thermocas score-cohort \
  --cohort data/derived/gse322563_native_differential.yaml \
  --catalog data/derived/catalog_hg19_chr5_6_10.jsonl \
  --output data/derived/scored_gse322563_native_differential.jsonl

uv run thermocas benchmark \
  --scored data/derived/scored_gse322563_native_differential.jsonl \
  --positives data/derived/epic_v2_positives/positives_roth_validated.txt \
  --cohort-name GSE322563-native-validated-V2.5 \
  --score-field p_therapeutic_selectivity --top-k 20 \
  --no-enforce-holdout \
  --output examples/gse322563_native_roth_labels/bench_validated_differential.jsonl

# Annotate the top 20 with gene / CGI / repeat / DNase context + Markdown companion
uv run python scripts/annotate_top_hits.py \
  --scored data/derived/scored_gse322563_native_differential.jsonl \
  --top-k 20 \
  --rmsk data/raw/ucsc/rmsk.txt.gz \
  --dnase data/raw/ucsc/wgEncodeRegDnaseClusteredV3.txt.gz \
  --positives data/derived/epic_v2_positives/positives_roth_wide.txt \
  --output examples/gse322563_native_roth_labels/top20_annotated_v25.tsv \
  --markdown examples/gse322563_native_roth_labels/top20_annotated_v25.md
```

All tables and figures in this manuscript are reproducible from the committed artifacts alone — no re-scoring required.

---

## Appendix B · Notes on this manuscript file

`MANUSCRIPT.md` is the Bioinformatics-submission-shaped version. The audit-trail memo that preserves the full development history (V2 misspecification narrative, pre-registration block, decision-table iterations, and the tag-stability discussion) remains as `PAPER.md` in the same repository. Both files point at the same underlying code and benchmark artifacts. Where a discussion point is elided in this manuscript for length, `PAPER.md` carries the longer form.
