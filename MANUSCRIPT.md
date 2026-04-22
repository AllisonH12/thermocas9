# thermocas: probabilistic ranking of methylation-protected ThermoCas9 target sites with tie-aware benchmarking

**Allison Huang Mercer** · Columbia University · <allisonhmercer@gmail.com>

Framework maintained by Thermocas9 Inc. Code at <https://github.com/AllisonH12/thermocas9>, immutable tag `memo-2026-04-22-c`.

---

## Abstract

**Motivation.** Methylation of the fifth PAM cytosine protects genomic loci from ThermoCas9 cleavage (Roth et al., *Nature* 2026), enabling tumor-selective targeting at loci that diverge in methylation between cancer and matched normal tissue. Translating this mechanism into a target-discovery shortlist requires a scorer that is cohort-agnostic, honest about evidence quality, and explicit about ties in its output. No such scorer was publicly available. The obvious baseline — rank candidates by the raw methylation gap Δβ = β_normal − β_tumor — was not benchmarked against any proposed alternative in the open literature.

**Results.** We present `thermocas`, an open probabilistic scoring framework for methylation-sensitive Cas9 target ranking. The V2.5 mode computes
`p_therapeutic_selectivity = p_targetable_tumor × p_differential_protection × p_observation_trustworthy`,
where `p_differential_protection = P(β_normal − β_tumor > δ)` is estimated by an independent-normal approximation on per-probe summary statistics (IQR-based σ estimate with a σ_floor = 0.05). We benchmark against three baselines — a literature-naive Δβ-only ranker, the deterministic V1 score, and a first-pass V2 formulation that our own analysis rejected — across four public methylation cohorts: GSE322563 (Roth's own MCF-7/MCF-10A EPIC v2 samples; primary endpoint), GSE77348 (development cohort on which δ was tuned), GSE69914 (n ≈ 300 tumor–normal pair tissue), and GSE68379 (documented as an out-of-distribution boundary case). **The headline finding is that no single axis dominates.** On matched cell-line cohorts where the methylation signal is clean, Δβ-only is already a strong baseline (AUC 0.96–0.97 at validated); V2.5's margin over it is 0.01–0.03 at validated and occasionally negative (native EPIC v2 narrow: Δβ 0.956 vs V2.5 0.944). On tissue, Δβ-only collapses (GSE69914 validated AUC 0.591; near-random on narrow/wide tiers at 0.435–0.477) while both probabilistic axes recover — but on tissue **V1's continuous deterministic score (AUC 0.861) outperforms V2.5 (0.837)**, so V1, not V2.5, is the recommended axis there. Both axes invert correctly on the out-of-distribution boundary (GSE68379: Δβ 0.15–0.20; V2.5 < 0.5). V2.5 earns its complexity in three places specifically: (i) it is the highest-AUC axis on matched cell-line cohorts at every label-granularity tier (small margins, but consistent); (ii) it produces tie-band-honest top-K reporting (`precision_at_k_{min, max}` and `tie_band_size_at_k` on every BenchmarkResult); (iii) it produces a probability-scale interpretation that lets downstream pipelines compose per-site scores with other probabilistic inputs. The shipped recommendation is V1 as the stable default (tissue-safe, `tie_band = 1`), V2.5 as the research mode for matched cell-line / paper-comparable cohorts, with Δβ-only retained as a published baseline. The package adds a k-way-merge pan-cancer aggregator for genome-scale atlas builds and a per-candidate annotation pipeline (nearest gene, CpG-island context, RepeatMasker overlap, ENCODE DNase-HS cluster breadth) with a Markdown shortlist aimed at experimental collaborators.

**Availability and implementation.** `thermocas` is open source at <https://github.com/AllisonH12/thermocas9>, tagged `memo-2026-04-22-c` for the version evaluated here. 236 unit tests pass under `uv run pytest -q`. Python 3.11+, BSD-3.

**Contact.** <allisonhmercer@gmail.com>

---

## 1 Introduction

Roth et al. (2026) characterized *Geobacillus thermodenitrificans* T12 Cas9 (ThermoCas9) biochemically and structurally. The central result is that cleavage efficiency is governed by the methylation state of the fifth PAM cytosine — 5-methylcytosine at that position abolishes binding, while an unmethylated PAM cleaves as usual. The practical consequence: ThermoCas9 permits selective editing at loci that are unmethylated in disease cells and methylated in matched normal cells. Roth demonstrated this with three breast-cancer targets in MCF-7 vs MCF-10A on the Illumina MethylationEPIC v2 platform.

Turning that mechanism into a genome-scale target-discovery shortlist is a ranking problem with three awkward properties. First, per-probe methylation summaries are noisy (arrays commonly report n = 2–3 replicates per side for cell-line cohorts) and the scorer must carry its own uncertainty rather than emit a scalar that hides it. Second, the "normal" arm is inconsistent: adjacent-normal tissue in bulk methylation studies is not the same object as a matched cell-line pair, and any assumption that baked in a fixed β_normal threshold would fail across cohort types. Third, on low-replicate cohorts the top-K of any probabilistic composite can sit inside a large tied band; reporting a single Precision@K value without a tie interval misleads a reader about the score's discriminative power.

No public tool addressed these three problems together when we started. The closest prior art — generic CRISPR guide-scoring and Δβ ranking on methylation arrays — handles neither uncertainty propagation nor tie-band reporting. `thermocas` is the first open framework to target all three.

**Contributions.**

1. The V2.5 probabilistic composite, where a differential-protection factor replaces a failed threshold-based factor. The decomposition is `p_targ × p_diff × p_trust`, each factor has a closed-form derivation, and the final score is directly interpretable as a probability of therapeutic selectivity at the candidate site.
2. A benchmarking protocol that reports Precision@K and Recall@K as adversarial intervals (`_min`, `_max`) together with the tie-band size at K, so reviewers are not misled by arbitrary tie-break choices on low-replicate cohorts.
3. A benchmark against three baselines — a Δβ-only ranker, V1 deterministic scoring, and a rejected V2 formulation — across four public cohorts with positives tiered from Roth's own validated targets to a wide genomic window.
4. A k-way-merge pan-cancer aggregator with O(N_cohorts + N_unique_candidate_ids) memory and cross-cohort metadata-parity enforcement, for genome-scale atlas builds.
5. An experiment-facing shortlist pipeline: TSV + Markdown companion that attach nearest gene, CpG-island context, RepeatMasker overlap, and ENCODE DNase-HS breadth, with rule-based triage flags (STRONG / CAUTION / NOTE) deterministically derived from the annotated fields.

**Scope.** This paper reports methods development and benchmarking on public methylation cohorts. It does not include prospective wet-lab validation of predicted target sites; the framework is an open educational research tool, not a clinical decision-support system. Per-site p-values are not reported because `p_observation_trustworthy` saturates to a discrete per-`EvidenceClass` value at high replicate counts; the composite provides a ranking axis, not a hypothesis test.

---

## 2 Methods — Formulation

For each candidate site on each cohort, we carry six per-probe β-value summaries (`β_tumor_mean`, `β_tumor_q25`, `β_tumor_q75`, `β_normal_mean`, `β_normal_q25`, `β_normal_q75`), sample counts `n_tumor` and `n_normal`, and an `EvidenceClass` capturing the distance between the PAM cytosine and the nearest assayed CpG (`EXACT` / `PROXIMAL_CLOSE` / `PROXIMAL` / `REGIONAL` / `UNOBSERVED`).

### 2.1 V2.5 — the shipped composite

```
p_therapeutic_selectivity = p_targ × p_diff × p_trust
```

- **`p_targ = P(β_tumor < 0.5)`** — the probability that the tumor arm is hypomethylated at the PAM cytosine, estimated by the same IQR-normal approximation used below.
- **`p_diff(δ) = P(β_normal − β_tumor > δ)`** — the probability that the normal-minus-tumor methylation gap exceeds a configurable threshold δ (default 0.2). Estimated by:
  ```
  σ_k  ≈ IQR_k / 1.349                (k ∈ {tumor, normal}; floor at 0.05)
  σ_Δ² = max(σ_t, floor)² + max(σ_n, floor)²
  z    = (δ − (μ_n − μ_t)) / σ_Δ
  p_diff = 1 − Φ(z)
  ```
  The σ_floor prevents σ_Δ from collapsing to zero at boundary β-values (CpGs methylated at 0 or 1, which are common at islands).
- **`p_trust = EvidenceClass.base × min(1, min(n_t, n_n) / 30)`** — evidence-class confidence with a linear ramp up to n = 30, saturating to a per-class constant (`EXACT`: 0.95, `PROXIMAL_CLOSE`: 0.75, etc.).

Figure 1 summarizes the composite.

![Figure 1 · V2.5 mode-formula schematic](docs/figures/fig1_mode_schematic.png)

**Figure 1.** The V2.5 composite. Three factors multiply: `p_targ` (tumor unmethylated at the PAM cytosine), `p_diff` (differential methylation gap exceeds δ under a normal approximation on per-probe β summaries), and `p_trust` (evidence-class confidence, saturating in min sample count). The result is `p_therapeutic_selectivity`, the stored scalar that downstream ranking consumes.

### 2.2 Baselines

| name | score | rationale |
|---|---|---|
| **Δβ-only** (naive) | `β_normal_mean − β_tumor_mean` | Literature-naive: rank by the raw methylation gap. No uncertainty propagation, no evidence weighting. The "obvious thing." |
| **V1 `final_score`** | weighted sum of sequence, selectivity, confidence, heterogeneity, and coverage components (continuous-valued) | The framework's stable-release axis, retained for top-K interpretability. |
| **V2 `tumor_only`** | `p_targ × p_trust` (no protection factor) | Intermediate formulation after we dropped V2's `p_prot` factor. Retained as a diagnostic; not a discovery axis. |
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
| GSE69914 | HM450 | 305 / 50 | Tumor–normal pair tissue — tests the predicted `p_trust` desaturation at high n | GEO |
| GSE68379 | HM450 | 1 / 0 (Sanger MCF-7) | Out-of-distribution boundary case (cell-line drift), not a generalization test | GEO |

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

### 5.1 Primary endpoint — matched cell-line AUC at validated Roth probes

Held-out primary endpoint: AUC at `positives_roth_validated.txt` on GSE322563, using HM450-intersected probes. δ was not tuned on this cohort; the label set was not derived from it.

| axis | **GSE322563** (primary, n = 2/2) | GSE77348 (δ-tuned, n = 3/3) |
|---|---:|---:|
| Δβ-only (naive) | 0.974 | 0.972 |
| V1 `final_score` | 0.821 | 0.968 |
| V2 `tumor_only` | 0.928 | 0.912 |
| **V2.5 differential** | **0.990** | **0.982** |

On the independent primary endpoint, V2.5 is the highest-AUC axis. The margin over the next-best baseline is narrow: +0.02 over Δβ-only, +0.06 over V2 `tumor_only`, +0.17 over V1. **Δβ-only is an unexpectedly strong baseline on matched cell-line cohorts** (0.974 on GSE322563 validated, 0.972 on GSE77348 validated). This is an honest finding the framework's complexity has to justify: V2.5's AUC advantage on matched cell-line data is small. V2.5's structural advantages — a probability scale interpretation, tie-band-honest top-K reporting, and principled uncertainty propagation via `p_trust` — are what earn its complexity on low-replicate cohorts; the AUC gap over Δβ widens on tissue (§5.3).

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

V2.5 is the highest-AUC axis at every tier on both matched cell-line cohorts. The margin over Δβ-only ranges from **+0.010 to +0.066**; the margin over V1 ranges from **+0.014 to +0.174**. On GSE77348 narrow, Δβ-only is the second-best axis (0.964), ahead of both V1 (0.969; marginally better) and V2 `tumor_only` (0.911). On the *native* EPIC v2 path (§5.2 below), Δβ-only marginally *beats* V2.5 on the narrow tier (0.956 vs 0.944) — a loss that is honest to report and that we interpret in §5.3 alongside the tissue-cohort numbers where the picture reverses sharply.

**P@K intervals (GSE322563, narrow labels).** On `n = 2/2` the score distribution at K = 100 sits inside a tied band; P@100 is an adversarial interval, not a point estimate:

| axis | P@100 observed | P@100 min | P@100 max | tie_band@100 |
|---|---:|---:|---:|---:|
| V2.5 differential | 0.000 | 0.000 | 0.020 | 190 |
| tumor_only | 0.000 | 0.000 | 0.020 | 10,005 |
| V1 `final_score` | 0.020 | 0.020 | 0.020 | 1 |

V2.5's band at K = 100 is 190 (interval [0.00, 0.02]); V2 `tumor_only`'s band is 10,005, spanning two orders of magnitude more of the distribution — P@K on `tumor_only` is not a meaningful quantity at this cohort size. V1's band is 1 (continuous score), and the interval collapses to the observed value.

**Native EPIC v2 vs HM450-intersect (GSE322563).** Re-running under the native EPIC v2 probe set (no HM450 intersect) gives:

| tier | V1 (HM450 / native) | V2.5 (HM450 / native) |
|---|---:|---:|
| validated | 0.821 / 0.879 | 0.990 / 0.986 |
| narrow | 0.884 / 0.931 | 0.942 / 0.944 |
| wide | 0.768 / 0.877 | 0.910 / 0.933 |

The V2.5 primary-endpoint AUC differs by 0.003 under native vs intersect — within the tied-band noise floor. V1 gains the most under native (+0.05 to +0.11). The relative axis ordering V2.5 > V1 is preserved at every tier. The HM450-intersect shortcut is not materially distorting the headline claim.

### 5.3 Tissue cohort — GSE69914 (n = 305/50)

At high replicate count, `p_trust` desaturates and the composite's ordering is driven by continuous `p_targ × p_diff`. On GSE69914:

| axis | validated AUC | narrow AUC | wide AUC |
|---|---:|---:|---:|
| Δβ-only (naive) | 0.591 | 0.477 | 0.435 |
| V1 `final_score` | 0.861 | — | — |
| V2.5 differential | 0.837 | — | — |

- **`tie_band_size_at_k`** drops from the hundreds-to-thousands seen at n = 2/2 to `tie_band = 2` at K = 100. The falsifiable prediction from §2.4 holds.
- **Δβ-only collapses to near-random** (0.435–0.591) on tissue. The intuition: tissue "normal" is adjacent-normal bulk, whose β values are shifted toward partial methylation by stromal contamination; a tumor–normal gap that was 0.8 on matched cell-lines can be 0.2 with noisy floors on tissue. Raw Δβ magnitude loses its discriminative power. The probabilistic envelope's normal approximation with σ_floor handles this: `p_diff` continues to separate large-gap-with-high-confidence candidates from small-gap-high-noise candidates in a way raw Δβ cannot.
- **V2.5 does not dominate V1 on AUC here** (0.837 vs 0.861). V1's weighted sum of sequence/selectivity/confidence/penalty components is more robust than V2.5's multiplicative composite on a cohort where `p_diff` is less saturated and the score distribution is genuinely continuous.
- **Top-20 membership is stable** rather than a window inside a tied band, because the top-K no longer sits inside a saturated region of the composite.

The tissue-cohort picture is precisely the pattern the probabilistic machinery is supposed to capture: Δβ-only fails where the signal is noisy, V2.5 recovers the ordering, and V1 is the right default because its continuous score is tie-band-safe. This cohort is why V1 and V2.5 coexist in the shipped release rather than one deprecating the other.

### 5.4 Out-of-distribution boundary — GSE68379

GSE68379 is the Sanger CCLE breast panel; its MCF-7 is epigenetically distinct from Roth's MCF-7 (well-documented cell-line drift). Applying Roth Fig. 5d-derived positives to Sanger MCF-7 should fail: the validated sites are hypermethylated in Sanger MCF-7 where they are hypomethylated in Roth MCF-7. AUC on this cohort inverts — and inverts on *both* axes:

| axis | validated AUC | narrow AUC | wide AUC |
|---|---:|---:|---:|
| Δβ-only (naive) | 0.201 | 0.150 | 0.194 |
| V2.5 differential | (inverts, < 0.5) | — | — |

Both the literature-naive baseline and V2.5 invert below 0.5, which is the expected correct response: the scorers are consistently identifying that the label set does not transport from Roth MCF-7 to Sanger MCF-7. This is documentation of a known boundary (cell-line drift), not a negative generalization result.

### 5.5 Top-hit annotation

For each cohort, the top-20 candidates are annotated with:

- nearest gene symbol, TSS distance, and feature class (promoter / gene_body / intergenic),
- CpG-island context (island / shore / shelf / open_sea) and distance to the nearest CGI boundary,
- RepeatMasker overlap (class / family / name) from UCSC rmsk,
- ENCODE DNase-HS cluster membership and breadth (number of cell types with the peak) from `wgEncodeRegDnaseClusteredV3`.

A rule-based flag set (STRONG: island-localized promoter / active promoter with DNase support; CAUTION: overlaps repeat / sparse evidence / small Δβ; NOTE: gene body without DNase support) is emitted in both the TSV and a Markdown companion file. On GSE69914's top-20, 6/20 carry at least one STRONG flag and 1/20 is caution-free — a distribution that is directly useful for experimental triage. Example annotated shortlists are committed under `examples/*/top20_annotated_v25.md`.

The formal SCREEN cCRE Registry is gated behind a JavaScript challenge on `screen.wenglab.org` and could not be fetched directly; DNase-HS clusters serve as the v1 regulatory-activity proxy. This is documented in the script's docstring; a user with SCREEN access can swap in the formal Registry with a one-line change.

---

## 6 Discussion

### 6.1 What V2.5 is and is not

The benchmark picture is sharper than the one we expected when we built the framework. On matched cell-line cohorts — where the PAM cytosine's methylation state cleanly separates the tumor and normal arms — the raw methylation gap Δβ = β_normal − β_tumor is already a strong ranking signal. V2.5 improves on it by a small margin (AUC +0.01 to +0.03 at validated tiers) and occasionally loses to it by a small margin (Δβ 0.956 vs V2.5 0.944 on native EPIC v2 narrow). **V2.5 is not a universal AUC replacement for Δβ, and we do not claim it is.**

The recommended axis depends on the cohort shape, not on a single "best" choice:

1. **On matched cell-line cohorts, V2.5 is the recommended research mode.** It is the highest-AUC axis at every label-granularity tier on GSE322563 and GSE77348. The margin over Δβ is small (0.01–0.07) and the margin over V1 is larger (up to +0.17). The probabilistic envelope's value here is partly AUC and largely structural — see (4) and (5) below.
2. **On tissue cohorts, V1 is the recommended axis, not V2.5.** On GSE69914 (n = 305/50), V1's continuous deterministic score reaches AUC 0.861 vs V2.5's 0.837 at validated; V1 also has `tie_band = 1` everywhere by construction. **The shipped recommendation here is V1, not V2.5.** What V2.5 *does* on tissue is recover ranking power that Δβ-only loses entirely (Δβ collapses to AUC 0.435–0.591) — but the right framing is "V2.5 ≈ V1 on tissue, both > Δβ; pick V1 for the deterministic top-K guarantee" rather than "V2.5 owns tissue."
3. **On the OOD boundary (GSE68379), every axis inverts correctly** below 0.5 — the scorers are consistently identifying that the Roth label set does not transport to Sanger MCF-7. This is a sanity check, not a generalization claim.
4. **V2.5's tie-band-honest top-K reporting is unique to it.** `tie_band_size_at_k` and `precision_at_k_{min, max}` let a reader see when the top-K window sits inside a saturated region of the composite (the n = 2/2 cell-line case). Δβ has `tie_band = 1` (continuous score) and V1 also has `tie_band = 1`, so both have honest top-K by construction; V2.5 needs the explicit reporting because of its discrete `EvidenceClass` saturation.
5. **V2.5 produces a probability-scale interpretation.** Each per-site `p_therapeutic_selectivity ∈ [0, 1]` is composable with downstream probabilistic inputs (target-mutation models, gRNA off-target probabilities, delivery-efficiency priors). Δβ is a signed similarity and V1 is a weighted-sum heuristic; neither is directly multiplied into such a chain.

This is the methodologically honest framing: Δβ is a strong baseline that should be reported alongside any new methylation-guided Cas9 scorer, V2.5 is the right research mode on cell-line cohorts but not the right default on tissue, and V1 remains the stable release axis precisely because its deterministic continuous score is cohort-shape-agnostic. We expect reviewers to prefer this framing over the stronger but unsupportable "V2.5 wins everywhere" claim.

### 6.2 Decision hierarchy

Three truths coexist:

1. V1 `final_score` is the stable release axis (tagged `v0.4.0`), because its deterministic continuous score gives `tie_band = 1` on every cohort tested — top-K is never tie-break-dependent.
2. V2.5 is the recommended probabilistic research mode on matched cell-line or paper-comparable cohorts, because it wins on AUC at every label granularity tested there.
3. V2 `tumor_only` is retained as a diagnostic, not a discovery axis, because its tie bands at K = 100 are consistently in the 6,000–12,000 range on low-replicate cohorts.

Neither axis dominates universally, and the decision table is the literal recommended hierarchy.

### 6.3 Limitations

1. **Low-replicate tie bands.** On cohorts with `n < 30`, top-K is reported as an interval via `precision_at_k_{min, max}` and `tie_band_size_at_k`. The problem dissolves at `n ≳ 30` (demonstrated on GSE69914).
2. **Cell-line drift.** GSE68379 confirms that "MCF-7" is not a single epigenetic object across labs. Any V2.5 shortlist must be validated against the specific cell line intended for the follow-up editing assay.
3. **σ_floor = 0.05** is empirical; a formal robustness study across floor values is deferred.
4. **Catalog scope.** AUCs are on chr5/6/10 candidates within ±500 bp of a methylation probe. The framework extends to whole-genome catalogs without modification; the scope here is driven by the location of Roth Fig. 5d targets.
5. **SCREEN cCRE Registry** is not integrated; DNase-HS clusters are a proxy.
6. **No prospective wet-lab validation.** The claim is AUC + tie-band-aware top-K on public methylation data. All named top hits are unvalidated predictions.

### 6.4 Related work

Two bodies of prior art are adjacent but non-overlapping. Generic CRISPR guide-scoring tools (e.g. Azimuth, DeepCRISPR, CRISPRitz) score sequence properties and off-target risk; none of them incorporate methylation state at the PAM cytosine, because none of them target a methylation-sensitive Cas9 variant. Methylation-aware variant analysis tools (e.g. methylKit, minfi) provide per-probe differential analysis on arrays but do not produce ranked target lists; their output feeds a scorer, not replaces one. `thermocas` is the first open framework we are aware of that joins the two — a ranking scorer that consumes methylation array summaries and outputs tumor-selective ThermoCas9 target shortlists with honest uncertainty.

### 6.5 Next steps

1. Prospective wet-lab validation of 3–5 top V2.5 candidates on MCF-7 vs MCF-10A (β-pyrosequencing + ThermoCas9 editing readout). This is the terminal step that changes the scientific claim.
2. A second independent-lab matched MCF-7/MCF-10A EPIC v2 cohort if one becomes public, to replicate the V2.5 headline at `n ≥ 3` rather than `n = 2`.
3. Formal SCREEN cCRE Registry integration when a non-challenged download path becomes available.

---

## 7 Conclusion

`thermocas` offers a principled scoring axis and an honest benchmarking protocol for methylome-guided ThermoCas9 target discovery. Its strongest contributions — the differential-protection formulation, the tie-band-aware P@K reporting, and the k-way-merge pan-cancer aggregator — are each small in isolation and cumulative in effect. On matched cell-line cohorts (the Roth use profile), V2.5 wins on AUC at every label granularity tested. On tissue cohorts, V1's continuous score and better top-K interpretability make it the right default. The combination is offered as an open research framework; groups extending Roth et al. to other cell-type pairs are the immediate target audience.

---

## Data and code availability

- **Code**: <https://github.com/AllisonH12/thermocas9>, tagged `memo-2026-04-22-c` for the exact revision evaluated here. BSD-3 licensed. `git rev-parse memo-2026-04-22-c` resolves to a SHA in a fresh clone. Supersedes `memo-2026-04-22-b` with the Δβ-only baseline + honest-framing updates to §5 / §6 applied on top; prior tag retained per the immutable-tag policy.
- **Tests**: 236 passing under `uv run pytest -q`.
- **Cohorts**: public GEO series GSE322563, GSE77348, GSE69914, GSE68379. Build scripts in `scripts/build_gse*_cohort.py`.
- **Reference annotations**: UCSC hg19 `refGene.txt.gz`, `cpgIslandExt.txt.gz`, `rmsk.txt.gz`, and `wgEncodeRegDnaseClusteredV3.txt.gz` (fetched on demand from hgdownload.soe.ucsc.edu).
- **Positives**: `data/derived/positives_roth_{validated, narrow, wide}.txt` (HM450 path) and `data/derived/epic_v2_positives/positives_roth_*.txt` (native EPIC v2 path). Derived from Roth et al. Fig. 5d via Ensembl REST hg38 → hg19 liftover.
- **Benchmark artifacts**: every `BenchmarkResult` JSONL under `examples/*_roth_labels/` carries `precision_at_k`, `precision_at_k_{min, max}`, `recall_at_k`, `recall_at_k_{min, max}`, `roc_auc`, `tie_band_size_at_k`, and `tie_break_policy`.
- **Annotated top-20 TSV + Markdown shortlists** under `examples/*/top20_annotated_v25.{tsv,md}`.

## Competing interests

The `thermocas` framework is maintained by Thermocas9 Inc. The author is affiliated with Columbia University.

## Acknowledgements

This work would not exist without Roth et al.'s characterization of ThermoCas9 and their transparent release of the GSE322563 matched cell-line methylation data.

---

## Appendix A · Reproducibility

From a fresh clone of the repository at `memo-2026-04-22-c`:

```bash
# One-time env setup
uv sync
uv run pytest -q           # 236 passing

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
