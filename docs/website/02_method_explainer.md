# Ranking methylation-protected ThermoCas9 targets — what `p_targ × gap × p_trust` actually does

*A plain-language walkthrough of the `thermocas` scoring framework. Audience: bioinformaticians and Cas-engineering researchers who want the gist without reading a 100-page memo.*

---

## The biology, in one paragraph

Roth et al. (*Nature* 2026) showed that ThermoCas9 — a thermophilic Cas9 from *Geobacillus thermodenitrificans* T12 — refuses to cleave when the fifth cytosine of its PAM is methylated. An unmethylated PAM cleaves as usual; a 5-methylcytosine at that position abolishes binding. The therapeutic implication is that loci which are unmethylated in tumor cells but methylated in matched normal cells become *selectively editable* — the editor only cuts the disease side. Roth demonstrated this with three breast-cancer targets in MCF-7 vs MCF-10A.

Turning that mechanism into a genome-scale shortlist is a ranking problem. Which loci, across millions of candidates, look like the best bets?

## Why a generic guide-scorer does not work

Off-the-shelf CRISPR guide-scoring tools (Azimuth, DeepCRISPR, CRISPRitz) score sequence properties and off-target risk. They have no notion of methylation state at the PAM cytosine, because they were not designed for a methylation-sensitive Cas9 variant. Methylation-array differential-analysis tools (methylKit, minfi, limma) produce per-probe statistics but not ranked target lists; their output feeds a scorer, it does not replace one.

Three problems, specifically, made it worth building a new scorer:

1. **Replicate counts are tiny.** Methylation cohorts often ship with `n = 2–3` samples per side. Any scorer must carry its own uncertainty rather than emit a scalar that hides it.
2. **The "normal" arm is inconsistent.** Adjacent-normal tissue in a bulk methylation study is not the same object as a matched cell-line pair. A scorer that bakes in a fixed β_normal threshold will fail across cohort types.
3. **Top-K can sit inside a tied band.** On low-replicate cohorts the top of any composite can collapse into hundreds or thousands of equally-scored candidates. Reporting a single Precision@K without that interval misleads the reader.

## The compositional skeleton

`thermocas` answers all three with the same multiplicative form:

```
p_therapeutic_selectivity  =  p_targ  ×  (gap factor)  ×  p_trust
```

Each factor maps to one independently-replaceable question.

**`p_targ` — *is the tumor side actually targetable?*** This is the probability that the tumor-arm methylation at the PAM cytosine is *low enough* to permit cleavage. It is computed by integrating a method-of-moments Beta posterior on the per-probe summary β, with a piecewise-linear fallback when the moment match degenerates. A site cannot be a candidate if the editor cannot bind to it.

**Gap factor — *is the tumor-vs-normal contrast meaningful?*** Two instances ship in the current tag:

- **V2.5-diff** (`tumor_plus_differential_protection`) — `p_diff = P(β_normal − β_tumor > δ)` under an independent-normal approximation on the per-probe β summaries with `δ = 0.20` and `σ_floor = 0.05`. The original V2.5 composite, retained for backward-compatibility and AUC parity on cell-line cohorts.
- **V2.5-sigmoid** (`tumor_plus_gap_sigmoid`) — `sigmoid((β_n − β_t − δ) / σ_fixed)` with `δ = 0.20` and `σ_fixed ≈ 0.0707 ≈ √2 · σ_floor`. A fixed-bandwidth logistic on the same gap signal. The worked / default whole-genome prioritization axis selected in the §5.2.2 frozen WG panel.

The point is not which gap factor "wins" — the point is that the slot is *replaceable*. A future axis (an exact difference-of-Betas, an SE-on-mean variant, a mixture-aware version) can be plugged in without rewriting the scorer.

**`p_trust` — *is the per-probe observation trustworthy?*** Discrete `EvidenceClass` (EXACT, PROXIMAL_CLOSE, PROXIMAL, REGIONAL, UNOBSERVED) times a sample-count ramp `min(1, min(n_t, n_n) / 30)`. A site whose methylation is nominated by a regional proxy on a 2-sample cohort cannot dominate the top-K just because its `p_targ × gap` happens to peak.

The multiplicative form is a **gating-style ranking heuristic**: a candidate is penalized if any required component is weak. An additive score over the same signals would let a strong component compensate for a failed gate, which is less aligned with a selectivity screen. The product is *not* a calibrated joint probability — `p_targ`, the gap factor, and `p_trust` are correlated on real catalogs (EXACT-class records cluster at extreme β values that also drive `p_targ`). Read the score as a probability-scale ranking axis throughout, not as "0.9 means 90% of 0.9-scored sites edit successfully."

## Tie-band-aware benchmarking — the actual durable contribution

This is the part of the paper most likely to outlive any specific gap-factor choice.

On a `n = 2/2` matched cell-line cohort with a probabilistic composite, the score distribution at K = 100 routinely sits inside a tied band. V2.5-diff at the whole-genome scale produces tied bands of **421 to 1,493 records** at K = 100 across the four evaluated cohort paths — a single Precision@100 number for that score is meaningless without an interval.

Every `BenchmarkResult` row that `thermocas` produces carries:

- `precision_at_k` (the natural-order observed value);
- `precision_at_k_{min, max}` (the adversarial interval over the tied band);
- `recall_at_k` and `recall_at_k_{min, max}`;
- `tie_band_size_at_k` (so the reader sees the denominator);
- a mid-rank Mann–Whitney `roc_auc` (no tie-break ambiguity);
- the `tie_break_policy` used.

Switching from V2.5-diff to V2.5-sigmoid does not change matched-cell-line AUC by more than 0.002, but it collapses every WG `tie_band@100` to **1**. That is a strict improvement in top-K usability on the same biological signal — and it is only visible because the benchmark contract reports tie bands explicitly.

## What the benchmark says

Four public methylation cohorts, scored on three positives tiers (`validated` from the Roth Fig. 5d targets; `narrow` ±50 bp; `wide` ±500 bp).

- **GSE322563 (Roth MCF-7/MCF-10A)** is the primary endpoint. On both the HM450-intersect and native EPIC v2 ingest paths, V2.5-diff and V2.5-sigmoid place all three Roth Fig. 5d positives in the upper 0.06–4.6% of millions of WG candidates.
- A **1,000,000-draw random-triple null** on the n = 3 primary endpoint places V2.5-diff at one-sided `p ≈ 4 × 10⁻⁶ / 1.5 × 10⁻⁵ / 2.8 × 10⁻⁵` across the three matched-cell-line cohort paths — strictly tighter than raw Δβ-only's `8.6 × 10⁻⁵ / 3.1 × 10⁻⁴ / 1.1 × 10⁻⁴` on the same paths. This is an against-random check, not a paired superiority test (the paired sign-flip null floors at `p = 0.125` with `n_pos = 3`).
- The **headline tissue stress-test** is on GSE69914 (n = 305 primary tumor + 50 healthy-donor breast tissue, unpaired) under transported Roth labels. A probe-level limma-style moderated-`t` baseline is competitive on cell lines (AUC 0.959–0.991) but **collapses to AUC 0.573 on bulk tissue**. V2.5-sigmoid recovers tissue to **AUC 0.862** on the same probe-level inputs — a +0.29 swing that confirms the `p_targ × p_trust` *wrapping* is the load-bearing contribution, not the gap factor's specific shape. The cohort-level tissue gain is per-positive heterogeneous: GATA3 carries the matched-pool signal, ESR1 is matched-near-random under V2.5-sigmoid.

## Honest scope

A few things the paper is explicit about and that we want a website reader to see:

- The validated positives set is **n = 3**. AUC at this scale is a *rank-lift summary*, not an inferential discovery-performance estimate.
- V2.5-sigmoid was selected by post-repair sensitivity and whole-genome stress testing on the same benchmark family. Its prospective ranking utility remains to be validated on newly generated labels or wet-lab follow-up. The detailed external-validation design lives in `docs/notes/external_validation_instruction.md`.
- On matched cell lines, **raw Δβ-only is already strong** (AUC 0.96–0.97). V2.5 variants add small but consistent margins (+0.010 to +0.080 across the nine tier × path rows). The *transferable* contribution is the tie-band-aware benchmark contract and whole-genome tied-band diagnosis, not a claim that any one gap factor is universally best.

## Where to go next

- **Read the long version**: [`PAPER.pdf`](https://github.com/AllisonH12/thermocas9/blob/paper-5-10j/PAPER.pdf) (`paper-5-10j`).
- **Read the journal-shaped version**: [`MANUSCRIPT.pdf`](https://github.com/AllisonH12/thermocas9/blob/memo-2026-04-22-bw/MANUSCRIPT.pdf) (`memo-2026-04-22-bw`).
- **Reproduce a row of the AUC table on a laptop**: see the [reproducibility tutorial](./06_reproducibility_tutorial.md).
- **Discuss a wet-lab collaboration**: <allisonhmercer@gmail.com>.

---

*Educational research framework. Not peer-reviewed. Not a clinical decision-support system. Cites Roth et al., Nature 2026, [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).*
