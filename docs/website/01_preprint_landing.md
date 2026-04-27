# `thermocas` — methylome-guided ranking of ThermoCas9 target sites

**Status.** Open educational research framework. Not peer-reviewed. Not a clinical decision-support system. Hypothesis generation only — no claim of prospective wet-lab validation absent the editing study described in [§6.3 / §6.5 / `external_validation_instruction.md`](https://github.com/AllisonH12/thermocas9/blob/memo-2026-04-22-bu/docs/notes/external_validation_instruction.md).

---

## Preprint

**Title.** *Compositional probability-scale scoring and tie-band-aware benchmarking for methylome-guided ThermoCas9 target-site ranking.*

**Author.** Allison Huang, Columbia University. <allisonhmercer@gmail.com>

**Date.** 2026-04-22.

**Versions.**

- **`PAPER.pdf`** — long-form technical memo (audit-trail revision, includes V2 → V2.5 development history and Roth System B HEK293T/HCT116 appendix). Cite tag `paper-5-10j`. [Download PDF](https://github.com/AllisonH12/thermocas9/blob/paper-5-10j/PAPER.pdf).
- **`MANUSCRIPT.pdf`** — *Bioinformatics*-submission-shaped short version. Cite tag `memo-2026-04-22-bu`. [Download PDF](https://github.com/AllisonH12/thermocas9/blob/memo-2026-04-22-bu/MANUSCRIPT.pdf).
- **bioRxiv DOI**: *to be added on posting*.
- **Zenodo DOI** (citable archive of the tagged source): *to be minted on first GitHub release*.

The two PDFs share the same code, benchmark artifacts, and figures; they differ in length, framing, and target audience.

## What the paper is about

Roth et al. (*Nature* 2026, [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z)) showed that methylation of the fifth PAM cytosine protects genomic loci from ThermoCas9 cleavage, opening a route to tumor-selective editing at loci whose methylation state diverges between cancer and matched normal tissue.

`thermocas` turns that mechanism into a target-prioritization ranking. The contribution is split into two parts:

1. **A compositional probability-scale scoring skeleton** — `p_targ × (gap factor) × p_trust`. Each of the three factors is independently replaceable; the choice of gap factor is an instance, not a load-bearing methodological claim. Two gap-factor instances ship in this tag (V2.5-diff and V2.5-sigmoid), with V2.5-sigmoid as the worked default.
2. **A tie-band-aware benchmarking contract** — every `BenchmarkResult` row reports `precision_at_k_{min, max}`, `tie_band_size_at_k`, and a mid-rank Mann–Whitney AUC, so a reader cannot be misled by arbitrary tie-break choices on cohorts where many candidates carry identical scores.

## What the paper is *not*

- It does **not** include prospective wet-lab editing validation. The closest analytical evidence is a 1,000,000-draw random-triple permutation null on the three Roth Fig. 5d positives.
- It does **not** report per-site *p*-values; `p_observation_trustworthy` saturates by `EvidenceClass` rather than continuously. The composite is a ranking axis, not a calibrated hypothesis test.
- It does **not** make a uniform-superiority claim across all restricted universes. In the GSE69914 `EXACT + PROXIMAL_CLOSE` tissue subset (where ESR1 is the only evaluable positive), V2.5-sigmoid trails V2.5-diff, raw Δβ-only, and the limma-style baseline.

## Code, data, and tests

- **Repository**: <https://github.com/AllisonH12/thermocas9> · BSD-3 · Python 3.11+
- **Tests**: 245 passing under `uv run pytest -q`
- **Cohorts**: public GEO series GSE322563, GSE77348, GSE69914, GSE68379. Build scripts produce the per-probe summary TSVs deterministically from the GEO supplementary files.
- **Verifier**: `scripts/verify_manuscript_claims.py` cross-checks selected numerical claims against the committed bench JSONLs and `probabilistic.py` constants. Must pass before any dated memo tag is cut.
- **Reproducibility tutorial**: see [the laptop walkthrough](./06_reproducibility_tutorial.md) on this site.

## How to cite

Until the bioRxiv DOI is live, cite the immutable git tag:

> Huang, A. (2026). *Compositional probability-scale scoring and tie-band-aware benchmarking for methylome-guided ThermoCas9 target-site ranking.* Technical memo, version `paper-5-10j`. <https://github.com/AllisonH12/thermocas9/tree/paper-5-10j>

For the *Bioinformatics*-shaped short version, cite `memo-2026-04-22-bu` instead.

## Contact

Methods or reproducibility questions: <allisonhmercer@gmail.com>. For wet-lab collaboration on the prospective external-validation study described in `docs/notes/external_validation_instruction.md`, the same address.
