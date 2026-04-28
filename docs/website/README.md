# Website-ready drafts for thermocas9.com

This directory holds publication-ready prose drafts intended to be lifted
to <https://thermocas9.com>. They are version-controlled here so the
canonical text stays in lockstep with the framework tags.

| File | Purpose | Target audience | Effort to publish |
|---|---|---|---|
| `01_preprint_landing.md` | Landing page for the bioRxiv preprint + repo + Zenodo DOI | Bioinformatics / methylation / Cas-engineering readers arriving from a citation or a tweet | Drop-in: copy + paste, swap in the bioRxiv DOI when minted |
| `02_method_explainer.md` | Plain-language methods walkthrough — what the compositional skeleton is, why tie-band-aware benchmarking matters | Researchers who will not read a 100-page memo but want to know what we built | Light edit: pick the headline figure and finalize the disclaimer block |
| `03_atlas.md` | Per-cohort target-shortlist atlas — fig4 per-positive WG-rank dot-plot + per-cohort top-100 tables linked to the underlying TSVs / JSON | Visitors who want to see *what the framework actually outputs* on the four publication cohorts, with the ESR1 reversal honest-disclosure visible inline | Light edit: wire `atlas_top100.json` into an interactive Astro / Observable Plot table component |
| `06_reproducibility_tutorial.md` | Step-by-step `uv run` walkthrough that reproduces the GSE322563-native primary endpoint end-to-end on a laptop | Methods-paper readers who want to verify a row of the AUC table or score a custom cohort | Drop-in: code blocks are copied from the verified Reproducibility Appendix |

Atlas data files (consumed by `03_atlas.md`):

| File | Size | Schema |
|---|---:|---|
| `atlas/per_positive_wg_percentile.json` | ~5 KB | 4 cohorts × 3 positives × 2 axes; rank, denominator, percentile |
| `atlas/atlas_top100.json` | ~200 KB | 4 cohorts × 100 candidates each, 15 columns (rank, candidate_id, chrom, position, strand, PAM, score, Δβ, p_trust, gene, TSS distance, feature class, CGI context, is_positive) |

Boundary rules that apply to every page on this directory:

1. The framework is **educational research**, not clinical decision support. Every page must preserve that boundary either inline or via a footer block.
2. Any numerical claim is sourced from the immutable tag (`paper-5-10j` for `PAPER.md`-rooted claims, `memo-2026-04-22-bw` for `MANUSCRIPT.md`-rooted claims). When the tag advances, re-check the numbers in these drafts before re-publishing.
3. Roth et al. *Nature* 2026 (DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z)) is cited on every page that mentions the underlying biology.
