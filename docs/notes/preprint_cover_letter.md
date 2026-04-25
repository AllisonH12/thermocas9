# Preprint cover letter — draft

**Status.** Draft. Edit and customize for the target venue before sending.
**Submission ladder.**
1. **bioRxiv** first (no cover letter required; the body below doubles as the preprint pitch). Free, citable, journal-compatible for most biomedical venues.
2. **Bioinformatics (Oxford)** as the primary peer-reviewed target — scope is explicit about new computational methods compared to state of the art on real biological data, which fits.
3. **Bioinformatics Advances** or **BMC Bioinformatics** as the acceptance-optimized backup (explicitly welcoming to benchmarked methods/software contributions without prospective-validation expectations).
4. **The CRISPR Journal** only if a field-specific audience is preferred over broader computational-methods readership — narrower but higher topical familiarity.
5. **Nature Methods** only if the smallest-possible prospective validation is added first (see variant section at the bottom). Without that, expect desk reject.
**Manuscript reference.** `PAPER.md` (or `MANUSCRIPT.md` for the Bioinformatics-shaped version) in <https://github.com/AllisonH12/thermocas9> at tag `memo-2026-04-22-bn`.

---

## Universal cover-letter body (works for a generalist methods venue)

Dear Editors,

We are submitting our manuscript, *"Compositional probability-scale scoring and tie-band-aware benchmarking for methylome-guided ThermoCas9 target-site ranking"*, for consideration.

**What the paper contributes.** The recent description of ThermoCas9 (Roth et al., *Nature* 2026) established that methylation of the fifth PAM cytosine protects genomic loci from cleavage, opening a route to tumor-selective targeting at loci whose methylation state diverges between cancer and matched normal tissue. Translating that mechanism into a target-prioritization ranking (hypothesis generation; not a discovery claim absent prospective wet-lab validation) requires a scorer that (a) is cohort-agnostic, (b) is honest about evidence quality, and (c) reports ties explicitly rather than pretending to rank-order indistinguishable candidates. No such scorer was publicly available when we started this work.

The paper's core contribution is a **compositional probability-scale scoring skeleton** paired with a **tie-band-aware benchmarking contract** specialized to the methylation-sensitive-Cas9 ranking problem. The skeleton is
`probability-scale selectivity score = p_targetable_tumor × (gap factor) × p_observation_trustworthy`,
where each factor is independently replaceable — the specific choice of gap factor is an instance, not a load-bearing methodological claim. Two gap-factor instances ship as selectable `probabilistic_mode` enum values in this tag: **V2.5-diff** (`tumor_plus_differential_protection`, the original V2.5 composite with `p_diff = P(β_normal − β_tumor > δ)`) and **V2.5-sigmoid** (`tumor_plus_gap_sigmoid`, a fixed-bandwidth logistic). **V2.5-sigmoid is the recommended probabilistic prioritization axis**, chosen on the basis of a frozen whole-genome panel (§5.2.2) that shows it matches V2.5-diff's AUC on matched cell-line cohorts within 0.002 but strictly beats V2.5-diff's top-K usability (`tie_band@100 = 1` vs 421–1,493 under the WG denominator) and improves tissue AUC by +0.05 to +0.08. V2.5-diff is retained for backward compatibility and AUC parity with earlier scored JSONLs. A probe-level limma-style moderated-t DMR baseline (Smyth 2004 empirical-Bayes variance prior, sample-level β matrices, pure-Python implementation with R `limma::lmFit + eBayes` parity reported in §5.8) is included alongside and collapses to near-random on tissue (0.573) where V2.5-sigmoid reaches 0.862 — confirming that the `p_targ × p_trust` wrapping, not the gap-factor shape, is the load-bearing contribution. (The tissue stress test is under transported Roth labels, not tissue-specific ThermoCas9 validation; GATA3 drives the matched-pool signal while ESR1 is matched-near-random under V2.5-sigmoid in the §5.9 feature-matched controls.) The scorer is benchmarked on four public methylation cohorts (GSE322563 matched cell lines as the independent primary endpoint on both HM450-intersect and native EPIC v2 ingest paths, GSE77348 as the development cohort on which the differential threshold δ was tuned, GSE69914 as 305 sporadic breast tumors + 50 healthy-donor breast tissue samples — unpaired by design with adjacent-normal arms excluded by the cohort build script, GSE68379 as the out-of-distribution boundary) against a 3-tier positives set derived from the Roth paper's validated targets. We also ship:

- a `tumor_only` analysis-only path that avoids the β_normal prior when the normal arm is unreliable,
- tie-band-aware Precision@K and Recall@K reporting (P@K_min, P@K_max, tie_band_size_at_k) so a reader cannot be misled by arbitrary tie-break choices,
- an annotation pipeline that attaches nearest gene, CpG-island context, RepeatMasker overlap, and ENCODE DNase-HS cluster breadth to each shortlisted candidate, plus a Markdown companion to the TSV aimed at experimental collaborators,
- a streaming k-way-merge pan-cancer aggregator with cross-cohort metadata-parity enforcement, so the framework scales to genome-scale atlas builds without loading every cohort into memory.

Benchmark `BenchmarkResult` JSONLs, positives lists, annotated top-20 TSV + Markdown shortlists, figures, and the test suite are all committed at the immutable tag `memo-2026-04-22-bn`. The large per-cohort scored-candidate JSONLs (`data/derived/scored_*.jsonl`, tens of millions of records each) are gitignored but fully reproducible from the committed build scripts and cohort YAMLs — the exact `uv run` invocations are listed in the manuscript's reproducibility appendix. 245 unit tests pass.

**What the paper does NOT claim.** We want to be explicit about scope. This is a *methods and benchmarking* paper on public data. It does not include prospective wet-lab validation of new target sites; the `thermocas` framework is an open educational research tool, not a clinical decision-support system. Per-site p-values are not reported because `p_observation_trustworthy` saturates by `EvidenceClass`, not continuously; what `p_therapeutic_selectivity` provides is a defensible *ranking axis*, not a hypothesis test.

**Why we think this venue.** The paper's core contribution is methodological — a probability-scale scoring formulation, a benchmarking protocol with honest tie-band reporting, and an open framework with audit-level test coverage — rather than a biological discovery. The target-prioritization problem it addresses is of immediate interest to the small but growing methylation-sensitive-CRISPR community catalyzed by the Roth paper. We believe your readership is well-placed to evaluate the rigor of both the formulation and the benchmarking.

**Competing interests.** The author declares no competing interests.

**Authorship & data availability.** Allison Huang is the sole author and corresponding author. All cohort data are from public GEO accessions (GSE322563, GSE77348, GSE69914, GSE68379). Benchmarks (`BenchmarkResult` JSONLs), positives lists, annotated top-20 TSV + Markdown shortlists, figures, scripts, and the test suite are committed at the tagged revision; the large per-cohort scored-candidate JSONLs (`data/derived/scored_*.jsonl`, tens of millions of records each) are gitignored but fully reproducible from the committed build scripts.

**Suggested reviewers** *(fill in per venue)*. Reviewers who would be able to evaluate this work well fall into three non-overlapping areas: methylation array bioinformatics, CRISPR target-ranking / gRNA-scoring benchmarking methodology, and probabilistic model validation with class-imbalanced ground truth.

We appreciate your consideration.

Sincerely,
Allison Huang
Columbia University
<allisonhmercer@gmail.com>

---

## *Nature Methods* variant — additional framing if pursuing that ceiling

Append before "Suggested reviewers":

> **On validation.** We recognize that a *Nature Methods* reader will reasonably ask for prospective wet-lab validation of predicted target sites in a cell line not used for scorer development. We are explicit in the manuscript (§6) that this work does *not* include such validation and is offered as a methods and benchmarking contribution on public data. If the editors believe validation is a precondition for consideration at this venue, we are preparing the smallest-possible prospective cohort (3 predicted-tumor-hypo / 2 predicted-normal-methylated guides in MCF-7 vs MCF-10A, β-pyrosequencing readout) as a companion submission; we would welcome editorial guidance on whether to hold the present paper until that reads out or submit as-is with the validation as a follow-up.

---

## Notes for the sender

- **If submitting to bioRxiv only**: no cover letter is required. Use the first four paragraphs of the body as the bioRxiv abstract prompt and paste them into the summary field. The claim-narrowing paragraph ("What the paper does NOT claim") is especially useful here — bioRxiv readers see method papers from anonymous-ish accounts all the time, and explicit scope-limiting reads as professional, not defensive.
- **Timing**: the Roth follow-up note (`roth_followup_2026-04-22.md`) carries more weight if it goes *after* the preprint is up, so readers of the follow-up can cite a Google Scholar result rather than a GitHub tag. Sequence: preprint first, Roth note second, same week.
- **Tag-reference check before submission**: `git rev-parse memo-2026-04-22-bn` must still resolve, and both `PAPER.md` and `MANUSCRIPT.md` citations must still read `memo-2026-04-22-bn`. Do not let development on `main` implicitly move the reference.
- **Reviewer suggestions**: the three fields listed under "Suggested reviewers" are generic by design — fill in specific names per venue from recent program committees / corresponding authors in each area.
- **Revision checkpoint**: if this letter sits more than a week before submission, re-verify (a) the tag, (b) the test count (`uv run pytest -q`), and (c) the set of committed benchmark artifacts under `examples/`. Update the numbers if they've moved.
