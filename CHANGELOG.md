# Changelog

All notable changes to the `thermocas` framework. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased] — submission-freeze post-`memo-2026-04-22`

> Submission-freeze cycle for the bioRxiv / *Bioinformatics* preprint.
> Stable release is still `v0.4.0`; cite `memo-2026-04-22-bd` for the
> current submission-freeze state. See `git log memo-2026-04-22..main`
> for the full per-commit history.

### Added

- **Δβ-only baseline** (`naive_selectivity = β_normal − β_tumor`)
  benchmarked across all cohort paths × positives tiers. Reported
  alongside V1, V2 `tumor_only`, and V2.5 differential in every
  results table. V2.5 beats Δβ-only on every matched-cell-line row
  (+0.010 to +0.080); on tissue, V2.5 (0.711–0.773) far exceeds
  Δβ-only (0.435–0.591) — where the probabilistic envelope earns
  its complexity.
- **`MANUSCRIPT.md`** — Bioinformatics-submission-shaped sibling
  to `PAPER.md`. Structured abstract, leads with the
  "no single axis dominates" framing.
- **Annotation pipeline**: `scripts/annotate_top_hits.py` gains
  `--rmsk` / `--dnase` / `--markdown` flags producing per-candidate
  triage cards with rule-based flags for experimental collaborators.
- **Pan-cancer aggregator**: `aggregate_streaming(...)` — k-way
  merge over pre-sorted cohort iterables with
  `O(N_unique_candidate_ids)` memory. `thermocas aggregate
  --streaming` flag added.
- **Submission-time guards**: `scripts/verify_manuscript_claims.py`
  cross-checks numeric / universal claims in MANUSCRIPT.md,
  PAPER.md, and README.md against committed bench JSONLs and
  `probabilistic.py` constants. Must pass before any dated memo
  tag is cut.
- **Deterministic PDFs**: `scripts/render_paper_pdf.sh` takes a
  source path and sources `**Date.**` from the MD; exports
  `SOURCE_DATE_EPOCH` so Typst embeds deterministic
  `CreationDate` / `ModDate`. Byte-identical re-renders verified
  within `pandoc 3.9` + `typst 0.14.2`.
- **Figures**: `fig2_auc_bars` (4-axis × 4-cohort-path × 3-tier
  cross-cohort AUC) and `fig3_topgene_heatmap` (5-column top-20
  gene-presence; cross-cell-line intersection *CELF2*, *XPNPEP1*;
  tissue cohort has zero gene overlap with any cell-line shortlist).

### Shipped recommendation

- **V1 as the overall package default** (backward compatibility +
  `tie_band = 1` by construction regardless of cohort shape).
- **V2.5-sigmoid (`tumor_plus_gap_sigmoid`) as the recommended
  probabilistic prioritization axis** (hypothesis generation; not
  a discovery claim absent prospective wet-lab validation) across
  every tested non-boundary cohort shape (matched cell-line at
  n = 2/2 to 3/3 and primary tissue at n ≳ 30/side). Ships in this tag as a first-class
  `probabilistic_mode` enum value with `sigma_fixed` cohort-YAML
  field + iff-semantics validators. PAPER.md §5.2.2 WG panel:
  matches V2.5-diff on matched cell-line AUC within 0.002 but
  strictly beats V2.5-diff on top-K usability (`tie_band@100 = 1`
  vs 421–1,493 under WG denominator), and improves tissue AUC by
  +0.05 to +0.08.
- **V2.5-diff (`tumor_plus_differential_protection`) retained for
  backward compatibility and AUC parity** with pre-ag scored
  JSONLs. Still a selectable mode, but no longer the recommended
  prioritization axis on any tested regime.
- **Δβ-only / Δβ_z / V1 / V2 `tumor_only` / limma-style moderated-t** retained
  as published baselines in §5.1 / §5.2.2.

---

## [Unreleased — earlier] — V2.5 experimental (on `main`, not tagged)

> Not a release. V2.5 ships as an opt-in mode on `main`. The stable default
> remains v0.4.0 (`tumor_only`). V1 `final_score` is still the recommended
> target-discovery ranking. See `V2_5_REVIEW.md` for the full self-review
> including the rank-tie finding that argued against tagging.

### Added

- **V2.5 differential-protection probabilistic mode** (opt-in).
  `probabilistic_mode: tumor_plus_differential_protection` on the cohort
  YAML selects a new composite `p_sel = p_targ × p_diff × p_trust`, where
  `p_diff = P(β_normal − β_tumor > δ)` is computed via an independent-normal
  approximation on tumor/normal summary statistics. Replaces the
  threshold-based `p_protected_normal` (which encoded a static
  "normal-is-methylated-above-0.5" assumption and was empirically
  anti-predictive on bulk comparators) with a differential factor that
  makes no absolute claim about the normal side — it asks only whether
  the normal-vs-tumor gap exceeds δ.
- **`differential_delta`** (default `0.2`) added to `CohortConfig`. Only
  consulted when the mode is `tumor_plus_differential_protection`.
- **`ProbabilisticScore.p_differential_protection`** and
  **`ProbabilisticScore.differential_delta`** fields — populated iff the
  mode is differential; a Pydantic validator enforces the iff both ways
  so malformed records cannot round-trip through JSONL.
- **`p_differential_protection(obs, delta, sigma_floor=0.05)`** public
  function in `thermocas.probabilistic`; exported from the top-level
  package.
- **`scripts/v2_5_differential_offline.py`** — the offline δ sweep that
  motivated the mode. Reads `scored_surrogate_tumor_only.jsonl` and
  reports AUC + P@100 for δ ∈ {0.2, 0.3, 0.4, 0.5} against both positive
  sets; no package-code changes.

### Benchmark label repair + cross-cohort empirical readout

**Label repair.** Roth Fig. 5d names three exact validated MCF-7/MCF-10A
target sites (`EGFLAM T11`, `ESR1 T17`, `GATA3 T18`) in hg38. These were
lifted to hg19 via Ensembl REST (`/map/human/GRCh38/.../GRCh37`) and
used to build three new positives files:

- `positives_roth_validated.txt` — one candidate_id per target (n=3)
- `positives_roth_narrow.txt` — NNNNCGA within ±50 bp (n=28)
- `positives_roth_wide.txt` — NNNNCGA within ±500 bp (n=142)

Our per-probe β values at the three target positions reproduce Roth's
Fig. 5d β values exactly (EGFLAM 0.01/0.49, ESR1 0.07/0.94,
GATA3 0.02/0.31) — cross-check confirming liftover + EPIC-v2 → HM450
intersect are both correct.

The earlier `positives_tight.txt` / `positives.txt` gene-symbol sets
(970 / 1687 entries) are retained as auxiliary gene-universe labels but
are no longer the authoritative supervision target: diagnostic on
GSE322563 showed only 23% of "Roth-gene" probes have `β_n − β_t > 0.2`
on the actual Roth samples.

**Cross-cohort matrix under repaired labels.** AUC
(`tumor_only` → `differential` → `V1`):

| cohort | regime | `n`/side | validated | narrow | wide | V2.5 tie_band |
|---|---|:---:|---|---|---|---:|
| GSE322563 | Roth cell lines | 2/2 | 0.928 / **0.990** / 0.821 | 0.886 / **0.942** / 0.884 | 0.871 / **0.910** / 0.768 | 190 |
| GSE77348 | MCF-7/MCF-10A surrogate | 3/3 | 0.912 / **0.982** / 0.968 | 0.911 / **0.983** / 0.969 | 0.887 / **0.949** / 0.931 | 299 |
| GSE69914 | primary tissue | 305/50 | **0.803** / 0.773 / 0.660 | **0.843** / 0.711 / 0.539 | **0.874** / 0.726 / 0.435 | **2** |

Bold = best AUC in that row. `tumor_only` tie_band is 6,540–11,848
everywhere. `V1` tie_band is 1 everywhere.

**Three findings:**

1. V2.5 wins AUC on matched cell-line cohorts by clear margins. Label
   repair lifted V2.5 on GSE322563 from 0.694 → 0.990 at validated;
   V1 from 0.541 → 0.821. Both were hurt by the old noisy labels; V2.5
   more so.
2. On primary tissue, `tumor_only` has higher AUC but unusable top-K
   (6,540-record tied band at K=100). V2.5 is second-best AUC and the
   only probabilistic axis whose top-K is interpretable on tissue.
3. V2.5's tie_band behaves correctly with `n`: 190 at n=2 → 299 at n=3
   → 2 at n=305/50, confirming the `p_trust`-saturation prediction from
   V2_5_REVIEW §3.

**Doc framing.** V2.5 is now the recommended probabilistic mode for
matched cell-line / paper-comparable biology, with low-`n` tie-band
caveats. `tumor_only` stays analysis-only. V1 stays available as a
continuous-score fallback for top-K stability across all cohort types.
V2.5 is **not** promoted to unconditional default; the recommendation
stays cohort-type-dependent.

GSE68379 (Sanger GDSC 52-breast-line panel × GSE69914 external healthy
normal, n=52/50) has also been run and is documented as an
**out-of-distribution boundary case** in `V2_5_REVIEW.md` §8.1 — not
a fourth generalization cohort. The three Roth-validated target probes
are methylated in Sanger's MCF-7 (β 0.63–0.92) and unmethylated in
Roth's MCF-7 (β 0.01–0.07), so the transferred Roth label set is not
biologically valid on this cohort; inversion of low-tumor-methylation
axes (V2.5, tumor_only) is the expected behavior under label–biology
mismatch and does not bear on V2.5's cohort-level claims. The
recommended-axis table is unchanged.

### Tests

- 197 → 203. Behavior-level contract tests cover mode dispatch, δ
  configurability, strict monotonicity in δ, composite formula equals
  `p_targ × p_diff × p_trust` exactly, audit-field consistency, YAML
  round-trip, `p_diff = 0.5` at the breakpoint, `ProbabilisticScore`
  composite-consistency validator (P2 fix), and deterministic tie-break
  + `tie_band_size_at_k` reporting in `evaluate_ranking` (P1 fix).

### Reporting / annotation

- `scripts/annotate_top_hits.py` gains optional `--rmsk` and `--dnase`
  flags. Output adds `in_repeat`, `repeat_class`, `repeat_family`,
  `repeat_name`, `in_dnase_cluster`, `dnase_source_count`. The DNase-HS
  clusters from `wgEncodeRegDnaseClusteredV3` are a stand-in for the
  ENCODE SCREEN cCRE Registry (which is gated behind a JS challenge on
  `screen.wenglab.org` and cannot be fetched directly).

### Pan-cancer aggregation

- New `aggregate_streaming(...)` (and `thermocas aggregate --streaming`)
  performs a k-way merge over pre-sorted cohort JSONLs. Candidate-side
  memory grows in `N_unique_candidate_ids` (a small `seen_cid → metadata`
  map, ~100 B per entry) rather than multiplying by `N_cohorts ×
  sizeof(ScoredCandidate)` as in the in-memory path. Both paths share
  the per-candidate emit logic; under valid input they produce
  byte-identical output.
- **Parity fix**: both `aggregate(...)` and `aggregate_streaming(...)`
  now raise `ValueError` on an intra-cohort duplicate `candidate_id`
  (previously: the in-memory path silently overwrote with the later
  record, while the streaming path rejected — a contract divergence
  under malformed input).
- **Behavior change**: in-memory `aggregate(...)` emission order now
  breaks ties on `(chrom, pos, family)` by `candidate_id` ascending
  (previously: dict-insertion order, i.e. cohort-iteration-order
  dependent). Strict improvement in determinism; only affects edge
  cases where two distinct `candidate_id`s share the same
  `(chrom, pos, family)` tuple.

---

## [0.4.0] — 2026-04-20 (V2.4 — `probabilistic_mode` config surface)

### Added

- **`probabilistic_mode` on `CohortConfig`** — explicit scoring policy
  for the V2 composite, with two initial values:
  - `tumor_only` (new default): `p_sel = p_targ × p_trust`
  - `tumor_plus_normal_protection` (opt-in, preserved): `p_sel = p_targ × p_prot × p_trust`
  The mode is recorded on every emitted `ProbabilisticScore` for
  downstream auditability. No silent behavior change on YAMLs that didn't
  specify a mode: they pick up the new `tumor_only` default, which
  matches the V3.1 published guidance.
- **Phase 5b factor ablation** (`scripts/phase5_v2_ablation.py`) —
  localized V2's rank-metric failure to `p_protected_normal` rather than
  to the multiplicative composition or to `p_trust`, by computing AUC for
  each factor, pair, and full composite on the MCF-7/MCF-10A surrogate.
  Results (loose AUC): `p_targ` 0.683, `v1_final` 0.657, `naive` 0.608,
  `p_targ × p_prot` 0.503, `p_prot` alone 0.384 (inverted). Archived
  at `data/derived/phase5_v2_ablation.txt`.

### Changed

- **Default probabilistic composite is `tumor_only`.** Cohorts that
  expect the old `p_targ × p_prot × p_trust` behavior must now declare
  `probabilistic_mode: tumor_plus_normal_protection` explicitly. The
  rationale is documented in README §"V2 scoring modes" and in the
  Phase 5b ablation report: `p_protected_normal` encodes
  `P(β_normal > 0.5)`, which is anti-predictive on cohorts where the
  normal comparator doesn't systematically methylate target promoters.

### Tests

- 182 → 197. New tests cover the mode-dispatch contract, that the mode
  field round-trips on every `ProbabilisticScore`, and that cohort YAMLs
  reject unknown mode strings at load time.

## [0.3.1] — 2026-04-20 (V3.1 — real-data validation + hardening)

### Added

- **Real-data benchmark.** Full pipeline against TCGA-BRCA (793 tumor + 97
  normal HM450 files, ~12 GB), catalog filtered to within 500 bp of HM450
  probes on chr5+chr6+chr10 (~3 M candidates), benchmarked against 1687
  gene-targeted positives (ESR1, GATA3, EGFLAM, VEGFA). LumA-only sub-cohort
  (285 files) via PAM50 subtyping from UCSC Xena. End-to-end results:
  - V1 `final_score` AUC 0.635 (bulk) / 0.656 (LumA)
  - V2 `p_therapeutic_selectivity` AUC 0.553–0.564 (underperforms)
  - `naive_selectivity` baseline AUC 0.597–0.629
- **`naive_selectivity` score axis** in `evaluate_ranking`. Ablation baseline
  = `β_normal − β_tumor`; lets you check whether the framework's machinery
  adds value over a 3-line predictor.
- **`--probe-annotation` + `--probe-window-bp` flags on `build-catalog`.**
  Drops ~95% of would-be UNOBSERVED candidates for genome-scale builds.
- **GDCBackend HTTP retries** (3 attempts, 1 s + 2 s backoff, 300 s timeout
  per call). Recovered from a multi-hour-kill on the first uncapped BRCA pull.
- **LumA summary script** (`scripts/build_luma_summary.py`) as analysis glue —
  builds a subtype-restricted summary TSV from cached GDC files.
- **Phase 4 diagnostics script** (`scripts/phase4_diagnostics.py`) —
  distribution analysis, top-N inspection, positives-set ablation, Roth
  coordinate-system verification.
- **Known-limitations section in README** explicitly documenting:
  - V2 `P(protected_normal)` encodes a "normal is methylated at target"
    assumption that breaks for adjacent-normal bulk benchmarks.
  - `aggregate()` memory cost is `O(N_candidates × N_cohorts)` (not streaming).
  - Positives-selection is coordinate-system-sensitive; Roth uses hg38.

### Fixed

Reviewer cycles over this release surfaced and fixed **16 real bugs** beyond
the initial V3 release. Grouped by theme:

- **GDC workflow integrity:** `LocalSummaryBackend` introduced so
  `gdc-fetch` outputs feed `score-cohort`; `GDCBackend` no longer claims
  the `MethylationBackend` interface it couldn't satisfy.
- **Benchmark contract:** `evaluate_ranking` actually filters by
  `held_out_chromosomes`; `missing_score_policy` parameter (default
  `rank_last`) so V2/V3 axis benchmarks don't silently shrink the
  evaluation set.
- **PAM model validation:** exhaustive ACGT⁴⁻⁸ enumeration rejects
  mixed-width regexes; duplicate family names rejected at load time.
- **Evidence-class defaults:** `exact_bp` tightened 1 → 0 (matching the
  "same CpG" contract); all config models switched to `extra='forbid'`
  so YAML typos error instead of silently defaulting.
- **Input integrity:** `read_beta_matrix` rejects duplicate probe_id
  rows AND duplicate sample IDs in the header; `_load_probes` rejects
  duplicate probe_ids in the annotation TSV; `split_by_subtype` refactored
  to reuse `read_beta_matrix` so it inherits all three checks.
- **Spacer extraction:** minus-strand candidates handled correctly
  (critical-C index is `len(seq)−1−fwd_idx` after reverse-complement,
  not `min(ccp, LOCAL_CONTEXT_HALF_WIDTH)`).
- **Summarization quartiles:** `method="inclusive"` + hard clamp to
  observed range so small subtype splits don't get extrapolated values.
- **Pan-cancer:** `normal_risk_max` renamed to `normal_protection_max`
  (it was computing best-case protection, not worst-case risk);
  `normal_protection_min` added for the honest worst-case view;
  `aggregate()` rejects metadata mismatch when same `candidate_id` maps
  to different `(chrom, pos, family)` across cohorts.
- **Atomic writes:** `write_jsonl_atomic` inserts `.tmp` *before* the
  final suffix so `.gz` outputs actually round-trip through gzip.
- **CLI hardening:** `inspect` streams instead of slurping; rejects
  `--top 0` at argparse time; `.gz` JSONL artifacts readable; Python 3.14
  argparse compatibility (escaped `%` in help strings).
- **Catalog streaming:** `stream_catalog` dropped its `seen_ids` set —
  candidate IDs are deterministic from coordinates and PAM family names
  are now validated unique, so the defensive state was growing with
  candidate count and violating the streaming contract.

### Tests

- 144 → 182 (38 new regression tests, 0.9 s total runtime).

## [0.3.0] — 2026-04-20 (V3)

### Added

- **Beta-distribution CDF** for the probabilistic decomposition. Method-of-moments
  fit from (mean, IQR-derived σ); regularized incomplete beta `I_x(α, β)`
  computed via Lentz's continued fraction (stdlib only). The V2 piecewise-linear
  estimator is preserved as the fallback for degenerate cases.
- **gRNA spacer scoring** (`thermocas.grna`). New `SpacerScore` model captures
  GC content, Wallace-rule Tm, longest mononucleotide run, and a 4-bp hairpin
  proxy. Optional via `score_candidate(..., compute_spacer=True)` and
  `thermocas score-cohort --spacer`.
- **Cross-validation benchmark harness** (`thermocas.benchmark`). `split_by_chrom`
  for held-out-chromosome splits; `evaluate_ranking` for top-K precision /
  recall and ROC-AUC (Mann-Whitney U formulation). New `BenchmarkResult` model
  + `thermocas benchmark` subcommand. Score axis selectable via `--score-field`
  (deterministic, probabilistic, or spacer-based).
- **`thermocas inspect`** subcommand: pretty-print summary stats for any JSONL
  artifact (catalog / scored / aggregate / benchmark) with auto-detection.
- `git` history initialized; CI workflow (`.github/workflows/ci.yml`);
  `CONTRIBUTING.md`.

### Fixed

- `PamFamily` validator now uses **exhaustive ACGT enumeration** (lengths 1–8)
  to determine motif width instead of a hard-coded probe corpus, so legitimate
  custom PAMs like `CCCCCGA` validate. Mixed-width regexes (e.g.
  `A|[ACGT][ACGT][ACGT][ACGT]CG[AG]`) are now rejected up front.

### Tests

- 144 → 147 (Beta math, gRNA scoring, benchmark metrics, CLI inspect).
- All assertions remain green in 0.75s.

## [0.2.0] — 2026-04-20 (V2)

### Added

- **Probabilistic targetability scoring** (`thermocas.probabilistic`):
  `P(targetable_tumor) × P(protected_normal) × P(observation_trustworthy)` →
  `ProbabilisticScore`. Optional via `score_candidate(..., compute_probabilistic=True)`
  and `thermocas score-cohort --probabilistic`.
- **`LocalArrayBackend.split_by_subtype()`** factory: splits a tumor matrix
  using a `sample_id, subtype` TSV into per-subtype backends. CLI exposes this
  as `score-cohort --sample-subtypes`, fanning out to one JSONL per subtype
  with `cohort_name = '<cohort>::<subtype>'`.
- **Real `GDCBackend`** over stdlib `urllib` with on-disk file cache. New
  `thermocas gdc-fetch` subcommand exports LocalArrayBackend-compatible
  summary TSVs from a TCGA cohort.

### Fixed

- `PamFamily` rejects malformed configs (`regex='A', critical_c_offset=4`)
  at YAML load instead of silently producing out-of-bounds coordinates.
- `ScoredCandidate` validator now also asserts
  `probabilistic.cohort_name == observation.cohort_name`.
- `_piecewise_linear_cdf` mean-only fallback is a sharp step at the mean;
  no longer fabricates a smooth linear CDF from a lone point estimate.

### Tests

- 28 → 132. Adds Beta math regressions, V2 backend mocks, subtype splitting,
  CLI flag coverage, model invariant regressions.

## [0.1.0] — 2026-04-19 (V1)

### Added

- Initial framework scaffolding.
- Pydantic v2 data model: `CandidateSite`, `MethylationObservation`,
  `EvidenceClass`, `ScoreComponents`, `ScoredCandidate`, `PanCancerAggregate`,
  `CohortConfig`, `Penalties`, `EvidenceThresholds`.
- Deterministic scoring schema:
  `final_score = sequence × selectivity × confidence − penalties`,
  with selectivity rewarding both mean and quantile separation.
- PAM-model loader + matcher with overlap-safe lookahead regex.
- Probe-to-site evidence classification (exact / proximal_close / proximal /
  regional / unobserved).
- Genome-wide PAM catalog builder (FASTA streamer, gzip-aware,
  strand-correct local context).
- `MethylationBackend` ABC + `LocalArrayBackend` over TSVs.
- Cohort adapter (`score_cohort`).
- Pan-cancer aggregator (recurrence, exclusivity, normal-tissue risk,
  + `top_recurrent` / `top_exclusive` helpers).
- CLI: `build-catalog`, `score-cohort`, `aggregate`.
- Tests for every public function (28 in initial release).
