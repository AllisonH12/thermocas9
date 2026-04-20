# Changelog

All notable changes to the `thermocas` framework. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
