# Changelog

All notable changes to the `thermocas` framework. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
