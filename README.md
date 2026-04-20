# thermocas — methylome-guided ThermoCas9 target-site discovery framework

A modular, pan-cohort ranking system for identifying genomic loci where ThermoCas9
is **likely to cut tumor DNA but not matched normal DNA**, based on local PAM-cytosine
methylation. Built around the mechanism reported in
Roth et al., *Nature* (2026), DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

## The central abstraction

```text
candidate_site
  = sequence-compatible PAM
  + inferred methylation state at the critical PAM cytosine
  + tumor / normal contrast
  + confidence of observation
```

The framework operates in four layers:

1. **Site enumeration** — genome-wide catalog of ThermoCas9-compatible PAM instances.
2. **Methylation evidence mapping** — assign each candidate an evidence class (exact / proximal / regional / unobserved) by relating it to assayed CpG probes.
3. **Cohort-specific scoring** — compute targetability per cancer type from tumor-vs-normal methylation summaries.
4. **Cross-cohort aggregation** — pan-cancer recurrence, exclusivity, and risk metrics.

## Why this design

The reported ThermoCas9 mechanism is **PAM-centric**: the fifth-position cytosine of the
`5'-NNNNCNR-3'` PAM family — including `5'-NNNNCGA-3'` and `5'-NNNNCCA-3'` — is the
methylation-sensing residue. Methylation of that cytosine (5mCpG or 5mCpC) blocks DNA
binding; protospacer methylation has little effect.

GDC-harmonized methylation data are reported as beta values for **probe-associated
CpG sites**, not for every cytosine in the genome. The framework therefore models
methylation evidence as a first-class object with explicit confidence classes.

## Repository layout

```text
thermocas-framework/
├── pyproject.toml
├── README.md
├── config/
│   ├── pam_model.yaml             # ThermoCas9 PAM families
│   └── cohorts/
│       └── brca_example.yaml      # one cohort spec
├── src/thermocas/
│   ├── models.py                  # Pydantic data model
│   ├── scoring.py                 # candidate score schema
│   ├── pam_model.py               # PAM regex matcher
│   ├── evidence.py                # probe-to-site evidence classifier
│   ├── config.py                  # YAML config loaders
│   └── cli.py                     # `thermocas` entry point
├── tests/                         # pytest suite
└── examples/
    └── synthetic_pipeline.py      # end-to-end demo on synthetic data
```

## Status

| Layer                                  | State                                                       |
|----------------------------------------|-------------------------------------------------------------|
| Data model + scoring schema            | V1 — implemented                                            |
| PAM-model loader + matcher             | V1 — overlap-safe; V3 — exhaustive ACGT⁴…⁸ width validation |
| Probe-to-site evidence classification  | V1 — implemented                                            |
| Cohort + PAM YAML configs              | V1 — implemented                                            |
| Genome-wide PAM catalog builder        | V1 — FASTA streamer, gzip-aware, strand-correct context     |
| Cohort adapter                         | V1 — implemented                                            |
| Pan-cancer aggregator                  | V1 — recurrence / exclusivity / risk                        |
| CLI: `build-catalog` / `score-cohort` / `aggregate` | V1 — implemented                               |
| Probabilistic score (P(targ) × P(prot) × P(trust))  | V2 — implemented                              |
| Subtype-aware backend factory (`split_by_subtype`)  | V2 — implemented                              |
| Real `GDCBackend` with HTTP + on-disk cache         | V2 — implemented (stdlib urllib)              |
| CLI: `--probabilistic`, `--sample-subtypes`, `gdc-fetch` | V2 — implemented                         |
| **Beta-distribution CDF** (regularized incomplete beta via Lentz CF) | **V3 — implemented**          |
| **gRNA spacer scoring** (GC, Tm, runs, hairpin)     | **V3 — implemented**                          |
| **Cross-validation benchmark harness** (P@K, R@K, ROC-AUC)  | **V3 — implemented**                  |
| **CLI: `--spacer`, `benchmark`, `inspect`**         | **V3 — implemented**                          |
| Repo readiness (git, CHANGELOG, CI workflow)        | **V3 — implemented**                          |
| End-to-end synthetic + CLI pipelines + tests | **144 tests** (in-process + CLI + mocked GDC + Beta math) |

## Quickstart

```bash
# Create env (uv recommended)
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"

# Run tests
pytest

# In-process synthetic demo (no CLI, all-Python)
python examples/synthetic_pipeline.py

# Full V1 CLI pipeline: FASTA → catalog → score 2 cohorts → pan-cancer atlas
bash examples/v1_pipeline.sh
```

## CLI

```bash
# 1. Scan a reference for ThermoCas9-compatible PAMs
thermocas build-catalog \
    --reference data/hg38.fa.gz \
    --pam-model config/pam_model.yaml \
    --output catalog.jsonl

# 2. Score the catalog against one cohort's methylation data
thermocas score-cohort \
    --catalog catalog.jsonl \
    --pam-model config/pam_model.yaml \
    --cohort config/cohorts/brca_example.yaml \
    --backend local \
    --probe-annotation tcga_brca_probes.tsv \
    --tumor-beta tcga_brca_tumor.tsv \
    --normal-beta tcga_brca_normal.tsv \
    --output scored.brca.jsonl

# 3. Aggregate scored cohorts into a pan-cancer atlas
thermocas aggregate \
    --scored BRCA=scored.brca.jsonl LUAD=scored.luad.jsonl COAD=scored.coad.jsonl \
    --output panatlas.jsonl
```

Each stage takes files in and writes JSONL out, so any stage is restartable.

## V2 features

```bash
# Add a probabilistic decomposition to each ScoredCandidate:
#   P(targetable_in_tumor) × P(protected_in_normal) × P(observation_trustworthy)
thermocas score-cohort ... --probabilistic --output scored.jsonl

# Fan a single tumor matrix out across PAM50 subtypes:
#   sample_subtypes.tsv has columns sample_id, subtype.
# Outputs scored.LumA.jsonl, scored.LumB.jsonl, ... with cohort_name = 'BRCA::LumA' etc.
thermocas score-cohort ... --sample-subtypes sample_subtypes.tsv --output scored.jsonl

# Pull a TCGA cohort from the live GDC API and export
# LocalArrayBackend-compatible summary TSVs (cached to disk for re-runs):
thermocas gdc-fetch \
    --project TCGA-BRCA \
    --platform HM450 \
    --sample-type both \
    --cache-dir   results/gdc_cache \
    --output-dir  results/gdc_brca
```

The probabilistic decomposition (V3) uses a **method-of-moments Beta(α, β)
fit** from `(mean, IQR/1.349)`, with the regularized incomplete beta
`I_x(α, β)` computed in pure stdlib via Lentz's continued fraction. When the
moments admit no Beta (σ² ≥ μ(1−μ)) or only one quantile is present, the
estimator falls back to the V2 piecewise-linear CDF; with only `mean` it
collapses to a sharp step. Trust scales with `EvidenceClass` and saturates
with sample count.

## V3 features

```bash
# Score the protospacer (gRNA) alongside the methylation-driven ranking
thermocas score-cohort ... --spacer --output scored.jsonl

# Cross-validation benchmark — held-out chromosomes + a positives list
thermocas benchmark \
    --scored scored.jsonl \
    --positives expected_targets.txt \
    --cohort-name BRCA \
    --top-k 20 \
    --score-field final_score \
    --held-out-chromosomes chrX chr22 \
    --output benchmark.jsonl

# Score by V2 probabilistic composite instead:
thermocas benchmark ... --score-field p_therapeutic_selectivity ...
# Or by V3 spacer quality:
thermocas benchmark ... --score-field spacer_final_score ...

# Quick summary of any JSONL artifact
thermocas inspect catalog.jsonl
thermocas inspect scored.brca.jsonl --top 20
thermocas inspect panatlas.jsonl
thermocas inspect benchmark.jsonl
```

## Data model — one paragraph

A `CandidateSite` is a typed record describing a single ThermoCas9-compatible
PAM instance in a reference genome. A `MethylationObservation` summarizes
tumor-vs-normal methylation in a cohort with an `EvidenceClass` indicating how
directly the observation maps onto the candidate's critical PAM cytosine. A
`ScoredCandidate` carries both, plus a `ScoreComponents` breakdown and a
`final_score`. All objects are Pydantic v2 models — JSON / YAML safe and validated.

## Scoring — one equation

```text
final_score
  = sequence_score
    × selectivity_score
    × confidence_score
    − heterogeneity_penalty
    − low_coverage_penalty
```

with

```text
selectivity_score
  = max(0, normal_mean - tumor_mean)
  + max(0, normal_q25  - tumor_q75)
```

quantile separation rewards clean class separation, not just shifts in mean.

## Citation

If you use this framework, please cite:

> Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.
> Molecular basis for methylation-sensitive editing by Cas9.
> *Nature* (2026). DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

## License

MIT.
