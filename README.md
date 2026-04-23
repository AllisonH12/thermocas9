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
| End-to-end synthetic + CLI pipelines + tests        | **V3 — 144 tests**                            |
| **Real-data validation** (TCGA-BRCA, 3 M candidates, 890 samples) | **V3.1 — implemented; AUC 0.66 on gene-targeted Roth positives** |
| **Input integrity hardening** (16 reviewer-cycle fixes)           | **V3.1 — implemented**                        |
| **Known-limitations documented** (V2 biology assumption, pan-cancer memory) | **V3.1 — implemented**              |
| **V2 formulation review** (additive log-odds / cohort-typed prior) | deferred pending cell-line benchmark       |
| End-to-end tests + real-data pipeline                             | **182 tests** (in-process + CLI + mocked GDC + Beta math + live GDC smoke) |

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

# Pull a TCGA cohort from the live GDC API. gdc-fetch exports per-probe
# summary TSVs (probe_id, n, mean, q25, q75) that LocalSummaryBackend reads
# directly — no raw beta matrices are materialized. Raw files are cached on
# disk for re-runs.
thermocas gdc-fetch \
    --project TCGA-BRCA \
    --platform HM450 \
    --sample-type both \
    --probe-annotation path/to/hm450_probes.tsv \
    --cache-dir   results/gdc_cache \
    --output-dir  results/gdc_brca

# Then score against the exported cohort using --backend summary:
thermocas score-cohort \
    --catalog catalog.jsonl \
    --pam-model config/pam_model.yaml \
    --cohort config/cohorts/brca_example.yaml \
    --backend summary \
    --probe-annotation results/gdc_brca/probes.tsv \
    --tumor-summary    results/gdc_brca/tumor_summary.tsv \
    --normal-summary   results/gdc_brca/normal_summary.tsv \
    --output scored.brca.jsonl
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

## Recommended score axis

Axis choice depends on cohort type. We have benchmarks on three cohorts
under a repaired positives list derived from the three Roth et al.
Fig. 5d validated target coordinates (`EGFLAM T11`, `ESR1 T17`,
`GATA3 T18`), lifted hg38 → hg19 and cross-checked against our per-probe
β values:

| cohort type | example | recommended axis | why |
|---|---|---|---|
| **matched cell-line / paper-comparable** | GSE322563 (Roth actual), GSE77348 (δ-tuning dev cohort) | **V2.5 `tumor_plus_differential_protection`** | Highest-AUC discovery axis at every label granularity (validated / narrow / wide). +0.01 to +0.08 over Δβ-only; +0.01 to +0.17 over V1 on the Roth-validated label set. |
| **primary tumor tissue** | GSE69914 (n=305 / 50) | **V2.5 `tumor_plus_differential_protection`** | Highest-AUC discovery axis on tissue too: +0.113 validated, +0.172 narrow, +0.291 wide over V1. V2's `tumor_only` has the highest raw AUC (0.803–0.874) but `tie_band = 6,540` at K=100 excludes it from discovery use (retained as a diagnostic). V1 collapses on wide (AUC 0.435 — at chance). V2.5's tie-band is 2 at K=100; top-K is effectively deterministic. |
| **any cohort, top-K stability priority / backward compat** | — | V1 `final_score` | V1's continuous-valued deterministic score has `tie_band = 1` at every K regardless of cohort shape. This is the stable-release default (tag `v0.4.0`) for backward compatibility; not the AUC leader anywhere. |

```bash
# cell-line / paper-comparable cohort:
thermocas benchmark ... --score-field p_therapeutic_selectivity ...
#                         (with probabilistic_mode: tumor_plus_differential_protection)

# when you specifically need a clean top-K on any cohort:
thermocas benchmark ... --score-field final_score ...
```

V2.5 is **not** the unconditional stable-release default — V1 `final_score`
is, for backward compatibility and its deterministic `tie_band = 1` top-K
guarantee. V2.5 is the recommended *research* mode across every cohort
shape we tested (matched cell-line *and* tissue), subject to low-`n`
tie-band caveats documented in the next section. `tumor_only` stays
framed as analysis-only. See `V2_5_REVIEW.md` §8 for the full 3-cohort
× 3-label-set × 3-mode matrix.

## Limitations and caveats

### V2 scoring modes (`probabilistic_mode`)

The cohort YAML's `probabilistic_mode` controls which factors enter the
probabilistic composite `p_therapeutic_selectivity`:

| mode | composite | suitable for |
|---|---|---|
| `tumor_only` (default) | `p_targ × p_trust` | **analysis-only** — competitive AUC on tissue, but tie_band at K=100 ranges **5,271–14,914 across the five cohort paths tested** (GSE68379 5,271; GSE69914 6,540; GSE322563 HM450 10,005; GSE77348 11,848; GSE322563 native 14,914). Do NOT use for target-list generation. |
| `tumor_plus_normal_protection` | `p_targ × p_prot × p_trust` | opt-in; known anti-predictive on TCGA-BRCA bulk and inverted on the MCF-7/MCF-10A surrogate. Retained for audit. |
| `tumor_plus_differential_protection` (V2.5) | `p_targ × p_diff × p_trust` where `p_diff = P(β_n − β_t > δ)` | **experimental**; requires `differential_delta` (default 0.2). Highest-AUC discovery axis at every non-boundary cohort × tier combination tested (12/12 rows across **GSE322563 HM450, GSE322563 native EPIC v2, GSE77348, GSE69914**); tie_band scales correctly with `n` (190 at n=2 → 2 at n=305/50). Not the raw-AUC leader on tissue (tumor_only is), but the only probabilistic axis whose top-K is usable there. |

Both modes emit `p_targetable_tumor`, `p_protected_normal`, and
`p_observation_trustworthy` for auditability — the `mode` field on
`ProbabilisticScore` records which policy was in effect. The ablation below
drove the decision.

**Phase 5 ablation result on the MCF-7 vs MCF-10A surrogate (GSE77348,
NOT the Roth samples):**

| score axis | AUC (loose) | AUC (tight) | P@100 (tight) |
|---|---|---|---|
| `v2_tumor_only` (p_targ × p_trust) | **0.733** | **0.770** | **0** |
| `p_targ` alone | 0.683 | 0.717 | — |
| `v1_final` | 0.657 | 0.628 | > 0 (Roth-pattern top-10) |
| `naive_selectivity` (β_normal − β_tumor) | 0.608 | 0.571 | — |
| `p_targ × p_prot` (ablation) | 0.503 | 0.488 | — |
| `p_prot` alone | 0.384 | 0.343 (inverted) | — |

`p_protected_normal` is anti-predictive. Its multiplicative composition
with `p_targ` destroys the signal. `tumor_only` mode avoids this, but
its top-100 collapses onto candidates where both tumor AND normal are
low-methylated (always-unmethylated loci, not cancer-selective) — AUC
and P@100 diverge. The shipped recommendation today (see the mode
table above and PAPER.md §6.1):

- `V1 final_score` is the **stable-release default** — deterministic,
  `tie_band = 1` on every cohort tested, kept for backward
  compatibility and top-K determinism.
- `V2.5 tumor_plus_differential_protection` is the **recommended
  probabilistic ranking axis** on every non-boundary cohort shape
  tested (matched cell-line at n = 2/2 and 3/3 and primary tissue at
  n = 305/50). On n = 2/2 cell-line cohorts the visible top-K
  should be read as a top tied candidate class rather than a ranked
  shortlist (see §6.1).
- `V2 tumor_only` is **analysis-only** — competitive AUC on tissue
  but `tie_band@K=100` is in the thousands; not a discovery axis.

A differential-based `p_prot` (`P(β_normal − β_tumor > δ)` rather than
the threshold-based `P(β_normal > 0.5)`) is the natural next step for
restoring the missing selectivity signal without reintroducing the
static-threshold assumption. **Wired in as V2.5**:
`probabilistic_mode = "tumor_plus_differential_protection"` with
`differential_delta = 0.2` (default).

#### Label-repair note (important)

An earlier `positives_tight.txt` benchmark set that selected Roth-gene
HM450 probes by gene-symbol membership was **noisy supervision** on these
cohorts. Diagnostic on GSE322563 showed only 23% of "Roth-gene"
probes had `β_n − β_t > 0.2` on the actual Roth samples — most were
gene-universe members, not validated target loci. Roth Fig. 5d names
three exact validated target sites in hg38 (`EGFLAM T11`, `ESR1 T17`,
`GATA3 T18`) which we lifted to hg19 and used to rebuild the positives
files:

- `data/derived/positives_roth_validated.txt` — 3 candidate_ids, one per Roth target
- `data/derived/positives_roth_narrow.txt` — 28 candidate_ids, NNNNCGA within ±50 bp
- `data/derived/positives_roth_wide.txt` — 142 candidate_ids, NNNNCGA within ±500 bp

Our per-probe β values at those three sites reproduce Roth's Fig. 5d
β values exactly (EGFLAM 0.01/0.49, ESR1 0.07/0.94, GATA3 0.02/0.31),
confirming the hg38 → hg19 liftover and the EPIC-v2 → HM450 probe
intersect are both correct. Use these files as the primary supervision
target; the older `positives_tight.txt` / `positives.txt` remain in the
repo as auxiliary gene-universe labels but are not the authoritative
benchmark.

#### Cross-cohort matrix under repaired labels

AUC (`tumor_only` → `differential` → `V1`) across all three cohorts:

| cohort | regime | `n` | validated AUC | narrow AUC | wide AUC | V2.5 tie_band |
|---|---|:---:|---|---|---|---:|
| **GSE322563** | Roth cell lines | 2/2 | 0.928 / **0.990** / 0.821 | 0.886 / **0.942** / 0.884 | 0.871 / **0.910** / 0.768 | 190 |
| **GSE77348** | MCF-7/MCF-10A surrogate | 3/3 | 0.912 / **0.982** / 0.968 | 0.911 / **0.983** / 0.969 | 0.887 / **0.949** / 0.931 | 299 |
| **GSE69914** | primary tissue | 305/50 | **0.803** / 0.773 / 0.660 | **0.843** / 0.711 / 0.539 | **0.874** / 0.726 / 0.435 | **2** |

Bold = best AUC in that row. `tumor_only` tie_band is 6,540–11,848 across
all three cohorts — its top-K is never usable. `V1` tie_band is 1 on
every cohort.

Three findings:

1. **V2.5 wins AUC on cell-line cohorts by clear margins.** The label
   repair lifted V2.5 on GSE322563 from 0.694 → 0.990 (validated) and on
   GSE77348 from 0.721 → 0.982. V1 improved similarly (0.541 → 0.821 on
   GSE322563) — so both were being penalized by noise in the old
   supervision target, but V2.5 penalized more.
2. **On primary tissue, `tumor_only` has higher AUC but unusable top-K.**
   V2.5 is the second-best AUC axis there and the only probabilistic
   axis whose top-K is interpretable.
3. **V2.5's tie_band behaves correctly with `n`.** 190 at n=2 → 299 at
   n=3 → 2 at n=305/50. The `n`-saturation prediction from V2_5_REVIEW §3
   holds exactly.

**Bottom line.** V2.5 fixes the biological mis-specification of
`p_protected_normal` (threshold → margin) and is the strongest
probabilistic axis on matched cell-line cohorts, but it is not an
unconditional default. V1 stays available as a continuous-score fallback
with guaranteed top-K stability. `tumor_only` stays analysis-only. V2.5
remains **experimental-on-main**.

A separate GSE68379 run (Sanger GDSC breast panel × external normal)
produced inverted AUC at the Roth-validated positives because Sanger's
MCF-7 is methylated at exactly the sites where Roth's MCF-7 is
unmethylated — a label-portability boundary case, not a V2.5
generalization failure. Documented in `V2_5_REVIEW.md` §8.1; does not
change the recommended-axis table.

See [`V2_5_REVIEW.md`](V2_5_REVIEW.md) §8 for the full cross-cohort
matrix and decision record.

### The V2 probabilistic `P(protected_normal)` factor encodes an assumption

V2 probabilistic scoring decomposes therapeutic selectivity as

    P(therapeutic_selectivity) = P(targetable_tumor)
                               × P(protected_normal)
                               × P(observation_trustworthy)

and operationalizes each factor from per-probe methylation summaries. The
**"protected_normal" factor asks whether the PAM-site cytosine is methylated
in the normal comparator** — specifically `P(β_normal > 0.5)` under a Beta fit
from `(mean, q25, q75)`. That question *assumes the target is methylated in
normal tissue*.

This assumption holds in the Roth et al. (Nature 2026) validation scheme,
where MCF-10A is methylated at `ESR1` / `GATA3` promoters (the genes are
silenced in non-tumorigenic mammary epithelium). **The assumption does NOT
hold for TCGA-BRCA adjacent-normal bulk**, where most gene promoters of
breast-expressed genes are already unmethylated. For that evaluation target,
V2's `P(protected_normal)` systematically zeroes out on biologically-correct
positives and the multiplicative composition underperforms the simpler
`β_normal − β_tumor` baseline — on a Phase 4 TCGA-BRCA run against 1687
gene-targeted positives, V1 `final_score` hit AUC 0.656 while V2
`p_therapeutic_selectivity` hit 0.553 (worse than the naive differential at
0.629).

**Practical guidance:**

| Normal comparator behavior | Recommended score axis |
|---|---|
| Normal tissue methylates the target (cell-line-specific: MCF-10A, HMEC, etc.) | V2 `p_therapeutic_selectivity` |
| Normal tissue does not methylate the target (adjacent-normal bulk, expressed genes) | V1 `final_score` or `naive_selectivity` |
| Unknown or mixed | V1 `final_score` |

The code is mathematically correct; the V2 scoring axis is specialized to
cell-line-style matched comparisons. Reformulating V2 (additive log-odds,
weighted geometric mean, or a cohort-typed `normal_is_methylated` flag) is
deferred until the cell-line-to-cell-line benchmark lands.

### Pan-cancer aggregation: two paths

`thermocas aggregate` groups ScoredCandidates across cohorts and emits in
sorted (chrom, pos, family, candidate_id) order. Two paths share the same
per-candidate emit contract and produce byte-identical output under valid
input:

- **Default (in-memory)**: buffers the full candidate × cohort matrix
  before the first record is yielded — peak cost
  `O(N_unique_candidate_ids × N_cohorts × sizeof(ScoredCandidate))`.
  Fine for small cohorts; for 3M candidates × 10 cohorts this is ~60 GB
  and not tractable.
- **`--streaming`**: k-way merge over pre-sorted cohort JSONLs.
  Candidate-side memory grows in `N_unique_candidate_ids` only (a small
  `seen_cid → metadata` map, ~100 B per entry, used to enforce
  cross-cohort metadata parity) rather than multiplying by
  `N_cohorts × sizeof(ScoredCandidate)`. Requires each input to be
  sorted by `(chrom, critical_c_pos, pam_family, candidate_id)`.

### Positives selection is coordinate-system-sensitive

Phase 3/4 benchmarks in this repo used the HM450 manifest (hg19) throughout,
matched methylation by `probe_id` (build-agnostic), and selected positives by
gene-symbol intersection (build-robust). The Roth paper reports target
coordinates in hg38, so any benchmark that uses raw Roth `(chrom, pos)`
tuples *directly* needs liftOver to hg19 first — or a full catalog rebuild
against hg38.

## Reproducibility

Manuscript / memo PDFs (`MANUSCRIPT.pdf`, `PAPER.pdf`) at
**`memo-2026-04-22-k` and later** tags are byte-identical from a
fresh clone given the same toolchain version pair (verified on
`pandoc 3.9` + `typst 0.14.2`). The render helper
(`scripts/render_paper_pdf.sh`) sources the date from the source
MD's `**Date.**` line and exports `SOURCE_DATE_EPOCH` so Typst
embeds a deterministic `CreationDate` / `ModDate`. A different
`typst` minor version will produce a different PDF from the same
source — that is renderer drift across versions, not pipeline
non-determinism. Earlier tags in the `memo-2026-04-22-*` series
do not have this guarantee: `-i` and earlier had a render-script
bug that leaked title-block metadata into the rendered PDF body,
and `-j` still injected wall-clock `CreationDate` / `ModDate` so
re-renders diff at the metadata level even when content is
identical. See the PAPER.md tag ledger for the per-tag history.

Run `scripts/verify_manuscript_claims.py` before cutting any new
`memo-*` tag; it cross-checks numeric / universal claims in
`MANUSCRIPT.md`, `PAPER.md`, and `README.md` against the committed
bench JSONLs and source-code constants. The script's docstring
enumerates exactly what it checks and three known coverage gaps.

## Citation

If you use this framework, please cite:

> Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.
> Molecular basis for methylation-sensitive editing by Cas9.
> *Nature* (2026). DOI [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).

## License

MIT.
