# Reproducibility tutorial — score and benchmark a methylation cohort on a laptop

*Audience: methods-paper readers who want to verify a row of the AUC table or score a custom cohort. Time: ~15 minutes after the raw inputs are downloaded.*

This walkthrough reproduces the **GSE322563 native EPIC v2** primary-endpoint row of the paper end-to-end, starting from a fresh clone. It is the same pipeline used to produce every `BenchmarkResult` row in the paper.

---

## Prerequisites

- macOS or Linux. Tested on Apple Silicon (Darwin 25.x).
- [`uv`](https://docs.astral.sh/uv/) ≥ 0.4 — `curl -LsSf https://astral.sh/uv/install.sh | sh`. The framework uses `uv` for reproducible Python environments and inline-script dependency management.
- ~5 GB of disk for raw GEO + UCSC inputs.
- A working internet connection for the one-time `curl` block below.

No R, no LaTeX, no Conda required.

## Step 1 — clone the repository at the immutable tag

```bash
git clone https://github.com/AllisonH12/thermocas9.git
cd thermocas9
git checkout memo-2026-04-22-bx        # MANUSCRIPT-shaped tag

# One-time environment setup — installs Python 3.11 + all locked deps
uv sync

# Sanity check — should print "245 passed"
uv run pytest -q
```

If the test count is anything other than 245, stop and open an issue — that means the lock file or the source tree is out of step with the tag.

## Step 2 — pull the raw GEO and UCSC inputs

These files are ~3 GB compressed. They live under `data/raw/` and are gitignored.

```bash
# GEO supplementary — beta matrix + phenotypes
mkdir -p data/raw/gse322563 data/raw/epic_v2 data/raw/ucsc data/raw/hg19

curl -L -o data/raw/gse322563/GSE322563_beta_matrix_EPIC_v2.txt.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE322nnn/GSE322563/suppl/GSE322563_beta_matrix_EPIC_v2.txt.gz"

# EPIC v2 platform annotation (Illumina GPL33022)
curl -L -o data/raw/epic_v2/GPL33022_family.soft.gz \
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL33nnn/GPL33022/soft/GPL33022_family.soft.gz"

# UCSC hg19 reference annotations (gene model + CpG islands + repeats + DNase-HS)
for f in refGene.txt.gz cpgIslandExt.txt.gz rmsk.txt.gz wgEncodeRegDnaseClusteredV3.txt.gz; do
  curl -L -o "data/raw/ucsc/$f" "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/$f"
done

# UCSC hg19 reference FASTA — only chr5/6/10 (the chromosomes hosting all
# three Roth Fig. 5d positives: ESR1 on chr6, EGFLAM on chr5, GATA3 on chr10).
# The catalog-build step (Step 5) needs a single concatenated FASTA.
for chrom in chr5 chr6 chr10; do
  curl -L -o "data/raw/hg19/${chrom}.fa.gz" \
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/${chrom}.fa.gz"
  gunzip -f "data/raw/hg19/${chrom}.fa.gz"
done
cat data/raw/hg19/chr5.fa data/raw/hg19/chr6.fa data/raw/hg19/chr10.fa \
  > data/raw/hg19/hg19_chr5_6_10.fa
```

## Step 3 — build the EPIC v2 probe annotation

This lifts the EPIC v2 SOFT annotation onto hg19 and restricts to chr5/6/10 (the chromosomes hosting all three Roth Fig. 5d positives — *ESR1* on chr6, *EGFLAM* on chr5, *GATA3* on chr10). Takes ~2 minutes; needs `pyliftover` (already locked into the `uv` environment).

```bash
uv run python scripts/build_epic_v2_probes.py \
  --soft data/raw/epic_v2/GPL33022_family.soft.gz \
  --output data/raw/probes_hg19_epic_v2.tsv
```

## Step 4 — build the cohort summary

Collapses the per-sample β matrix into per-probe `(β_mean, β_iqr, n)` summaries for the tumor and normal arms. ~30 seconds.

```bash
uv run python scripts/build_gse322563_native_epic_v2_cohort.py \
  --beta-matrix data/raw/gse322563/GSE322563_beta_matrix_EPIC_v2.txt.gz \
  --epic-v2-probes data/raw/probes_hg19_epic_v2.tsv \
  --output-dir data/derived/gse322563_native_epic_v2_cohort/
```

## Step 5 — build the candidate catalog

The candidate catalog (~3M PAM-bearing loci across chr5/6/10, restricted to within 500 bp of an EPIC v2 probe) is **not** distributed via git — it is gitignored under `data/derived/catalog_*.jsonl` and rebuilt deterministically from the reference FASTA + the shipped PAM model + the EPIC v2 probe annotation. ~3–5 minutes on a modern laptop.

```bash
uv run thermocas build-catalog \
  --reference data/raw/hg19/hg19_chr5_6_10.fa \
  --pam-model config/pam_model.yaml \
  --probe-annotation data/raw/probes_hg19_epic_v2.tsv \
  --output data/derived/catalog_hg19_chr5_6_10_epic_v2.jsonl
```

The same `--pam-model config/pam_model.yaml` value must be passed to `score-cohort` in the next step — the catalog header records the model SHA and the scorer refuses to run against a catalog produced by a different PAM model.

## Step 6 — score the cohort

This walks the catalog, joins each candidate to its nearest probe, and applies V2.5-diff under the cohort YAML. ~90 seconds on a modern laptop.

```bash
uv run thermocas score-cohort \
  --catalog data/derived/catalog_hg19_chr5_6_10_epic_v2.jsonl \
  --cohort data/derived/gse322563_native_differential.yaml \
  --pam-model config/pam_model.yaml \
  --backend summary \
  --probe-annotation data/raw/probes_hg19_epic_v2.tsv \
  --tumor-summary data/derived/gse322563_native_epic_v2_cohort/tumor_summary.tsv \
  --normal-summary data/derived/gse322563_native_epic_v2_cohort/normal_summary.tsv \
  --probabilistic \
  --output data/derived/scored_gse322563_native_differential.jsonl
```

Two flags are easy to miss and both matter:

- `--backend summary` selects the per-probe-summary backend (the format `build_gse322563_native_epic_v2_cohort.py` produced in Step 4). The default `local` backend expects raw per-sample β matrices and will reject the summary TSVs.
- `--probabilistic` is what populates the `p_therapeutic_selectivity` field that Step 7's `--score-field` reads. Without it, the scored JSONL only carries the deterministic V1 `final_score` and the V2.5 benchmark row will not reproduce.

To score under **V2.5-sigmoid** instead of V2.5-diff, edit the cohort YAML in place — there is no separately-shipped sigmoid YAML for this cohort path. Change two fields in `data/derived/gse322563_native_differential.yaml`:

```yaml
probabilistic_mode: tumor_plus_gap_sigmoid   # was: tumor_plus_differential_protection
sigma_fixed: 0.0707                          # add this line; ≈ √2 · σ_floor
```

The cohort loader's iff-semantics validators enforce that `sigma_fixed` is present iff `probabilistic_mode == tumor_plus_gap_sigmoid`. For the multi-bandwidth WG sweep referenced in `MANUSCRIPT.md §5.6` and `PAPER.md §5.2.2`, see `scripts/sigmoid_bandwidth_sweep.py` and `scripts/sigmoid_delta_sigma_wg_sweep.py`, which re-score programmatically rather than via materialized YAMLs.

## Step 7 — benchmark against the Roth validated positives

```bash
uv run thermocas benchmark \
  --scored data/derived/scored_gse322563_native_differential.jsonl \
  --positives data/derived/epic_v2_positives/positives_roth_validated.txt \
  --cohort-name GSE322563-native-validated-V2.5 \
  --score-field p_therapeutic_selectivity --top-k 20 \
  --no-enforce-holdout \
  --output examples/gse322563_native_roth_labels/bench_validated_differential.jsonl
```

Inspect the result with `cat` or `jq`. The `roc_auc`, `precision_at_k`, `precision_at_k_{min, max}`, `recall_at_k`, `recall_at_k_{min, max}`, `tie_band_size_at_k`, and `tie_break_policy` fields are exactly what `MANUSCRIPT.md §5` reports. The `roc_auc` value should match the paper to four decimal places.

## Step 8 — annotate the top 20 with biological context

Joins each top-K hit to nearest gene (refGene), CpG-island context, RepeatMasker overlap, and ENCODE DNase-HS cluster breadth. Emits both a TSV (machine-readable) and a Markdown shortlist (collaborator-readable).

```bash
uv run python scripts/annotate_top_hits.py \
  --scored data/derived/scored_gse322563_native_differential.jsonl \
  --top-k 20 \
  --rmsk data/raw/ucsc/rmsk.txt.gz \
  --dnase data/raw/ucsc/wgEncodeRegDnaseClusteredV3.txt.gz \
  --positives data/derived/epic_v2_positives/positives_roth_wide.txt \
  --output examples/gse322563_native_roth_labels/top20_annotated_v25.tsv \
  --markdown examples/gse322563_native_roth_labels/top20_annotated_v25.md
```

The committed Markdown file under `examples/gse322563_native_roth_labels/top20_annotated_v25.md` should be byte-identical to the one you just produced.

## What you have just reproduced

- One row of the §5.1 / §5.2 primary-endpoint AUC table.
- One row of the §5.5 top-hit annotation table.
- The full set of tie-band-aware benchmark fields the paper relies on.

## Going further

- **Score V2.5-sigmoid instead**: edit `data/derived/gse322563_native_differential.yaml` per the in-place mode swap shown in Step 6 (set `probabilistic_mode: tumor_plus_gap_sigmoid` and add `sigma_fixed: 0.0707`) and re-run from Step 6.
- **Score a different cohort**: build scripts for GSE77348 / GSE69914 / GSE68379 sit beside the GSE322563 build script under `scripts/build_gse*_cohort.py`.
- **Score a custom cohort**: write a cohort YAML against your own per-probe summary TSV. The schema is documented in `src/thermocas/cohort.py` and validated on load.
- **Run the whole-genome panel**: re-run Steps 5–7 with the WG catalog (`--reference` against the full hg19 FASTA — concatenate all autosomes + chrX + chrY rather than only chr5/6/10 — and write to `data/derived/catalog_hg19_wg_epic_v2.jsonl`). The same `data/derived/gse322563_native_differential.yaml` cohort YAML is reused; only the `--catalog` flag changes. Plan for ~10 minutes per (mode × cohort) at the WG scale; the panel SHA256s in `MANUSCRIPT.md §5.6` give you a deterministic sanity check on the candidate universe.

## When something does not match

If a `BenchmarkResult` value diverges from the paper:

1. Confirm `git rev-parse memo-2026-04-22-bx` resolves to the same SHA as the cloned tip.
2. Run `uv run python scripts/verify_manuscript_claims.py`. A green run means the bench artifacts and constants the paper cites are intact.
3. If the verifier is green but your numbers differ, the divergence is almost always either (a) a different `--score-field` flag, (b) a different positives-list tier, or (c) a stale `data/raw/` download. Re-pull from Step 2 and re-run.

---

*Educational research framework. Not peer-reviewed. Not a clinical decision-support system. Cites Roth et al., Nature 2026, [10.1038/s41586-026-10384-z](https://doi.org/10.1038/s41586-026-10384-z).*
