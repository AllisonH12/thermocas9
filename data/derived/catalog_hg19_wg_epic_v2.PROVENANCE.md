# catalog_hg19_wg_epic_v2.jsonl — provenance

Whole-genome probe-window ThermoCas9 target-site catalog using the
**native EPIC v2** probe annotation (GPL33022 lifted hg38 → hg19,
all chromosomes). Companion to `catalog_hg19_wg.jsonl` (HM450 path)
for the GSE322563 native ingest in PAPER.md §5.2.2 / §5.2.3.

## Build commands

```
# Step 1 — extend the EPIC v2 probe annotation to whole-genome:
uv run scripts/build_epic_v2_probes.py \
    --soft data/raw/epic_v2/GPL33022_family.soft.gz \
    --output data/raw/probes_hg19_epic_v2_wg.tsv \
    --chroms all

# Step 2 — build the WG EPIC v2 catalog:
thermocas build-catalog \
    --reference data/raw/hg19/hg19_wg.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg19_epic_v2_wg.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg19_wg_epic_v2.jsonl
```

## Summary

- **Probe annotation**: 929,801 EPIC v2 probes lifted hg38 → hg19
  successfully (out of 937,691 in GPL33022; 339 liftover failures /
  cross-chromosome lifts skipped).
- **Candidate PAM families**: `NNNNCGA`, `NNNNCCA`.
- **Probe-window filter**: ±500 bp.
- **Candidates emitted**: **35,380,431** (1.79× the HM450 catalog
  19.79M, reflecting EPIC v2's denser probe coverage).
- **File size**: 11 GB (JSONL).
- **Chromosome coverage**: chr1-22 + chrX + chrY.

## Input checksums

See `catalog_hg19_wg_epic_v2.PROVENANCE.sha256` for SHA256 of:
  - `catalog_hg19_wg_epic_v2.jsonl` (the catalog)
  - `probes_hg19_epic_v2_wg.tsv` (the WG probe annotation)

The source `hg19_wg.fa` and `config/pam_model.yaml` checksums are
in `catalog_hg19_wg.PROVENANCE.sha256` (same source files used for
both whole-genome catalogs).

## Downstream

The catalog is gitignored (11 GB) but is fully reproducible from
the build commands above. The canonical citation for this catalog
is the SHA256 `39df8f0f…` paired with `git rev-parse` of the
enclosing tag.
