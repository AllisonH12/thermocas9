# catalog_hg38_wg_epic_v2.jsonl — provenance

Whole-genome probe-window ThermoCas9 target-site catalog using the
**native EPIC v2** probe annotation in **hg38 coordinates** (the platform-
native build — the GPL33022 SOFT source is hg38-native, so this catalog
has no liftover step). Companion to `catalog_hg19_wg_epic_v2.jsonl`
(the hg19-lifted variant) for the atlas Roadmap 24 hg38 release line.

## Build commands

```
# Step 1 — extend the EPIC v2 probe annotation to whole-genome hg38:
uv run scripts/build_epic_v2_probes.py \
    --soft data/raw/epic_v2/GPL33022_family.soft.gz \
    --output data/raw/probes_hg38_epic_v2_wg.tsv \
    --target-build hg38 \
    --chroms all

# Step 2 — build the WG hg38 EPIC v2 catalog:
thermocas build-catalog \
    --reference data/raw/hg38/hg38_wg.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg38_epic_v2_wg.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg38_wg_epic_v2.jsonl
```

## Summary

- **Probe annotation**: 930,140 native EPIC v2 probes (chr1-22 + chrX + chrY);
  out of 937,691 in GPL33022. Drops are the 7,551 probes on alt-loci /
  chrM / contigs outside the standard 24-chromosome set. **No liftover
  failures** — the platform is hg38-native, so MAPINFO is used as-is.
- **Candidate PAM families**: `NNNNCGA`, `NNNNCCA`.
- **Probe-window filter**: ±500 bp.
- **Candidates emitted**: **35,406,213** (vs hg19-lifted EPIC-v2's
  35,380,431; gain of 25,782, +0.07%, reflects avoiding the hg19
  liftover's 339-probe drop).
- **File size**: 11 GB (JSONL).
- **Chromosome coverage**: chr1-22 + chrX + chrY (24 sequences).
- **Build wallclock**: ~13 min (single-threaded, on the same developer
  laptop that built the hg19 catalog in ~12 min).

## Input checksums

- `catalog_hg38_wg_epic_v2.jsonl` SHA256:
  `b57abd429f53f292deae6cbf88f1c566d56c199055597ebc73ff069495f0110a`

The PAM model (`config/pam_model.yaml`, sha `ff987b34…`) is the same
file used for the hg19 catalogs — PAM definitions are genome-build
independent. The hg38 reference FASTA is `data/raw/hg38/hg38_wg.fa`
(concatenation of UCSC per-chromosome FASTAs for chr1-22, X, Y;
md5 `a7f2de9a443ba3d8a55e8dabc4f05299`, ~2.9 GB).

## Downstream

The catalog is gitignored (11 GB) but is fully reproducible from
the build commands above. The canonical citation for this catalog
revision is the SHA256 `b57abd42…` paired with `git rev-parse` of
the enclosing tag.

## Cross-build comparison (hg19 vs hg38 EPIC-v2 catalogs)

| Field | hg19 EPIC-v2 | hg38 EPIC-v2 |
|---|---:|---:|
| probes in scope | 929,801 | 930,140 |
| coordinate origin | lifted hg38 → hg19 | hg38-native (no lift) |
| candidates emitted | 35,380,431 | 35,406,213 |
| catalog SHA256 | `39df8f0f…` | `b57abd42…` |
| genome build field | hg19 | hg38 |
