# catalog_hg19_wg.jsonl — provenance

Whole-genome probe-window ThermoCas9 target-site catalog. Frozen
before scoring per the audit discipline: the catalog's SHA256 below
is what every downstream scored JSONL under `scored_*_wg_*.jsonl`
was computed against.

## Build command

```
thermocas build-catalog \
    --reference data/raw/hg19/hg19_wg.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg19.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg19_wg.jsonl
```

Wallclock: ~12 min 36 s on one developer laptop (single-threaded
Python; the `stream_catalog` path is used for the whole-genome
build so peak RSS stays under one FASTA chromosome's worth of
decoded sequence).

## Summary

- **Candidate PAM families**: `NNNNCGA`, `NNNNCCA` (§3.1 / §7 methods).
- **Probe-window filter**: ±500 bp from any HM450 probe.
- **Candidates emitted**: **19,787,820**.
- **File size**: 5.9 GB (JSONL, one record per candidate).
- **Chromosome coverage**: chr1-22 + chrX + chrY (24 sequences).

## Input checksums

See `catalog_hg19_wg.PROVENANCE.sha256` for SHA256 of:
  - the output catalog JSONL
  - the source `hg19_wg.fa` (concatenation of UCSC per-chromosome
    FASTAs for chr1-22, X, Y)
  - `probes_hg19.tsv` (HM450 probe annotation, 485,513 probes)
  - `config/pam_model.yaml` (ThermoCas9 PAM definitions)

Any scored JSONL computed against this catalog should cite the
catalog SHA256 at scoring time, so the denominator (negative
universe) used for that AUC is recoverable from the tag alone.

## Downstream

The catalog is gitignored (5.9 GB exceeds reasonable repo size)
but is fully reproducible from the committed sources above. The
canonical citation for this catalog revision is the catalog SHA256
`d20661c5…` paired with `git rev-parse` of the enclosing tag.
