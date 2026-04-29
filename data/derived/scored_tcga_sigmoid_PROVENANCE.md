# scored_tcga_*_wg_sigmoid.jsonl — provenance (V2.5-sigmoid)

V2.5-sigmoid (`tumor_plus_gap_sigmoid`) WG scored JSONLs for the
TCGA-COAD and TCGA-LUAD cohorts. These unblock atlas roadmap item 11
(V2.5-sigmoid WG release) for two of the four cohorts the v0.3 atlas
release line will eventually consume; BRCA and LIHC sigmoid scoring
pending separate runs.

## Build commands

```text
thermocas score-cohort \
    --catalog data/derived/catalog_hg19_wg.jsonl \
    --cohort config/cohorts/tcga_coad_v0_sigmoid.yaml \
    --pam-model config/pam_model.yaml \
    --backend summary \
    --probe-annotation data/derived/tcga_coad_summary/probes.tsv \
    --tumor-summary  data/derived/tcga_coad_summary/tumor_summary.tsv \
    --normal-summary data/derived/tcga_coad_summary/normal_summary.tsv \
    --probabilistic \
    --output data/derived/scored_tcga_coad_wg_sigmoid.jsonl

thermocas score-cohort \
    --catalog data/derived/catalog_hg19_wg.jsonl \
    --cohort config/cohorts/tcga_luad_v0_sigmoid.yaml \
    --pam-model config/pam_model.yaml \
    --backend summary \
    --probe-annotation data/derived/tcga_luad_summary/probes.tsv \
    --tumor-summary  data/derived/tcga_luad_summary/tumor_summary.tsv \
    --normal-summary data/derived/tcga_luad_summary/normal_summary.tsv \
    --probabilistic \
    --output data/derived/scored_tcga_luad_wg_sigmoid.jsonl
```

## Output artifacts

| Artifact | Rows | Size | SHA256 |
|---|---:|---:|---|
| `scored_tcga_coad_wg_sigmoid.jsonl` | 19,787,820 | 25,488,678,160 B (~25.5 GB) | `edeb39138dc031a3e75bef0a712ff7024aa889e2af87ea7975e7733f50e21047` |
| `scored_tcga_luad_wg_sigmoid.jsonl` | 19,787,820 | 25,480,786,506 B (~25.5 GB) | `29887c372bc09f85dd476e0044ac17ed79cab4f05948ce8f5abb5fbdf1b27f04` |

Row count matches `catalog_hg19_wg.jsonl` exactly (19,787,820 rows; sha256 `d20661c5d5fc0c42491d9a94ef9485f69dbf071ea9a69f43a887d6f82b1357dc`). One row per HM450 hg19 WG ThermoCas9 PAM catalog site.

## Scoring policy

- `probabilistic_mode: tumor_plus_gap_sigmoid`
- `differential_delta: 0.2`
- `sigma_fixed: 0.07071067811865475` (= √2 × σ_floor; PAPER.md §5.2.1 bandwidth-robust default)
- `evidence_thresholds`, `penalties`: same as the V2.5-diff cohort configs (zero diff)

The sigmoid mode replaces the V2.5-diff threshold-based `p_differential_protection` with a smooth `p_gap_sigmoid` over the tumor−normal gap distribution. PAPER.md §5.2.2 documents the tie-band reduction at WG K=100 vs V2.5-diff.

## Atlas consumption

These artifacts are referenced by SHA256 from the future atlas v0.3 release configs. Until BRCA + LIHC sigmoid JSONLs exist, the atlas's v0.3 line stays pending — see atlas-side `docs/v25_sigmoid_ingestion_contract.md`.

## Pending sibling artifacts

- `scored_brca_wg_sigmoid.jsonl` — needs to score the three BRCA cohorts (GSE322563 cell-line, GSE77348 cell-line, GSE69914 tissue) under sigmoid mode. Existing differential JSONLs at `scored_gse322563_wg_differential.jsonl` etc. The work is mechanical from the existing per-series cohort YAMLs.
- `scored_tcga_lihc_wg_sigmoid.jsonl` — needs the TCGA-LIHC tumor GDC fetch + cohort summary build first (see `tcga_lihc_v0.yaml`).

## Build wall-time

~11 minutes per cohort on this hardware (24 GB box, Apple Silicon). Both COAD + LUAD ran in parallel. CPU is the bottleneck (probabilistic_score's beta-CDF computation per row × 19.8M rows × probabilistic-component evaluations).
