# scored_*_wg_sigmoid.jsonl — provenance (V2.5-sigmoid)

V2.5-sigmoid (`tumor_plus_gap_sigmoid`) WG scored JSONLs for **all four**
cohorts the atlas v0.3 line will consume: TCGA-COAD, TCGA-LUAD, BRCA
(three sub-cohorts: GSE322563 + GSE69914 + surrogate), and TCGA-LIHC.
Together they unblock atlas roadmap item 11.

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
| `scored_gse322563_wg_sigmoid.jsonl` (BRCA cell-line, n=2/2) | 19,787,820 | 25,218,545,759 B (~25.2 GB) | `e7f229b12713d269b41b2a0998138697e6dabfbcc74530f86159a3c5fb609b89` |
| `scored_gse69914_wg_sigmoid.jsonl` (BRCA tissue, n=305/50) | 19,787,820 | 25,746,487,172 B (~25.7 GB) | `fbcb8e7d40f98a40e3ec6f10b9cfd88aa8bdfddeece5a5c0931dabce171bb9ac` |
| `scored_surrogate_wg_sigmoid.jsonl` (BRCA cell-line surrogate, n=3/3) | 19,787,820 | 25,463,904,450 B (~25.5 GB) | `4c450df82384a1e630304f49301a25ffca79bff04885a15d5410e0306ca5459f` |
| `scored_tcga_lihc_wg_sigmoid.jsonl` (TCGA-LIHC, n=377/50) | TBD | TBD | TBD (scoring in progress) |

Row count matches `catalog_hg19_wg.jsonl` exactly (19,787,820 rows; sha256 `d20661c5d5fc0c42491d9a94ef9485f69dbf071ea9a69f43a887d6f82b1357dc`). One row per HM450 hg19 WG ThermoCas9 PAM catalog site.

## Scoring policy

- `probabilistic_mode: tumor_plus_gap_sigmoid`
- `differential_delta: 0.2`
- `sigma_fixed: 0.07071067811865475` (= √2 × σ_floor; PAPER.md §5.2.1 bandwidth-robust default)
- `evidence_thresholds`, `penalties`: same as the V2.5-diff cohort configs (zero diff)

The sigmoid mode replaces the V2.5-diff threshold-based `p_differential_protection` with a smooth `p_gap_sigmoid` over the tumor−normal gap distribution. PAPER.md §5.2.2 documents the tie-band reduction at WG K=100 vs V2.5-diff.

## Atlas consumption

These artifacts are referenced by SHA256 from the future atlas v0.3 release configs. Until BRCA + LIHC sigmoid JSONLs exist, the atlas's v0.3 line stays pending — see atlas-side `docs/v25_sigmoid_ingestion_contract.md`.

## Pending

- `scored_tcga_lihc_wg_sigmoid.jsonl` — scoring in progress (LIHC tumor GDC fetch + summaries done 2026-04-29; both differential and sigmoid scoring kicked off). SHA256 + size will land here when scoring completes.

## Build wall-time

~11 minutes per cohort on this hardware (24 GB box, Apple Silicon). Both COAD + LUAD ran in parallel. CPU is the bottleneck (probabilistic_score's beta-CDF computation per row × 19.8M rows × probabilistic-component evaluations).
