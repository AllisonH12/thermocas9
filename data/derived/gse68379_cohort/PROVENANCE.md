# GSE68379 breast cell-line tumor cohort

- β matrix: data/raw/gse68379/GSE68379_Matrix.processed.txt.gz
- tumor samples: 52 breast cell lines (primary_site=breast per GSE68379 metadata)
- probes with summaries: 485,512
- platform: HM450 (GPL13534)

## Normal arm note

GSE68379 has no in-study normal. Pair the emitted `tumor_summary.tsv`
with an external healthy-normal `normal_summary.tsv` (typically
`data/derived/gse69914_cohort/normal_summary.tsv`, n=50 healthy
donor breast on the same HM450 platform).

Cross-series composition is a known caveat — Sanger (GSE68379) and
the external normal cohort were processed by different labs with
different scanners. Treat AUC and tie-band behavior as the primary
signal; interpret absolute β values with caution.
