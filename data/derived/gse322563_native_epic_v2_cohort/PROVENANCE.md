# GSE322563 native EPIC v2 cohort build

- β matrix: data/raw/gse322563/GSE322563_beta_matrix_EPIC_v2.txt.gz
- tumor arm:  MCF71, MCF72  (n=2)
- normal arm: MCF10A1, MCF10A2  (n=2)
- probe universe: data/raw/probes_hg19_epic_v2.tsv (EPIC v2, hg19-lifted)
- rows streamed from β matrix: 937,690
- rows kept (probe in native universe): 147,928
- probe IDs are kept in canonical EPIC v2 form with `_BC##` / `_TC##` /
  `_TO##` / `_BO##` suffixes preserved (no HM450 intersect).
- probes with summaries: tumor=147,928, normal=147,928
