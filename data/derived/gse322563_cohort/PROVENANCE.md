# GSE322563 cohort build — Roth et al. Nature 2026 actual samples

- β matrix: data/raw/gse322563/GSE322563_beta_matrix_EPIC_v2.txt.gz
- tumor arm (MCF-7):   MCF71, MCF72   (n=2)
- normal arm (MCF-10A): MCF10A1, MCF10A2  (n=2)
- platform: Infinium MethylationEPIC v2.0 (GPL33022)
- intersect policy: canonical probe_id match after stripping `_BCxx` suffix
- EPIC-v2 probe rows read:           937,690
- HM450 probe universe:              485,512
- probes retained (EPIC-v2 ∩ HM450): 391,635  (80.7% of HM450)
- Roth-gene probe retention:
    - EGFLAM: 40/44 (91%)
    - EMX1: 27/33 (82%)
    - ESR1: 66/69 (96%)
    - GATA3: 60/72 (83%)
    - PRDX4: 9/10 (90%)
    - VEGFA: 18/20 (90%)

## Caveats

- **n=2 per side** means `p_trust` saturates at 0.0633 for EXACT-evidence
  candidates. Tie-band at the top of any p_trust-multiplied composite is
  expected to be large; treat P@K as secondary to AUC and top-list biology.
- EPIC v2 drops some HM450 probes; retention reported above.
- This is the paper-comparable biological test, NOT a high-n rank-metric
  validation.

## Citation

Roth M.O., Shu Y., Zhao Y., Trasanidou D., Hoffman R.D., et al.
Molecular basis for methylation-sensitive editing by Cas9.
Nature (2026). DOI 10.1038/s41586-026-10384-z.
GEO accession: GSE322563 (per supplementary reporting summary).
Main paper's data-availability statement prints 'GSE32256' (typo).
