# Roth HEK293T/HCT116 System B Transport Addendum

This addendum is committed after `prereg-transport` and before any System B
scoring. It records the transport outcome and the secondary-backend check for
the VEGFA controls.

## Tag B Transport Outcome

`data/derived/roth_hek_hct_transport.tsv` resolves the ENCODE RRBS transport
gate frozen in `prereg-coords`.

Endpoint 1 status:

| target | cell-line status | consequence |
|---|---|---|
| VEGFA T3 | low coverage in HEK293T and HCT116 | excluded from the selectivity claim |
| VEGFA T9 | low coverage in HEK293T and HCT116 | excluded from the selectivity claim |
| EMX1 T4 | confirmed in HEK293T and HCT116 | retained |
| PRDX4 T5 | confirmed in HEK293T and HCT116 | retained |

Endpoint 2 status:

| target | HEK293T status | consequence |
|---|---|---|
| T4, T10, T11, T12 | confirmed unmethylated/editable | retained |
| T14, T15 | confirmed methylated/not editable | retained |
| T3, T9 | low coverage | excluded |
| T5 | confirmed methylated/not editable in HEK direction | retained |
| T13 | ambiguous | reported as boundary, excluded from binary AUC |

## Secondary-Backend Check

The primary discriminator in the original Endpoint 1 design was the rank
separation between the direction-specific selective positive and VEGFA T9.
Because T9 does not transport under ENCODE RRBS, a secondary backend was checked
before scoring.

Checked sources:

| source | cells | finding | decision |
|---|---|---|---|
| ENCODE Experiment search, WGBS assay | HCT116, HEK293, HEK293T | no WGBS experiments found for these cell-line terms; ENCODE has RRBS/DNAme-array here | not usable as a T9 rescue |
| GEO GSE60106 / GSM1465024 WGBS | HCT116 | VEGFA-control window has sparse one-sample coverage; T9 nearest row has coverage 4 and T3 nearest row has coverage 7 | fails the frozen coverage standard |
| GEO GSE224406 BSseq_Ctrl | HEK293T | VEGFA T9 nearest covered cytosine is 25 bp away with beta 0.818 at coverage 11, conflicting with Roth's unmethylated T9 label; T3 remains below coverage 10 at the nearest rows | not transport-confirming |
| EPIC v2 probe catalog | HCT116 available publicly; no matched HEK293T/HCT116 pair found | nearest EPIC v2 probes to T9/T3 are 231-543 bp away from the PAM cytosines, not the PAM-cytosine transport check frozen for RRBS | not usable as a T9 rescue |

The secondary-backend check does not change the pre-registered ENCODE RRBS
transport thresholds. No new methylation backend is added for VEGFA T3/T9 before
scoring.

## Scoring Consequence

The original Endpoint 1 primary selectivity claim is not evaluable on the public
ENCODE RRBS backend because VEGFA T9 is not transport-confirmed. In particular:

- do not report selective-positive-vs-T9 separation as a System B success claim;
- do not report T9 top-100/top-1000 absence as a V2.5 success claim;
- do not interpret a T9 result from Roth BSS alone as an independent-backend
  validation result.

Any tag-C scoring is therefore a transport-confirmed-subset analysis only:

- Endpoint 1 may report EMX1 T4 and PRDX4 T5 directionality against the
  transport-confirmed opposite-direction target;
- Endpoint 1 must be described as a diagnostic subset, not the pre-registered
  T9 discriminator;
- Endpoint 2 may report binary editability AUC only on transport-confirmed
  HEK293T calls, excluding T3, T9, and ambiguous T13.

This downgrade is frozen before scored JSONLs or benchmark tables are produced.

## Source URLs

- ENCODE WGBS data standard:
  https://www.encodeproject.org/data-standards/wgbs/
- GEO HCT116 WGBS sample:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1465024
- GEO HEK293T WGBS control series:
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE224406
