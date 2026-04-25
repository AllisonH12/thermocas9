# Roth HEK293T/HCT116 System B Pre-Registration

This note freezes the planned System B validation before ENCODE RRBS transport
resolution, scoring, or manuscript interpretation. The benchmark tests whether
the frozen ThermoCas9 prioritization pipeline respects methylation polarity on
Roth et al.'s independent HEK293T/HCT116 cell-line system.

## Source Targets

Primary labels are in `data/positives/positives_roth_hek_hct_v0.tsv`.

Endpoint 1 is the selectivity direction-flip panel:

| target | role |
|---|---|
| EMX1 T4 | HEK293T-selective positive |
| PRDX4 T5 | HCT116-selective positive |
| VEGFA T9 | editable but non-selective control |
| VEGFA T3 | protected/off control in both cell lines |

Endpoint 2 is the HEK293T editability panel:

| target | role |
|---|---|
| T9, T10, T11, T12, T4 in HEK direction | unmethylated/editable |
| T3, T13, T14, T15, T5 in HEK direction | methylated/not editable |

Endpoint 1 is reported as a pattern-of-ranks table, not as an AUC headline.
Endpoint 2 may be summarized as binary editability AUC because it has more
binary calls.

## Direction Convention

Use target/protected-comparator wording in prose. The scoring code still uses
`tumor` and `normal` field names.

| direction | scoring tumor side | scoring normal side | selective positive |
|---|---|---|---|
| HEK target / HCT protected | HEK293T | HCT116 | EMX1 T4 |
| HCT target / HEK protected | HCT116 | HEK293T | PRDX4 T5 |

Do not pool directions into one Endpoint 1 AUC.

## Sign Pattern

| Axis | Direction-specific selective positive | VEGFA T9 editable/non-selective | VEGFA T3 both-protected | Opposite-direction T4/T5 |
|---|---|---|---|---|
| `tumor_only` | High | High | Low | Low if target-side methylated |
| `p_targ * p_prot` | High | Not high | Low | Low |
| V2.5-diff | High | Not high | Low | Low |
| V2.5-sigmoid | High | Not high | Low | Low |
| delta-beta-only | High | Not high | Not high | Low / direction-flipped |
| limma-style moderated-t | High | Not high | Not high | Low / direction-flipped |

The `p_targ * p_prot` row is an intended-regime positive control for the
deprecated static-comparator-methylation axis. It can succeed here without
weakening the Appendix A deprecation argument, because this curated panel
satisfies the assumption that the protected comparator methylates the
selective target.

## Success Criteria

Primary selectivity claim: under V2.5-diff and V2.5-sigmoid, in each
transport-confirmed direction, the direction-specific selective positive ranks
above VEGFA T9, above the opposite-direction differential target, and above
VEGFA T3.

T9 falsification criterion: VEGFA T9 must not enter the top-100 or top-1000
under V2.5-sigmoid against the tag-A-frozen denominator. T9 ranking high under
`tumor_only` is expected and diagnostic of the V2 failure mode.

Comparator coherence: delta-beta-only and limma-style moderated-t should rank
the selective positive high and T9 not high. T9-high under either comparator is
investigated as a backend or label-transport issue before interpretation.

V2.5 variant tracking: V2.5-diff and V2.5-sigmoid should give similar
directionality calls on this low-replicate cell-line RRBS cohort. Large
divergence triggers backend inspection and is not interpreted as a winner.

## Transport Check

Resolve `transport_status` before scoring, using independent ENCODE RRBS data.
Roth BSS/editing labels remain the validation labels; ENCODE RRBS is used only
to decide whether the label transports to the public methylation backend.

Thresholds:

- methylated: beta > 0.7
- unmethylated: beta < 0.3
- ambiguous: 0.3 <= beta <= 0.7
- coverage: >= 10 reads at the PAM cytosine or nearest assayed CpG within the
  pre-registered distance bins, in >= 2 replicates per cell line

Transport statuses:

- `confirmed`: independent beta agrees with Roth's methylation state
- `ambiguous`: independent beta is in [0.3, 0.7]
- `low_coverage`: coverage rule fails
- `non_transportable`: independent beta disagrees with Roth's state

Only `confirmed` sites enter the Section 5.10 selectivity claim. Other sites
are reported as transport-boundary rows before scoring interpretation.

Coverage is a measurement gate and provenance field. It does not upgrade
genomic `EvidenceClass`, which remains a distance class from assayed CpG to
critical PAM cytosine.

## Frozen Parameters

- differential delta: 0.2
- V2.5-sigmoid sigma_fixed: sqrt(2) * 0.05 ~= 0.0707106781
- V2.5-diff sigma_floor: 0.05
- p_trust ramp_n: 30
- EvidenceClass trust weights: package defaults in `src/thermocas/probabilistic.py`
- unmethylated threshold for `p_targ`: 0.30
- methylated threshold for `p_prot`: 0.50

No delta, sigma, threshold, catalog, or EvidenceClass retuning is permitted
after transport resolution.

## Catalog Rule

Primary denominator is `data/derived/catalog_hg19_wg.jsonl` with SHA256
`d20661c5d5fc0c42491d9a94ef9485f69dbf071ea9a69f43a887d6f82b1357dc`.

Exact resolution requires matching candidate_id, hg19 critical-C coordinate,
strand, 23 bp spacer, and PAM family. T4, T5, and T10-T15 resolve exactly.
VEGFA T3 and T9 are present in hg19 but outside the frozen HM450 probe-window
catalog, so tag A must freeze an augmented denominator:
`catalog_hg19_wg_plus_roth_hek_hct.jsonl`.

The augmented catalog is the frozen WG catalog plus named additions for the
missing Roth candidates. It must receive a new SHA256 at tag A. Do not switch
denominators after scoring.

## Git Tags

| tag suffix | contents |
|---|---|
| `prereg-coords` | labels, cohort YAMLs, thresholds, sign-pattern note, transport stub; no transport table, scores, or interpretation |
| `prereg-transport` | resolved transport table only; no scores |
| `prereg-scored` | scored JSONLs and benchmark TSVs only; no prose interpretation |
| `paper-5-10` | manuscript interpretation referencing the scored artifacts |
