# External Validation Instruction

This note defines the external-validation study needed for a future
stronger ThermoCas9 methods claim. "External validation" here means
**new editing-outcome labels generated after freezing the scorer**. A new
public methylation cohort, or another retrospective benchmark against
known Roth labels, can harden transport and ranking assumptions, but it
does not by itself validate prospective editing selectivity.

The validation claim should be:

> Using a frozen V2.5-sigmoid scorer, we prospectively selected
> ThermoCas9-compatible PAM sites from methylation data only, then tested
> whether the selected sites showed target-cell editing with
> comparator-cell sparing.

## What Counts As External

| Level | What it tests | Counts as external validation? |
|---|---|---|
| New public methylation cohort only | Whether methylation ranks transfer | Weak; not editing validation |
| New public methylation + known Roth labels | Label transport / stress test | Useful, not decisive |
| Frozen V2.5-sigmoid panel + targeted bisulfite only | Whether nominated PAM methylation states transport | Good pre-editing check |
| Frozen panel + genomic-DNA cleavage assay | Methylation-sensitive cleavage outside live-cell context | Strong mechanism check |
| Frozen panel + live-cell editing in both cell types | Actual selective editing | Best external validation |
| Independent lab repeats live-cell panel | Reproducibility and generalizability | Best journal-level validation |

The HEK293T/HCT116 System B public-data attempt belongs in the useful
diagnostic category: it was biologically independent and correctly
pre-registered, but the decisive editable/non-selective VEGFA T9 control
failed public-backend transport because of low coverage. It improves the
audit trail, not the headline selectivity claim.

## Recommended Design

Run a prospective, non-MCF-7 cell-system panel.

1. Pick a new target/protected cell pair.
2. Freeze the scorer, hyperparameters, candidate universe, methylation
   backend, inclusion rules, and analysis code.
3. Score all eligible candidate PAM sites using methylation data only.
4. Select a stratified wet-lab panel before seeing editing outcomes.
5. Confirm PAM methylation by targeted bisulfite amplicon sequencing.
6. Edit both cell types with ThermoCas9 or CE-ThermoCas9.
7. Measure editing outcomes and test whether V2.5-sigmoid enriches for
   selective editing: high editing in target cells, low editing in the
   protected comparator.

## Cell-System Selection

Choose a system satisfying all four requirements:

1. Not MCF-7 / MCF-10A.
2. Methylation data are good enough to score genome-wide or near
   genome-wide.
3. Both cell types are practical for delivery and editing.
4. The target/protected contrast is biologically meaningful.

Candidate system types:

| System type | Pros | Cons |
|---|---|---|
| Cancer cell line vs normal-like line | Closest to therapeutic framing | Matched methylomes may be hard |
| Isogenic methylation-perturbed pair | Clean causal methylation contrast | Less disease-realistic |
| Two unrelated cell lines with strong contrast | Fast and easy | More vulnerable to drift |
| Patient-derived organoid vs matched normal organoid | Strongest biology | Expensive and slow |
| HEK293T/HCT116 with newly generated BSS/editing | Builds on Roth System B | Not disease/normal; T9 must be newly measured |

HEK293T/HCT116 is a reasonable conceptual template because Roth used
ENCODE RRBS and targeted bisulfite sequencing to nominate methylation
states before editing validation. For a real external validation,
however, the missing methylation and editing outcomes must be generated
prospectively under the frozen panel rules.

## Freeze Before Selection

Before choosing any wet-lab targets, freeze:

- code tag;
- V2.5-sigmoid mode: `tumor_plus_gap_sigmoid`;
- `differential_delta = 0.2`;
- `sigma_fixed = sqrt(2) * 0.05 ~= 0.0707`;
- `p_trust` weights and EvidenceClass rules;
- candidate universe and SHA256;
- methylation backend and cohort build script;
- inclusion/exclusion rules;
- top-K and matched-negative reporting rules;
- primary endpoint and QC exclusions;
- analysis script that will consume editing outcomes.

Do not retune after seeing targeted bisulfite or editing results unless
the retuned analysis is explicitly labeled exploratory.

## Stratified Panel

Do not test only the top V2.5-sigmoid sites. A top-only panel can show
whether selected candidates work, but it cannot show why V2.5-sigmoid is
better than `tumor_only` or raw delta-beta.

A journal-strength panel should contain roughly 24-36 assayable sites
after methylation transport filtering. Start with 40-60 nominated sites
because some will fail targeted bisulfite coverage or methylation-state
confirmation.

Suggested 30-site edit panel:

| Candidate class | Count | Purpose |
|---|---:|---|
| High V2.5-sigmoid | 8-10 | Primary predicted selective hits |
| High `tumor_only`, low gap | 5-6 | T9-like editable-but-nonselective decoys |
| High delta-beta-only, lower V2.5 | 4-5 | Tests whether `p_targ` / `p_trust` add value |
| V2.5-diff/sigmoid discordant | 3-4 | Tests gap-factor choice |
| Target-side methylated negatives | 4-5 | Tests `p_targ` specificity |
| Feature-matched random negatives | 4-5 | Denominator controls |

A minimal pilot can use 12-16 sites, but that should be described as
assay feasibility, not definitive external validation.

## PAM Methylation Confirmation

Before editing, perform targeted bisulfite amplicon sequencing at every
selected PAM cytosine in both cell types. This separates methylation
backend transport from editing biology.

Pre-register methylation states:

| State | Suggested rule |
|---|---|
| Target-cell unmethylated | methylation fraction < 0.30 |
| Comparator protected | comparator - target > 0.20, preferably comparator > 0.50 or > 0.70 |
| Nonselective editable decoy | both cell types < 0.30 |
| Protected/off negative | target cell > 0.70 |
| Ambiguous | 0.30-0.70 or low coverage |

If public methylation predicts a selective PAM but targeted bisulfite
does not confirm it, classify the site as a methylation-transport
failure. Do not count it as an editing failure.

## Editing Assay

Use the same editing context the scorer is intended to support. If the
claim is live-cell prioritization, perform live-cell editing.

Recommended assay:

- ThermoCas9 or CE-ThermoCas9 delivered as RNP or mRNA;
- both target and comparator cell types tested for every guide;
- at least 3 biological replicates;
- indels quantified by amplicon NGS;
- non-targeting controls;
- delivery/transfection controls;
- known ThermoCas9-positive control if available;
- optional methylation-insensitive nuclease or nearby SpyCas9 guide as a
  delivery/chromatin-accessibility control.

A genomic-DNA cleavage assay can be used as a lower-cost mechanism check
before live-cell editing. It should not be presented as a substitute for
cellular selectivity.

## Primary Endpoint

Do not make raw editing rate the primary endpoint. The scorer is for
selective editing.

For each candidate:

```text
E_target      = editing percentage in target cell
E_comparator  = editing percentage in protected comparator
selectivity   = E_target - E_comparator
selectivity_ratio = (E_target + epsilon) / (E_comparator + epsilon)
```

Predefine a binary validated selective hit. Example:

```text
PAM methylation transport confirmed
and E_target >= 5 percentage points above no-guide background
and E_comparator <= 1-2%
and E_target / E_comparator >= 3
```

Thresholds may be adjusted to assay sensitivity, but they must be frozen
before editing outcomes are unblinded.

## Expected Success Pattern

| Candidate class | Expected result |
|---|---|
| High V2.5-sigmoid | Enriched for selective editing |
| High `tumor_only`, low gap | Edits target but also edits comparator, or fails selectivity |
| High delta-beta-only, lower V2.5 | Mixed comparison group |
| Target-side methylated negatives | Low target-cell editing |
| Feature-matched random negatives | Low selective-hit rate |

The key comparison is not "does V2.5-sigmoid find editable sites?" The
key comparison is:

> V2.5-sigmoid enriches for target-cell editing with comparator sparing,
> while `tumor_only` enriches for targetability but admits
> nonselective, always-unmethylated loci.

That directly validates the V2 to V2.5 rationale.

## Blinding And Preregistration

For credibility:

1. Freeze scoring output and panel-selection rules.
2. Generate the wet-lab panel from frozen rules.
3. Randomize candidate order.
4. Blind the lab to candidate class and score.
5. Run targeted bisulfite, editing, and amplicon sequencing.
6. Lock QC exclusions before decoding candidate class.
7. Analyze with the pre-written script.

This is the difference between testing interesting targets and
externally validating a ranking method.

## Reporting

For each validation system, report:

- candidate ID, gene, PAM, strand, and EvidenceClass;
- V2.5-sigmoid score/rank;
- V2.5-diff, `tumor_only`, delta-beta-only, and limma-style ranks;
- targeted bisulfite methylation in both cell types;
- editing percentage in both cell types and replicates;
- selectivity difference and selectivity ratio;
- binary selective-hit label;
- QC status and exclusion reason, if any.

Primary summaries:

- top-K selective-hit rate;
- ROC-AUC and PR-AUC for selective-hit labels;
- rank correlation with continuous selectivity;
- matched-negative empirical p-values;
- paired comparisons against delta-beta-only and `tumor_only`;
- candidate-level table with methylation and editing outcomes.

For small panels, the candidate-level table is more convincing than a
single AUC.

## Two-Phase Execution Plan

### Phase 1: Transport And Assay Feasibility

Use 40-60 nominated sites from one new cell pair.

1. Score the frozen candidate universe.
2. Select the stratified nomination panel.
3. Run targeted bisulfite at every PAM in both cell types.
4. Exclude ambiguous, low-coverage, or non-transported sites using
   pre-registered rules.
5. Freeze the final 24-36 site edit panel.

This phase answers whether the methylation backend can nominate testable
differential PAMs.

### Phase 2: Live-Cell Editing

Test the retained panel in both cell types.

Primary analysis:

> Among methylation-confirmed candidates, does V2.5-sigmoid rank
> validated selective hits above decoys and matched negatives?

This is the external validation.

## Deliverables

Before wet-lab work:

- frozen code tag;
- frozen candidate universe SHA256;
- scored candidate file;
- panel-selection TSV;
- preregistration note with endpoint and QC rules;
- analysis script with an empty outcome template.

After targeted bisulfite:

- PAM methylation table;
- transport-status table;
- frozen retained edit panel.

After editing:

- raw amplicon NGS QC outputs;
- per-candidate editing outcome table;
- primary analysis report;
- exploratory analysis report, if any, clearly separated.

## Bottom Line

The right external validation is prospective wet-lab validation of a
frozen, stratified V2.5-sigmoid panel in a non-MCF system. The decisive
endpoint is selective editing, not methylation transfer alone and not
target-cell editability alone. Public methylation-only datasets can
harden the memo, but they do not close the validation gap.
