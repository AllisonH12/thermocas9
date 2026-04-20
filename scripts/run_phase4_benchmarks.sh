#!/usr/bin/env bash
# Phase 4 — full-pipeline benchmark on the uncapped TCGA-BRCA cohort.
# Runs after gdc-fetch has finished writing data/raw/gdc_brca/cohort_full/.
#
# Outputs (all under data/derived/):
#   scored_brca_full.jsonl                — 3M scored candidates against full BRCA
#   scored_brca_luma_full.jsonl           — same, against LumA-only sub-cohort
#   bench_full_<axis>.jsonl               — full-cohort benchmarks per score axis
#   bench_lumafull_<axis>.jsonl           — LumA-cohort benchmarks per score axis

set -euo pipefail

REPO=/Volumes/Data/code/thermocas-framework
PATH=/usr/bin:/bin:$PATH

if ! command -v thermocas >/dev/null 2>&1; then
    if [[ -x "$REPO/.venv/bin/thermocas" ]]; then
        export PATH="$REPO/.venv/bin:$PATH"
    fi
fi

COHORT_FULL="$REPO/data/raw/gdc_brca/cohort_full"

echo "=== Stage A: full BRCA cohort ==="
thermocas score-cohort \
    --catalog       "$REPO/data/derived/catalog_hg19_chr5_6_10.jsonl" \
    --pam-model     "$REPO/config/pam_model.yaml" \
    --cohort        "$REPO/data/derived/brca_real.yaml" \
    --backend       summary \
    --probe-annotation "$COHORT_FULL/probes.tsv" \
    --tumor-summary    "$COHORT_FULL/tumor_summary.tsv" \
    --normal-summary   "$COHORT_FULL/normal_summary.tsv" \
    --probabilistic \
    --output        "$REPO/data/derived/scored_brca_full.jsonl"

for axis in final_score p_therapeutic_selectivity naive_selectivity; do
    thermocas benchmark \
        --scored      "$REPO/data/derived/scored_brca_full.jsonl" \
        --positives   "$REPO/data/derived/positives.txt" \
        --cohort-name "BRCA-full / $axis" \
        --top-k 100 \
        --score-field "$axis" \
        --no-enforce-holdout \
        --output      "$REPO/data/derived/bench_full_${axis}.jsonl" 2>&1 | tail -1
done

echo
echo "=== Stage B: LumA-only sub-cohort (full) ==="
mkdir -p "$REPO/data/derived/luma_cohort_full"
.venv/bin/python "$REPO/scripts/build_luma_summary.py" \
    --cache-dir       "$REPO/data/raw/gdc_brca/cache" \
    --file-id-map     "$REPO/data/derived/brca_tumor_file_id_to_submitter.tsv" \
    --luma-submitters "$REPO/data/derived/luma_submitter_ids.txt" \
    --output          "$REPO/data/derived/luma_cohort_full/tumor_summary.tsv"
cp "$COHORT_FULL/normal_summary.tsv" "$REPO/data/derived/luma_cohort_full/normal_summary.tsv"
cp "$COHORT_FULL/probes.tsv" "$REPO/data/derived/luma_cohort_full/probes.tsv"

thermocas score-cohort \
    --catalog       "$REPO/data/derived/catalog_hg19_chr5_6_10.jsonl" \
    --pam-model     "$REPO/config/pam_model.yaml" \
    --cohort        "$REPO/data/derived/brca_luma.yaml" \
    --backend       summary \
    --probe-annotation "$REPO/data/derived/luma_cohort_full/probes.tsv" \
    --tumor-summary    "$REPO/data/derived/luma_cohort_full/tumor_summary.tsv" \
    --normal-summary   "$REPO/data/derived/luma_cohort_full/normal_summary.tsv" \
    --probabilistic \
    --output        "$REPO/data/derived/scored_brca_luma_full.jsonl"

for axis in final_score p_therapeutic_selectivity naive_selectivity; do
    thermocas benchmark \
        --scored      "$REPO/data/derived/scored_brca_luma_full.jsonl" \
        --positives   "$REPO/data/derived/positives.txt" \
        --cohort-name "BRCA-LumA-full / $axis" \
        --top-k 100 \
        --score-field "$axis" \
        --no-enforce-holdout \
        --output      "$REPO/data/derived/bench_lumafull_${axis}.jsonl" 2>&1 | tail -1
done

echo
echo "=== Combined comparison ==="
python3 -c "
import json, glob
print(f'{\"cohort\":<25}  {\"axis\":<28}  {\"n_pos\":>5}  {\"AUC\":>6}  {\"P@100\":>6}')
print('-' * 80)
for fn in sorted(glob.glob('$REPO/data/derived/bench_*.jsonl')):
    with open(fn) as f:
        for ln in f:
            r = json.loads(ln)
            n = fn.split('/')[-1].replace('bench_', '').replace('.jsonl', '')
            cohort, _, axis = r['cohort_name'].partition(' / ')
            auc = f'{r[\"roc_auc\"]:.3f}' if r.get('roc_auc') is not None else 'n/a'
            pk = f'{r[\"precision_at_k\"]:.3f}' if r.get('precision_at_k') is not None else 'n/a'
            print(f'{cohort:<25}  {axis:<28}  {r[\"n_positives\"]:>5}  {auc:>6}  {pk:>6}')
"
