#!/usr/bin/env bash
# V2 demo — exercises the V2 features added on top of the V1 CLI:
#   * --probabilistic: emit ProbabilisticScore alongside the deterministic score
#   * --sample-subtypes: fan out a single tumor matrix into per-subtype JSONLs
# Builds on the same synthetic FASTA + probe layout as v1_pipeline.sh.
#
# Skipped here (network-dependent, runs against the live GDC API):
#   thermocas gdc-fetch --project TCGA-BRCA \
#       --cache-dir results/gdc_cache \
#       --output-dir results/gdc_brca

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${REPO_ROOT}/results/v2_demo"
mkdir -p "$WORK"
cd "$REPO_ROOT"

if ! command -v thermocas >/dev/null 2>&1; then
    if [[ -x "$REPO_ROOT/.venv/bin/thermocas" ]]; then
        export PATH="$REPO_ROOT/.venv/bin:$PATH"
    else
        echo "thermocas CLI not found — install the package first." >&2
        exit 1
    fi
fi

PAM="$REPO_ROOT/config/pam_model.yaml"

echo "→ writing synthetic inputs into $WORK"

# Reuse a 1-chromosome FASTA with one BRCA-relevant locus.
cat > "$WORK/ref.fa" <<'EOF'
>chr1
AAAAAAAAAAACGTCGAGGGGGGGGGGGGGGGGAAAACCAGGGGAAAACCAGGGGGGGGGGAAAATCGAA
EOF

cat > "$WORK/probes.tsv" <<'EOF'
probe_id	chrom	pos
cg001	chr1	14
cg002	chr1	40
cg003	chr1	52
EOF

# Tumor matrix with 6 samples spanning two subtypes (LumA + LumB).
# cg001 is BRCA-selective in LumA but not LumB.
cat > "$WORK/tumor.tsv" <<'EOF'
probe_id	T1	T2	T3	T4	T5	T6
cg001	0.04	0.06	0.05	0.85	0.82	0.88
cg002	0.50	0.55	0.52	0.48	0.51	0.49
cg003	0.92	0.95	0.91	0.93	0.94	0.96
EOF

cat > "$WORK/normal.tsv" <<'EOF'
probe_id	n1	n2	n3
cg001	0.86	0.90	0.88
cg002	0.50	0.52	0.49
cg003	0.20	0.22	0.18
EOF

# T1, T2, T3 are LumA (BRCA-selective at cg001); T4, T5, T6 are LumB (not selective).
cat > "$WORK/subtypes.tsv" <<'EOF'
sample_id	subtype
T1	LumA
T2	LumA
T3	LumA
T4	LumB
T5	LumB
T6	LumB
EOF

cat > "$WORK/cohort.yaml" <<'EOF'
cohort:
  name: BRCA
  tumor_dataset: TCGA-BRCA
  normal_dataset: TCGA-BRCA-normal
  platform: HM450
  min_samples_tumor: 3
  min_samples_normal: 3
EOF

echo
echo "→ Stage 1 — build-catalog"
thermocas build-catalog \
    --reference "$WORK/ref.fa" \
    --pam-model "$PAM" \
    --output "$WORK/catalog.jsonl"

echo
echo "→ Stage 2a — score-cohort with --probabilistic (single cohort)"
thermocas score-cohort \
    --catalog "$WORK/catalog.jsonl" \
    --pam-model "$PAM" \
    --cohort "$WORK/cohort.yaml" \
    --backend local \
    --probe-annotation "$WORK/probes.tsv" \
    --tumor-beta "$WORK/tumor.tsv" \
    --normal-beta "$WORK/normal.tsv" \
    --probabilistic \
    --output "$WORK/scored.brca.jsonl"

echo
echo "→ Stage 2b — score-cohort --sample-subtypes (fan out into per-subtype JSONLs)"
thermocas score-cohort \
    --catalog "$WORK/catalog.jsonl" \
    --pam-model "$PAM" \
    --cohort "$WORK/cohort.yaml" \
    --backend local \
    --probe-annotation "$WORK/probes.tsv" \
    --tumor-beta "$WORK/tumor.tsv" \
    --normal-beta "$WORK/normal.tsv" \
    --sample-subtypes "$WORK/subtypes.tsv" \
    --probabilistic \
    --output "$WORK/scored.brca.jsonl"

echo
echo "→ Stage 3 — aggregate (LumA vs LumB shows subtype-specific addressability)"
thermocas aggregate \
    --scored "BRCA-LumA=$WORK/scored.brca.LumA.jsonl" "BRCA-LumB=$WORK/scored.brca.LumB.jsonl" \
    --high-score-threshold 0.20 \
    --output "$WORK/panatlas_subtype.jsonl"

echo
echo "→ top 5 by pan_cancer_score (subtype atlas):"
python3 -c "
import json
rows = [json.loads(l) for l in open('$WORK/panatlas_subtype.jsonl')]
rows.sort(key=lambda r: r['pan_cancer_score'], reverse=True)
print(f'{\"candidate_id\":<48}  {\"obs\":>4}  {\"high\":>5}  {\"pan\":>6}  {\"recur\":>6}  {\"excl\":>6}  cohort_scores')
for r in rows[:5]:
    cs = ', '.join(f'{k}={v:.2f}' for k, v in sorted(r['cohort_scores'].items()))
    print(f'{r[\"candidate_id\"]:<48}  {r[\"n_cohorts_observed\"]:>4}  {r[\"n_cohorts_high_score\"]:>5}  '
          f'{r[\"pan_cancer_score\"]:>6.3f}  {r[\"recurrence\"]:>6.2f}  {r[\"exclusivity\"]:>6.3f}  {cs}')
"

echo
echo "→ probabilistic score samples (single-cohort run, top 3):"
python3 -c "
import json
rows = [json.loads(l) for l in open('$WORK/scored.brca.jsonl')]
rows.sort(key=lambda r: (r.get('probabilistic') or {}).get('p_targetable_tumor', 0)
                       * (r.get('probabilistic') or {}).get('p_protected_normal', 0)
                       * (r.get('probabilistic') or {}).get('p_observation_trustworthy', 0),
          reverse=True)
print(f'{\"candidate_id\":<48}  {\"P(targ)\":>8}  {\"P(prot)\":>8}  {\"P(trust)\":>9}  {\"P(sel)\":>7}')
for r in rows[:3]:
    p = r.get('probabilistic')
    if not p:
        continue
    psel = p['p_targetable_tumor'] * p['p_protected_normal'] * p['p_observation_trustworthy']
    print(f'{r[\"candidate\"][\"candidate_id\"]:<48}  '
          f'{p[\"p_targetable_tumor\"]:>8.3f}  {p[\"p_protected_normal\"]:>8.3f}  '
          f'{p[\"p_observation_trustworthy\"]:>9.3f}  {psel:>7.3f}')
"
echo
echo "✓ V2 pipeline complete. Artifacts under $WORK/"
