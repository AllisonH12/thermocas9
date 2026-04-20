#!/usr/bin/env bash
# End-to-end V1 demo: synthetic FASTA → catalog → score 2 cohorts → pan-cancer atlas.
# Runs in ~1 second on a laptop. Verifies every CLI subcommand works on real files.
#
# Usage:  bash examples/v1_pipeline.sh
#         (run from the repo root, with the venv activated:
#            source .venv/bin/activate)

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORK="${REPO_ROOT}/results/v1_demo"
mkdir -p "$WORK"
cd "$REPO_ROOT"

# Use the in-repo venv's `thermocas` if PATH doesn't already have one.
if ! command -v thermocas >/dev/null 2>&1; then
    if [[ -x "$REPO_ROOT/.venv/bin/thermocas" ]]; then
        export PATH="$REPO_ROOT/.venv/bin:$PATH"
    else
        echo "thermocas CLI not found — install the package first:" >&2
        echo "  uv venv && source .venv/bin/activate && uv pip install -e ." >&2
        exit 1
    fi
fi

echo "→ writing synthetic inputs into $WORK"

# 1. Synthetic reference: two chromosomes, each with several embedded
#    NNNNCGA / NNNNCCA windows so the catalog is non-trivial.
cat > "$WORK/ref.fa" <<'EOF'
>chr1
AAAAAAAAAAACGTCGAGGGGGGGGGGGGGGGGAAAACCAGGGGAAAACCAGGGGGGGGGGAAAATCGAA
>chr2
TTTTACGTCGAGGGGGGGGGGGGGAAAACCAGGAAAACCAGGGGAAAATCCAGGGGAAAACCAGGGGAAAA
EOF

# 2. Probe annotation — covers some but not all candidate positions, exercising
#    every EvidenceClass.
cat > "$WORK/probes.tsv" <<'EOF'
probe_id	chrom	pos
cg001	chr1	14
cg002	chr1	40
cg003	chr1	52
cg004	chr2	8
cg005	chr2	36
cg006	chr2	68
EOF

# 3. Tumor + normal beta matrices.
#    cg001: BRCA cancer-low (selective)            cg004: LUAD cancer-low (selective)
#    cg002: BRCA equal (no selectivity)            cg005: LUAD equal
#    cg003: BRCA cancer-high (wrong direction)     cg006: LUAD cancer-low
cat > "$WORK/brca_tumor.tsv" <<'EOF'
probe_id	t1	t2	t3	t4	t5
cg001	0.04	0.06	0.05	0.08	0.07
cg002	0.50	0.55	0.52	0.48	0.51
cg003	0.92	0.95	0.91	0.93	0.94
cg004	0.55	0.50	0.52	0.53	0.51
cg005	0.55	0.50	0.52	0.53	0.51
cg006	0.55	0.50	0.52	0.53	0.51
EOF

cat > "$WORK/brca_normal.tsv" <<'EOF'
probe_id	n1	n2	n3
cg001	0.86	0.90	0.88
cg002	0.50	0.52	0.49
cg003	0.20	0.22	0.18
cg004	0.55	0.50	0.52
cg005	0.55	0.50	0.52
cg006	0.55	0.50	0.52
EOF

cat > "$WORK/luad_tumor.tsv" <<'EOF'
probe_id	t1	t2	t3	t4	t5
cg001	0.55	0.50	0.52	0.53	0.51
cg002	0.55	0.50	0.52	0.53	0.51
cg003	0.55	0.50	0.52	0.53	0.51
cg004	0.05	0.07	0.04	0.06	0.08
cg005	0.50	0.52	0.48	0.51	0.49
cg006	0.06	0.05	0.07	0.04	0.08
EOF

cat > "$WORK/luad_normal.tsv" <<'EOF'
probe_id	n1	n2	n3
cg001	0.55	0.50	0.52
cg002	0.55	0.50	0.52
cg003	0.55	0.50	0.52
cg004	0.84	0.88	0.90
cg005	0.55	0.50	0.52
cg006	0.86	0.92	0.89
EOF

# 4. Cohort YAMLs (smaller min_samples than production so tiny demo data passes).
cat > "$WORK/brca.yaml" <<'EOF'
cohort:
  name: BRCA
  tumor_dataset: TCGA-BRCA
  normal_dataset: TCGA-BRCA-normal
  platform: HM450
  min_samples_tumor: 5
  min_samples_normal: 3
EOF

cat > "$WORK/luad.yaml" <<'EOF'
cohort:
  name: LUAD
  tumor_dataset: TCGA-LUAD
  normal_dataset: TCGA-LUAD-normal
  platform: HM450
  min_samples_tumor: 5
  min_samples_normal: 3
EOF

PAM="$REPO_ROOT/config/pam_model.yaml"

echo
echo "→ Stage 1 — build-catalog"
thermocas build-catalog \
    --reference "$WORK/ref.fa" \
    --pam-model "$PAM" \
    --output "$WORK/catalog.jsonl"

echo
echo "→ Stage 2 — score-cohort BRCA"
thermocas score-cohort \
    --catalog "$WORK/catalog.jsonl" \
    --pam-model "$PAM" \
    --cohort "$WORK/brca.yaml" \
    --backend local \
    --probe-annotation "$WORK/probes.tsv" \
    --tumor-beta "$WORK/brca_tumor.tsv" \
    --normal-beta "$WORK/brca_normal.tsv" \
    --output "$WORK/scored.brca.jsonl"

echo
echo "→ Stage 2 — score-cohort LUAD"
thermocas score-cohort \
    --catalog "$WORK/catalog.jsonl" \
    --pam-model "$PAM" \
    --cohort "$WORK/luad.yaml" \
    --backend local \
    --probe-annotation "$WORK/probes.tsv" \
    --tumor-beta "$WORK/luad_tumor.tsv" \
    --normal-beta "$WORK/luad_normal.tsv" \
    --output "$WORK/scored.luad.jsonl"

echo
echo "→ Stage 3 — aggregate (pan-cancer atlas)"
thermocas aggregate \
    --scored "BRCA=$WORK/scored.brca.jsonl" "LUAD=$WORK/scored.luad.jsonl" \
    --high-score-threshold 0.20 \
    --output "$WORK/panatlas.jsonl"

echo
echo "→ top 5 pan-cancer aggregates by pan_cancer_score:"
python3 -c "
import json
rows = [json.loads(l) for l in open('$WORK/panatlas.jsonl')]
rows.sort(key=lambda r: r['pan_cancer_score'], reverse=True)
print(f'{\"candidate_id\":<48}  {\"obs\":>4}  {\"high\":>5}  {\"pan\":>6}  {\"recur\":>6}  {\"excl\":>6}  cohort_scores')
for r in rows[:5]:
    cs = ', '.join(f'{k}={v:.2f}' for k, v in sorted(r['cohort_scores'].items()))
    print(f'{r[\"candidate_id\"]:<48}  {r[\"n_cohorts_observed\"]:>4}  {r[\"n_cohorts_high_score\"]:>5}  '
          f'{r[\"pan_cancer_score\"]:>6.3f}  {r[\"recurrence\"]:>6.2f}  {r[\"exclusivity\"]:>6.3f}  {cs}')
"
echo
echo "✓ V1 pipeline complete. Artifacts under $WORK/"
