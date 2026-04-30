# hg38 prerequisites plan (atlas Roadmap 24 / framework feature branch)

**Status (2026-04-30):** plan-only, no artifacts built. Companion to atlas
`docs/phase4e_hg38_track_design.md` (Roadmap 24). This doc is the framework-
side counterpart, enumerating exactly which artifacts the framework needs to
produce before atlas can run an hg38 build.

**Branch posture:** this doc lives on `hg38-prep`. Stay off `main` until
owner approves merge — hg38 catalog generation is iterative and a feature
branch is cleaner per atlas design-doc decision point #5.

## Snapshot of current framework state

- Zero hg38 footprint: `find . -name '*hg38*'` → empty.
- All catalogs are hg19-aligned: `catalog_hg19_wg.jsonl` (HM450,
  19,787,820 candidates), `catalog_hg19_wg_epic_v2.jsonl` (EPIC-v2,
  35,380,431 candidates).
- `src/thermocas/catalog.py` is genome-build-agnostic at the column
  level: it accepts an arbitrary `--reference` FASTA and `--probe-
  annotation` TSV and emits candidates in the FASTA's coordinate
  system. No hg19 string anywhere in `src/thermocas/`.
- `scripts/build_epic_v2_probes.py` lifts hg38 → hg19 via pyliftover;
  the GPL33022 SOFT source is hg38-native. A hg38-native EPIC-v2 probe
  TSV is the same script with the liftover step skipped (one-flag
  additive change OR a sibling script).

## Atlas-side prereqs (verified 2026-04-30)

All ✅ per atlas `docs/phase4e_hg38_track_design.md`. The atlas can ingest
hg38 artifacts as soon as the framework produces them; no atlas code change
required.

## Framework-side prerequisites (the long pole)

| # | Artifact | Cost | Status |
|---|---|---|---|
| 1 | `data/raw/hg38/hg38_wg.fa` | hours (~3 GB DL) | not started |
| 2 | `data/raw/probes_hg38_hm450.tsv` | days (liftover OR Bioconductor fetch) | not started |
| 3 | `data/raw/probes_hg38_epic_v2_wg.tsv` | hours (script tweak — skip liftover) | not started |
| 4 | `data/raw/ucsc_hg38/{rmsk,refGene,cpgIslandExt,wgEncodeRegDnaseClusteredV3}.txt.gz` | hours | not started |
| 5 | `data/derived/catalog_hg38_wg.jsonl` (HM450) | hours run-time | not started |
| 6 | `data/derived/catalog_hg38_wg_epic_v2.jsonl` (EPIC-v2) | hours run-time | not started |
| 7 | One per-cancer cohort summary on hg38 coords | weeks (re-process) OR fast (liftover) | not started |

Phase 4e-1 (atlas minimum-viable hg38) needs items 1, 4, one of {2, 3},
one of {5, 6}, and #7 for the chosen cancer/platform.

## Concrete fetch commands (DO NOT RUN without owner approval)

The owner explicitly asked: do not download large tracks without confirmation.
These are written down so the next session can execute on green-light.

```bash
# Item 1 — hg38 reference FASTA (~900 MB compressed)
mkdir -p data/raw/hg38
cd data/raw/hg38
for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do
    curl -O "http://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/chr${c}.fa.gz"
done
gunzip *.fa.gz
cat chr{1..22}.fa chrX.fa chrY.fa > hg38_wg.fa

# Item 4 — UCSC hg38 tracks (mirrors hg19 layout under data/raw/ucsc/)
mkdir -p data/raw/ucsc_hg38
cd data/raw/ucsc_hg38
curl -O http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz
curl -O http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/refGene.txt.gz
curl -O http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cpgIslandExt.txt.gz
curl -O http://hgdownload.cse.ucsc.edu/goldenpath/hg38/encDCC/wgEncodeRegDnaseClusteredV3.txt.gz
```

## Probe-annotation paths

### EPIC-v2 (hg38-native; cleanest path)

The SOFT file is hg38-native. The current `build_epic_v2_probes.py`
lifts hg38 → hg19; for hg38 output we skip the lift. Two options:

**Option A — extend the existing script (additive flag, recommended):**

```python
p.add_argument(
    "--target-build", default="hg19", choices=("hg19", "hg38"),
    help="Output coordinate system. hg19 (default) lifts via pyliftover; "
         "hg38 emits MAPINFO unchanged (SOFT is hg38-native)."
)
# In the loop:
if args.target_build == "hg38":
    new_chrom, new_pos1 = chrom, pos_hg38
else:
    hits = lo.convert_coordinate(chrom, pos_hg38 - 1)
    if not hits:
        n_failed += 1
        continue
    new_chrom, new_pos0, *_ = hits[0]
    if new_chrom not in CHROMS:
        n_failed += 1
        continue
    new_pos1 = new_pos0 + 1
```

Run command (when ready):

```bash
uv run scripts/build_epic_v2_probes.py \
    --soft data/raw/epic_v2/GPL33022_family.soft.gz \
    --output data/raw/probes_hg38_epic_v2_wg.tsv \
    --target-build hg38 \
    --chroms all
```

Expected output: ~937,691 probes (up to 8K more than hg19 because there
are no liftover failures; SOFT file is the source of truth at hg38).

**Option B — sibling script.** Cleaner if `build_epic_v2_probes.py` is
hg19-frozen by audit policy. Defer; ask owner.

### HM450 (hg38; two paths)

The current `data/raw/probes_hg19.tsv` is the HM450 manifest's hg19
coordinate set (485,513 probes). For hg38:

**Path A — pyliftover from hg19** (fast, ~99% retention):

```bash
uv run scripts/liftover_probes.py \
    --input data/raw/probes_hg19.tsv \
    --output data/raw/probes_hg38_hm450.tsv \
    --source hg19 --target hg38
```

(`liftover_probes.py` does not yet exist — would be a small new script.)

**Path B — Bioconductor `IlluminaHumanMethylation450k.db` hg38 build**
(authoritative, slower setup):

```r
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Switch to hg38 annotation package and write probe_id, chrom, pos TSV.
```

Recommend Path A for Phase 4e-1 (consistent with cohort-summary
liftover decision), Path B for Phase 4e-2.

## Catalog build commands (when prereqs land)

```bash
# HM450 hg38
thermocas build-catalog \
    --reference data/raw/hg38/hg38_wg.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg38_hm450.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg38_wg.jsonl

# EPIC-v2 hg38
thermocas build-catalog \
    --reference data/raw/hg38/hg38_wg.fa \
    --pam-model config/pam_model.yaml \
    --probe-annotation data/raw/probes_hg38_epic_v2_wg.tsv \
    --probe-window-bp 500 \
    --output data/derived/catalog_hg38_wg_epic_v2.jsonl
```

PAM model is platform-independent (`config/pam_model.yaml`,
sha `ff987b34…`); same model used for hg19 and hg38 catalogs. No
re-derivation needed.

Expected counts:

- HM450 hg38: ~19.8M candidates (probe denominator unchanged from hg19;
  number of PAMs in 500 bp windows shifts very slightly with coordinate
  remapping — should match hg19 to within 0.1-0.5%).
- EPIC-v2 hg38: ~35.4M candidates (probe denominator slightly higher
  than hg19 — 937,691 vs 929,801 because no liftover failures —
  matching small uplift expected).

## Cohort-summary path

Phase 4e-1 owner-recommended path: liftover existing hg19 cohort summary
to hg38 (fast, ~0.3-0.5% probe drop). Phase 4e-2: re-process IDAT files
against hg38-mapped probe annotations (clean, weeks per cancer).

The cohort-summary parquet has a `chrom` + `pos` column per probe; a
liftover utility could be ~50 lines.

## What this branch ships

This branch should land:

- This doc (you're reading it).
- The `--target-build hg38` extension to `scripts/build_epic_v2_probes.py`
  (Option A). Additive only; default behaviour unchanged.
- A `scripts/liftover_probes.py` for HM450 hg19 → hg38 (small new script,
  pyliftover-based).

It should NOT land:

- Downloaded artifacts (FASTA, UCSC tracks).
- Generated catalogs.
- Any change to `main`-branch behaviour.

## Open owner decisions (mirroring atlas design doc)

1. Phase 4e-1 cancer choice (BRCA recommended for cross-build comparison
   with v0.4).
2. Phase 4e-1 platform: HM450 (smaller catalog, broader compat) vs
   EPIC-v2 (cleanest cross-build given GSE322563 native EPIC-v2 was
   already lifted hg38→hg19).
3. Cohort-summary strategy: liftover (4e-1) vs re-process (4e-2).
4. Merge posture: stay on `hg38-prep` until catalog build validates,
   then merge to `main`; or rebase periodically.

These all gate the actual build work; the prep work above is owner-
decision-orthogonal (ships in any of the variants).
