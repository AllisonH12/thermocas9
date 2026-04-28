#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Atlas top-100 builder — produces the per-cohort top-100 shortlists for the website.

Wraps `scripts/annotate_top_hits.py` over the four publication cohort
paths and emits:
  - examples/<cohort>_roth_labels/top100_atlas.{tsv,md}  per cohort
  - docs/website/atlas/atlas_top100.json                  combined slim
    JSON for the website (one row per (cohort, rank, candidate)) with
    enough columns for a searchable table without the per-row
    quantitative columns the paper covers in §5.5.

Skips RepeatMasker and DNase clustering on purpose: they cost ~500 MB
of UCSC table loads per invocation and the atlas table doesn't depend
on them. If a v2 atlas needs them, pass --rmsk / --dnase through.

Run:
  uv run python scripts/build_atlas_top100.py
"""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent

# (cohort_label, cohort_short, scored_jsonl, positives_path)
# scored_jsonl is the chr5/6/10 V2.5-diff scoring path; positives_path
# is the corresponding Roth Fig. 5d "wide" positive ID list.
COHORTS = [
    (
        "GSE322563 HM450",
        "gse322563_roth_labels",
        REPO / "data/derived/scored_gse322563_differential.jsonl",
        REPO / "data/derived/positives_roth_wide.txt",
    ),
    (
        "GSE322563 native EPIC v2",
        "gse322563_native_roth_labels",
        REPO / "data/derived/scored_gse322563_native_differential.jsonl",
        REPO / "data/derived/epic_v2_positives/positives_roth_wide.txt",
    ),
    (
        "GSE77348",
        "gse77348_roth_labels",
        REPO / "data/derived/scored_surrogate_differential.jsonl",
        REPO / "data/derived/positives_roth_wide.txt",
    ),
    (
        "GSE69914 tissue",
        "gse69914_roth_labels",
        REPO / "data/derived/scored_gse69914_differential.jsonl",
        REPO / "data/derived/positives_roth_wide.txt",
    ),
]

REFGENE = REPO / "data/raw/ucsc/refGene.txt.gz"
CPG_ISLANDS = REPO / "data/raw/ucsc/cpgIslandExt.txt.gz"

# Slim JSON column subset — keep the table visually scannable on the
# website without dumping the full 30+ paper columns.
JSON_COLUMNS = [
    "rank",
    "candidate_id",
    "chrom",
    "critical_c_pos",
    "strand",
    "pam_family",
    "pam",
    "score",
    "delta_beta",
    "p_trust",
    "nearest_gene",
    "tss_distance_bp",
    "feature_class",
    "cpg_island_context",
    "is_positive",
]


def run_annotate(scored: Path, positives: Path, out_tsv: Path,
                 out_md: Path) -> None:
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "uv", "run", "python", "scripts/annotate_top_hits.py",
        "--scored", str(scored),
        "--score-field", "p_therapeutic_selectivity",
        "--top-k", "100",
        "--refgene", str(REFGENE),
        "--cpg-islands", str(CPG_ISLANDS),
        "--positives", str(positives),
        "--output", str(out_tsv),
        "--markdown", str(out_md),
    ]
    print(f"  $ {' '.join(cmd)}", flush=True)
    subprocess.run(cmd, check=True, cwd=REPO)


def read_top100_tsv(tsv: Path) -> list[dict]:
    """Read the annotated TSV and project onto the slim JSON_COLUMNS schema."""
    out = []
    with tsv.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            slim = {}
            for col in JSON_COLUMNS:
                v = row.get(col, "")
                if col in {"rank", "critical_c_pos", "tss_distance_bp"}:
                    slim[col] = int(v) if v not in ("", "-") else None
                elif col in {"score", "delta_beta", "p_trust"}:
                    try:
                        slim[col] = float(v)
                    except (TypeError, ValueError):
                        slim[col] = None
                elif col == "is_positive":
                    slim[col] = v.strip() not in ("", "0", "FALSE", "False",
                                                  "false", "-")
                else:
                    slim[col] = v
            out.append(slim)
    return out


def main() -> None:
    if not REFGENE.exists() or not CPG_ISLANDS.exists():
        print(
            f"ERROR: missing UCSC inputs at {REFGENE} or {CPG_ISLANDS}.\n"
            "Pull them per the reproducibility tutorial Step 2.",
            file=sys.stderr,
        )
        sys.exit(1)

    combined = {
        "schema_version": 1,
        "description": (
            "Top-100 shortlist per publication cohort, scored under V2.5-diff "
            "(p_therapeutic_selectivity), annotated with nearest gene, "
            "feature class, and CpG-island context. Source pipeline: "
            "scripts/annotate_top_hits.py at top-K=100. RepeatMasker and "
            "DNase clusters intentionally omitted from the slim schema; the "
            "per-cohort TSV/MD under examples/<cohort>_roth_labels/ has the "
            "full 30+ column form."
        ),
        "score_axis": "V2.5-diff (p_therapeutic_selectivity)",
        "cohorts": [],
    }

    for label, short, scored, positives in COHORTS:
        if not scored.exists():
            print(f"  SKIP {label}: {scored} not on disk", file=sys.stderr)
            continue

        out_dir = REPO / "examples" / short
        out_tsv = out_dir / "top100_atlas.tsv"
        out_md = out_dir / "top100_atlas.md"

        print(f"[{label}]", flush=True)
        run_annotate(scored, positives, out_tsv, out_md)
        rows = read_top100_tsv(out_tsv)
        print(f"  -> {len(rows)} rows", flush=True)

        combined["cohorts"].append({
            "cohort": label,
            "short": short,
            "n_rows": len(rows),
            "tsv_relpath": str(out_tsv.relative_to(REPO)),
            "rows": rows,
        })

    out_json = REPO / "docs/website/atlas/atlas_top100.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(combined, indent=2) + "\n")
    print(f"wrote {out_json.relative_to(REPO)}")


if __name__ == "__main__":
    main()
