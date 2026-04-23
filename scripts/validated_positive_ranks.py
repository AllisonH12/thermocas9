#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Per-cohort rank table for the n = 3 Roth-validated positives.

Reads the `data/derived/scored_*.jsonl` files in repo, filters to the
three `positives_roth_validated.txt` candidate_ids, and emits the
table of per-positive ranks/percentiles/scores/β values for every
cohort × scoring axis the manuscript reports on.

Output: `examples/validated_positive_ranks.tsv` plus a Markdown
companion `examples/validated_positive_ranks.md` for the §5.1 table.

Design choices:
- We compute rank by re-sorting the JSONL on the requested score axis
  in memory (memory cost ≈ 200–400 MB per file at ~3M rows; acceptable
  on a developer laptop). Streaming a partial sort would be faster
  but obscures the tie-break logic that the benchmark uses.
- Tie-break order is `(-score, candidate_id)` ascending, matching the
  benchmark's `tie_break_policy = "candidate_id_ascending"`. Reported
  rank is 1-based; reported percentile rank is `100 × (1 − rank/N)`
  so higher = better.
- For the V1 final_score axis we use `final_score`. For V2.5 we use
  `probabilistic.p_therapeutic_selectivity`. For tumor_only we re-read
  the tumor_only-mode JSONLs (separate file). For Δβ-only we compute
  `β_normal_mean − β_tumor_mean` from the observation block.
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "validated_positive_ranks.tsv"
OUT_MD = REPO / "examples" / "validated_positive_ranks.md"

# (cohort_label, axis_label, scored_path, score_field)
# score_field is one of:
#   "final_score"    → top-level d["final_score"]
#   "p_therapeutic"  → d["probabilistic"]["p_therapeutic_selectivity"]
#   "delta_beta"     → β_normal_mean − β_tumor_mean (from observation)
RUNS = [
    ("GSE322563 HM450",      "V1 final_score", "scored_gse322563_differential.jsonl",            "final_score"),
    ("GSE322563 HM450",      "V2 tumor_only",  "scored_gse322563_tumor_only.jsonl",              "p_therapeutic"),
    ("GSE322563 HM450",      "V2.5 diff",      "scored_gse322563_differential.jsonl",            "p_therapeutic"),
    ("GSE322563 HM450",      "Δβ-only",        "scored_gse322563_differential.jsonl",            "delta_beta"),
    ("GSE322563 native v2",  "V1 final_score", "scored_gse322563_native_differential.jsonl",     "final_score"),
    ("GSE322563 native v2",  "V2 tumor_only",  "scored_gse322563_native_tumor_only.jsonl",       "p_therapeutic"),
    ("GSE322563 native v2",  "V2.5 diff",      "scored_gse322563_native_differential.jsonl",     "p_therapeutic"),
    ("GSE322563 native v2",  "Δβ-only",        "scored_gse322563_native_differential.jsonl",     "delta_beta"),
    ("GSE77348",             "V1 final_score", "scored_surrogate_differential.jsonl",            "final_score"),
    ("GSE77348",             "V2 tumor_only",  "scored_surrogate_tumor_only.jsonl",              "p_therapeutic"),
    ("GSE77348",             "V2.5 diff",      "scored_surrogate_differential.jsonl",            "p_therapeutic"),
    ("GSE77348",             "Δβ-only",        "scored_surrogate_differential.jsonl",            "delta_beta"),
    ("GSE69914",             "V1 final_score", "scored_gse69914_differential.jsonl",             "final_score"),
    ("GSE69914",             "V2 tumor_only",  "scored_gse69914_tumor_only.jsonl",               "p_therapeutic"),
    ("GSE69914",             "V2.5 diff",      "scored_gse69914_differential.jsonl",             "p_therapeutic"),
    ("GSE69914",             "Δβ-only",        "scored_gse69914_differential.jsonl",             "delta_beta"),
]

POS_GENE = {
    "chr5:38258943+:NNNNCGA":  "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA":  "GATA3",
}


@dataclass
class Row:
    cohort: str
    axis: str
    candidate_id: str
    gene: str
    rank: int
    n_total: int
    percentile: float
    score: float
    beta_tumor: float | None
    beta_normal: float | None
    delta_beta: float | None


def load_positives() -> list[str]:
    return [line.strip() for line in POSITIVES.read_text().splitlines() if line.strip()]


def score_of(record: dict, field: str) -> float | None:
    """Extract the score for one record under the requested axis.

    Returns None if the record cannot be scored on this axis (for
    Δβ-only, that means either β is missing — those records sort to
    the bottom of the ranking deterministically).
    """
    if field == "final_score":
        return record.get("final_score")
    if field == "p_therapeutic":
        return record.get("probabilistic", {}).get("p_therapeutic_selectivity")
    if field == "delta_beta":
        obs = record.get("observation", {})
        bt, bn = obs.get("beta_tumor_mean"), obs.get("beta_normal_mean")
        if bt is None or bn is None:
            return None
        return bn - bt
    raise ValueError(f"unknown score field: {field}")


def rank_one_run(scored_path: Path, score_field: str, positives: list[str]) -> dict[str, tuple[int, int, float, float | None, float | None]]:
    """One pass over scored_path, returning {candidate_id: (rank, n_total, score, β_t, β_n)} for each positive."""
    # Load: list of (score_for_sort, candidate_id, β_t, β_n). For
    # records with score=None on this axis we coerce to -inf so they
    # sort to the bottom under the same tie-break policy.
    rows: list[tuple[float, str, float | None, float | None]] = []
    pos_set = set(positives)
    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            cid = d["candidate"]["candidate_id"]
            s = score_of(d, score_field)
            obs = d.get("observation", {})
            bt, bn = obs.get("beta_tumor_mean"), obs.get("beta_normal_mean")
            # Negate score for ascending sort = descending by score; missing → bottom.
            sort_key = -s if s is not None else float("inf")
            rows.append((sort_key, cid, bt, bn))
            if cid in pos_set:
                pass  # remember? we'll find them after the sort
    rows.sort(key=lambda r: (r[0], r[1]))
    n_total = len(rows)

    found: dict[str, tuple[int, int, float, float | None, float | None]] = {}
    for rank_idx, (sort_key, cid, bt, bn) in enumerate(rows, start=1):
        if cid in pos_set:
            score = -sort_key if sort_key != float("inf") else None
            found[cid] = (rank_idx, n_total, score, bt, bn)
            if len(found) == len(pos_set):
                break
    return found


def main() -> None:
    positives = load_positives()
    print(f"validated positives: {positives}")

    out_rows: list[Row] = []
    for cohort, axis, fname, field in RUNS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP missing: {path}")
            continue
        print(f"scoring {path.name} on {field} ...")
        ranks = rank_one_run(path, field, positives)
        for cid in positives:
            if cid not in ranks:
                continue
            rank, n_total, score, bt, bn = ranks[cid]
            db = (bn - bt) if (bt is not None and bn is not None) else None
            out_rows.append(Row(
                cohort=cohort,
                axis=axis,
                candidate_id=cid,
                gene=POS_GENE[cid],
                rank=rank,
                n_total=n_total,
                percentile=100.0 * (1.0 - rank / n_total),
                score=score if score is not None else float("nan"),
                beta_tumor=bt,
                beta_normal=bn,
                delta_beta=db,
            ))

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\taxis\tgene\tcandidate_id\trank\tn_total\tpercentile\tscore\tbeta_tumor\tbeta_normal\tdelta_beta\n")
        for r in out_rows:
            fh.write(
                f"{r.cohort}\t{r.axis}\t{r.gene}\t{r.candidate_id}\t{r.rank}\t{r.n_total}\t"
                f"{r.percentile:.4f}\t{r.score:.6g}\t"
                f"{'' if r.beta_tumor is None else f'{r.beta_tumor:.4f}'}\t"
                f"{'' if r.beta_normal is None else f'{r.beta_normal:.4f}'}\t"
                f"{'' if r.delta_beta is None else f'{r.delta_beta:.4f}'}\n"
            )
    print(f"wrote {OUT_TSV} ({len(out_rows)} rows)")

    # Markdown companion: one sub-table per cohort, axis × gene grid of ranks.
    by_cohort: dict[str, list[Row]] = defaultdict(list)
    for r in out_rows:
        by_cohort[r.cohort].append(r)

    md = ["# Per-positive ranks of the three Roth Fig. 5d validated targets",
          "",
          "Generated by `scripts/validated_positive_ranks.py` from the per-",
          "cohort scored JSONLs in `data/derived/`. One row per",
          "(cohort × scoring axis × validated positive). Rank is 1-based",
          "under the benchmark's `(-score, candidate_id)` ascending tie-break",
          "(`tie_break_policy = candidate_id_ascending`). Percentile rank is",
          "`100 × (1 − rank/N)` (higher = better).",
          ""]
    for cohort in ["GSE322563 HM450", "GSE322563 native v2", "GSE77348", "GSE69914"]:
        md.append(f"## {cohort}")
        md.append("")
        md.append("| axis | gene | rank | N | percentile | score | β_t | β_n | Δβ |")
        md.append("|---|---|---:|---:|---:|---:|---:|---:|---:|")
        for r in by_cohort.get(cohort, []):
            md.append(
                f"| {r.axis} | *{r.gene}* | {r.rank:,} | {r.n_total:,} | "
                f"{r.percentile:.3f}% | {r.score:.4g} | "
                f"{'' if r.beta_tumor is None else f'{r.beta_tumor:.3f}'} | "
                f"{'' if r.beta_normal is None else f'{r.beta_normal:.3f}'} | "
                f"{'' if r.delta_beta is None else f'{r.delta_beta:.3f}'} |"
            )
        md.append("")
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
