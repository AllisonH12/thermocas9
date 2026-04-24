"""Evidence-class distributions for full candidate catalogs and
top-K windows per cohort.

Addresses M1 from the third-party review (memo-2026-04-22-as): the
manuscript reports `EvidenceClass` (EXACT / PROXIMAL_CLOSE / PROXIMAL /
REGIONAL / UNOBSERVED) as a load-bearing component of the `p_trust`
factor, but does not report the actual distribution of evidence
classes for either the full chr5/6/10 catalog or the top-K windows
the score axes select. This script fills that gap.

For each of the four primary cohorts:

  * Full chr5/6/10 catalog: count of candidates by EvidenceClass.
  * Top-100 by V2.5-diff `p_therapeutic_selectivity` (from the
    committed scored JSONL): count of candidates by EvidenceClass.
  * Per-class median `p_targ`, median `p_diff`, median `p_trust`
    summary in the top-100 window — exposes whether one class is
    driving the saturation or whether the mix is balanced.

Outputs:

  examples/evidence_class_distribution.tsv
  examples/evidence_class_distribution.md

Run from the repo root:

    uv run python scripts/evidence_class_distribution.py

The scored JSONLs at `data/derived/scored_<cohort>_differential.jsonl`
are gitignored (multi-GB per cohort) but reproducible — see the
PAPER.md reproducibility appendix for the full `uv run` invocations.
"""

from __future__ import annotations

import heapq
import json
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

REPO = Path(__file__).resolve().parent.parent
EXAMPLES = REPO / "examples"
DERIVED = REPO / "data" / "derived"

# (cohort_label, scored_jsonl_path).
COHORTS = [
    ("GSE322563 HM450",      DERIVED / "scored_gse322563_differential.jsonl"),
    ("GSE322563 native v2",  DERIVED / "scored_gse322563_native_differential.jsonl"),
    ("GSE77348",             DERIVED / "scored_surrogate_differential.jsonl"),
    ("GSE69914 (tissue)",    DERIVED / "scored_gse69914_differential.jsonl"),
]
EVIDENCE_CLASSES = ["exact", "proximal_close", "proximal", "regional", "unobserved"]
TOP_K = 100


def _stream_scored(jsonl: Path):
    with jsonl.open() as f:
        for line in f:
            yield json.loads(line)


def _full_catalog_counts(jsonl: Path) -> Counter:
    """Count EvidenceClass across every candidate in the scored JSONL."""
    counts: Counter[str] = Counter()
    for rec in _stream_scored(jsonl):
        ev = rec["observation"]["evidence_class"]
        counts[ev] += 1
    return counts


def _top_k_summary(jsonl: Path, k: int = TOP_K) -> tuple[Counter, dict]:
    """Top-K by V2.5-diff `p_therapeutic_selectivity` plus per-class
    (n, median p_targ, median p_diff, median p_trust) inside the
    top-K window."""
    heap: list[tuple[float, int, dict]] = []  # min-heap of (score, idx, rec)
    idx = 0
    for rec in _stream_scored(jsonl):
        prob = rec.get("probabilistic")
        if not prob:
            continue
        score = prob.get("p_therapeutic_selectivity")
        if score is None:
            continue
        item = (score, idx, rec)
        if len(heap) < k:
            heapq.heappush(heap, item)
        elif score > heap[0][0]:
            heapq.heapreplace(heap, item)
        idx += 1

    counts: Counter[str] = Counter()
    raw: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {"p_targ": [], "p_diff": [], "p_trust": []}
    )
    for _, _, rec in heap:
        ev = rec["observation"]["evidence_class"]
        prob = rec["probabilistic"]
        counts[ev] += 1
        raw[ev]["p_targ"].append(prob.get("p_targetable_tumor", 0.0))
        p_diff = prob.get("p_differential_protection", prob.get("p_protected_normal", 0.0))
        raw[ev]["p_diff"].append(p_diff)
        raw[ev]["p_trust"].append(prob.get("p_observation_trustworthy", 0.0))

    summary: dict[str, dict] = {}
    for ev, d in raw.items():
        summary[ev] = {
            "n": len(d["p_targ"]),
            "p_targ": median(d["p_targ"]) if d["p_targ"] else float("nan"),
            "p_diff": median(d["p_diff"]) if d["p_diff"] else float("nan"),
            "p_trust": median(d["p_trust"]) if d["p_trust"] else float("nan"),
        }
    return counts, summary


def _format_full_table(per_cohort: dict[str, Counter]) -> str:
    """Markdown table of EvidenceClass counts (and within-cohort %)
    across the full chr5/6/10 catalogs of all cohorts."""
    header = "| EvidenceClass | " + " | ".join(per_cohort.keys()) + " |"
    sep = "|---|" + ":---|" * len(per_cohort)
    lines = [header, sep]
    totals = {c: sum(counts.values()) for c, counts in per_cohort.items()}
    for ev in EVIDENCE_CLASSES:
        row = [f"`{ev}`"]
        for cohort in per_cohort:
            n = per_cohort[cohort].get(ev, 0)
            pct = (100.0 * n / totals[cohort]) if totals[cohort] else 0.0
            row.append(f"{n:,} ({pct:.1f}%)")
        lines.append("| " + " | ".join(row) + " |")
    lines.append("| **total** | " + " | ".join(f"{totals[c]:,}" for c in per_cohort) + " |")
    return "\n".join(lines)


def _format_topk_table(per_cohort: dict[str, Counter]) -> str:
    """Markdown table of EvidenceClass counts (and within-cohort %)
    inside the top-100 window of each cohort."""
    header = "| EvidenceClass | " + " | ".join(per_cohort.keys()) + " |"
    sep = "|---|" + ":---|" * len(per_cohort)
    lines = [header, sep]
    for ev in EVIDENCE_CLASSES:
        row = [f"`{ev}`"]
        for cohort in per_cohort:
            n = per_cohort[cohort].get(ev, 0)
            row.append(f"{n}/100" if n else "—")
        lines.append("| " + " | ".join(row) + " |")
    return "\n".join(lines)


def _format_topk_medians(per_cohort: dict[str, dict]) -> str:
    """Markdown table of per-class (n, median p_targ, median p_diff,
    median p_trust) inside top-100 per cohort."""
    lines = [
        "| cohort | EvidenceClass | n in top-100 | median p_targ | median p_diff | median p_trust |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for cohort, by_class in per_cohort.items():
        for ev in EVIDENCE_CLASSES:
            s = by_class.get(ev)
            if not s or s["n"] == 0:
                continue
            lines.append(
                f"| {cohort} | `{ev}` | {s['n']} | "
                f"{s['p_targ']:.3f} | {s['p_diff']:.3f} | {s['p_trust']:.3f} |"
            )
    return "\n".join(lines)


def main() -> int:
    full_counts: dict[str, Counter] = {}
    topk_counts: dict[str, Counter] = {}
    topk_medians: dict[str, dict] = {}

    for cohort, jsonl in COHORTS:
        if not jsonl.exists():
            print(f"  SKIP {cohort} — missing {jsonl.relative_to(REPO)}")
            continue
        print(f"  scanning {cohort}: {jsonl.name}")
        full_counts[cohort] = _full_catalog_counts(jsonl)
        topk_counts[cohort], topk_medians[cohort] = _top_k_summary(jsonl)

    # TSV emission (machine-readable companion to the .md).
    tsv_rows = ["cohort\tscope\tevidence_class\tcount\tpct"]
    for cohort, counts in full_counts.items():
        total = sum(counts.values())
        for ev in EVIDENCE_CLASSES:
            n = counts.get(ev, 0)
            pct = (100.0 * n / total) if total else 0.0
            tsv_rows.append(f"{cohort}\tfull_chr5_6_10\t{ev}\t{n}\t{pct:.4f}")
    for cohort, counts in topk_counts.items():
        for ev in EVIDENCE_CLASSES:
            n = counts.get(ev, 0)
            pct = (100.0 * n / TOP_K) if TOP_K else 0.0
            tsv_rows.append(f"{cohort}\ttop_100\t{ev}\t{n}\t{pct:.4f}")
    (EXAMPLES / "evidence_class_distribution.tsv").write_text("\n".join(tsv_rows) + "\n")

    # Markdown summary.
    md = [
        "# EvidenceClass distribution — full catalogs vs top-100 windows",
        "",
        "Generated by `scripts/evidence_class_distribution.py`. Addresses",
        "M1 from the third-party review (memo-2026-04-22-as): the",
        "`EvidenceClass` distance bin is the `p_trust` factor's discrete",
        "saturation lever, but the manuscript did not previously report",
        "what fraction of full-catalog candidates fall into each bin or",
        "what the bin distribution looks like inside the top-100 window",
        "the score axes select.",
        "",
        "## Full chr5/6/10 catalog distribution",
        "",
        "Counts (and within-cohort percentages) by `EvidenceClass` across",
        "every candidate in the per-cohort scored JSONL on chr5/6/10",
        "(±500 bp probe scope).",
        "",
        _format_full_table(full_counts),
        "",
        "## Top-100 window distribution (V2.5-diff axis)",
        "",
        "Top-100 candidates by `p_therapeutic_selectivity` under V2.5-diff",
        "(the score field stored in the committed scored JSONLs). On",
        "matched cell-line cohorts the top-100 sits inside the saturated",
        "tied band documented in §3.5 + §5.1 — its membership is the",
        "documented tie-break policy's choice, not the scorer's",
        "discrimination — but the EvidenceClass mix is informative",
        "regardless of whether the tied-band records are individually",
        "rank-stable.",
        "",
        _format_topk_table(topk_counts),
        "",
        "## Per-class score medians inside top-100 windows",
        "",
        "Median `p_targ`, `p_diff`, and `p_trust` for the candidates in",
        "each cohort × evidence-class cell of the top-100 window. Lets",
        "the reader see whether the top-100's saturation is driven by a",
        "single class (e.g. EXACT records hitting `p_trust ≈ 0.95`) or",
        "a balanced mix.",
        "",
        _format_topk_medians(topk_medians),
        "",
        "## Notes",
        "",
        "- V2.5-sigmoid top-100 EvidenceClass distribution is",
        "  structurally identical to V2.5-diff's at this filtering",
        "  step: both axes share the same `p_targ × ... × p_trust`",
        "  outer factors and only differ in the gap-factor slot.",
        "  V2.5-sigmoid breaks V2.5-diff's whole-genome tied bands",
        "  (§5.2.2 + Figure 2b) by replacing the discrete-step `p_diff`",
        "  with a smooth sigmoid; the EvidenceClass mix it draws from",
        "  is the same.",
        "- Reproduce by running the manuscript-appendix `uv run` chain to",
        "  rebuild each cohort's scored JSONL, then `uv run python",
        "  scripts/evidence_class_distribution.py`.",
        "",
    ]
    (EXAMPLES / "evidence_class_distribution.md").write_text("\n".join(md))
    print(f"\nwrote {EXAMPLES / 'evidence_class_distribution.tsv'}")
    print(f"wrote {EXAMPLES / 'evidence_class_distribution.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
