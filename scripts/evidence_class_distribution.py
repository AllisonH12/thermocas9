"""Evidence-class distributions for full candidate catalogs and
top-K windows per cohort.

Addresses M1 from the third-party review (memo-2026-04-22-as): the
manuscript reports `EvidenceClass` (EXACT / PROXIMAL_CLOSE / PROXIMAL /
REGIONAL / UNOBSERVED) as a load-bearing component of the `p_trust`
factor, but does not report the actual distribution of evidence
classes for either the full chr5/6/10 catalog or the top-K windows
the score axes select. This script fills that gap.

For each of the four evaluated cohort paths (GSE322563 HM450 = primary
endpoint; GSE322563 native EPIC v2 = parallel ingest path; GSE77348 =
δ-development cohort; GSE69914 = high-`n` tissue cohort):

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
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median

# V2.5-sigmoid hyperparameters (defaults from src/thermocas/probabilistic.py).
V25_SIGMOID_DELTA = 0.2
V25_SIGMOID_SIGMA_FIXED = math.sqrt(2) * 0.05  # ≈ 0.0707, matches DEFAULT_GAP_SIGMOID_SIGMA_FIXED


def _v25_sigmoid_score(rec: dict) -> float | None:
    """Recompute V2.5-sigmoid p_therapeutic_selectivity directly from
    the V2.5-diff scored JSONL fields. Lets the EvidenceClass + p_trust
    audit be performed against the V2.5-sigmoid axis without rerunning
    the full scoring pass."""
    prob = rec.get("probabilistic")
    obs = rec.get("observation")
    if not prob or not obs:
        return None
    p_t = prob.get("p_targetable_tumor")
    p_r = prob.get("p_observation_trustworthy")
    if p_t is None or p_r is None:
        return None
    mu_t = obs.get("beta_tumor_mean")
    mu_n = obs.get("beta_normal_mean")
    if mu_t is None or mu_n is None:
        return 0.0  # mirrors p_gap_sigmoid's "no fabrication from one-sided data"
    x = ((mu_n - mu_t) - V25_SIGMOID_DELTA) / V25_SIGMOID_SIGMA_FIXED
    if x >= 0:
        e = math.exp(-x)
        p_gap = 1.0 / (1.0 + e)
    else:
        e = math.exp(x)
        p_gap = e / (1.0 + e)
    return p_t * p_gap * p_r

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


def _top_k_summary(
    jsonl: Path, k: int = TOP_K, axis: str = "v25_diff"
) -> tuple[Counter, dict]:
    """Top-K by `axis` plus per-class (n, median p_targ, median
    gap-factor, median p_trust) inside the top-K window.

    `axis` ∈ {"v25_diff", "v25_sigmoid"}. V2.5-diff uses the stored
    `p_therapeutic_selectivity`; V2.5-sigmoid recomputes the score from
    the same JSONL via `_v25_sigmoid_score`. Both share `p_targ` and
    `p_trust`; the gap factor differs.
    """
    heap: list[tuple[float, int, dict]] = []  # min-heap of (score, idx, rec)
    idx = 0
    for rec in _stream_scored(jsonl):
        prob = rec.get("probabilistic")
        if not prob:
            continue
        if axis == "v25_diff":
            score = prob.get("p_therapeutic_selectivity")
        elif axis == "v25_sigmoid":
            score = _v25_sigmoid_score(rec)
        else:
            raise ValueError(f"unknown axis: {axis!r}")
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
        lambda: {"p_targ": [], "p_gap": [], "p_trust": []}
    )
    for _, _, rec in heap:
        ev = rec["observation"]["evidence_class"]
        prob = rec["probabilistic"]
        counts[ev] += 1
        raw[ev]["p_targ"].append(prob.get("p_targetable_tumor", 0.0))
        if axis == "v25_diff":
            p_gap = prob.get(
                "p_differential_protection", prob.get("p_protected_normal", 0.0)
            )
        else:
            obs = rec["observation"]
            mu_t = obs.get("beta_tumor_mean")
            mu_n = obs.get("beta_normal_mean")
            if mu_t is None or mu_n is None:
                p_gap = 0.0
            else:
                x = (
                    (mu_n - mu_t) - V25_SIGMOID_DELTA
                ) / V25_SIGMOID_SIGMA_FIXED
                p_gap = (
                    1.0 / (1.0 + math.exp(-x))
                    if x >= 0
                    else math.exp(x) / (1.0 + math.exp(x))
                )
        raw[ev]["p_gap"].append(p_gap)
        raw[ev]["p_trust"].append(prob.get("p_observation_trustworthy", 0.0))

    summary: dict[str, dict] = {}
    for ev, d in raw.items():
        summary[ev] = {
            "n": len(d["p_targ"]),
            "p_targ": median(d["p_targ"]) if d["p_targ"] else float("nan"),
            "p_gap": median(d["p_gap"]) if d["p_gap"] else float("nan"),
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


def _format_topk_medians(per_cohort: dict[str, dict], gap_label: str) -> str:
    """Markdown table of per-class (n, median p_targ, median gap-factor,
    median p_trust) inside top-100 per cohort."""
    lines = [
        f"| cohort | EvidenceClass | n in top-100 | median p_targ | median {gap_label} | median p_trust |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for cohort, by_class in per_cohort.items():
        for ev in EVIDENCE_CLASSES:
            s = by_class.get(ev)
            if not s or s["n"] == 0:
                continue
            lines.append(
                f"| {cohort} | `{ev}` | {s['n']} | "
                f"{s['p_targ']:.3f} | {s['p_gap']:.3f} | {s['p_trust']:.3f} |"
            )
    return "\n".join(lines)


def main() -> int:
    full_counts: dict[str, Counter] = {}
    topk_counts_diff: dict[str, Counter] = {}
    topk_medians_diff: dict[str, dict] = {}
    topk_counts_sigmoid: dict[str, Counter] = {}
    topk_medians_sigmoid: dict[str, dict] = {}

    for cohort, jsonl in COHORTS:
        if not jsonl.exists():
            print(f"  SKIP {cohort} — missing {jsonl.relative_to(REPO)}")
            continue
        print(f"  scanning {cohort}: {jsonl.name}")
        full_counts[cohort] = _full_catalog_counts(jsonl)
        topk_counts_diff[cohort], topk_medians_diff[cohort] = _top_k_summary(
            jsonl, axis="v25_diff"
        )
        topk_counts_sigmoid[cohort], topk_medians_sigmoid[cohort] = _top_k_summary(
            jsonl, axis="v25_sigmoid"
        )

    # TSV emission (machine-readable companion to the .md).
    tsv_rows = ["cohort\tscope\taxis\tevidence_class\tcount\tpct"]
    for cohort, counts in full_counts.items():
        total = sum(counts.values())
        for ev in EVIDENCE_CLASSES:
            n = counts.get(ev, 0)
            pct = (100.0 * n / total) if total else 0.0
            tsv_rows.append(
                f"{cohort}\tfull_chr5_6_10\t-\t{ev}\t{n}\t{pct:.4f}"
            )
    for axis_name, topk in [
        ("v25_diff", topk_counts_diff),
        ("v25_sigmoid", topk_counts_sigmoid),
    ]:
        for cohort, counts in topk.items():
            for ev in EVIDENCE_CLASSES:
                n = counts.get(ev, 0)
                pct = (100.0 * n / TOP_K) if TOP_K else 0.0
                tsv_rows.append(
                    f"{cohort}\ttop_100\t{axis_name}\t{ev}\t{n}\t{pct:.4f}"
                )
    (EXAMPLES / "evidence_class_distribution.tsv").write_text(
        "\n".join(tsv_rows) + "\n"
    )

    # Markdown summary.
    md = [
        "# EvidenceClass distribution — full catalogs vs top-100 windows",
        "",
        "Generated by `scripts/evidence_class_distribution.py`. Addresses",
        "M1 from the third-party review (memo-2026-04-22-as) and F-2 from",
        "the follow-up review of -aw: the `EvidenceClass` distance bin is",
        "the `p_trust` factor's discrete saturation lever, but the",
        "manuscript did not previously report what fraction of full-catalog",
        "candidates fall into each bin or what the bin distribution looks",
        "like inside the top-100 window the score axes select. Top-100 is",
        "now reported separately under both **V2.5-diff** (the stored",
        "`p_therapeutic_selectivity` in the scored JSONLs) and",
        "**V2.5-sigmoid** (recomputed in-script from the same JSONL via",
        f"`p_gap_sigmoid` with δ = {V25_SIGMOID_DELTA},",
        f"σ_fixed = √2 × σ_floor ≈ {V25_SIGMOID_SIGMA_FIXED:.4f}) so the",
        "audit matches the recommended axis directly rather than assuming",
        "the two top-100 mixes are identical.",
        "",
        "## Full chr5/6/10 catalog distribution",
        "",
        "Counts (and within-cohort percentages) by `EvidenceClass` across",
        "every candidate in the per-cohort scored JSONL on chr5/6/10",
        "(±500 bp probe scope).",
        "",
        _format_full_table(full_counts),
        "",
        "## Top-100 window distribution — V2.5-diff (`tumor_plus_differential_protection`)",
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
        _format_topk_table(topk_counts_diff),
        "",
        "## Top-100 window distribution — V2.5-sigmoid (`tumor_plus_gap_sigmoid`)",
        "",
        "Top-100 candidates by V2.5-sigmoid `p_targ × p_gap_sigmoid ×",
        f"p_trust` (recomputed in-script with δ = {V25_SIGMOID_DELTA},",
        f"σ_fixed ≈ {V25_SIGMOID_SIGMA_FIXED:.4f} from the same scored",
        "JSONL fields). V2.5-sigmoid is the recommended probabilistic",
        "prioritization axis (PAPER.md §5.2.2 / §6.1); the table below",
        "is the direct audit of its top-100 EvidenceClass mix rather than",
        "an inferred one.",
        "",
        _format_topk_table(topk_counts_sigmoid),
        "",
        "## Per-class score medians inside top-100 windows",
        "",
        "Median `p_targ`, gap-factor, and `p_trust` for the candidates in",
        "each cohort × evidence-class cell of the top-100 window. Lets",
        "the reader see whether the top-100's saturation is driven by a",
        "single class (e.g. EXACT records hitting `p_trust ≈ 0.95`) or",
        "a balanced mix.",
        "",
        "### V2.5-diff top-100",
        "",
        _format_topk_medians(topk_medians_diff, "p_diff"),
        "",
        "### V2.5-sigmoid top-100",
        "",
        _format_topk_medians(topk_medians_sigmoid, "p_gap_sigmoid"),
        "",
        "## Notes",
        "",
        "- Both V2.5 variants share the outer `p_targ × ... × p_trust`",
        "  factors; only the gap-factor slot differs. The V2.5-sigmoid",
        "  top-100 EvidenceClass mix is reported directly above rather",
        "  than assumed to match V2.5-diff's. Compare the two top-100",
        "  tables side-by-side to see how the smooth-sigmoid response",
        "  reshapes the cohort-specific mix vs the discrete-step `p_diff`.",
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
