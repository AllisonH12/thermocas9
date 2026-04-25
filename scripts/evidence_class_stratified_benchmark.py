#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""EvidenceClass-controlled benchmark for validated Roth positives.

Reviewer question: does V2.5-sigmoid rank validated candidates because
of methylation logic, or because `p_trust` strongly rewards candidates
near assayed probes?

This script benchmarks four axes within EvidenceClass-controlled
universes:

  * exact_only
  * exact_plus_proximal_close
  * each EvidenceClass separately

For each cohort/path, each catalog scope (WG and chr5/6/10), each
universe, and each axis, it reports:

  * n_total, n_pos, n_neg
  * validated-label AUC (mid-rank/tie-aware)
  * per-positive ranks and percentiles
  * tie_band@100 and P@100 [min, max]

The implementation streams the scored JSONLs twice and does not retain
the full ranking in memory. V2.5-sigmoid and Delta-beta-only are
recomputed from the committed V2.5-diff JSONL fields; limma-style scores
are joined from the per-probe moderated-t TSVs.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import json
import math
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"
DERIVED = REPO / "data" / "derived"

CHR_SUBSET = {"chr5", "chr6", "chr10"}
EVIDENCE_CLASSES = ("exact", "proximal_close", "proximal", "regional", "unobserved")
K = 100

DEFAULT_DELTA = 0.2
DEFAULT_SIGMA_FIXED = math.sqrt(2) * 0.05

POS_GENE = {
    "chr5:38258943+:NNNNCGA": "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA": "GATA3",
}
GENE_ORDER = ("ESR1", "EGFLAM", "GATA3")

AXES = ("delta_beta_only", "v25_diff", "v25_sigmoid", "limma_style")
SIGMOID_ONLY = ("v25_sigmoid",)
SCOPES = ("wg", "chr5_6_10")


def strip_probe_suffix(probe_id: str) -> str:
    if "_" in probe_id:
        return probe_id.split("_", 1)[0]
    return probe_id


def load_limma(path: Path) -> dict[str, float]:
    """Return {canonical_probe_id: t_mod}.

    EPIC v2 probe suffixes collapse to their canonical probe id. When
    duplicates collapse, keep the most extreme absolute t-statistic,
    matching the existing limma candidate-mapping script.
    """
    out: dict[str, float] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            probe_id = row["probe_id"]
            try:
                t_mod = float(row["t_mod"])
            except ValueError:
                t_mod = 0.0
            if not math.isfinite(t_mod):
                t_mod = 0.0
            key = strip_probe_suffix(probe_id)
            old = out.get(key)
            if old is None or abs(t_mod) > abs(old):
                out[key] = t_mod
    return out


def sigmoid(x: float) -> float:
    if x >= 0:
        e = math.exp(-x)
        return 1.0 / (1.0 + e)
    e = math.exp(x)
    return e / (1.0 + e)


def reverse_lex_key(s: str) -> str:
    """Invert ASCII lexicographic order for heap comparisons.

    Top-K ordering is score desc, candidate_id asc. A min-heap keeps the
    worst retained item at the root. For equal scores, the worse item is
    the lexicographically larger candidate_id; inverting characters lets
    tuple ordering model that without a custom comparator.
    """
    return "".join(chr(255 - ord(ch)) for ch in s)


def rank_key(score: float, candidate_id: str) -> tuple[float, str]:
    return (score, reverse_lex_key(candidate_id))


def universe_axes_for(evidence_class: str) -> tuple[tuple[str, tuple[str, ...]], ...]:
    """Return benchmark universes and axes relevant for one record.

    Full four-axis comparisons are only material where the three validated
    positives can appear: exact-only and exact+proximal-close. Per-class
    non-EXACT rows are retained for auditability, but they have n_pos = 0
    under the exact Roth validated labels, so V2.5-sigmoid alone is enough
    to show that dependence without spending full-axis work on every record.
    """
    out: list[tuple[str, tuple[str, ...]]] = [(f"class_{evidence_class}", SIGMOID_ONLY)]
    if evidence_class == "exact":
        out[0] = ("class_exact", AXES)
        out.append(("exact_only", AXES))
        out.append(("exact_plus_proximal_close", AXES))
    elif evidence_class == "proximal_close":
        out.append(("exact_plus_proximal_close", AXES))
    return tuple(out)


def v25_sigmoid_score(rec: dict) -> float:
    obs = rec["observation"]
    prob = rec["probabilistic"]
    mu_t = obs.get("beta_tumor_mean")
    mu_n = obs.get("beta_normal_mean")

    if mu_t is None or mu_n is None:
        p_gap_sigmoid = 0.0
    else:
        delta_beta = mu_n - mu_t
        p_gap_sigmoid = sigmoid((delta_beta - DEFAULT_DELTA) / DEFAULT_SIGMA_FIXED)

    p_targ = prob.get("p_targetable_tumor") or 0.0
    p_trust = prob.get("p_observation_trustworthy") or 0.0
    return p_targ * p_gap_sigmoid * p_trust


def score_record(rec: dict, limma_t: dict[str, float]) -> dict[str, float]:
    obs = rec["observation"]
    prob = rec["probabilistic"]
    mu_t = obs.get("beta_tumor_mean")
    mu_n = obs.get("beta_normal_mean")

    if mu_t is None or mu_n is None:
        delta_beta = float("-inf")
    else:
        delta_beta = mu_n - mu_t

    probe_id = obs.get("probe_id")
    if probe_id is None:
        limma_score = 0.0
    else:
        t_mod = limma_t.get(probe_id)
        if t_mod is None:
            t_mod = limma_t.get(strip_probe_suffix(probe_id))
        limma_score = -t_mod if t_mod is not None else 0.0

    return {
        "delta_beta_only": delta_beta,
        "v25_diff": prob.get("p_therapeutic_selectivity") or 0.0,
        "v25_sigmoid": v25_sigmoid_score(rec),
        "limma_style": limma_score,
    }


@dataclass
class RunningMetric:
    n_total: int = 0
    positive_scores: dict[str, float] = field(default_factory=dict)
    top_heap: list[tuple[tuple[float, str], str, float]] = field(default_factory=list)

    # Filled after pass 1.
    cutoff_score: float | None = None
    top_ids: set[str] = field(default_factory=set)

    # Filled in pass 2.
    n_neg: int = 0
    auc_wins: float = 0.0
    rank_counts: dict[str, int] = field(default_factory=lambda: defaultdict(int))
    tie_band: int = 0
    band_pos: int = 0
    above_cutoff: int = 0
    above_cutoff_pos: int = 0

    def observe_first_pass(self, candidate_id: str, score: float, is_positive: bool) -> None:
        self.n_total += 1
        if is_positive:
            self.positive_scores[candidate_id] = score
        key = rank_key(score, candidate_id)
        item = (key, candidate_id, score)
        if len(self.top_heap) < K:
            heapq.heappush(self.top_heap, item)
        elif key > self.top_heap[0][0]:
            heapq.heapreplace(self.top_heap, item)

    def finalize_first_pass(self) -> None:
        if not self.top_heap:
            return
        ordered = sorted(self.top_heap, key=lambda item: item[0], reverse=True)
        self.cutoff_score = ordered[min(K, len(ordered)) - 1][2]
        self.top_ids = {candidate_id for _, candidate_id, _ in ordered[:K]}

    def observe_second_pass(self, candidate_id: str, score: float, is_positive: bool) -> None:
        if not is_positive:
            self.n_neg += 1
            for pos_score in self.positive_scores.values():
                if pos_score > score:
                    self.auc_wins += 1.0
                elif pos_score == score:
                    self.auc_wins += 0.5

        for pos_id, pos_score in self.positive_scores.items():
            if score > pos_score or (score == pos_score and candidate_id < pos_id):
                self.rank_counts[pos_id] += 1

        if self.cutoff_score is None:
            return
        if score > self.cutoff_score:
            self.above_cutoff += 1
            if is_positive:
                self.above_cutoff_pos += 1
        elif score == self.cutoff_score:
            self.tie_band += 1
            if is_positive:
                self.band_pos += 1

    def row(self, cohort: str, scope: str, universe: str, axis: str) -> dict[str, str]:
        n_pos = len(self.positive_scores)
        if n_pos == 0 or self.n_neg == 0:
            auc = ""
        else:
            auc = f"{self.auc_wins / (n_pos * self.n_neg):.6f}"

        k_eff = min(K, self.n_total)
        if k_eff == 0 or self.cutoff_score is None:
            p_obs = p_min = p_max = ""
        else:
            observed_pos = sum(1 for cid in self.top_ids if cid in self.positive_scores)
            k_band = max(0, k_eff - self.above_cutoff)
            drop_slots = max(0, self.tie_band - k_band)
            min_pos = self.above_cutoff_pos + max(0, self.band_pos - drop_slots)
            max_pos = self.above_cutoff_pos + min(self.band_pos, k_band)
            p_obs = f"{observed_pos / k_eff:.4f}"
            p_min = f"{min_pos / k_eff:.4f}"
            p_max = f"{max_pos / k_eff:.4f}"

        ranks: dict[str, str] = {}
        percentiles: dict[str, str] = {}
        for cid, gene in POS_GENE.items():
            if cid in self.positive_scores:
                rank = self.rank_counts[cid] + 1
                ranks[gene] = str(rank)
                percentiles[gene] = f"{100.0 * (1.0 - ((rank - 1) / self.n_total)):.4f}"
            else:
                ranks[gene] = ""
                percentiles[gene] = ""

        return {
            "cohort": cohort,
            "scope": scope,
            "universe": universe,
            "axis": axis,
            "n_total": str(self.n_total),
            "n_pos": str(n_pos),
            "n_neg": str(self.n_neg),
            "auc": auc,
            "tie_band@100": str(self.tie_band if self.cutoff_score is not None else ""),
            "P@100_obs": p_obs,
            "P@100_min": p_min,
            "P@100_max": p_max,
            "ESR1_rank": ranks["ESR1"],
            "EGFLAM_rank": ranks["EGFLAM"],
            "GATA3_rank": ranks["GATA3"],
            "ESR1_percentile": percentiles["ESR1"],
            "EGFLAM_percentile": percentiles["EGFLAM"],
            "GATA3_percentile": percentiles["GATA3"],
        }


def make_metrics() -> dict[tuple[str, str, str], RunningMetric]:
    metrics: dict[tuple[str, str, str], RunningMetric] = {}
    for scope in SCOPES:
        for universe in ("exact_only", "exact_plus_proximal_close", "class_exact"):
            for axis in AXES:
                metrics[(scope, universe, axis)] = RunningMetric()
        for ev in ("proximal_close", "proximal", "regional", "unobserved"):
            metrics[(scope, f"class_{ev}", "v25_sigmoid")] = RunningMetric()
    return metrics


def stream_jsonl(path: Path):
    with path.open() as fh:
        for line in fh:
            yield json.loads(line)


def scopes_for(chrom: str) -> tuple[str, ...]:
    if chrom in CHR_SUBSET:
        return ("wg", "chr5_6_10")
    return ("wg",)


def first_pass(
    scored_path: Path,
    positives: set[str],
    limma_t: dict[str, float],
    metrics: dict[tuple[str, str, str], RunningMetric],
) -> None:
    for rec in stream_jsonl(scored_path):
        candidate_id = rec["candidate"]["candidate_id"]
        evidence_class = rec["observation"]["evidence_class"]
        is_positive = candidate_id in positives
        full_scores: dict[str, float] | None = None
        sigmoid_score: float | None = None
        for scope in scopes_for(rec["candidate"]["chrom"]):
            for universe, axes in universe_axes_for(evidence_class):
                if axes == SIGMOID_ONLY:
                    if sigmoid_score is None:
                        sigmoid_score = v25_sigmoid_score(rec)
                    axis_scores = {"v25_sigmoid": sigmoid_score}
                else:
                    if full_scores is None:
                        full_scores = score_record(rec, limma_t)
                    axis_scores = {axis: full_scores[axis] for axis in axes}
                for axis, score in axis_scores.items():
                    metrics[(scope, universe, axis)].observe_first_pass(
                        candidate_id, score, is_positive
                    )
    for metric in metrics.values():
        metric.finalize_first_pass()


def second_pass(
    scored_path: Path,
    positives: set[str],
    limma_t: dict[str, float],
    metrics: dict[tuple[str, str, str], RunningMetric],
) -> None:
    for rec in stream_jsonl(scored_path):
        candidate_id = rec["candidate"]["candidate_id"]
        evidence_class = rec["observation"]["evidence_class"]
        is_positive = candidate_id in positives
        full_scores: dict[str, float] | None = None
        sigmoid_score: float | None = None
        for scope in scopes_for(rec["candidate"]["chrom"]):
            for universe, axes in universe_axes_for(evidence_class):
                if axes == SIGMOID_ONLY:
                    if sigmoid_score is None:
                        sigmoid_score = v25_sigmoid_score(rec)
                    axis_scores = {"v25_sigmoid": sigmoid_score}
                else:
                    if full_scores is None:
                        full_scores = score_record(rec, limma_t)
                    axis_scores = {axis: full_scores[axis] for axis in axes}
                for axis, score in axis_scores.items():
                    metrics[(scope, universe, axis)].observe_second_pass(
                        candidate_id, score, is_positive
                    )


def run_cohort(label: str, scored_path: Path, limma_path: Path, positives: set[str]) -> list[dict[str, str]]:
    scored_path = scored_path if scored_path.is_absolute() else REPO / scored_path
    limma_path = limma_path if limma_path.is_absolute() else REPO / limma_path
    print(f"cohort: {label}")
    print(f"  scored: {scored_path.relative_to(REPO)}")
    print(f"  limma : {limma_path.relative_to(REPO)}")
    limma_t = load_limma(limma_path)
    print(f"  limma probes: {len(limma_t):,}")
    metrics = make_metrics()
    print("  pass 1/2")
    first_pass(scored_path, positives, limma_t, metrics)
    print("  pass 2/2")
    second_pass(scored_path, positives, limma_t, metrics)

    rows: list[dict[str, str]] = []
    for (scope, universe, axis), metric in sorted(metrics.items()):
        if metric.n_total == 0:
            continue
        rows.append(metric.row(label, scope, universe, axis))
    return rows


def write_tsv(rows: list[dict[str, str]], path: Path) -> None:
    columns = [
        "cohort",
        "scope",
        "universe",
        "axis",
        "n_total",
        "n_pos",
        "n_neg",
        "auc",
        "tie_band@100",
        "P@100_obs",
        "P@100_min",
        "P@100_max",
        "ESR1_rank",
        "EGFLAM_rank",
        "GATA3_rank",
        "ESR1_percentile",
        "EGFLAM_percentile",
        "GATA3_percentile",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(col, "") for col in columns) + "\n")


def fmt_auc(row: dict[str, str]) -> str:
    return f"{float(row['auc']):.3f}" if row.get("auc") else "—"


def fmt_int(row: dict[str, str], key: str) -> str:
    return f"{int(row[key]):,}" if row.get(key) else "—"


def fmt_pct(row: dict[str, str], key: str) -> str:
    return f"{float(row[key]):.2f}" if row.get(key) else "—"


def positive_genes(row: dict[str, str]) -> str:
    genes = [gene for gene in GENE_ORDER if row.get(f"{gene}_rank")]
    return ", ".join(genes) if genes else "none"


def lookup(rows: list[dict[str, str]], cohort: str, scope: str, universe: str, axis: str) -> dict[str, str] | None:
    for row in rows:
        if (
            row["cohort"] == cohort
            and row["scope"] == scope
            and row["universe"] == universe
            and row["axis"] == axis
        ):
            return row
    return None


def write_md(rows: list[dict[str, str]], path: Path) -> None:
    cohorts = list(dict.fromkeys(row["cohort"] for row in rows))
    md = [
        "# EvidenceClass-stratified benchmark",
        "",
        "Generated by `scripts/evidence_class_stratified_benchmark.py`.",
        "Scores are evaluated against the three Roth Fig. 5d validated",
        "positives under EvidenceClass-controlled candidate universes.",
        "V2.5-sigmoid and Delta-beta-only are recomputed from the committed",
        "V2.5-diff scored JSONLs; limma-style scores are joined from the",
        "probe-level moderated-t TSVs. This is a validated-label rank-lift",
        "analysis, not prospective target validation.",
        "",
        "The full-WG panel in §5.2.2 remains the all-candidate comparison.",
        "This artifact controls the candidate universe by EvidenceClass; it",
        "therefore answers a different question: whether the recommended axis",
        "still behaves well once probe-distance confidence is held fixed or",
        "restricted to high-confidence candidates.",
        "",
        "## EXACT-only vs high-confidence universe (WG)",
        "",
        "In the current candidate mapping, none of the three exact Roth",
        "validated candidate IDs is `EXACT` EvidenceClass: ESR1 is",
        "`PROXIMAL_CLOSE`, GATA3 is `PROXIMAL`, and EGFLAM is `REGIONAL`",
        "on the reported cohort paths. Therefore, EXACT-only has `n_pos = 0`",
        "and is a denominator audit rather than an evaluable rank-lift",
        "endpoint; EXACT+PROXIMAL_CLOSE evaluates ESR1 only.",
        "",
        "| cohort | universe | n_pos | positives in universe | Delta-beta AUC | V2.5-diff AUC | V2.5-sigmoid AUC | limma-style AUC | V2.5-sigmoid tie@100 | V2.5-sigmoid ESR1 %ile |",
        "|---|---|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for cohort in cohorts:
        for universe in ("exact_only", "exact_plus_proximal_close"):
            row_by_axis = {
                axis: lookup(rows, cohort, "wg", universe, axis) for axis in AXES
            }
            sig = row_by_axis["v25_sigmoid"]
            if sig is None:
                continue
            md.append(
                f"| {cohort} | `{universe}` | "
                f"{sig['n_pos']} | {positive_genes(sig)} | "
                f"{fmt_auc(row_by_axis['delta_beta_only'])} | "
                f"{fmt_auc(row_by_axis['v25_diff'])} | "
                f"{fmt_auc(sig)} | "
                f"{fmt_auc(row_by_axis['limma_style'])} | "
                f"{fmt_int(sig, 'tie_band@100')} | "
                f"{fmt_pct(sig, 'ESR1_percentile')} |"
            )

    md.extend(
        [
            "",
            "## Within-class V2.5-sigmoid AUC (WG)",
            "",
            "The three validated positives occupy three different non-EXACT",
            "EvidenceClass bins in the candidate mapping: ESR1 in",
            "`PROXIMAL_CLOSE`, GATA3 in `PROXIMAL`, and EGFLAM in `REGIONAL`.",
            "The per-class rows below are therefore single-positive rank-lift",
            "checks, not independent three-positive replications.",
            "",
            "| cohort | universe | n_total | n_pos | V2.5-sigmoid AUC | tie@100 | ESR1 %ile | EGFLAM %ile | GATA3 %ile |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for cohort in cohorts:
        for universe in [f"class_{ev}" for ev in EVIDENCE_CLASSES]:
            row = lookup(rows, cohort, "wg", universe, "v25_sigmoid")
            if row is None:
                continue
            md.append(
                f"| {cohort} | `{universe}` | {int(row['n_total']):,} | {row['n_pos']} | "
                f"{fmt_auc(row)} | {fmt_int(row, 'tie_band@100')} | "
                f"{fmt_pct(row, 'ESR1_percentile')} | "
                f"{fmt_pct(row, 'EGFLAM_percentile')} | "
                f"{fmt_pct(row, 'GATA3_percentile')} |"
            )

    md.extend(
        [
            "",
            "## Reproduce",
            "",
            "```bash",
            "uv run python scripts/evidence_class_stratified_benchmark.py",
            "```",
            "",
            "This reads the gitignored multi-GB WG scored JSONLs listed in the",
            "PAPER.md reproducibility appendix. The TSV contains the full",
            "cohort × scope × universe × axis grid, including chr5/6/10 rows.",
            "",
        ]
    )
    path.write_text("\n".join(md))


def default_specs() -> list[tuple[str, Path, Path]]:
    return [
        (
            "GSE322563 HM450",
            DERIVED / "scored_gse322563_wg_differential.jsonl",
            DERIVED / "limma_gse322563_probes.tsv",
        ),
        (
            "GSE322563 native v2",
            DERIVED / "scored_gse322563_native_wg_differential.jsonl",
            DERIVED / "limma_gse322563_probes.tsv",
        ),
        (
            "GSE77348",
            DERIVED / "scored_surrogate_wg_differential.jsonl",
            DERIVED / "limma_gse77348_probes.tsv",
        ),
        (
            "GSE69914 (tissue)",
            DERIVED / "scored_gse69914_wg_differential.jsonl",
            DERIVED / "limma_gse69914_probes.tsv",
        ),
    ]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--cohort",
        action="append",
        metavar="LABEL:SCORED_JSONL:LIMMA_TSV",
        help="Repeatable. Omit to run the four default WG cohort paths.",
    )
    ap.add_argument(
        "--positives",
        type=Path,
        default=DERIVED / "positives_roth_validated.txt",
    )
    ap.add_argument(
        "--output-prefix",
        type=Path,
        default=EXAMPLES / "evidence_class_stratified_benchmark",
    )
    args = ap.parse_args()

    positives = {
        line.strip()
        for line in args.positives.read_text().splitlines()
        if line.strip()
    }
    print(f"positives: {len(positives)}")

    if args.cohort:
        specs = []
        for item in args.cohort:
            label, scored, limma = item.split(":", 2)
            specs.append((label, Path(scored), Path(limma)))
    else:
        specs = default_specs()

    rows: list[dict[str, str]] = []
    for label, scored_path, limma_path in specs:
        if not scored_path.exists():
            raise FileNotFoundError(scored_path)
        if not limma_path.exists():
            raise FileNotFoundError(limma_path)
        rows.extend(run_cohort(label, scored_path, limma_path, positives))

    write_tsv(rows, args.output_prefix.with_suffix(".tsv"))
    write_md(rows, args.output_prefix.with_suffix(".md"))
    print(f"wrote {args.output_prefix.with_suffix('.tsv')}")
    print(f"wrote {args.output_prefix.with_suffix('.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
