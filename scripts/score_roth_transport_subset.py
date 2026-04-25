"""Score the transport-confirmed Roth HEK/HCT System B subset.

This script is intentionally narrower than `thermocas score-cohort`: the RRBS
transport check can assign different nearest covered cytosines to the two sides
of a candidate, while `LocalSummaryBackend` assumes one shared probe_id per
candidate. For the tag-C subset, score only the pre-registered Roth targets
with a side-specific nearest-covered-CpG policy matching the transport script.

Outputs:
  - data/derived/roth_hek_hct_subset_scores.tsv
  - data/derived/roth_hek_hct_subset_benchmark.tsv

No manuscript prose is emitted here.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from dataclasses import dataclass
from pathlib import Path

import check_methylation_transport as transport

from thermocas.config import load_cohort_config
from thermocas.evidence import classify_evidence
from thermocas.models import EvidenceClass, MethylationObservation
from thermocas.probabilistic import (
    p_differential_protection,
    p_gap_sigmoid,
    p_observation_trustworthy,
    p_protected_normal,
    p_targetable_tumor,
    probabilistic_score,
)

DEFAULT_HEK_COHORT = Path("config/cohorts/hek_target_hct_protected.yaml")
DEFAULT_HCT_COHORT = Path("config/cohorts/hct_target_hek_protected.yaml")
DEFAULT_LABELS = Path("data/positives/positives_roth_hek_hct_v0.tsv")
DEFAULT_TRANSPORT = Path("data/derived/roth_hek_hct_transport.tsv")
DEFAULT_OUTPUT_SCORES = Path("data/derived/roth_hek_hct_subset_scores.tsv")
DEFAULT_OUTPUT_BENCHMARK = Path("data/derived/roth_hek_hct_subset_benchmark.tsv")

AXES = (
    "p_targ_only",
    "tumor_only",
    "p_targ_x_p_prot",
    "v25_diff",
    "v25_sigmoid",
    "delta_beta_only",
)


@dataclass(frozen=True)
class SideSummary:
    beta_mean: float
    beta_q25: float
    beta_q75: float
    n: int
    min_distance_bp: int
    max_distance_bp: int
    replicates: str


@dataclass(frozen=True)
class Direction:
    label: str
    target_cell_line: str
    comparator_cell_line: str
    cohort_yaml: Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--transport", type=Path, default=DEFAULT_TRANSPORT)
    parser.add_argument("--hek-cohort", type=Path, default=DEFAULT_HEK_COHORT)
    parser.add_argument("--hct-cohort", type=Path, default=DEFAULT_HCT_COHORT)
    parser.add_argument("--output-scores", type=Path, default=DEFAULT_OUTPUT_SCORES)
    parser.add_argument("--output-benchmark", type=Path, default=DEFAULT_OUTPUT_BENCHMARK)
    parser.add_argument("--min-reads", type=int, default=10)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--max-nearest-cpg-distance-bp", type=int, default=50)
    parser.add_argument("--encode-dir", type=Path, default=transport.DEFAULT_ENCODE_DIR)
    return parser


def read_label_rows(path: Path) -> list[dict[str, str]]:
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def read_transport_status(path: Path) -> dict[tuple[str, str], str]:
    out: dict[tuple[str, str], str] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[(row["target_id"], row["cell_line"])] = row["transport_status"]
    return out


def summarize_hits(hits: list[transport.BedHit]) -> SideSummary | None:
    if not hits:
        return None
    total_cov = sum(hit.coverage for hit in hits)
    mean = sum(hit.beta * hit.coverage for hit in hits) / total_cov

    betas = [hit.beta for hit in hits]
    if len(betas) < 3:
        q25 = q75 = mean
    else:
        q25, _q50, q75 = statistics.quantiles(betas, n=4, method="inclusive")
        obs_min, obs_max = min(betas), max(betas)
        q25 = max(obs_min, min(obs_max, q25))
        q75 = max(obs_min, min(obs_max, q75))
        q25 = min(q25, mean)
        q75 = max(q75, mean)

    return SideSummary(
        beta_mean=mean,
        beta_q25=q25,
        beta_q75=q75,
        n=len(hits),
        min_distance_bp=min(hit.distance_bp for hit in hits),
        max_distance_bp=max(hit.distance_bp for hit in hits),
        replicates=transport.format_replicates(hits),
    )


def build_side_summaries(
    targets: list[transport.Target],
    *,
    min_reads: int,
    min_replicates: int,
    max_nearest_cpg_distance_bp: int,
    encode_dir: Path,
) -> dict[tuple[str, str], SideSummary | None]:
    sources = transport.parse_sources([], encode_dir)
    hits_by_source = {
        source: transport.nearest_hits_for_source(
            source,
            targets,
            max_nearest_cpg_distance_bp,
            min_reads,
        )
        for source in sources
    }

    out: dict[tuple[str, str], SideSummary | None] = {}
    for target in targets:
        for cell_line in ("HEK293T", "HCT116"):
            hits = [
                hit
                for source in sources
                if source.cell_line == cell_line
                for hit in [hits_by_source[source][target.target_id]]
                if hit is not None
            ]
            out[(target.target_id, cell_line)] = (
                summarize_hits(hits) if len(hits) >= min_replicates else None
            )
    return out


def make_observation(
    *,
    candidate_id: str,
    cohort_name: str,
    target_summary: SideSummary | None,
    comparator_summary: SideSummary | None,
    thresholds,
) -> MethylationObservation:
    if target_summary is None or comparator_summary is None:
        return MethylationObservation(
            candidate_id=candidate_id,
            cohort_name=cohort_name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )

    distance = max(target_summary.max_distance_bp, comparator_summary.max_distance_bp)
    evidence_class = classify_evidence(distance, thresholds)
    if evidence_class == EvidenceClass.UNOBSERVED:
        return MethylationObservation(
            candidate_id=candidate_id,
            cohort_name=cohort_name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )

    return MethylationObservation(
        candidate_id=candidate_id,
        cohort_name=cohort_name,
        evidence_class=evidence_class,
        evidence_distance_bp=distance,
        probe_id=(
            f"rrbs:target_maxdist={target_summary.max_distance_bp};"
            f"comparator_maxdist={comparator_summary.max_distance_bp}"
        ),
        beta_tumor_mean=target_summary.beta_mean,
        beta_tumor_q25=target_summary.beta_q25,
        beta_tumor_q75=target_summary.beta_q75,
        n_samples_tumor=target_summary.n,
        beta_normal_mean=comparator_summary.beta_mean,
        beta_normal_q25=comparator_summary.beta_q25,
        beta_normal_q75=comparator_summary.beta_q75,
        n_samples_normal=comparator_summary.n,
    )


def target_only_observation(
    *,
    candidate_id: str,
    cohort_name: str,
    target_summary: SideSummary | None,
    thresholds,
) -> MethylationObservation:
    if target_summary is None:
        return MethylationObservation(
            candidate_id=candidate_id,
            cohort_name=cohort_name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )
    evidence_class = classify_evidence(target_summary.max_distance_bp, thresholds)
    if evidence_class == EvidenceClass.UNOBSERVED:
        return MethylationObservation(
            candidate_id=candidate_id,
            cohort_name=cohort_name,
            evidence_class=EvidenceClass.UNOBSERVED,
        )
    return MethylationObservation(
        candidate_id=candidate_id,
        cohort_name=cohort_name,
        evidence_class=evidence_class,
        evidence_distance_bp=target_summary.max_distance_bp,
        probe_id=f"rrbs:target_maxdist={target_summary.max_distance_bp}",
        beta_tumor_mean=target_summary.beta_mean,
        beta_tumor_q25=target_summary.beta_q25,
        beta_tumor_q75=target_summary.beta_q75,
        n_samples_tumor=target_summary.n,
    )


def axis_scores(
    obs: MethylationObservation,
    target_obs: MethylationObservation,
    sigma_fixed: float | None,
) -> dict[str, float]:
    p_targ = p_targetable_tumor(target_obs)
    p_prot = p_protected_normal(obs)
    p_trust = p_observation_trustworthy(obs)
    delta = (
        0.0
        if obs.beta_tumor_mean is None or obs.beta_normal_mean is None
        else obs.beta_normal_mean - obs.beta_tumor_mean
    )
    sigma = 0.07071067811865477 if sigma_fixed is None else sigma_fixed
    return {
        "p_targ_only": p_targ,
        "tumor_only": probabilistic_score(obs, mode="tumor_only").p_therapeutic_selectivity,
        "p_targ_x_p_prot": p_targetable_tumor(obs) * p_prot * p_trust,
        "v25_diff": p_targetable_tumor(obs) * p_differential_protection(obs) * p_trust,
        "v25_sigmoid": p_targetable_tumor(obs) * p_gap_sigmoid(obs, sigma_fixed=sigma) * p_trust,
        "delta_beta_only": delta,
    }


def auc_midrank(rows: list[tuple[float, bool]]) -> float:
    n_pos = sum(1 for _score, is_pos in rows if is_pos)
    n_neg = len(rows) - n_pos
    if n_pos == 0 or n_neg == 0:
        return math.nan

    # Ascending ranks; higher scores are better, so this is standard ROC AUC.
    ordered = sorted(enumerate(rows), key=lambda x: x[1][0])
    rank_sum = 0.0
    i = 0
    while i < len(ordered):
        j = i + 1
        while j < len(ordered) and ordered[j][1][0] == ordered[i][1][0]:
            j += 1
        midrank = (i + 1 + j) / 2.0
        for _idx, (_score, is_pos) in ordered[i:j]:
            if is_pos:
                rank_sum += midrank
        i = j
    return (rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def fmt_float(value: float | None) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and math.isnan(value):
        return "nan"
    return f"{value:.6f}"


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    label_rows = read_label_rows(args.labels)
    targets = transport.read_targets(args.labels)
    statuses = read_transport_status(args.transport)
    summaries = build_side_summaries(
        targets,
        min_reads=args.min_reads,
        min_replicates=args.min_replicates,
        max_nearest_cpg_distance_bp=args.max_nearest_cpg_distance_bp,
        encode_dir=args.encode_dir,
    )

    directions = [
        Direction("HEK_target_HCT_protected", "HEK293T", "HCT116", args.hek_cohort),
        Direction("HCT_target_HEK_protected", "HCT116", "HEK293T", args.hct_cohort),
    ]
    cohorts = {d.label: load_cohort_config(d.cohort_yaml) for d in directions}

    score_rows: list[dict[str, str]] = []
    score_lookup: dict[tuple[str, str], dict[str, float]] = {}
    for direction in directions:
        cohort = cohorts[direction.label]
        for row in label_rows:
            target_id = row["target_id"]
            target_summary = summaries[(target_id, direction.target_cell_line)]
            comparator_summary = summaries[(target_id, direction.comparator_cell_line)]
            obs = make_observation(
                candidate_id=row["candidate_id"],
                cohort_name=cohort.name,
                target_summary=target_summary,
                comparator_summary=comparator_summary,
                thresholds=cohort.evidence_thresholds,
            )
            target_obs = target_only_observation(
                candidate_id=row["candidate_id"],
                cohort_name=cohort.name,
                target_summary=target_summary,
                thresholds=cohort.evidence_thresholds,
            )
            scores = axis_scores(obs, target_obs, cohort.sigma_fixed)
            score_lookup[(direction.label, target_id)] = scores
            target_status = statuses.get(
                (target_id, direction.target_cell_line), "not_applicable"
            )
            comparator_status = statuses.get(
                (target_id, direction.comparator_cell_line), "not_applicable"
            )
            selectivity_distance = (
                "" if obs.evidence_distance_bp is None else str(obs.evidence_distance_bp)
            )
            target_distance = (
                ""
                if target_obs.evidence_distance_bp is None
                else str(target_obs.evidence_distance_bp)
            )
            beta_target = None if target_summary is None else target_summary.beta_mean
            beta_comparator = (
                None if comparator_summary is None else comparator_summary.beta_mean
            )

            score_rows.append(
                {
                    "direction": direction.label,
                    "target_id": target_id,
                    "gene": row["gene"],
                    "candidate_id": row["candidate_id"],
                    "endpoint_membership": row["endpoint_membership"],
                    "editable_in_HEK293T": row["editable_in_HEK293T"],
                    "editable_in_HCT116": row["editable_in_HCT116"],
                    "target_cell_line": direction.target_cell_line,
                    "comparator_cell_line": direction.comparator_cell_line,
                    "target_transport_status": target_status,
                    "comparator_transport_status": comparator_status,
                    "selectivity_evidence_class": obs.evidence_class.value,
                    "selectivity_evidence_distance_bp": selectivity_distance,
                    "target_evidence_class": target_obs.evidence_class.value,
                    "target_evidence_distance_bp": target_distance,
                    "beta_target_mean": fmt_float(beta_target),
                    "beta_comparator_mean": fmt_float(beta_comparator),
                    "n_target": "" if target_summary is None else str(target_summary.n),
                    "n_comparator": "" if comparator_summary is None else str(comparator_summary.n),
                    "target_distance_minmax": distance_pair(target_summary),
                    "comparator_distance_minmax": distance_pair(comparator_summary),
                    **{axis: fmt_float(scores[axis]) for axis in AXES},
                }
            )

    args.output_scores.parent.mkdir(parents=True, exist_ok=True)
    with args.output_scores.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(score_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(score_rows)

    benchmark_rows = build_benchmark_rows(label_rows, statuses, score_lookup)
    with args.output_benchmark.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(benchmark_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(benchmark_rows)

    print(f"wrote {args.output_scores}")
    print(f"wrote {args.output_benchmark}")
    return 0


def distance_pair(summary: SideSummary | None) -> str:
    if summary is None:
        return ""
    return f"{summary.min_distance_bp}-{summary.max_distance_bp}"


def build_benchmark_rows(
    label_rows: list[dict[str, str]],
    statuses: dict[tuple[str, str], str],
    score_lookup: dict[tuple[str, str], dict[str, float]],
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []

    endpoint1_specs = [
        ("HEK_target_HCT_protected", "T4", "T5"),
        ("HCT_target_HEK_protected", "T5", "T4"),
    ]
    for direction, positive_target, negative_target in endpoint1_specs:
        for axis in AXES:
            scored = [
                (score_lookup[(direction, positive_target)][axis], True),
                (score_lookup[(direction, negative_target)][axis], False),
            ]
            rows.append(
                {
                    "endpoint": "E1_transport_confirmed_directionality",
                    "direction": direction,
                    "axis": axis,
                    "n_pos": "1",
                    "n_neg": "1",
                    "auc": fmt_float(auc_midrank(scored)),
                    "positive_targets": positive_target,
                    "negative_targets": negative_target,
                    "note": "T9 discriminator excluded by prereg-transport-addendum",
                }
            )

    # Endpoint 2 is an HEK293T editability diagnostic. Use only targets whose
    # HEK293T transport row is confirmed, excluding low-coverage T3/T9 and
    # ambiguous T13 per the addendum.
    hek_confirmed = [
        row
        for row in label_rows
        if "E2" in row["endpoint_membership"]
        and statuses.get((row["target_id"], "HEK293T")) == "confirmed"
        and row["editable_in_HEK293T"] in {"true", "false"}
    ]
    positives = [row["target_id"] for row in hek_confirmed if row["editable_in_HEK293T"] == "true"]
    negatives = [row["target_id"] for row in hek_confirmed if row["editable_in_HEK293T"] == "false"]
    for axis in ("p_targ_only",):
        scored = [
            (score_lookup[("HEK_target_HCT_protected", target_id)][axis], target_id in positives)
            for target_id in positives + negatives
        ]
        rows.append(
            {
                "endpoint": "E2_HEK293T_editability_confirmed",
                "direction": "HEK_target_HCT_protected",
                "axis": axis,
                "n_pos": str(len(positives)),
                "n_neg": str(len(negatives)),
                "auc": fmt_float(auc_midrank(scored)),
                "positive_targets": ",".join(positives),
                "negative_targets": ",".join(negatives),
                "note": "confirmed HEK293T transport only; T3/T9/T13 excluded",
            }
        )

    return rows


if __name__ == "__main__":
    raise SystemExit(main())
