"""`thermocas` CLI.

Subcommands:

    thermocas build-catalog --reference REF.fa --pam-model PAM.yaml --output catalog.jsonl
    thermocas score-cohort  --catalog catalog.jsonl --cohort cohort.yaml \
                            --backend local --probe-annotation A.tsv \
                            --tumor-beta T.tsv --normal-beta N.tsv \
                            --output scored.jsonl
    thermocas aggregate     --scored COHORT1=scored1.jsonl COHORT2=scored2.jsonl ... \
                            --output panatlas.jsonl

Designed so each stage is restartable: every input/output is a file, every
output is JSONL, no hidden global state.
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Iterable
from pathlib import Path

from thermocas import __version__
from thermocas.catalog import stream_catalog
from thermocas.cohort import score_cohort
from thermocas.config import load_cohort_config
from thermocas.io import read_jsonl, write_jsonl_atomic
from thermocas.methylation_backend import (
    LocalArrayBackend,
    LocalSummaryBackend,
    MethylationBackend,
)
from thermocas.models import CandidateSite, CohortConfig, ScoredCandidate
from thermocas.pam_model import PamModel
from thermocas.pan_cancer import aggregate


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.command is None:
        parser.print_help()
        return 0
    try:
        return args.func(args)
    except Exception as e:
        print(f"thermocas: error: {e}", file=sys.stderr)
        return 1


def _positive_int(value: str) -> int:
    """argparse `type=` for an integer ≥ 1."""
    try:
        v = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected integer, got {value!r}") from None
    if v < 1:
        raise argparse.ArgumentTypeError(f"expected integer >= 1, got {v}")
    return v


def _nonneg_int(value: str) -> int:
    """argparse `type=` for an integer ≥ 0."""
    try:
        v = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"expected integer, got {value!r}") from None
    if v < 0:
        raise argparse.ArgumentTypeError(f"expected integer >= 0, got {v}")
    return v


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="thermocas",
        description="Methylome-guided ThermoCas9 target-site discovery framework.",
    )
    p.add_argument("-V", "--version", action="version", version=f"thermocas {__version__}")
    sub = p.add_subparsers(dest="command", metavar="<command>")

    # build-catalog
    bc = sub.add_parser(
        "build-catalog",
        help="Scan a reference FASTA for ThermoCas9-compatible PAM sites.",
    )
    bc.add_argument("--reference", required=True, type=Path,
                    help="Reference genome FASTA (.fa, .fa.gz, .fasta, .fasta.gz)")
    bc.add_argument("--pam-model", required=True, type=Path,
                    help="YAML PAM model (e.g. config/pam_model.yaml)")
    bc.add_argument("--probe-annotation", type=Path,
                    help="Optional probe annotation TSV; when provided, only candidates "
                         "within --probe-window-bp of an assayed probe are kept "
                         "(drops most would-be UNOBSERVED candidates)")
    bc.add_argument("--probe-window-bp", type=int, default=500,
                    help="Maximum bp distance from a probe to keep a candidate "
                         "(default 500, matching the V1 'regional' EvidenceClass cap)")
    bc.add_argument("--output", required=True, type=Path,
                    help="Output JSONL of CandidateSite records")
    bc.set_defaults(func=_cmd_build_catalog)

    # score-cohort
    sc = sub.add_parser(
        "score-cohort",
        help="Score a candidate catalog against one cohort's methylation data.",
    )
    sc.add_argument("--catalog", required=True, type=Path,
                    help="JSONL of CandidateSite records (from build-catalog)")
    sc.add_argument("--cohort", required=True, type=Path,
                    help="YAML cohort config (e.g. config/cohorts/brca_example.yaml)")
    sc.add_argument("--pam-model", required=True, type=Path,
                    help="YAML PAM model — must match the one used to build the catalog")
    sc.add_argument("--backend", choices=["local", "summary", "gdc"], default="local",
                    help="local = raw beta matrices; summary = per-probe summary TSVs "
                         "(the format gdc-fetch produces); gdc = unsupported live mode")
    sc.add_argument("--probe-annotation", type=Path,
                    help="(local + summary backends) TSV with probe_id, chrom, pos columns")
    sc.add_argument("--tumor-beta", type=Path,
                    help="(local backend) TSV beta matrix for the tumor cohort")
    sc.add_argument("--normal-beta", type=Path,
                    help="(local backend) TSV beta matrix for the normal cohort")
    sc.add_argument("--tumor-summary", type=Path,
                    help="(summary backend) TSV with probe_id, n, mean, q25, q75 for tumor")
    sc.add_argument("--normal-summary", type=Path,
                    help="(summary backend) TSV with probe_id, n, mean, q25, q75 for normal")
    sc.add_argument("--sample-subtypes", type=Path,
                    help=("(local backend, V2) optional TSV with sample_id, subtype "
                          "columns. When provided, the tumor matrix is split into "
                          "per-subtype submatrices and one output JSONL is written "
                          "per subtype, with cohort_name = '<cohort>::<subtype>'."))
    sc.add_argument("--probabilistic", action="store_true",
                    help="V2 — also compute ProbabilisticScore for each candidate")
    sc.add_argument("--spacer", action="store_true",
                    help="V3 — also compute SpacerScore (gRNA design-quality "
                         "heuristics on the 20-nt protospacer)")
    sc.add_argument("--output", required=True, type=Path,
                    help=("Output JSONL of ScoredCandidate records. With --sample-subtypes, "
                          "this is treated as a *prefix* and one file is written per subtype: "
                          "<output>.<subtype>.jsonl"))
    sc.set_defaults(func=_cmd_score_cohort)

    # inspect
    ip = sub.add_parser(
        "inspect",
        help="V3 — pretty-print summary stats for any JSONL artifact (catalog / scored / aggregate).",
    )
    ip.add_argument("file", type=Path, help="JSONL file produced by another subcommand")
    ip.add_argument("--top", type=_positive_int, default=10,
                    help="Show top-N records by score (default: 10, must be >=1)")
    ip.add_argument("--head", type=_nonneg_int, default=0,
                    help="Show first N raw records instead of top-N by score")
    ip.set_defaults(func=_cmd_inspect)

    # benchmark
    bn = sub.add_parser(
        "benchmark",
        help="V3 — evaluate ranking quality of a scored cohort against a positives list.",
    )
    bn.add_argument("--scored", required=True, type=Path,
                    help="JSONL of ScoredCandidate (from score-cohort)")
    bn.add_argument("--positives", required=True, type=Path,
                    help="Text file with one candidate_id per line (ground-truth positives)")
    bn.add_argument("--cohort-name", required=True,
                    help="Label written into the BenchmarkResult.cohort_name field")
    bn.add_argument("--top-k", type=int, default=10)
    bn.add_argument("--score-field", default="final_score",
                    choices=["final_score", "p_therapeutic_selectivity",
                             "spacer_final_score", "naive_selectivity"],
                    help="Which scalar to rank by. naive_selectivity = β_normal − β_tumor "
                         "(baseline ablation: if this matches final_score's AUC, the "
                         "framework's machinery isn't adding value).")
    bn.add_argument("--missing-score-policy", default="rank_last",
                    choices=["rank_last", "drop", "error"],
                    help="How to handle candidates lacking the requested sub-score "
                         "(default: rank_last — they count toward n_total but never "
                         "outrank a candidate that has the sub-score)")
    bn.add_argument("--held-out-chromosomes", nargs="*", default=[],
                    help="Chromosomes excluded from training; only candidates on "
                         "these chroms are evaluated unless --no-enforce-holdout")
    bn.add_argument("--no-enforce-holdout", action="store_true",
                    help="Record held_out_chromosomes as metadata but evaluate the "
                         "full scored set (not the cross-validation default)")
    bn.add_argument("--output", required=True, type=Path,
                    help="Output JSONL containing one BenchmarkResult record")
    bn.set_defaults(func=_cmd_benchmark)

    # gdc-fetch
    gd = sub.add_parser(
        "gdc-fetch",
        help="V2 — download a TCGA cohort from GDC and export "
             "LocalSummaryBackend-compatible per-probe summary TSVs.",
    )
    gd.add_argument("--project", required=True,
                    help="GDC project_id, e.g. TCGA-BRCA")
    gd.add_argument("--platform", default="HM450", choices=["HM450", "EPIC"])
    gd.add_argument("--sample-type", default="both",
                    choices=["both", "tumor", "normal"],
                    help="Which sample side to fetch")
    gd.add_argument("--cache-dir", required=True, type=Path,
                    help="Directory for raw GDC file cache (re-used across runs)")
    gd.add_argument("--output-dir", required=True, type=Path,
                    help="Output directory for summary TSVs")
    gd.add_argument("--probe-annotation", type=Path,
                    help="Optional probe-annotation TSV (probe_id, chrom, pos) to copy "
                         "into the output dir as probes.tsv. Required for the downstream "
                         "score-cohort --backend summary step.")
    gd.add_argument("--max-files", type=int,
                    help="Cap files per side for testing or quota control "
                         "(default: download all matching files)")
    gd.set_defaults(func=_cmd_gdc_fetch)

    # aggregate
    agg = sub.add_parser(
        "aggregate",
        help="Aggregate per-cohort scored candidates into a pan-cancer atlas.",
    )
    agg.add_argument(
        "--scored",
        required=True,
        nargs="+",
        metavar="COHORT=PATH",
        help="One or more COHORT_NAME=path/to/scored.jsonl pairs",
    )
    agg.add_argument("--high-score-threshold", type=float, default=0.30,
                     help="Cohort-level final_score threshold for 'addressable' (default 0.30)")
    agg.add_argument("--output", required=True, type=Path,
                     help="Output JSONL of PanCancerAggregate records")
    agg.set_defaults(func=_cmd_aggregate)

    return p


# ---------- command implementations ----------


def _cmd_build_catalog(args: argparse.Namespace) -> int:
    from thermocas.catalog import probe_window_filter

    pam_model = PamModel.from_yaml(args.pam_model)
    region_filter = None
    if args.probe_annotation is not None:
        region_filter = probe_window_filter(args.probe_annotation, args.probe_window_bp)
    candidates = stream_catalog(args.reference, pam_model, region_filter=region_filter)
    n = write_jsonl_atomic(args.output, candidates)
    suffix = (
        f" (probe-filtered, ±{args.probe_window_bp} bp)"
        if region_filter is not None else ""
    )
    print(f"build-catalog: wrote {n} candidates{suffix} → {args.output}")
    return 0


def _cmd_score_cohort(args: argparse.Namespace) -> int:
    pam_model = PamModel.from_yaml(args.pam_model)
    cohort = load_cohort_config(args.cohort)

    if args.sample_subtypes is not None:
        return _cmd_score_cohort_by_subtype(args, pam_model, cohort)

    backend = _build_backend(args)
    candidates = read_jsonl(args.catalog, CandidateSite)
    scored = score_cohort(
        candidates, backend, cohort, pam_model,
        compute_probabilistic=args.probabilistic,
        compute_spacer=args.spacer,
    )
    n = write_jsonl_atomic(args.output, scored)
    extras = []
    if args.probabilistic:
        extras.append("probabilistic")
    if args.spacer:
        extras.append("spacer")
    extras_str = f" (+ {', '.join(extras)})" if extras else ""
    print(f"score-cohort[{cohort.name}]: wrote {n} scored candidates{extras_str} → {args.output}")
    return 0


def _cmd_score_cohort_by_subtype(
    args: argparse.Namespace,
    pam_model: PamModel,
    cohort: CohortConfig,
) -> int:
    """V2 — fan out across subtypes, writing one JSONL per subtype."""

    if args.backend != "local":
        raise ValueError("--sample-subtypes is currently only supported for --backend local")
    for required, val in (
        ("--probe-annotation", args.probe_annotation),
        ("--tumor-beta", args.tumor_beta),
        ("--normal-beta", args.normal_beta),
    ):
        if val is None:
            raise ValueError(f"--sample-subtypes requires {required}")

    backends = LocalArrayBackend.split_by_subtype(
        probe_annotation=args.probe_annotation,
        tumor_beta=args.tumor_beta,
        normal_beta=args.normal_beta,
        sample_subtypes=args.sample_subtypes,
    )

    out_prefix = args.output
    total = 0
    extras_list: list[str] = []
    if args.probabilistic:
        extras_list.append("probabilistic")
    if args.spacer:
        extras_list.append("spacer")
    extras_str = f" (+ {', '.join(extras_list)})" if extras_list else ""
    for subtype, backend in sorted(backends.items()):
        sub_cohort = cohort.model_copy(update={"name": f"{cohort.name}::{subtype}"})
        candidates = read_jsonl(args.catalog, CandidateSite)
        scored = score_cohort(
            candidates, backend, sub_cohort, pam_model,
            compute_probabilistic=args.probabilistic,
            compute_spacer=args.spacer,
        )
        sub_out = out_prefix.with_name(f"{out_prefix.stem}.{subtype}{out_prefix.suffix}")
        n = write_jsonl_atomic(sub_out, scored)
        total += n
        print(f"score-cohort[{sub_cohort.name}]: wrote {n} scored candidates{extras_str} → {sub_out}")
    print(f"score-cohort[{cohort.name}]: {total} candidates across {len(backends)} subtype(s)")
    return 0


def _cmd_inspect(args: argparse.Namespace) -> int:
    """V3 — print a quick summary of any JSONL artifact, auto-detecting the record type.

    Streams the file: histograms run incrementally, top-N uses a heap so memory
    stays O(N + top_k) regardless of artifact size. This matches the framework's
    streaming design — `inspect` works on a 100 GB catalog without OOM.

    Detection peeks at the first record's keys: if it has a `candidate` field
    we treat it as ScoredCandidate; if it has `n_cohorts_observed` we treat
    it as PanCancerAggregate; if it has `precision_at_k` it's BenchmarkResult;
    otherwise it's CandidateSite.
    """

    import heapq
    import json

    from thermocas.io import _open_text  # noqa: PLC0415

    # First pass: detect the record type from the first non-empty line.
    with _open_text(args.file, "rt") as f:
        first_line = ""
        for ln in f:
            if ln.strip():
                first_line = ln.strip()
                break
    if not first_line:
        print(f"{args.file}: empty")
        return 0
    first = json.loads(first_line)
    if "candidate" in first and "components" in first:
        kind = "ScoredCandidate"
    elif "n_cohorts_observed" in first:
        kind = "PanCancerAggregate"
    elif "precision_at_k" in first:
        kind = "BenchmarkResult"
    elif "candidate_id" in first and "pam_family" in first:
        kind = "CandidateSite"
    else:
        kind = "unknown"

    # Second pass: stream. For --head, stop after N lines.
    if args.head > 0:
        # Don't even count the rest; just print the first N raw records.
        with _open_text(args.file, "rt") as f:
            count = 0
            for ln in f:
                if not ln.strip():
                    continue
                if count >= args.head:
                    break
                print(ln.rstrip())
                count += 1
        return 0

    n_total = 0
    if kind == "ScoredCandidate":
        # Match `evaluate_ranking` exactly: final_score descending,
        # candidate_id ascending inside ties. `heapq.nsmallest` preserves the
        # same O(N log K) behavior while relying on Python's native string
        # ordering, which correctly handles prefix relations like `a` < `ab`.
        def _iter_scored():
            nonlocal n_total
            with _open_text(args.file, "rt") as f:
                for ln in f:
                    if not ln.strip():
                        continue
                    n_total += 1
                    r = json.loads(ln)
                    yield {
                        "score": r.get("final_score", 0.0),
                        "candidate_id": r["candidate"]["candidate_id"],
                        "evidence": r["observation"]["evidence_class"],
                    }

        top = heapq.nsmallest(
            args.top,
            _iter_scored(),
            key=lambda x: (-x["score"], x["candidate_id"]),
        )
        print(f"{args.file}: {n_total} records ({kind})")
        print(f"  top {len(top)} by final_score:")
        print(f"  {'rank':>4}  {'candidate_id':<48}  {'evidence':<14}  {'score':>8}")
        for i, row in enumerate(top, 1):
            print(
                f"  {i:>4}  {row['candidate_id']:<48}  "
                f"{row['evidence']:<14}  {row['score']:>8.3f}"
            )

    elif kind == "PanCancerAggregate":
        heap: list[tuple[float, str, int, float]] = []
        with _open_text(args.file, "rt") as f:
            for ln in f:
                if not ln.strip():
                    continue
                n_total += 1
                r = json.loads(ln)
                key = (
                    r.get("pan_cancer_score", 0.0),
                    r["candidate_id"],
                    r["n_cohorts_observed"],
                    r["recurrence"],
                )
                if len(heap) < args.top:
                    heapq.heappush(heap, key)
                elif key > heap[0]:
                    heapq.heapreplace(heap, key)
        top = sorted(heap, reverse=True)
        print(f"{args.file}: {n_total} records ({kind})")
        print(f"  top {len(top)} by pan_cancer_score:")
        print(f"  {'rank':>4}  {'candidate_id':<48}  {'obs':>4}  {'pan':>6}  {'recur':>6}")
        for i, (pan, cid, obs, recur) in enumerate(top, 1):
            print(f"  {i:>4}  {cid:<48}  {obs:>4}  {pan:>6.3f}  {recur:>6.2f}")

    elif kind == "BenchmarkResult":
        # BenchmarkResult files are tiny (one record per cohort run); just print all.
        print(f"{args.file}: ({kind})")
        with _open_text(args.file, "rt") as f:
            for ln in f:
                if not ln.strip():
                    continue
                n_total += 1
                r = json.loads(ln)
                print(f"  cohort={r['cohort_name']}  n={r['n_total']}  pos={r['n_positives']}  "
                      f"P@{r['top_k']}={_fmt(r.get('precision_at_k'))}  "
                      f"R@{r['top_k']}={_fmt(r.get('recall_at_k'))}  "
                      f"AUC={_fmt(r.get('roc_auc'))}")

    elif kind == "CandidateSite":
        # Streaming histograms — O(distinct values) memory.
        by_chrom: dict[str, int] = {}
        by_family: dict[str, int] = {}
        with _open_text(args.file, "rt") as f:
            for ln in f:
                if not ln.strip():
                    continue
                n_total += 1
                r = json.loads(ln)
                by_chrom[r["chrom"]] = by_chrom.get(r["chrom"], 0) + 1
                by_family[r["pam_family"]] = by_family.get(r["pam_family"], 0) + 1
        print(f"{args.file}: {n_total} records ({kind})")
        print(f"  by chromosome: {dict(sorted(by_chrom.items()))}")
        print(f"  by PAM family: {dict(sorted(by_family.items()))}")
    else:
        # Unknown record type — just count.
        with _open_text(args.file, "rt") as f:
            n_total = sum(1 for ln in f if ln.strip())
        print(f"{args.file}: {n_total} records ({kind})")

    return 0


def _cmd_benchmark(args: argparse.Namespace) -> int:
    """V3 — evaluate ranking quality against a positives list."""

    from thermocas.benchmark import evaluate_ranking

    positives = {
        line.strip()
        for line in args.positives.read_text().splitlines()
        if line.strip() and not line.startswith("#")
    }
    if not positives:
        raise ValueError(f"{args.positives}: no positive candidate_ids found")

    scored = read_jsonl(args.scored, ScoredCandidate)
    result = evaluate_ranking(
        scored,
        positives=positives,
        cohort_name=args.cohort_name,
        top_k=args.top_k,
        held_out_chromosomes=list(args.held_out_chromosomes),
        enforce_holdout=not args.no_enforce_holdout,
        score_field=args.score_field,
        missing_score_policy=args.missing_score_policy,
    )
    write_jsonl_atomic(args.output, [result])
    print(
        f"benchmark[{args.cohort_name}]: n={result.n_total} "
        f"pos={result.n_positives} → "
        f"P@{result.top_k}={_fmt(result.precision_at_k)} "
        f"R@{result.top_k}={_fmt(result.recall_at_k)} "
        f"AUC={_fmt(result.roc_auc)} → {args.output}"
    )
    return 0


def _fmt(v: float | None) -> str:
    return f"{v:.3f}" if v is not None else "n/a"


def _cmd_gdc_fetch(args: argparse.Namespace) -> int:
    """V2 — download a cohort from GDC and export per-side summary TSVs."""

    import shutil

    args.output_dir.mkdir(parents=True, exist_ok=True)
    sides: list[tuple[str, str, str]] = []  # (label, sample_type, filename)
    if args.sample_type in ("tumor", "both"):
        sides.append(("tumor", "Primary Tumor", "tumor_summary.tsv"))
    if args.sample_type in ("normal", "both"):
        sides.append(("normal", "Solid Tissue Normal", "normal_summary.tsv"))

    from thermocas.methylation_backend import GDCBackend
    for label, sample_type, fname in sides:
        backend = GDCBackend(
            project_id=args.project,
            cache_dir=args.cache_dir,
            platform=args.platform,
            sample_type=sample_type,
            max_files=args.max_files,
        )
        out = args.output_dir / fname
        files = backend.list_files()
        cap = f" (capped to {args.max_files})" if args.max_files else ""
        print(f"gdc-fetch[{args.project}/{label}]: {len(files)} files matched{cap}")
        n = backend.export_summaries(out)
        print(f"gdc-fetch[{args.project}/{label}]: wrote {n} probe summaries → {out}")

    if args.probe_annotation is not None:
        dest = args.output_dir / "probes.tsv"
        shutil.copyfile(args.probe_annotation, dest)
        print(f"gdc-fetch[{args.project}]: copied probe annotation → {dest}")
    else:
        print(f"gdc-fetch[{args.project}]: no --probe-annotation supplied; "
              "downstream score-cohort will need a probes.tsv.")
    return 0


def _cmd_aggregate(args: argparse.Namespace) -> int:
    # Pass iterators through rather than `list(read_jsonl(...))` — the CLI
    # was preloading every cohort JSONL into memory before handing it to
    # `aggregate()`. `aggregate()` itself still has to group across cohorts
    # (requires O(N_candidates) memory for the sorted pan-cancer emission
    # — see its docstring), but the CLI no longer double-buffers.
    cohorts: dict[str, Iterable[ScoredCandidate]] = {}
    for spec in args.scored:
        if "=" not in spec:
            raise ValueError(
                f"--scored entries must be COHORT=path; got {spec!r}"
            )
        name, path = spec.split("=", 1)
        cohorts[name] = read_jsonl(Path(path), ScoredCandidate)

    aggregates = aggregate(cohorts, high_score_threshold=args.high_score_threshold)
    n = write_jsonl_atomic(args.output, aggregates)
    print(f"aggregate: {len(cohorts)} cohort(s) → {n} candidates → {args.output}")
    return 0


def _build_backend(args: argparse.Namespace) -> MethylationBackend:
    if args.backend == "local":
        missing = [
            name for name, val in (
                ("--probe-annotation", args.probe_annotation),
                ("--tumor-beta", args.tumor_beta),
                ("--normal-beta", args.normal_beta),
            ) if val is None
        ]
        if missing:
            raise ValueError(f"local backend requires: {', '.join(missing)}")
        return LocalArrayBackend(
            probe_annotation=args.probe_annotation,
            tumor_beta=args.tumor_beta,
            normal_beta=args.normal_beta,
        )
    if args.backend == "summary":
        missing = [
            name for name, val in (
                ("--probe-annotation", args.probe_annotation),
                ("--tumor-summary", args.tumor_summary),
                ("--normal-summary", args.normal_summary),
            ) if val is None
        ]
        if missing:
            raise ValueError(f"summary backend requires: {', '.join(missing)}")
        return LocalSummaryBackend(
            probe_annotation=args.probe_annotation,
            tumor_summary=args.tumor_summary,
            normal_summary=args.normal_summary,
        )
    if args.backend == "gdc":
        raise ValueError(
            "gdc backend is unsupported as a live MethylationBackend. Use "
            "`thermocas gdc-fetch` to materialize the cohort to disk, then "
            "rerun score-cohort with --backend summary."
        )
    raise ValueError(f"unknown backend: {args.backend}")


if __name__ == "__main__":
    raise SystemExit(main())
