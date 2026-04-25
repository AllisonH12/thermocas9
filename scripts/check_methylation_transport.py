"""Resolve Roth HEK/HCT target methylation transport before scoring.

This is an interface stub for the tag-A pre-registration. The implemented
script must compare Roth's BSS methylation labels in
`data/positives/positives_roth_hek_hct_v0.tsv` against independent ENCODE RRBS
measurements, then emit a transport-status TSV before any scored JSONL is
created.

Expected output columns:

    target_id
    gene
    cell_line
    expected_meth_state
    observed_beta
    coverage_replicates
    nearest_cpg_distance_bp
    transport_status
    reason

Allowed transport_status values:

    confirmed
    ambiguous
    low_coverage
    non_transportable
"""

from __future__ import annotations

import argparse
from pathlib import Path


DEFAULT_LABELS = Path("data/positives/positives_roth_hek_hct_v0.tsv")
DEFAULT_OUTPUT = Path("data/derived/roth_hek_hct_transport.tsv")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--methylated-beta-gt", type=float, default=0.7)
    parser.add_argument("--unmethylated-beta-lt", type=float, default=0.3)
    parser.add_argument("--min-reads", type=int, default=10)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--max-nearest-cpg-distance-bp", type=int, default=50)
    parser.add_argument(
        "--rrbs-bedmethyl",
        action="append",
        default=[],
        metavar="CELL_LINE=ACCESSION_OR_PATH",
        help=(
            "RRBS source for one replicate. Repeat for all sources, e.g. "
            "HEK293T=ENCFF001TMR. Fetch/materialization is not implemented yet."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    raise NotImplementedError(
        "RRBS transport parsing is not implemented yet. This stub freezes the "
        "CLI contract and output schema; implement it before tag prereg-transport. "
        f"labels={args.labels} output={args.output}"
    )


if __name__ == "__main__":
    raise SystemExit(main())
