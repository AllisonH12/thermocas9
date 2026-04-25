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
import csv
import gzip
import re
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path


ENCODE_ACCESSIONS: dict[str, tuple[str, ...]] = {
    "HEK293T": ("ENCFF001TMR", "ENCFF001TMQ", "ENCFF001TMS", "ENCFF001TMT"),
    "HCT116": ("ENCFF001TMM", "ENCFF001TMN"),
}


DEFAULT_LABELS = Path("data/positives/positives_roth_hek_hct_v0.tsv")
DEFAULT_OUTPUT = Path("data/derived/roth_hek_hct_transport.tsv")
DEFAULT_ENCODE_DIR = Path("data/raw/encode/rrbs")
BEDMETHYL_COLUMNS = 11


@dataclass(frozen=True)
class Source:
    cell_line: str
    accession_or_path: str
    path: Path


@dataclass(frozen=True)
class BedHit:
    accession: str
    pos0: int
    strand: str
    coverage: int
    beta: float
    distance_bp: int


@dataclass(frozen=True)
class Target:
    target_id: str
    gene: str
    chrom: str
    pos0: int
    strand: str
    meth_by_cell_line: dict[str, str]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--labels", type=Path, default=DEFAULT_LABELS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--methylated-beta-gt", type=float, default=0.7)
    parser.add_argument("--unmethylated-beta-lt", type=float, default=0.3)
    parser.add_argument("--min-reads", type=int, default=10)
    parser.add_argument("--min-replicates", type=int, default=2)
    parser.add_argument("--max-nearest-cpg-distance-bp", type=int, default=50)
    parser.add_argument("--encode-dir", type=Path, default=DEFAULT_ENCODE_DIR)
    parser.add_argument(
        "--rrbs-bedmethyl",
        action="append",
        default=[],
        metavar="CELL_LINE=ACCESSION_OR_PATH",
        help=(
            "RRBS source for one replicate. Repeat for all sources, e.g. "
            "HEK293T=ENCFF001TMR. Defaults to the Roth Methods ENCODE "
            "accessions when omitted."
        ),
    )
    return parser


def parse_sources(values: Iterable[str], encode_dir: Path) -> list[Source]:
    requested = list(values)
    if not requested:
        requested = [
            f"{cell_line}={accession}"
            for cell_line, accessions in ENCODE_ACCESSIONS.items()
            for accession in accessions
        ]

    sources: list[Source] = []
    for item in requested:
        try:
            cell_line, accession_or_path = item.split("=", 1)
        except ValueError as e:
            raise ValueError(f"bad --rrbs-bedmethyl {item!r}; expected CELL_LINE=ACCESSION_OR_PATH") from e
        cell_line = cell_line.strip()
        accession_or_path = accession_or_path.strip()
        if cell_line not in {"HEK293T", "HCT116"}:
            raise ValueError(f"unsupported cell line {cell_line!r}; expected HEK293T or HCT116")
        path = resolve_source(accession_or_path, encode_dir)
        sources.append(Source(cell_line=cell_line, accession_or_path=accession_or_path, path=path))
    return sources


def resolve_source(accession_or_path: str, encode_dir: Path) -> Path:
    path = Path(accession_or_path)
    if path.exists():
        return path

    if re.fullmatch(r"ENCFF[0-9A-Z]+", accession_or_path):
        encode_dir.mkdir(parents=True, exist_ok=True)
        out = encode_dir / f"{accession_or_path}.bed.gz"
        if not out.exists() or out.stat().st_size == 0:
            url = (
                "https://www.encodeproject.org/files/"
                f"{accession_or_path}/@@download/{accession_or_path}.bed.gz"
            )
            urllib.request.urlretrieve(url, out)
        return out

    raise FileNotFoundError(f"RRBS source {accession_or_path!r} is not a file or ENCODE accession")


def read_targets(path: Path) -> list[Target]:
    data_lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    targets: list[Target] = []
    for row in csv.DictReader(data_lines, delimiter="\t"):
        meth_by_cell_line = {
            "HEK293T": row["meth_state_HEK293T"],
            "HCT116": row["meth_state_HCT116"],
        }
        targets.append(
            Target(
                target_id=row["target_id"],
                gene=row["gene"],
                chrom=row["hg19_chrom"],
                pos0=int(row["hg19_pam_cytosine_pos"]) - 1,
                strand=row["hg19_strand"],
                meth_by_cell_line=meth_by_cell_line,
            )
        )
    return targets


def nearest_hits_for_source(
    source: Source,
    targets: list[Target],
    max_distance_bp: int,
    min_reads: int,
) -> dict[str, BedHit | None]:
    by_chrom: dict[str, list[Target]] = {}
    for target in targets:
        by_chrom.setdefault(target.chrom, []).append(target)

    best: dict[str, BedHit | None] = {target.target_id: None for target in targets}
    with gzip.open(source.path, "rt") as fh:
        for line in fh:
            if not line or line.startswith("track") or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < BEDMETHYL_COLUMNS:
                continue
            chrom = parts[0]
            chrom_targets = by_chrom.get(chrom)
            if not chrom_targets:
                continue
            try:
                pos0 = int(parts[1])
                strand = parts[5]
                coverage = int(parts[9])
                beta = float(parts[10]) / 100.0
            except (ValueError, IndexError):
                continue
            for target in chrom_targets:
                if coverage < min_reads:
                    continue
                distance = abs(pos0 - target.pos0)
                if distance > max_distance_bp:
                    continue
                current = best[target.target_id]
                # Prefer nearest. At equal distance, prefer same-strand evidence
                # because exact CpG PAMs should have a same-strand bedMethyl row.
                if current is None or _hit_key(distance, strand, target.strand) < _hit_key(
                    current.distance_bp, current.strand, target.strand
                ):
                    best[target.target_id] = BedHit(
                        accession=source.accession_or_path,
                        pos0=pos0,
                        strand=strand,
                        coverage=coverage,
                        beta=beta,
                        distance_bp=distance,
                    )
    return best


def _hit_key(distance: int, hit_strand: str, target_strand: str) -> tuple[int, int]:
    return (distance, 0 if hit_strand == target_strand else 1)


def classify_transport(
    expected_state: str,
    hits: list[BedHit],
    *,
    methylated_beta_gt: float,
    unmethylated_beta_lt: float,
    min_replicates: int,
) -> tuple[str, float | None, list[BedHit], str]:
    if expected_state == "na":
        return "not_applicable", None, [], "no Roth methylation state for this cell line"

    if len(hits) < min_replicates:
        return (
            "low_coverage",
            None,
            hits,
            f"{len(hits)} qualifying replicates after coverage gate; need {min_replicates}",
        )

    total_cov = sum(hit.coverage for hit in hits)
    observed_beta = sum(hit.beta * hit.coverage for hit in hits) / total_cov

    if unmethylated_beta_lt <= observed_beta <= methylated_beta_gt:
        return "ambiguous", observed_beta, hits, "observed beta is within the ambiguous interval"

    if expected_state == "methylated":
        if observed_beta > methylated_beta_gt:
            return "confirmed", observed_beta, hits, "observed beta agrees with Roth methylated label"
        return "non_transportable", observed_beta, hits, "observed beta is unmethylated but Roth label is methylated"

    if expected_state == "unmethylated":
        if observed_beta < unmethylated_beta_lt:
            return "confirmed", observed_beta, hits, "observed beta agrees with Roth unmethylated label"
        return "non_transportable", observed_beta, hits, "observed beta is methylated but Roth label is unmethylated"

    raise ValueError(f"unsupported methylation state {expected_state!r}")


def format_replicates(hits: list[BedHit]) -> str:
    if not hits:
        return "0"
    return ";".join(
        f"{hit.accession}:cov={hit.coverage}:beta={hit.beta:.4f}:dist={hit.distance_bp}:strand={hit.strand}"
        for hit in hits
    )


def write_transport_table(
    output: Path,
    targets: list[Target],
    sources: list[Source],
    *,
    methylated_beta_gt: float,
    unmethylated_beta_lt: float,
    min_reads: int,
    min_replicates: int,
    max_nearest_cpg_distance_bp: int,
) -> None:
    hits_by_source = {
        source: nearest_hits_for_source(source, targets, max_nearest_cpg_distance_bp, min_reads)
        for source in sources
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "target_id",
        "gene",
        "cell_line",
        "expected_meth_state",
        "observed_beta",
        "coverage_replicates",
        "nearest_cpg_distance_bp",
        "transport_status",
        "reason",
    ]
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for target in targets:
            for cell_line in ("HEK293T", "HCT116"):
                expected_state = target.meth_by_cell_line[cell_line]
                if expected_state == "na":
                    continue
                hits: list[BedHit] = []
                for source in sources:
                    if source.cell_line != cell_line:
                        continue
                    hit = hits_by_source[source][target.target_id]
                    if hit is not None:
                        hits.append(hit)
                status, observed_beta, qualifying, reason = classify_transport(
                    expected_state,
                    hits,
                    methylated_beta_gt=methylated_beta_gt,
                    unmethylated_beta_lt=unmethylated_beta_lt,
                    min_replicates=min_replicates,
                )
                nearest = min((hit.distance_bp for hit in qualifying), default="")
                writer.writerow(
                    {
                        "target_id": target.target_id,
                        "gene": target.gene,
                        "cell_line": cell_line,
                        "expected_meth_state": expected_state,
                        "observed_beta": "" if observed_beta is None else f"{observed_beta:.6f}",
                        "coverage_replicates": format_replicates(qualifying),
                        "nearest_cpg_distance_bp": nearest,
                        "transport_status": status,
                        "reason": reason,
                    }
                )


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    sources = parse_sources(args.rrbs_bedmethyl, args.encode_dir)
    targets = read_targets(args.labels)
    write_transport_table(
        args.output,
        targets,
        sources,
        methylated_beta_gt=args.methylated_beta_gt,
        unmethylated_beta_lt=args.unmethylated_beta_lt,
        min_reads=args.min_reads,
        min_replicates=args.min_replicates,
        max_nearest_cpg_distance_bp=args.max_nearest_cpg_distance_bp,
    )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
