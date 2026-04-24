"""Build a native EPIC v2 probe annotation in hg19 coordinates.

Source: GEO platform record GPL33022_family.soft.gz (Illumina Infinium
MethylationEPIC v2, ~937K probes, GRCh38 coordinates).

Output: data/raw/probes_hg19_epic_v2.tsv with columns
    probe_id    chrom    pos
where probe_id keeps the EPIC v2 beadchip-design suffix (`_BC##`,
`_TC##`, `_TO##`, `_BO##`) — the canonical native form. pos is 1-based
hg19 lifted from hg38 MAPINFO via pyliftover + UCSC's hg38ToHg19 chain
file (fetched on demand).

Restricts to chr5 / chr6 / chr10 — the working catalog scope. Probes
whose hg38 coordinate fails the liftover are skipped (no silent
fallback to hg38 coordinates).

Usage:

    pip install pyliftover            # one-time external dep
    python scripts/build_epic_v2_probes.py

Prereqs: data/raw/epic_v2/GPL33022_family.soft.gz present (fetched
separately from the GEO FTP).

This script exists specifically to replace the HM450-intersect
shortcut used in the earlier GSE322563 ingest (scripts/build_gse322563_cohort.py).
The native path keeps every EPIC v2 probe's canonical ID and
coordinate, removing the "80.7 % of HM450 probes retained" caveat
from the Roth-comparable pipeline (PAPER.md §4.4).
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

from pyliftover import LiftOver  # type: ignore

_DEFAULT_CHROMS = "chr5,chr6,chr10"


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--soft", required=True, type=Path,
                   help="GPL33022_family.soft.gz")
    p.add_argument("--output", required=True, type=Path,
                   help="Output TSV (probe_id, chrom, pos) in hg19")
    p.add_argument("--chain", default="hg38",
                   help="Source build name for pyliftover (hg38 → hg19)")
    p.add_argument(
        "--chroms", default=_DEFAULT_CHROMS,
        help=("Comma-separated chromosome filter (e.g. 'chr5,chr6,chr10') "
              "or 'all' for whole-genome (chr1..chr22, chrX, chrY). "
              f"Default: {_DEFAULT_CHROMS}"),
    )
    args = p.parse_args()
    if args.chroms.lower() == "all":
        CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    else:
        CHROMS = {c.strip() for c in args.chroms.split(",") if c.strip()}
    print(f"chromosome filter: {sorted(CHROMS)}", flush=True)

    print(f"loading {args.chain} → hg19 liftover chain...", flush=True)
    lo = LiftOver(args.chain, "hg19")  # downloads chain on first use

    print(f"streaming {args.soft.name}...", flush=True)
    n_read = 0
    n_scoped = 0
    n_mapped = 0
    n_failed = 0
    args.output.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(args.soft, "rt") as f, args.output.open("w") as out:
        out.write("probe_id\tchrom\tpos\n")

        # fast-forward to the data block
        for line in f:
            if line.startswith("!platform_table_begin"):
                break
        header = next(f).rstrip("\n").split("\t")
        try:
            id_col = header.index("ID")
            chr_col = header.index("CHR")
            pos_col = header.index("MAPINFO")
            build_col = header.index("Genome_Build")
        except ValueError as e:
            raise RuntimeError(
                f"GPL33022 SOFT table missing expected columns: {e}"
            ) from e

        for line in f:
            if line.startswith("!platform_table_end"):
                break
            parts = line.rstrip("\n").split("\t")
            n_read += 1
            probe_id = parts[id_col]
            chrom = parts[chr_col]
            if not chrom:
                continue
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            if chrom not in CHROMS:
                continue
            n_scoped += 1
            try:
                pos_hg38 = int(parts[pos_col])
            except (ValueError, IndexError):
                n_failed += 1
                continue
            build = parts[build_col] if build_col < len(parts) else ""
            if build not in ("", "GRCh38", "38"):
                # Unexpected build; skip rather than silently mishandle.
                n_failed += 1
                continue
            # pyliftover uses 0-based coords; SOFT MAPINFO is 1-based.
            hits = lo.convert_coordinate(chrom, pos_hg38 - 1)
            if not hits:
                n_failed += 1
                continue
            new_chrom, new_pos0, *_ = hits[0]
            if new_chrom not in CHROMS:
                # Cross-chromosome lift — reject.
                n_failed += 1
                continue
            new_pos1 = new_pos0 + 1
            out.write(f"{probe_id}\t{new_chrom}\t{new_pos1}\n")
            n_mapped += 1
            if n_read % 100_000 == 0:
                print(f"  read {n_read:,} probes...", flush=True)

    print(f"total probes in SOFT table:     {n_read:,}", flush=True)
    print(f"in scope ({sorted(CHROMS)}):    {n_scoped:,}", flush=True)
    print(f"lifted hg38 → hg19 successfully: {n_mapped:,}", flush=True)
    print(f"liftover failures / skipped:     {n_failed:,}", flush=True)
    print(f"→ {args.output}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
