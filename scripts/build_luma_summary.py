"""Build a LumA-only per-probe summary TSV from cached GDC files + a TCGA
PAM50 sample sheet.

Inputs:
    --cache-dir       directory of GDC-cached *.txt beta-value files
    --file-id-map     TSV with columns file_id, submitter_id (TCGA-XX-NNNN)
    --luma-submitters one TCGA-XX-NNNN per line (LumA case IDs)
    --output          output summary TSV: probe_id, n, mean, q25, q75

This script intentionally lives outside the framework — it's analysis glue.
The framework's GDCBackend doesn't yet support a "filter by submitter_id
allowlist" knob. If we do this kind of analysis often we should plumb it
into the framework proper.
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--file-id-map", required=True, type=Path)
    parser.add_argument("--luma-submitters", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    luma_submitters = {
        line.strip()
        for line in args.luma_submitters.read_text().splitlines()
        if line.strip()
    }
    print(f"loaded {len(luma_submitters)} LumA submitter IDs")

    file_to_submitter: dict[str, str] = {}
    with args.file_id_map.open() as f:
        next(f)  # skip header
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                file_to_submitter[parts[0]] = parts[1]
    print(f"loaded {len(file_to_submitter)} file_id → submitter mappings")

    luma_file_ids = {
        fid for fid, sub in file_to_submitter.items()
        if sub in luma_submitters
    }
    print(f"intersection: {len(luma_file_ids)} LumA methylation files")

    # Find which of those are actually cached
    cached_paths: list[Path] = []
    for fid in luma_file_ids:
        p = args.cache_dir / f"{fid}.txt"
        if p.exists() and p.stat().st_size > 0:
            cached_paths.append(p)
    print(f"of those, {len(cached_paths)} are present in {args.cache_dir}")

    if not cached_paths:
        print("Nothing to summarize; abort.")
        return 1

    # Aggregate per-probe across the cached LumA files
    per_probe: dict[str, list[float | None]] = {}
    for path in cached_paths:
        for row in path.read_text().splitlines():
            if not row or row.startswith("#"):
                continue
            parts = row.split("\t")
            if len(parts) < 2:
                continue
            probe_id = parts[0]
            if probe_id == "Composite Element REF":
                continue
            cell = parts[1].strip()
            if cell == "" or cell.upper() in {"NA", "NAN", "NULL"}:
                per_probe.setdefault(probe_id, []).append(None)
            else:
                try:
                    v = float(cell)
                except ValueError:
                    continue
                if 0.0 <= v <= 1.0:
                    per_probe.setdefault(probe_id, []).append(v)

    # Write summary TSV in the format LocalSummaryBackend consumes
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as g:
        g.write("probe_id\tn\tmean\tq25\tq75\n")
        for pid in sorted(per_probe):
            betas = [b for b in per_probe[pid] if b is not None]
            n = len(betas)
            if n == 0:
                g.write(f"{pid}\t0\tNA\tNA\tNA\n")
                continue
            mean = statistics.fmean(betas)
            if n < 2:
                q25 = q75 = mean
            else:
                qs = statistics.quantiles(betas, n=4, method="inclusive")
                q25, _q50, q75 = qs[0], qs[1], qs[2]
                obs_min, obs_max = min(betas), max(betas)
                q25 = max(obs_min, min(obs_max, q25))
                q75 = max(obs_min, min(obs_max, q75))
                mean = max(q25, min(q75, mean))
            g.write(f"{pid}\t{n}\t{mean:.6f}\t{q25:.6f}\t{q75:.6f}\n")
    print(f"wrote {len(per_probe)} probe summaries → {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
