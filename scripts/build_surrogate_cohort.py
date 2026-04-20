"""Build a SURROGATE MCF-7 vs MCF-10A cohort from GEO series GSE77348.

IMPORTANT: this is NOT the Roth et al. Nature 2026 paper's dataset (that
paper cites GSE32256, which resolves to a Paramecium series — upstream
data-availability error flagged to the authors). GSE77348 is an independent
2016 breast cell-line methylation series (DNMT3B paper, same HM450 platform)
used here only as a **structural surrogate** to answer the V2 question:
does `p_therapeutic_selectivity` recover sensible ranking when the normal
comparator is a cell line that actually methylates the target promoters?

Inputs: GEO signal-intensity files (Unmethylated + Methylated per sample).
Outputs: `tumor_summary.tsv` (MCF-7) + `normal_summary.tsv` (MCF-10A)
in the per-probe summary format LocalSummaryBackend consumes.

We restrict to the untreated replicates (`MCF7_nt_rep1/2/3` and
`MCF10A_nt_rep1/2/3`) to avoid mixing in treatment arms (B2/B3 DNMT3B
overexpression / 4OHT etc.) that change methylation.
"""

from __future__ import annotations

import argparse
import gzip
import re
import statistics
from pathlib import Path


def identify_untreated_columns(header: list[str], cell_prefix: str) -> list[int]:
    """Return column indices (into `header`) of Unmethylated-signal columns for
    the cell line's untreated replicates.

    Heuristic: column name starts with `{cell_prefix}-<digit>_` (e.g. `MCF7-1_`,
    `MCF7-2_`, `MCF10A-1_`). Treatment arms have different column-name
    conventions (`MCF7_B2-overexpressing_<slide>_<position>`).
    """

    # Case-insensitive "signal" — MCF-7 intensity file uses lowercase, MCF-10A uses capital
    pattern = re.compile(rf"^{re.escape(cell_prefix)}-\d+_\S+_Unmethylated signal$", re.IGNORECASE)
    return [i for i, h in enumerate(header) if pattern.match(h)]


def load_intensities(path: Path, cell_prefix: str) -> tuple[list[str], dict[str, list[float | None]]]:
    """Return (per-sample column IDs, probe_id → list[β per sample]).

    β = M / (M + U + 100), the standard HM450 offset-stabilized estimator.
    Returns None for samples where a signal is missing.
    """

    with gzip.open(path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        u_idxs = identify_untreated_columns(header, cell_prefix)
        if not u_idxs:
            raise ValueError(f"no untreated-replicate columns found for {cell_prefix!r}")
        # For each untreated U-column at idx i, the matching M-column is i+1
        # (Unmethylated, Methylated, DetectionPval triples).
        sample_ids = [header[i].split("_")[0] for i in u_idxs]
        betas_by_probe: dict[str, list[float | None]] = {}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            probe_id = parts[0]
            row_betas: list[float | None] = []
            for u_idx in u_idxs:
                try:
                    u = float(parts[u_idx])
                    m = float(parts[u_idx + 1])
                except (IndexError, ValueError):
                    row_betas.append(None)
                    continue
                denom = m + u + 100.0
                if denom <= 0:
                    row_betas.append(None)
                    continue
                beta = m / denom
                if 0.0 <= beta <= 1.0:
                    row_betas.append(beta)
                else:
                    row_betas.append(None)
            betas_by_probe[probe_id] = row_betas
    return sample_ids, betas_by_probe


def summarize(betas_by_probe: dict[str, list[float | None]]) -> dict[str, tuple[int, float | None, float | None, float | None]]:
    """Compute per-probe (n, mean, q25, q75) — same semantics as LocalSummaryBackend."""

    out: dict[str, tuple[int, float | None, float | None, float | None]] = {}
    for pid, betas in betas_by_probe.items():
        clean = [b for b in betas if b is not None]
        n = len(clean)
        if n == 0:
            out[pid] = (0, None, None, None)
            continue
        mean = statistics.fmean(clean)
        if n < 2:
            q25 = q75 = mean
        else:
            qs = statistics.quantiles(clean, n=4, method="inclusive")
            q25, _q50, q75 = qs[0], qs[1], qs[2]
            obs_min, obs_max = min(clean), max(clean)
            q25 = max(obs_min, min(obs_max, q25))
            q75 = max(obs_min, min(obs_max, q75))
            mean = max(q25, min(q75, mean))
        out[pid] = (n, mean, q25, q75)
    return out


def write_summary_tsv(path: Path, summary: dict[str, tuple[int, float | None, float | None, float | None]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as g:
        g.write("probe_id\tn\tmean\tq25\tq75\n")
        def na(x): return "NA" if x is None else f"{x:.6f}"
        for pid in sorted(summary):
            n, mean, q25, q75 = summary[pid]
            g.write(f"{pid}\t{n}\t{na(mean)}\t{na(q25)}\t{na(q75)}\n")
    return len(summary)


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--mcf7-intensities", required=True, type=Path)
    p.add_argument("--mcf10a-intensities", required=True, type=Path)
    p.add_argument("--probe-annotation", required=True, type=Path,
                   help="probes_hg19.tsv — copied through as probes.tsv")
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    import shutil

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for cell, intens_path, out_name in (
        ("MCF7",   args.mcf7_intensities,   "tumor_summary.tsv"),
        ("MCF10A", args.mcf10a_intensities, "normal_summary.tsv"),
    ):
        print(f"processing {cell} from {intens_path} ...")
        sample_ids, betas = load_intensities(intens_path, cell)
        print(f"  found {len(sample_ids)} untreated replicate columns: {sample_ids}")
        summary = summarize(betas)
        n_written = write_summary_tsv(args.output_dir / out_name, summary)
        n_with_mean = sum(1 for (_, m, _, _) in summary.values() if m is not None)
        print(f"  wrote {n_written} probe summaries ({n_with_mean} with valid mean) → "
              f"{args.output_dir / out_name}")

    dest_probes = args.output_dir / "probes.tsv"
    shutil.copyfile(args.probe_annotation, dest_probes)
    print(f"copied probe annotation → {dest_probes}")

    # Emit a provenance note alongside the cohort so the label never gets lost.
    prov = args.output_dir / "PROVENANCE.md"
    prov.write_text(
        "# Surrogate matched cell-line cohort\n\n"
        "**NOT** the Roth et al. Nature 2026 samples — their cited GEO "
        "accession (GSE32256) resolves to unrelated Paramecium data and has "
        "been flagged to the authors.\n\n"
        "This cohort is built from **GSE77348** (Ostler et al., 2016, \"DNMT3B "
        "gene isoforms revealed specific methylation...\"), HM450 platform:\n"
        "  - MCF-7  untreated replicates: 3 samples\n"
        "  - MCF-10A untreated replicates: 3 samples\n\n"
        "Used as a **structural surrogate** to test whether V2 "
        "p_therapeutic_selectivity recovers under a cell-line-style matched "
        "comparator where the normal actually methylates target promoters. "
        "Not a substitute for the Roth validation data; do not publish "
        "benchmark numbers without this provenance caveat.\n"
    )
    print(f"wrote provenance → {prov}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
