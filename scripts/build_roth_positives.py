"""Build Roth validated-target positives sets from Fig. 5d coordinates.

Fig. 5d (Roth et al., Nature 2026) reports three experimentally validated
MCF-7/MCF-10A target sites with `5'-NNNNCGA-3'` PAM:

    EGFLAM T11 (control)   chr5  38,258,842     MCF-10A β 0.49 / MCF-7 β 0.01
    ESR1    T17            chr6  151,690,043    MCF-10A β 0.94 / MCF-7 β 0.07
    GATA3   T18            chr10 8,045,425      MCF-10A β 0.31 / MCF-7 β 0.02

These are the **validated supervision targets** against which V2.5 should be
benchmarked — not the broad gene-probe set, which is mostly gene-universe
membership and includes many probes that don't show the Roth pattern on
actual Roth samples (see GSE322563 diagnostic: 23% of Roth-gene probes
have β_n − β_t > 0.2).

Three positives sets produced, increasingly strict → lax:

    positives_roth_validated.txt     single closest NNNNCGA candidate per target (n=3)
    positives_roth_narrow.txt        NNNNCGA candidates within ±50 bp per target
    positives_roth_wide.txt          NNNNCGA candidates within ±500 bp per target

The ±50 and ±500 windows capture Roth's "promoter/target-window" framing
without relying on a specific spacer offset convention. Only NNNNCGA PAM
candidates are included — NNNNCCA matches would be off-biology since all
three validated targets are in the CGA family.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

# Roth Fig. 5d validated target coordinates.
#
# The paper prints **hg38** positions (confirmed by cross-check: the nearest
# HM450 probe at the paper's hg38 coord for ESR1/GATA3 does not exist in
# hg19 within 500 bp, and the β values at the hg19-equivalent-by-position
# do not match Roth's reported β). Ensembl REST
# /map/human/GRCh38/X:Y-Y/GRCh37 was used to liftOver each target to hg19
# so the existing hg19 catalog can be re-used:
#
#     target        hg38                   hg19                   Δ bp
#     EGFLAM_T11    chr5:38,258,842     → chr5:38,258,944          +102
#     ESR1_T17      chr6:151,690,043    → chr6:152,011,178     +321,135
#     GATA3_T18     chr10:8,045,425     → chr10:8,087,388        +41,963
#
# hg19 coordinates (1-based, as Ensembl returns) — catalog stores 0-indexed.
ROTH_TARGETS_1BASED: list[tuple[str, str, int]] = [
    ("EGFLAM_T11", "chr5",  38_258_944),
    ("ESR1_T17",   "chr6",  152_011_178),
    ("GATA3_T18",  "chr10",  8_087_388),
]

PAM_FAMILY = "NNNNCGA"
NARROW_BP = 50
WIDE_BP = 500


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--catalog", required=True, type=Path,
                   help="catalog_hg19_chr5_6_10.jsonl")
    p.add_argument("--out-dir", required=True, type=Path,
                   help="Directory to write positives_roth_*.txt")
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    # Per-target bucketing. Buckets are populated during the single catalog pass.
    exact: dict[str, tuple[int, str]] = {}  # target_name → (|Δ|, cid)
    narrow: dict[str, list[str]] = {name: [] for name, _, _ in ROTH_TARGETS_1BASED}
    wide: dict[str, list[str]] = {name: [] for name, _, _ in ROTH_TARGETS_1BASED}

    with args.catalog.open() as f:
        for line in f:
            r = json.loads(line)
            if r["pam_family"] != PAM_FAMILY:
                continue
            chrom = r["chrom"]
            ccp_0 = r["critical_c_pos"]            # 0-indexed
            ccp_1 = ccp_0 + 1                      # 1-indexed, matches paper
            for name, tchrom, tpos_1 in ROTH_TARGETS_1BASED:
                if tchrom != chrom:
                    continue
                d = abs(ccp_1 - tpos_1)
                if d > WIDE_BP:
                    continue
                cid = r["candidate_id"]
                if d <= NARROW_BP:
                    narrow[name].append(cid)
                wide[name].append(cid)
                prev = exact.get(name)
                if prev is None or d < prev[0]:
                    exact[name] = (d, cid)

    # Emit sets + a small report.
    exact_ids = [cid for (_, cid) in exact.values()]
    (args.out_dir / "positives_roth_validated.txt").write_text(
        "\n".join(sorted(exact_ids)) + "\n"
    )
    narrow_all = sorted({cid for ids in narrow.values() for cid in ids})
    (args.out_dir / "positives_roth_narrow.txt").write_text(
        "\n".join(narrow_all) + "\n"
    )
    wide_all = sorted({cid for ids in wide.values() for cid in ids})
    (args.out_dir / "positives_roth_wide.txt").write_text(
        "\n".join(wide_all) + "\n"
    )

    print("Roth validated-target positives — built from Fig. 5d coordinates")
    print("-" * 68)
    print(f"{'target':<14} {'chrom':<6} {'target_1b':>12}  {'nearest_cid':<32}  {'Δbp':>4}")
    for name, tchrom, tpos_1 in ROTH_TARGETS_1BASED:
        if name in exact:
            d, cid = exact[name]
            print(f"{name:<14} {tchrom:<6} {tpos_1:>12}  {cid:<32}  {d:>4}")
        else:
            print(f"{name:<14} {tchrom:<6} {tpos_1:>12}  (no candidate within ±{WIDE_BP} bp)")
    print()
    print(f"positives_roth_validated.txt  (n={len(exact_ids)}): one closest NNNNCGA candidate per target")
    print(f"positives_roth_narrow.txt     (n={len(narrow_all)}): NNNNCGA within ±{NARROW_BP} bp")
    print(f"positives_roth_wide.txt       (n={len(wide_all)}): NNNNCGA within ±{WIDE_BP} bp")

    # Per-target window counts — important for P@K interpretability.
    print()
    print(f"{'target':<14}  {'narrow':>6}  {'wide':>6}")
    for name, _, _ in ROTH_TARGETS_1BASED:
        print(f"{name:<14}  {len(narrow[name]):>6}  {len(wide[name]):>6}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
