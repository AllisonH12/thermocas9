#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Map probe-level limma-eBayes t-stats onto candidates + compute AUC panel.

Per the user's instruction: **probe-level first, then candidate-mapped.**
limma operates at the probe level (§7 Methods). We then assign each
ThermoCas9 candidate the moderated t-statistic of its nearest-probe (the
same probe the V2.5 composite already binds to via the
`evidence_distance_bp` / `EvidenceClass` path). Ranking is then by
*signed* t (higher = stronger tumor-hypomethylation vs normal-methylation),
using `−t_mod` as the ranking score (since our limma group assignment
was tumor=1 / normal=0, so `Δ = β_t − β_n` and a *targetable* site has
Δ < 0 → t < 0 → −t > 0).

Output per cohort: AUC + tie_band@{20,100,500,1000} + P@100 interval
+ per-positive ranks under limma-t vs V2.5-sigmoid and V2.5-diff,
on both chr5/6/10 and WG subsets.

Inputs:
  --scored PATH          scored JSONL for the cohort (produces probe_id
                         per candidate under observation.probe_id)
  --limma PATH           limma_{cohort}_probes.tsv from run_limma_per_cohort.py
  --positives PATH       positives_roth_validated.txt
  --output-prefix PATH   writes PREFIX.md + PREFIX.tsv
"""

from __future__ import annotations

import argparse
import bisect
import csv
import json
from pathlib import Path

CHR_SUBSET = {"chr5", "chr6", "chr10"}

POS_GENE = {
    "chr5:38258943+:NNNNCGA":  "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA":  "GATA3",
}


def strip_suffix(pid: str) -> str:
    """EPIC v2 probe IDs carry `_BC##` / `_TC##` / `_TO##` / `_BO##` beadchip
    suffixes; strip to the canonical `cg*` / `ch*` / `nv*` id so HM450-style
    lookup works."""
    if "_" in pid:
        return pid.split("_")[0]
    return pid


def load_limma(path: Path) -> dict[str, float]:
    """{stripped_probe_id: t_mod}. Robust to NaN (→ 0.0). For EPIC v2 probes
    that collapse to the same HM450 key under suffix-stripping, keep the
    *most-extreme-|t|* value — the probe that would drive the ranking most.
    That matches the cohort-build behaviour (which also collapses EPIC v2
    suffixes when intersecting to HM450)."""

    out: dict[str, float] = {}
    with path.open() as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row["probe_id"]
            try:
                t = float(row["t_mod"])
            except ValueError:
                t = 0.0
            if not (t == t) or t == float("inf") or t == float("-inf"):
                t = 0.0
            key = strip_suffix(pid)
            existing = out.get(key)
            if existing is None or abs(t) > abs(existing):
                out[key] = t
    return out


def stream_candidates(scored_path: Path, limma_t: dict[str, float], positives: set[str]):
    """Yield (cid, chrom, in_chr_subset, probe_id, limma_score, v25_diff_score, is_positive)."""
    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            cid = d["candidate"]["candidate_id"]
            chrom = d["candidate"]["chrom"]
            obs = d["observation"]
            probe_id = obs.get("probe_id")
            v25 = d["probabilistic"]["p_therapeutic_selectivity"]
            # limma score: signed −t (so higher = more Cas9-selectable).
            if probe_id is None:
                score = 0.0
            else:
                # limma_t is keyed by stripped HM450-canonical id; if the
                # scored JSONL carries a raw EPIC v2 id we strip on lookup.
                t = limma_t.get(probe_id)
                if t is None:
                    t = limma_t.get(strip_suffix(probe_id))
                score = -t if t is not None else 0.0
            yield cid, chrom, chrom in CHR_SUBSET, probe_id, score, v25, (cid in positives)


def midrank_auc_and_positions(sorted_desc: list[tuple[float, str]], positives: set[str]) -> tuple[float, dict[str, int], dict[str, float]]:
    """Returns (auc, {cid: 1-based rank}, {cid: midrank})."""
    n_total = len(sorted_desc)
    n_pos = len(positives)
    n_neg = n_total - n_pos

    rank_of: dict[str, int] = {}
    for idx, (_, cid) in enumerate(sorted_desc, start=1):
        if cid in positives:
            rank_of[cid] = idx

    # Mid-rank AUC.
    asc_scores = [-k for k, _ in sorted_desc]; asc_scores.reverse()
    asc_cids = [c for _, c in sorted_desc]; asc_cids.reverse()
    midrank: dict[str, float] = {}
    pos_rem = set(positives)
    i = 0
    while i < n_total and pos_rem:
        s = asc_scores[i]; j = i
        while j < n_total and asc_scores[j] == s:
            j += 1
        mid = (i + 1 + j) / 2.0
        for k in range(i, j):
            cid = asc_cids[k]
            if cid in pos_rem:
                midrank[cid] = mid; pos_rem.remove(cid)
                if not pos_rem:
                    break
        i = j
    if n_pos == 0 or n_neg == 0:
        return float("nan"), rank_of, midrank
    sum_ranks = sum(midrank.values())
    auc = (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)
    return auc, rank_of, midrank


def tie_band_at_k(sorted_desc: list[tuple[float, str]], k: int) -> int:
    if k > len(sorted_desc):
        return 0
    cutoff = sorted_desc[k - 1][0]
    return sum(1 for s, _ in sorted_desc if s == cutoff)


def p_at_k_interval(sorted_desc: list[tuple[float, str]], positives: set[str], k: int) -> tuple[float, float, float, int]:
    if k > len(sorted_desc):
        return (float("nan"), float("nan"), float("nan"), 0)
    cutoff = sorted_desc[k - 1][0]
    tied = [(s, c) for s, c in sorted_desc if s == cutoff]
    tie_band = len(tied)
    top_k = sorted_desc[:k]
    committed_non_tied = [c for s, c in top_k if s != cutoff]
    committed_pos = sum(1 for c in committed_non_tied if c in positives)
    k_band = sum(1 for s, c in top_k if s == cutoff)
    band_pos = sum(1 for _, c in tied if c in positives)
    drop_slots = tie_band - k_band
    observed = sum(1 for _, c in top_k if c in positives)
    p_obs = observed / k
    p_min = (committed_pos + max(0, band_pos - drop_slots)) / k
    p_max = (committed_pos + min(band_pos, k_band)) / k
    return p_obs, p_min, p_max, tie_band


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--scored", required=True, type=Path)
    ap.add_argument("--limma", required=True, type=Path)
    ap.add_argument("--positives", required=True, type=Path)
    ap.add_argument("--output-prefix", required=True, type=Path)
    ap.add_argument("--cohort-label", required=True)
    args = ap.parse_args()

    positives = set(line.strip() for line in args.positives.read_text().splitlines() if line.strip())
    limma_t = load_limma(args.limma)
    print(f"loaded {len(limma_t):,} probe-level limma t-stats from {args.limma.name}")
    print(f"scoring candidates from {args.scored.name} against {args.cohort_label}")

    limma_wg: list[tuple[float, str]] = []
    limma_chr: list[tuple[float, str]] = []
    v25_wg: list[tuple[float, str]] = []
    v25_chr: list[tuple[float, str]] = []
    n_total_wg = 0; n_total_chr = 0
    n_unmatched_probe = 0
    for cid, chrom, in_chr, probe_id, limma_score, v25_score, _is_pos in stream_candidates(args.scored, limma_t, positives):
        n_total_wg += 1
        if in_chr:
            n_total_chr += 1
        if probe_id is not None and probe_id not in limma_t and strip_suffix(probe_id) not in limma_t:
            n_unmatched_probe += 1
        elif probe_id is None:
            n_unmatched_probe += 1
        limma_wg.append((-limma_score, cid))
        v25_wg.append((-v25_score, cid))
        if in_chr:
            limma_chr.append((-limma_score, cid))
            v25_chr.append((-v25_score, cid))
    print(f"  N (WG)         : {n_total_wg:,}")
    print(f"  N (chr5/6/10)  : {n_total_chr:,}")
    print(f"  probes with no limma t (joined as score=0): {n_unmatched_probe:,}")

    for lst in (limma_wg, limma_chr, v25_wg, v25_chr):
        lst.sort()

    rows = []
    for axis, wg_sort, chr_sort in (
        ("limma_t", limma_wg, limma_chr),
        ("v25_diff", v25_wg, v25_chr),
    ):
        for subset, sd, n_total in (("wg", wg_sort, n_total_wg), ("chr5_6_10", chr_sort, n_total_chr)):
            auc, rank_of, _ = midrank_auc_and_positions(sd, positives)
            p100_obs, p100_min, p100_max, _ = p_at_k_interval(sd, positives, 100)
            rows.append({
                "axis": axis, "subset": subset, "n_total": n_total,
                "auc": auc,
                "tb20": tie_band_at_k(sd, 20),
                "tb100": tie_band_at_k(sd, 100),
                "tb500": tie_band_at_k(sd, 500),
                "tb1000": tie_band_at_k(sd, 1000),
                "p100_obs": p100_obs, "p100_min": p100_min, "p100_max": p100_max,
                "ranks": {POS_GENE[c]: rank_of.get(c) for c in positives},
            })

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    tsv = args.output_prefix.with_suffix(".tsv")
    with tsv.open("w") as fh:
        fh.write("axis\tsubset\tn_total\tauc\ttie@20\ttie@100\ttie@500\ttie@1000"
                 "\tP@100_obs\tP@100_min\tP@100_max"
                 "\tESR1_rank\tEGFLAM_rank\tGATA3_rank\n")
        for r in rows:
            fh.write(
                f"{r['axis']}\t{r['subset']}\t{r['n_total']}\t{r['auc']:.6f}\t"
                f"{r['tb20']}\t{r['tb100']}\t{r['tb500']}\t{r['tb1000']}\t"
                f"{r['p100_obs']:.4f}\t{r['p100_min']:.4f}\t{r['p100_max']:.4f}\t"
                f"{r['ranks'].get('ESR1', '')}\t{r['ranks'].get('EGFLAM', '')}\t{r['ranks'].get('GATA3', '')}\n"
            )
    print(f"wrote {tsv}")

    md_lines = [f"# limma-eBayes vs V2.5-diff on {args.cohort_label}",
                "",
                f"Generated by `scripts/limma_candidate_ranking.py`. limma t-stats from",
                f"`{args.limma.name}`, per-candidate probe assignment from the V2.5",
                "`observation.probe_id` field (probe-level limma first, then candidate-",
                "mapped via nearest assigned probe — §7 Methods).",
                "",
                "## AUC + tie-band curve at validated labels (n_pos = 3)",
                "",
                "| axis | subset | AUC | tie@20 | tie@100 | tie@500 | tie@1000 | P@100 [min, max] |",
                "|---|---|---:|---:|---:|---:|---:|---|"]
    for r in rows:
        md_lines.append(
            f"| {r['axis']} | {r['subset']} | {r['auc']:.3f} "
            f"| {r['tb20']:,} | {r['tb100']:,} | {r['tb500']:,} | {r['tb1000']:,} "
            f"| [{r['p100_min']:.3f}, {r['p100_max']:.3f}] |"
        )

    md_lines.extend(["",
                     "## Per-positive ranks (1-based, `candidate_id` ascending tie-break)",
                     "",
                     "| axis | subset | ESR1 | EGFLAM | GATA3 |",
                     "|---|---|---:|---:|---:|"])
    for r in rows:
        pr = r["ranks"]
        md_lines.append(
            f"| {r['axis']} | {r['subset']} | "
            f"{pr.get('ESR1', '—') or '—'} | "
            f"{pr.get('EGFLAM', '—') or '—'} | "
            f"{pr.get('GATA3', '—') or '—'} |"
        )

    md = args.output_prefix.with_suffix(".md")
    md.write_text("\n".join(md_lines))
    print(f"wrote {md}")


if __name__ == "__main__":
    main()
