"""V2.5 OFFLINE experiment — differential-based p_protected_normal.

Replaces the threshold-based `p_prot = P(β_normal > 0.5)` with a differential-
based one:

    p_diff(δ) = P(β_normal − β_tumor > δ)

computed from tumor / normal summary statistics via a normal approximation:

    μ_Δ = μ_n − μ_t
    σ_Δ ≈ sqrt(σ_t² + σ_n²),   σ_k ≈ IQR_k / 1.349
    P(Δ > δ) = 1 − Φ((δ − μ_Δ) / σ_Δ)

This avoids the static-threshold assumption that broke `p_prot` on bulk
comparators (AUC 0.343 inverted) while keeping the factor cohort-agnostic.

This is an OFFLINE prototype — does NOT modify the installed package. If any
δ variant cleanly beats V1 final_score on BOTH AUC and P@100, wire it in as
a new `probabilistic_mode = "tumor_plus_differential_protection"` in V2.5.

Reads:
    data/derived/scored_surrogate_tumor_only.jsonl  (3M records)
    data/derived/positives.txt                      (loose, 1687)
    data/derived/positives_tight.txt                (tight promoter, 970)
    data/derived/probes_at_roth_genes.tsv           (probe → gene)

Writes:
    data/derived/v2_5_differential_offline.txt
"""

from __future__ import annotations

import heapq
import json
import math
from bisect import bisect_left, bisect_right
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DERIVED = REPO / "data" / "derived"
REPORT = DERIVED / "v2_5_differential_offline.txt"

SRC = DERIVED / "scored_surrogate_tumor_only.jsonl"

DELTAS = [0.2, 0.3, 0.4, 0.5]
TOPK = 100

SIGMA_FLOOR = 0.05  # avoid collapsing σ_Δ when IQR=0 at boundary betas


def phi(z: float) -> float:
    """Standard normal CDF via erf. Stdlib only."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def p_differential(
    mu_t: float, sigma_t: float,
    mu_n: float, sigma_n: float,
    delta: float,
) -> float:
    """P(β_n − β_t > δ) under independent-normal approximation."""
    sigma_sq = max(sigma_t, SIGMA_FLOOR) ** 2 + max(sigma_n, SIGMA_FLOOR) ** 2
    z = (delta - (mu_n - mu_t)) / math.sqrt(sigma_sq)
    return 1.0 - phi(z)


def fast_auc(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """Sort-based Mann-Whitney U normalized. Tie-aware."""
    if not pos_scores or not neg_scores:
        return None
    neg = sorted(neg_scores)
    n_neg = len(neg)
    wins = 0.0
    for p in pos_scores:
        strict = bisect_left(neg, p)
        ties = bisect_right(neg, p) - strict
        wins += strict + 0.5 * ties
    return wins / (len(pos_scores) * n_neg)


def _out(fh, *a, **kw):
    print(*a, **kw)
    print(*a, **kw, file=fh)


def main() -> int:
    positives_loose = set((DERIVED / "positives.txt").read_text().split())
    positives_tight = set((DERIVED / "positives_tight.txt").read_text().split())

    probe_to_gene: dict[str, str] = {}
    with (DERIVED / "probes_at_roth_genes.tsv").open() as f:
        next(f)
        for ln in f:
            pid, _chrom, _pos, gene = ln.rstrip().split("\t")
            probe_to_gene[pid] = gene

    # Axes:
    #   diff_{δ}          — p_diff(δ) alone
    #   diff_{δ}_x_trust  — p_diff(δ) × p_observation_trustworthy
    #   v1_final          — V1 final_score (reference)
    #   tumor_only        — p_targ × p_trust (current default V2 mode)
    axes: list[str] = []
    for d in DELTAS:
        axes.append(f"diff_{d:.1f}")
        axes.append(f"diff_{d:.1f}_x_trust")
    axes += ["v1_final", "tumor_only"]

    pos_loose: dict[str, list[float]] = {a: [] for a in axes}
    neg_loose: dict[str, list[float]] = {a: [] for a in axes}
    pos_tight: dict[str, list[float]] = {a: [] for a in axes}
    neg_tight: dict[str, list[float]] = {a: [] for a in axes}

    # bounded top-K heaps per axis
    tops: dict[str, list[tuple[float, str, dict]]] = {a: [] for a in axes}

    n_total = 0
    n_computable = 0
    n_skipped = 0
    with SRC.open() as f:
        for ln in f:
            n_total += 1
            r = json.loads(ln)
            cid = r["candidate"]["candidate_id"]
            obs = r["observation"]
            prob = r.get("probabilistic")
            if prob is None:
                n_skipped += 1
                continue

            mu_t = obs.get("beta_tumor_mean")
            mu_n = obs.get("beta_normal_mean")
            if mu_t is None or mu_n is None:
                n_skipped += 1
                continue

            q25_t = obs.get("beta_tumor_q25")
            q75_t = obs.get("beta_tumor_q75")
            q25_n = obs.get("beta_normal_q25")
            q75_n = obs.get("beta_normal_q75")

            sigma_t = (q75_t - q25_t) / 1.349 if (q25_t is not None and q75_t is not None) else SIGMA_FLOOR
            sigma_n = (q75_n - q25_n) / 1.349 if (q25_n is not None and q75_n is not None) else SIGMA_FLOOR

            p_trust = prob["p_observation_trustworthy"]
            p_targ = prob["p_targetable_tumor"]
            v1_final = r["final_score"]
            tumor_only_score = p_targ * p_trust

            scores: dict[str, float] = {
                "v1_final": v1_final,
                "tumor_only": tumor_only_score,
            }
            for d in DELTAS:
                pd = p_differential(mu_t, sigma_t, mu_n, sigma_n, d)
                scores[f"diff_{d:.1f}"] = pd
                scores[f"diff_{d:.1f}_x_trust"] = pd * p_trust

            n_computable += 1

            is_loose = cid in positives_loose
            is_tight = cid in positives_tight
            for a, s in scores.items():
                (pos_loose if is_loose else neg_loose)[a].append(s)
                (pos_tight if is_tight else neg_tight)[a].append(s)

                heap = tops[a]
                key = (s, cid, r)
                if len(heap) < TOPK:
                    heapq.heappush(heap, key)
                elif key > heap[0]:
                    heapq.heapreplace(heap, key)

    for a in axes:
        tops[a].sort(reverse=True)

    with REPORT.open("w") as fh:
        _out(fh, "V2.5 OFFLINE — differential-based p_protected_normal experiment")
        _out(fh, "=" * 80)
        _out(fh, f"input:       {SRC.name}")
        _out(fh, f"records:     {n_total} total, {n_computable} computable, {n_skipped} skipped")
        _out(fh, f"positives:   loose={len(positives_loose)}  tight={len(positives_tight)}")
        _out(fh, f"deltas:      {DELTAS}")
        _out(fh, f"sigma model: IQR/1.349, floor={SIGMA_FLOOR}")
        _out(fh, "")

        _out(fh, "AUC + P@100 table")
        _out(fh, "-" * 80)
        header = (f"{'axis':<22}  {'AUC loose':>10}  {'AUC tight':>10}  "
                  f"{'P@100 loose':>12}  {'P@100 tight':>12}")
        _out(fh, header)
        _out(fh, "-" * len(header))

        def _fmt(x: float | None) -> str:
            return f"{x:.3f}" if x is not None else "   --"

        for a in axes:
            auc_l = fast_auc(pos_loose[a], neg_loose[a])
            auc_t = fast_auc(pos_tight[a], neg_tight[a])
            top = tops[a]
            p100_l = sum(1 for (_, cid, _) in top if cid in positives_loose) / TOPK
            p100_t = sum(1 for (_, cid, _) in top if cid in positives_tight) / TOPK
            _out(fh, f"{a:<22}  {_fmt(auc_l):>10}  {_fmt(auc_t):>10}  "
                 f"{p100_l:>12.3f}  {p100_t:>12.3f}")

        _out(fh, "")
        _out(fh, "Reading guide:")
        _out(fh, "  * V1 final_score is the current-recommended discovery axis.")
        _out(fh, "  * tumor_only = current V2 default (p_targ × p_trust) — high AUC, P@100≈0.")
        _out(fh, "  * diff_δ restores a differential signal without the β_n>0.5 assumption.")
        _out(fh, "  * A candidate 'wins' only if it beats V1 on BOTH AUC and P@100.")
        _out(fh, "")

        # Top-10 per axis for biological sanity.
        def _print_top10(axis: str):
            _out(fh, "=" * 80)
            _out(fh, f"Top-10 by {axis}")
            _out(fh, "-" * 80)
            hdr = (f"{'rank':>4}  {'candidate_id':<28}  {'β_t':>5}  {'β_n':>5}  "
                   f"{'Δ':>6}  {'p_targ':>6}  {'p_trust':>7}  {'gene':<8}  pos?")
            _out(fh, hdr)
            for i, (_s, cid, r) in enumerate(tops[axis][:10], 1):
                obs = r["observation"]
                prob = r.get("probabilistic", {})
                bt = obs.get("beta_tumor_mean")
                bn = obs.get("beta_normal_mean")
                delta_v = (bn - bt) if (bt is not None and bn is not None) else None
                pid = obs.get("probe_id") or "-"
                gene = probe_to_gene.get(pid, "")
                mark = "LOOSE" if cid in positives_loose else "-"
                if cid in positives_tight:
                    mark += "+TIGHT"
                _out(fh, f"{i:>4}  {cid:<28}  "
                     f"{(f'{bt:.2f}' if bt is not None else '  --'):>5}  "
                     f"{(f'{bn:.2f}' if bn is not None else '  --'):>5}  "
                     f"{(f'{delta_v:+.2f}' if delta_v is not None else '   --'):>6}  "
                     f"{prob.get('p_targetable_tumor', 0):>6.3f}  "
                     f"{prob.get('p_observation_trustworthy', 0):>7.3f}  "
                     f"{gene:<8}  {mark}")
            _out(fh, "")

        _print_top10("v1_final")
        _print_top10("tumor_only")
        for d in DELTAS:
            _print_top10(f"diff_{d:.1f}_x_trust")

    print(f"report → {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
