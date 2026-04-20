"""Phase 5b — V2 factor ablation + factor-aware top-100 inspection.

Diagnostic sweep to localize WHY V2 `p_therapeutic_selectivity` underperforms
`naive_selectivity` on the matched MCF-7/MCF-10A surrogate cohort, given that
the "wrong normal comparator" rescue hypothesis (Phase 5 main) did not hold.

Near-uniform multiplicative factors (like `p_trust` at n=3) preserve rank, so
rank-metric failure must live in `p_targ`, `p_prot`, or the × composition
itself. This script computes AUC for each isolated factor + pair + full V2 +
a log-odds additive combiner, to decide whether the V2 rescue should start
with composition (V2.2) or trust (V2.1).

Reads:
    data/derived/scored_surrogate.jsonl
    data/derived/positives.txt           (loose, 1687)
    data/derived/positives_tight.txt     (tight promoter, 970)
    data/derived/probes_at_roth_genes.tsv (probe → gene lookup)

Writes:
    data/derived/phase5_v2_ablation.txt
"""

from __future__ import annotations

import json
import math
from bisect import bisect_left, bisect_right
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DERIVED = REPO / "data" / "derived"
REPORT = DERIVED / "phase5_v2_ablation.txt"


def fast_auc(pos_scores: list[float], neg_scores: list[float]) -> float | None:
    """Sort-based Mann-Whitney U normalized. O(n_pos log n_neg). Tie-aware."""

    if not pos_scores or not neg_scores:
        return None
    neg = sorted(neg_scores)
    n_neg = len(neg)
    wins = 0.0
    for p in pos_scores:
        strict = bisect_left(neg, p)               # neg < p
        ties   = bisect_right(neg, p) - strict     # neg == p
        wins += strict + 0.5 * ties
    return wins / (len(pos_scores) * n_neg)


def log_odds(p: float, eps: float = 1e-9) -> float:
    """logit(p) with clamping. NaN-safe on boundaries."""
    p = max(eps, min(1.0 - eps, p))
    return math.log(p / (1.0 - p))


def _out(fh, *args, **kwargs):
    print(*args, **kwargs)
    print(*args, **kwargs, file=fh)


def main() -> int:
    positives_loose = set((DERIVED / "positives.txt").read_text().split())
    positives_tight = set((DERIVED / "positives_tight.txt").read_text().split())

    probe_to_gene: dict[str, str] = {}
    with (DERIVED / "probes_at_roth_genes.tsv").open() as f:
        next(f)
        for ln in f:
            pid, chrom, pos, gene = ln.rstrip().split("\t")
            probe_to_gene[pid] = gene

    # ---------- pass 1: stream scored_surrogate and collect per-axis scores ----------

    # axes we want to score:
    #   v1_final        — framework V1 final_score
    #   v2_full         — p_targ × p_prot × p_trust (V2)
    #   v2_pair         — p_targ × p_prot (ablation: drop trust)
    #   p_targ          — P(targetable_tumor) alone
    #   p_prot          — P(protected_normal) alone
    #   naive           — β_normal - β_tumor (as before)
    #   additive_lo     — log_odds(p_targ) + log_odds(p_prot) (additive composition)
    #   weighted_lo     — same but weighted 1/1 with clamping to protect against extreme negatives
    #
    # For each axis, split into pos_scores/neg_scores for LOOSE and TIGHT positive sets.

    axes = ["v1_final", "v2_full", "v2_pair", "p_targ", "p_prot",
            "naive", "additive_lo"]
    pos_by_axis_loose: dict[str, list[float]] = {a: [] for a in axes}
    neg_by_axis_loose: dict[str, list[float]] = {a: [] for a in axes}
    pos_by_axis_tight: dict[str, list[float]] = {a: [] for a in axes}
    neg_by_axis_tight: dict[str, list[float]] = {a: [] for a in axes}

    # Keep top-100 by V1 and V2 via bounded heaps (score, cid, rec)
    import heapq
    TOPK = 100
    top_v1: list[tuple[float, str, dict]] = []
    top_v2: list[tuple[float, str, dict]] = []

    # For rank lookups, keep full lists sorted after the pass.
    all_v1_scores: list[tuple[float, str]] = []
    all_v2_scores: list[tuple[float, str]] = []

    n_total = 0
    with (DERIVED / "scored_surrogate.jsonl").open() as f:
        for ln in f:
            r = json.loads(ln)
            n_total += 1
            cid  = r["candidate"]["candidate_id"]
            obs  = r["observation"]
            prob = r.get("probabilistic")
            if prob is None:
                # V2-less record — only V1 + naive axes available.
                continue
            bt = obs.get("beta_tumor_mean")
            bn = obs.get("beta_normal_mean")

            v1_final = r["final_score"]
            p_t      = prob["p_targetable_tumor"]
            p_p      = prob["p_protected_normal"]
            p_r      = prob["p_observation_trustworthy"]
            v2_full  = p_t * p_p * p_r
            v2_pair  = p_t * p_p
            naive    = (bn - bt) if (bn is not None and bt is not None) else None
            add_lo   = log_odds(p_t) + log_odds(p_p)

            scores = {
                "v1_final":    v1_final,
                "v2_full":     v2_full,
                "v2_pair":     v2_pair,
                "p_targ":      p_t,
                "p_prot":      p_p,
                "naive":       naive if naive is not None else float("-inf"),
                "additive_lo": add_lo,
            }

            is_loose = cid in positives_loose
            is_tight = cid in positives_tight
            for a, s in scores.items():
                if s == float("-inf"):
                    continue
                (pos_by_axis_loose if is_loose else neg_by_axis_loose)[a].append(s)
                (pos_by_axis_tight if is_tight else neg_by_axis_tight)[a].append(s)

            # bounded heaps for top-100 on V1 and V2_full
            key_v1 = (v1_final, cid, r)
            if len(top_v1) < TOPK:
                heapq.heappush(top_v1, key_v1)
            elif key_v1 > top_v1[0]:
                heapq.heapreplace(top_v1, key_v1)

            key_v2 = (v2_full, cid, r)
            if len(top_v2) < TOPK:
                heapq.heappush(top_v2, key_v2)
            elif key_v2 > top_v2[0]:
                heapq.heapreplace(top_v2, key_v2)

            all_v1_scores.append((v1_final, cid))
            all_v2_scores.append((v2_full, cid))

    # ---------- report ----------

    top_v1.sort(reverse=True)
    top_v2.sort(reverse=True)

    # Rank lookups
    all_v1_scores.sort(reverse=True)
    all_v2_scores.sort(reverse=True)
    v1_rank_of: dict[str, int] = {cid: i + 1 for i, (_, cid) in enumerate(all_v1_scores)}
    v2_rank_of: dict[str, int] = {cid: i + 1 for i, (_, cid) in enumerate(all_v2_scores)}

    with REPORT.open("w") as fh:
        _out(fh, "Phase 5b — V2 factor ablation + factor-aware top-100 inspection")
        _out(fh, "=" * 80)
        _out(fh, f"scored records with V2 probabilistic: {n_total}")
        _out(fh, f"loose positives: {len(positives_loose)}  tight positives: {len(positives_tight)}")
        _out(fh, "")

        # ---------- ablation AUC table ----------

        _out(fh, "AUC table (higher = better; 0.5 = random; <0.5 = inverted)")
        _out(fh, "-" * 80)
        _out(fh, f"{'axis':<14}  {'loose AUC':>10}  {'tight AUC':>10}  "
             f"{'n_pos_loose':>11}  {'n_pos_tight':>11}")
        for a in axes:
            loose = fast_auc(pos_by_axis_loose[a], neg_by_axis_loose[a])
            tight = fast_auc(pos_by_axis_tight[a], neg_by_axis_tight[a])
            def f(x): return f"{x:.3f}" if x is not None else "   --"
            _out(fh, f"{a:<14}  {f(loose):>10}  {f(tight):>10}  "
                 f"{len(pos_by_axis_loose[a]):>11}  {len(pos_by_axis_tight[a]):>11}")

        _out(fh, "")
        _out(fh, "Diagnostic reading:")
        _out(fh, "  * If p_targ*p_prot ≈ v2_full → trust factor isn't the rank-metric problem.")
        _out(fh, "  * If additive_lo clearly beats the products → × composition is the problem.")
        _out(fh, "  * If p_targ alone ≈ p_targ*p_prot → p_prot adds nothing useful (bulk biology).")

        # ---------- top-100 by V1 ----------

        def _print_top(label: str, items: list[tuple[float, str, dict]]):
            _out(fh, "")
            _out(fh, "=" * 80)
            _out(fh, f"Top-{TOPK} by {label} on surrogate (showing top 20 + aggregates)")
            _out(fh, "=" * 80)
            header = (f"{'rank':>4}  {'candidate_id':<28}  {'β_t':>5}  {'β_n':>5}  "
                      f"{'Δ':>6}  {'p_targ':>6}  {'p_prot':>6}  {'v1r':>5}  {'v2r':>7}  "
                      f"{'probe':<11}  {'gene':<8}  pos?")
            _out(fh, header)
            _out(fh, "-" * len(header))
            for i, (_score, cid, r) in enumerate(items[:20], 1):
                obs  = r["observation"]
                prob = r.get("probabilistic", {})
                bt = obs.get("beta_tumor_mean");  bn = obs.get("beta_normal_mean")
                delta = (bn - bt) if (bt is not None and bn is not None) else None
                pid = obs.get("probe_id") or "-"
                gene = probe_to_gene.get(pid, "")
                mark = "LOOSE" if cid in positives_loose else "-"
                if cid in positives_tight: mark += "+TIGHT"
                _out(fh, f"{i:>4}  {cid:<28}  "
                     f"{(f'{bt:.2f}' if bt is not None else '  --'):>5}  "
                     f"{(f'{bn:.2f}' if bn is not None else '  --'):>5}  "
                     f"{(f'{delta:+.2f}' if delta is not None else '   --'):>6}  "
                     f"{prob.get('p_targetable_tumor', 0):>6.3f}  "
                     f"{prob.get('p_protected_normal', 0):>6.3f}  "
                     f"{v1_rank_of.get(cid, 0):>5}  {v2_rank_of.get(cid, 0):>7}  "
                     f"{pid:<11}  {gene:<8}  {mark}")

            # aggregates over full top-100
            in_loose = sum(1 for (_, cid, _) in items if cid in positives_loose)
            in_tight = sum(1 for (_, cid, _) in items if cid in positives_tight)
            genes = {}
            for (_, _, r) in items:
                g = probe_to_gene.get(r["observation"].get("probe_id") or "-", "")
                genes[g] = genes.get(g, 0) + 1
            _out(fh, "")
            _out(fh, f"Top-{TOPK} aggregates:")
            _out(fh, f"  positives (loose): {in_loose}/{TOPK}  (tight): {in_tight}/{TOPK}")
            _out(fh, f"  gene breakdown:    {dict(sorted(genes.items(), key=lambda t: -t[1]))}")

        _print_top("V1 final_score", top_v1)
        _print_top("V2 p_therapeutic_selectivity (full)", top_v2)

    print(f"report → {REPORT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
