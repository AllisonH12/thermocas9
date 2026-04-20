"""V2.4 top-100 inspection: V1 final_score vs tumor_only p_therapeutic_selectivity.

Per the earlier guidance: the more useful comparison is v1_final vs the
tumor-only probabilistic score, NOT vs the old (inverted) v2_full. Asks:
now that tumor_only mode gives AUC 0.73-0.77, is the *top* of its ranked
list biologically different from v1_final's top?
"""

from __future__ import annotations

import heapq
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DERIVED = REPO / "data" / "derived"


def main() -> int:
    positives_loose = set((DERIVED / "positives.txt").read_text().split())
    positives_tight = set((DERIVED / "positives_tight.txt").read_text().split())
    probe_to_gene: dict[str, str] = {}
    with (DERIVED / "probes_at_roth_genes.tsv").open() as f:
        next(f)
        for ln in f:
            pid, chrom, pos, gene = ln.rstrip().split("\t")
            probe_to_gene[pid] = gene

    K = 100
    top_v1: list[tuple[float, str, dict]] = []
    top_v24: list[tuple[float, str, dict]] = []

    # One pass through scored_surrogate_tumor_only.jsonl
    src = DERIVED / "scored_surrogate_tumor_only.jsonl"
    with src.open() as f:
        for ln in f:
            r = json.loads(ln)
            cid = r["candidate"]["candidate_id"]
            v1  = r["final_score"]
            prob = r.get("probabilistic")
            if prob is None:
                continue
            v24 = prob["p_therapeutic_selectivity"]

            for heap, score in ((top_v1, v1), (top_v24, v24)):
                key = (score, cid, r)
                if len(heap) < K:
                    heapq.heappush(heap, key)
                elif key > heap[0]:
                    heapq.heapreplace(heap, key)

    top_v1.sort(reverse=True)
    top_v24.sort(reverse=True)

    def aggregate(top: list[tuple[float, str, dict]], label: str) -> None:
        in_loose = sum(1 for (_, cid, _) in top if cid in positives_loose)
        in_tight = sum(1 for (_, cid, _) in top if cid in positives_tight)
        genes: dict[str, int] = {}
        chroms: dict[str, int] = {}
        for (_, _, r) in top:
            g = probe_to_gene.get(r["observation"].get("probe_id") or "-", "")
            genes[g] = genes.get(g, 0) + 1
            chroms[r["candidate"]["chrom"]] = chroms.get(r["candidate"]["chrom"], 0) + 1
        print(f"\n  {label} top-{K}:")
        print(f"    loose pos:   {in_loose}/{K}")
        print(f"    tight pos:   {in_tight}/{K}")
        print(f"    by chrom:    {dict(sorted(chroms.items()))}")
        print(f"    by gene:     {dict(sorted(genes.items(), key=lambda t: -t[1]))}")

    # identical-set analysis
    v1_set  = {cid for (_, cid, _) in top_v1}
    v24_set = {cid for (_, cid, _) in top_v24}
    overlap = v1_set & v24_set
    print(f"Top-{K} overlap between V1 and V2.4 tumor_only: {len(overlap)}")

    aggregate(top_v1,  "V1 final_score")
    aggregate(top_v24, "V2.4 p_sel (tumor_only)")

    # Show a few from each top-20 side-by-side.
    print("\nSide-by-side top-10 (rank, id, β_t, β_n, p_targ, p_trust, pos?):")
    print(f"\n  V1 final_score:")
    print(f"  {'rank':>4}  {'candidate_id':<28}  {'β_t':>5}  {'β_n':>5}  "
          f"{'p_targ':>6}  {'p_trust':>7}  pos?")
    for i, (s, cid, r) in enumerate(top_v1[:10], 1):
        obs = r["observation"]; prob = r.get("probabilistic", {})
        bt = obs.get("beta_tumor_mean"); bn = obs.get("beta_normal_mean")
        mark = "LOOSE" if cid in positives_loose else "-"
        if cid in positives_tight: mark += "+TIGHT"
        print(f"  {i:>4}  {cid:<28}  "
              f"{(f'{bt:.2f}' if bt is not None else '  --'):>5}  "
              f"{(f'{bn:.2f}' if bn is not None else '  --'):>5}  "
              f"{prob.get('p_targetable_tumor', 0):>6.3f}  "
              f"{prob.get('p_observation_trustworthy', 0):>7.3f}  {mark}")

    print(f"\n  V2.4 tumor_only p_sel:")
    print(f"  {'rank':>4}  {'candidate_id':<28}  {'β_t':>5}  {'β_n':>5}  "
          f"{'p_targ':>6}  {'p_trust':>7}  pos?")
    for i, (s, cid, r) in enumerate(top_v24[:10], 1):
        obs = r["observation"]; prob = r.get("probabilistic", {})
        bt = obs.get("beta_tumor_mean"); bn = obs.get("beta_normal_mean")
        mark = "LOOSE" if cid in positives_loose else "-"
        if cid in positives_tight: mark += "+TIGHT"
        print(f"  {i:>4}  {cid:<28}  "
              f"{(f'{bt:.2f}' if bt is not None else '  --'):>5}  "
              f"{(f'{bn:.2f}' if bn is not None else '  --'):>5}  "
              f"{prob.get('p_targetable_tumor', 0):>6.3f}  "
              f"{prob.get('p_observation_trustworthy', 0):>7.3f}  {mark}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
