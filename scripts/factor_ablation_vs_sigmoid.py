#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# ///
"""Factor ablation: does p_diff's σ-aware structure beat a fixed-bandwidth sigmoid?

The V2.5 composite is `p_targ × p_diff × p_trust` where
`p_diff = P(β_n − β_t > δ)` is a per-site normal-approximation
tail probability with a per-site σ_Δ derived from the IQRs
(floored at σ_floor = 0.05). The §3.5 binding-rate diagnostic shows
σ_floor dominates σ_Δ on essentially every observed record on low-n
matched cell-line cohorts (100.0% either-side, 99.5–99.9% both
sides). That raises a sharp reviewer-style question: if σ_floor is
determining σ_Δ almost everywhere on cell lines, is p_diff's
per-site structure actually adding anything over a bare sigmoid
centered at δ with fixed bandwidth?

Ablation axes (all use the same p_targ and p_trust as V2.5, same
δ = 0.2, so the test isolates the **gap factor** in the composite):

  - **full V2.5**                `p_targ × p_diff(σ_Δ) × p_trust`
  - **sigmoid ablation**         `p_targ × sigmoid((Δβ − δ) / σ_fixed) × p_trust`
                                  with σ_fixed = sqrt(2) × σ_floor ≈ 0.0707
                                  (the σ_Δ when σ_floor binds on both sides —
                                  the modal case on cell-line cohorts per §3.5)
  - **hard-threshold ablation**  `p_targ × 1[Δβ > δ] × p_trust`
                                  (a binary step; isolates whether smoothness
                                  matters at all)

Report: AUC at the n = 3 validated positives on every primary cohort
+ the GSE69914 tissue cohort. If V2.5 ≈ sigmoid on matched
cell-line cohorts, report that honestly: the per-site σ adaptation
is not what's doing the work there, and on low-n cohorts V2.5's
incremental value is structural (probability-scale composition +
tie-band reporting), not numerical.

Output: `examples/factor_ablation_vs_sigmoid.tsv` + `.md`.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCORED = REPO / "data" / "derived"
POSITIVES = SCORED / "positives_roth_validated.txt"
OUT_TSV = REPO / "examples" / "factor_ablation_vs_sigmoid.tsv"
OUT_MD = REPO / "examples" / "factor_ablation_vs_sigmoid.md"

_IQR_TO_STDEV = 1.349
DEFAULT_DELTA = 0.2
DEFAULT_SIGMA_FLOOR = 0.05
SIGMA_FIXED = math.sqrt(2) * DEFAULT_SIGMA_FLOOR  # ≈ 0.0707, σ_Δ when floor binds both sides

COHORTS = [
    ("GSE322563 HM450",     "scored_gse322563_differential.jsonl"),
    ("GSE322563 native v2", "scored_gse322563_native_differential.jsonl"),
    ("GSE77348",            "scored_surrogate_differential.jsonl"),
    ("GSE69914",            "scored_gse69914_differential.jsonl"),
]


def p_diff_full(mu_t, mu_n, q25_t, q75_t, q25_n, q75_n, delta=DEFAULT_DELTA, sigma_floor=DEFAULT_SIGMA_FLOOR):
    sigma_t = (q75_t - q25_t) / _IQR_TO_STDEV if (q25_t is not None and q75_t is not None) else 0.0
    sigma_n = (q75_n - q25_n) / _IQR_TO_STDEV if (q25_n is not None and q75_n is not None) else 0.0
    sigma_sq = max(sigma_t, sigma_floor) ** 2 + max(sigma_n, sigma_floor) ** 2
    z = (delta - (mu_n - mu_t)) / math.sqrt(sigma_sq)
    return 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))


def p_diff_sigmoid(mu_t, mu_n, delta=DEFAULT_DELTA, sigma_fixed=SIGMA_FIXED):
    x = (mu_n - mu_t - delta) / sigma_fixed
    # Numerically stable logistic
    if x >= 0:
        e = math.exp(-x)
        return 1.0 / (1.0 + e)
    e = math.exp(x)
    return e / (1.0 + e)


def p_diff_hard(mu_t, mu_n, delta=DEFAULT_DELTA):
    return 1.0 if (mu_n - mu_t) > delta else 0.0


def auc_midrank(scored_path: Path, positives: set[str], score_fn) -> float:
    """score_fn takes a dict record and returns a float (lower = bottom of sort)."""
    rows: list[tuple[float, str]] = []
    with scored_path.open() as fh:
        for line in fh:
            d = json.loads(line)
            cid = d["candidate"]["candidate_id"]
            s = score_fn(d)
            rows.append((-s, cid))  # descending sort by score
    n_total = len(rows)
    rows.sort()
    asc_scores = [-k for k, _ in rows]; asc_scores.reverse()
    asc_cids = [c for _, c in rows]; asc_cids.reverse()

    n_pos = len(positives); n_neg = n_total - n_pos
    midrank = {}
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
    sum_ranks = sum(midrank.values())
    return (sum_ranks - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def make_score_fn(gap_kind: str):
    """Return a function `dict -> float` for the full composite under the given gap factor."""
    def _f(d: dict) -> float:
        obs = d["observation"]; prob = d["probabilistic"]
        mu_t = obs.get("beta_tumor_mean"); mu_n = obs.get("beta_normal_mean")
        p_t = prob["p_targetable_tumor"]; p_r = prob["p_observation_trustworthy"]
        if mu_t is None or mu_n is None:
            return 0.0
        if gap_kind == "full":
            q25_t, q75_t = obs.get("beta_tumor_q25"), obs.get("beta_tumor_q75")
            q25_n, q75_n = obs.get("beta_normal_q25"), obs.get("beta_normal_q75")
            g = p_diff_full(mu_t, mu_n, q25_t, q75_t, q25_n, q75_n)
        elif gap_kind == "sigmoid":
            g = p_diff_sigmoid(mu_t, mu_n)
        elif gap_kind == "hard":
            g = p_diff_hard(mu_t, mu_n)
        else:
            raise ValueError(gap_kind)
        return p_t * g * p_r
    return _f


def main() -> None:
    positives = set(line.strip() for line in POSITIVES.read_text().splitlines() if line.strip())
    print(f"positives: {sorted(positives)}")
    print(f"σ_fixed = sqrt(2) * σ_floor = {SIGMA_FIXED:.4f}")

    rows = []
    table: dict[tuple[str, str], float] = {}
    for cohort, fname in COHORTS:
        path = SCORED / fname
        if not path.exists():
            print(f"SKIP {path}")
            continue
        for gap_kind in ("full", "sigmoid", "hard"):
            print(f"  {cohort} × {gap_kind} ...")
            auc = auc_midrank(path, positives, make_score_fn(gap_kind))
            table[(cohort, gap_kind)] = auc
            rows.append((cohort, gap_kind, auc))
            print(f"    AUC = {auc:.4f}")

    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as fh:
        fh.write("cohort\tgap_factor\tauc\n")
        for cohort, gap, auc in rows:
            fh.write(f"{cohort}\t{gap}\t{auc:.6f}\n")
    print(f"wrote {OUT_TSV}")

    md = ["# Factor ablation — p_diff vs fixed-bandwidth sigmoid vs hard threshold",
          "",
          "Generated by `scripts/factor_ablation_vs_sigmoid.py`.",
          "",
          "All three axes share the same `p_targ` and `p_trust` factors as V2.5",
          "and the same `δ = 0.2`. Only the *gap factor* varies:",
          "",
          "- **full V2.5** — `p_diff(σ_Δ)` with per-site σ_Δ derived from the IQRs (σ_floor = 0.05).",
          f"- **sigmoid** — `sigmoid((Δβ − δ) / σ_fixed)` with σ_fixed = √2 × σ_floor ≈ {SIGMA_FIXED:.4f}.",
          "  This is the σ_Δ V2.5 sees when σ_floor binds on both sides — the modal",
          "  case on cell-line cohorts per §3.5 (99.5–99.9% of records).",
          "- **hard** — `1[Δβ > δ]` (no smoothness). Isolates whether the smoothing matters.",
          "",
          "AUC at the n = 3 Roth-validated positives:",
          "",
          "| cohort | V2.5 (full p_diff) | sigmoid ablation | hard-threshold ablation |",
          "|---|---:|---:|---:|"]
    for cohort, _ in COHORTS:
        cells = []
        for gap in ("full", "sigmoid", "hard"):
            if (cohort, gap) in table:
                cells.append(f"{table[(cohort, gap)]:.3f}")
            else:
                cells.append("—")
        md.append(f"| **{cohort}** | " + " | ".join(cells) + " |")

    md.extend([
        "",
        "## Interpretation",
        "",
        "If (full V2.5) ≈ (sigmoid) on matched cell-line cohorts, the per-site",
        "σ_Δ adaptation in `p_diff` is not doing the work there — the floor is.",
        "The value of the V2.5 composite on low-n matched cell-line data is",
        "structural (bounded probability-scale composition + tie-band-honest",
        "top-K reporting + the p_targ / p_trust wrapping) rather than numerical",
        "improvement from p_diff's shape.",
        "",
        "On tissue (GSE69914), the σ_floor binding rate drops to 65% either-side",
        "and 43% both-sides (§3.5), so there is room for p_diff's per-site σ to",
        "do something the sigmoid cannot. The tissue cell in this table is the",
        "direct test.",
    ])
    OUT_MD.write_text("\n".join(md))
    print(f"wrote {OUT_MD}")


if __name__ == "__main__":
    main()
