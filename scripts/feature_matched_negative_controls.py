#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2.0"]
# ///
"""Feature-matched negative-universe controls for the three Roth
validated positives.

Reviewer ask (memo-2026-04-22-az / -ba follow-up): the WG AUC numbers
report the positives' rank against the *full* candidate universe.
That is a valid rank-lift summary but does not control for the
candidate features the scorer is allowed to use as discrimination
proxies (EvidenceClass / probe distance, PAM family, CpG context).
This script reports the positives' rank against a *feature-matched
negative universe* — candidates that share the positive's
distance/CpG/PAM profile — to answer "are the positives high versus
comparable candidates?" rather than "are they high versus everything?"

For each (cohort, scope, positive, axis) cell we report:

  * positive's score under the axis
  * size of the matched-negative pool (`n_matched_neg`)
  * positive's rank within the matched-null distribution (1 = best)
  * empirical p-value `(rank − 1) / n_matched_neg` — equivalently the
    fraction of matched negatives that score *strictly better* than the
    positive (lower is better; matches the convention of an
    against-random one-sided null)
  * positive's percentile within the matched-null

Matching keys (categorical, exact):

  - EvidenceClass (exact, proximal_close, proximal, regional, unobserved)
  - pam_family (NNNNCGA vs NNNNCCA)
  - is_cpg_pam (True / False — perfectly correlated with NNNNCGA but
    kept explicit so the matching rule is auditable)
  - chrom (same chromosome — eliminates chromosome-class confounders
    like centromere coverage, GC skew, etc.)

This is intentionally a strict matching set; loosening to
"any-chromosome-matched" or "same-chromosome-class" is documented
as future work in PAPER §6.3.

Outputs:
  examples/feature_matched_negative_controls.{tsv,md}

Usage:
  uv run python scripts/feature_matched_negative_controls.py
"""

from __future__ import annotations

import bisect
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"
DERIVED = REPO / "data" / "derived"

DEFAULT_DELTA = 0.2
DEFAULT_SIGMA_FIXED = math.sqrt(2) * 0.05  # ≈ 0.0707, matches src/thermocas/probabilistic.py

POS_GENE = {
    "chr5:38258943+:NNNNCGA": "EGFLAM",
    "chr6:152011177+:NNNNCGA": "ESR1",
    "chr10:8087387+:NNNNCGA": "GATA3",
}
GENE_ORDER = ("ESR1", "EGFLAM", "GATA3")
SCOPES = ("wg", "chr5_6_10")
CHR_SUBSET = {"chr5", "chr6", "chr10"}


def sigmoid(x: float) -> float:
    if x >= 0:
        e = math.exp(-x)
        return 1.0 / (1.0 + e)
    e = math.exp(x)
    return e / (1.0 + e)


@dataclass(frozen=True)
class Signature:
    evidence_class: str
    pam_family: str
    is_cpg_pam: bool
    chrom: str


def _signature(rec: dict) -> Signature:
    cand = rec["candidate"]
    obs = rec["observation"]
    return Signature(
        evidence_class=obs.get("evidence_class", "unobserved"),
        pam_family=cand.get("pam_family", ""),
        is_cpg_pam=bool(cand.get("is_cpg_pam", False)),
        chrom=cand.get("chrom", ""),
    )


def _scores(rec: dict, limma_t: dict[str, float]) -> dict[str, float]:
    """Compute every axis score from already-stored fields. Mirrors
    `scripts/evidence_class_stratified_benchmark.py`'s sign convention
    (limma is -t_mod, V2.5-sigmoid is recomputed in-script)."""
    obs = rec["observation"]
    prob = rec.get("probabilistic") or {}
    p_targ = prob.get("p_targetable_tumor") or 0.0
    p_trust = prob.get("p_observation_trustworthy") or 0.0
    p_diff = prob.get("p_therapeutic_selectivity") or 0.0  # V2.5-diff composite
    mu_t = obs.get("beta_tumor_mean")
    mu_n = obs.get("beta_normal_mean")
    if mu_t is None or mu_n is None:
        delta_beta = float("-inf")
        v25_sigmoid = 0.0
    else:
        delta_beta = mu_n - mu_t
        v25_sigmoid = p_targ * sigmoid((delta_beta - DEFAULT_DELTA) / DEFAULT_SIGMA_FIXED) * p_trust
    probe_id = obs.get("probe_id")
    if probe_id is None:
        limma_score = 0.0
    else:
        t = limma_t.get(probe_id)
        if t is None and "_" in probe_id:
            t = limma_t.get(probe_id.split("_", 1)[0])
        limma_score = -t if t is not None else 0.0
    return {
        "v25_sigmoid": v25_sigmoid,
        "v25_diff": p_diff,
        "delta_beta_only": delta_beta,
        "limma_style": limma_score,
    }


def _load_limma(path: Path) -> dict[str, float]:
    out: dict[str, float] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        pid_i = header.index("probe_id")
        t_i = header.index("t_mod")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            try:
                t = float(parts[t_i])
            except (ValueError, IndexError):
                continue
            if not math.isfinite(t):
                continue
            pid = parts[pid_i]
            base = pid.split("_", 1)[0] if "_" in pid else pid
            old = out.get(base)
            if old is None or abs(t) > abs(old):
                out[base] = t
    return out


def _scopes_for(chrom: str) -> tuple[str, ...]:
    return ("wg", "chr5_6_10") if chrom in CHR_SUBSET else ("wg",)


def _empirical_p_and_percentile(pos_score: float, neg_scores: list[float]) -> tuple[float, float, int]:
    """Return (empirical_p, percentile, rank).

    rank = 1 + (number of negatives with score > pos_score)
    empirical_p = (rank - 1) / n_matched_neg = fraction strictly better
    percentile = 100 * (n_matched_neg - rank + 1) / n_matched_neg, in
    [0, 100]; higher = better. ties contribute 0.5 to the strict-better
    count so the empirical p remains a defensible one-sided null.
    """
    if not neg_scores:
        return float("nan"), float("nan"), 0
    n = len(neg_scores)
    arr = sorted(neg_scores)
    # number of negatives strictly greater than the positive
    strict_better = n - bisect.bisect_right(arr, pos_score)
    # number tied exactly equal
    tied = bisect.bisect_right(arr, pos_score) - bisect.bisect_left(arr, pos_score)
    # 0.5 weight on ties so the empirical p is unbiased under the null.
    rank_float = strict_better + 0.5 * tied + 1.0
    rank = int(round(rank_float))
    p = (rank_float - 1.0) / n
    pct = 100.0 * (n - rank_float + 1.0) / n
    return p, pct, rank


def run_cohort(label: str, scored_path: Path, limma_path: Path, positives: set[str]):
    print(f"== {label}")
    print(f"  scored: {scored_path.relative_to(REPO)}")
    print(f"  limma : {limma_path.relative_to(REPO)}")
    limma_t = _load_limma(limma_path)
    print(f"  limma probes: {len(limma_t):,}")

    # Pass 1: find positives + their signatures + their scores per axis.
    pos_records: dict[str, tuple[Signature, dict[str, float]]] = {}
    with scored_path.open() as fh:
        for line in fh:
            rec = json.loads(line)
            cid = rec["candidate"]["candidate_id"]
            if cid in positives:
                pos_records[cid] = (_signature(rec), _scores(rec, limma_t))
    print(f"  positives found: {len(pos_records)} of {len(positives)}")

    # Pass 2: stream once, accumulate matched-neg score lists per
    # (positive, scope, axis). Memory cost is O(n_matched_neg × n_pos × n_scope × n_axis);
    # matched-neg pools are ~10⁴–10⁶ records, comfortably in memory.
    AXES = ("v25_sigmoid", "v25_diff", "delta_beta_only", "limma_style")
    matched: dict[tuple[str, str, str], list[float]] = defaultdict(list)
    with scored_path.open() as fh:
        for line in fh:
            rec = json.loads(line)
            cid = rec["candidate"]["candidate_id"]
            if cid in positives:
                continue
            sig = _signature(rec)
            scores: dict[str, float] | None = None
            chrom = rec["candidate"]["chrom"]
            for pos_id, (pos_sig, _) in pos_records.items():
                if sig != pos_sig:
                    continue
                if scores is None:
                    scores = _scores(rec, limma_t)
                for scope in _scopes_for(chrom):
                    if scope == "chr5_6_10" and chrom not in CHR_SUBSET:
                        continue
                    for axis in AXES:
                        matched[(pos_id, scope, axis)].append(scores[axis])

    # Build rows.
    rows = []
    for pos_id, (pos_sig, pos_scores) in pos_records.items():
        gene = POS_GENE[pos_id]
        for scope in SCOPES:
            for axis in AXES:
                neg_scores = matched.get((pos_id, scope, axis), [])
                pos_score = pos_scores[axis]
                p, pct, rank = _empirical_p_and_percentile(pos_score, neg_scores)
                rows.append(
                    {
                        "cohort": label,
                        "gene": gene,
                        "candidate_id": pos_id,
                        "evidence_class": pos_sig.evidence_class,
                        "pam_family": pos_sig.pam_family,
                        "is_cpg_pam": pos_sig.is_cpg_pam,
                        "chrom": pos_sig.chrom,
                        "scope": scope,
                        "axis": axis,
                        "pos_score": pos_score,
                        "n_matched_neg": len(neg_scores),
                        "rank": rank,
                        "empirical_p": p,
                        "percentile": pct,
                    }
                )
    return rows


def _format_md(rows: list[dict]) -> str:
    md = [
        "# Feature-matched negative-universe controls",
        "",
        "Generated by `scripts/feature_matched_negative_controls.py`.",
        "Each of the three Roth Fig. 5d validated positives is",
        "ranked against a *feature-matched* negative universe — the",
        "subset of WG candidates that share the positive's",
        "(`EvidenceClass`, `pam_family`, `is_cpg_pam`, `chrom`)",
        "signature exactly. Empirical *p*-value is the fraction of",
        "matched negatives that score strictly better than the",
        "positive (ties contribute 0.5).",
        "",
        "**This is a denominator-confounding audit**, not validation.",
        "It answers \"are the positives high versus comparable",
        "candidates?\" — not \"does the method discover new targets?\"",
        "Loosening the matching to any-chromosome-matched or",
        "same-chromosome-class is documented as future work",
        "(PAPER.md §6.3).",
        "",
        "## Per-positive matched signatures",
        "",
        "| gene | candidate_id | EvidenceClass | pam_family | is_cpg_pam | chrom |",
        "|---|---|---|---|---|---|",
    ]
    seen = set()
    for r in rows:
        key = r["candidate_id"]
        if key in seen:
            continue
        seen.add(key)
        md.append(
            f"| *{r['gene']}* | `{r['candidate_id']}` | "
            f"`{r['evidence_class']}` | `{r['pam_family']}` | "
            f"{r['is_cpg_pam']} | {r['chrom']} |"
        )

    for scope in SCOPES:
        md += [
            "",
            f"## Empirical *p*-value within matched negatives — scope `{scope}`",
            "",
            "| cohort | gene | axis | n_matched_neg | rank | empirical *p* | percentile |",
            "|---|---|---|---:|---:|---:|---:|",
        ]
        for r in rows:
            if r["scope"] != scope:
                continue
            if r["n_matched_neg"] == 0:
                p_str = "—"
                pct_str = "—"
                rank_str = "—"
            else:
                p_str = (
                    f"{r['empirical_p']:.4f}"
                    if r["empirical_p"] >= 1e-4
                    else f"{r['empirical_p']:.2e}"
                )
                pct_str = f"{r['percentile']:.2f}"
                rank_str = f"{r['rank']:,}"
            md.append(
                f"| {r['cohort']} | *{r['gene']}* | `{r['axis']}` | "
                f"{r['n_matched_neg']:,} | {rank_str} | {p_str} | {pct_str} |"
            )

    md += [
        "",
        "## Reproduce",
        "",
        "```bash",
        "uv run python scripts/feature_matched_negative_controls.py",
        "```",
        "",
        "Reads the gitignored multi-GB WG scored JSONLs listed in the",
        "PAPER.md reproducibility appendix and the per-cohort",
        "`limma_<cohort>_probes.tsv` per-probe moderated-t outputs.",
        "",
    ]
    return "\n".join(md)


def main() -> int:
    cohorts = [
        (
            "GSE322563 HM450",
            DERIVED / "scored_gse322563_wg_differential.jsonl",
            DERIVED / "limma_gse322563_probes.tsv",
        ),
        (
            "GSE322563 native v2",
            DERIVED / "scored_gse322563_native_wg_differential.jsonl",
            DERIVED / "limma_gse322563_probes.tsv",
        ),
        (
            "GSE77348",
            DERIVED / "scored_surrogate_wg_differential.jsonl",
            DERIVED / "limma_gse77348_probes.tsv",
        ),
        (
            "GSE69914 (tissue)",
            DERIVED / "scored_gse69914_wg_differential.jsonl",
            DERIVED / "limma_gse69914_probes.tsv",
        ),
    ]
    positives = set(POS_GENE.keys())

    rows: list[dict] = []
    for label, scored, limma in cohorts:
        if not scored.exists():
            print(f"SKIP {label}: missing {scored}")
            continue
        if not limma.exists():
            print(f"SKIP {label}: missing {limma}")
            continue
        rows.extend(run_cohort(label, scored, limma, positives))

    out_md = EXAMPLES / "feature_matched_negative_controls.md"
    out_tsv = EXAMPLES / "feature_matched_negative_controls.tsv"
    out_md.write_text(_format_md(rows))
    cols = [
        "cohort",
        "gene",
        "candidate_id",
        "evidence_class",
        "pam_family",
        "is_cpg_pam",
        "chrom",
        "scope",
        "axis",
        "pos_score",
        "n_matched_neg",
        "rank",
        "empirical_p",
        "percentile",
    ]
    with out_tsv.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            row_vals = []
            for c in cols:
                v = r[c]
                if isinstance(v, float):
                    row_vals.append(f"{v:.6g}" if math.isfinite(v) else "")
                elif isinstance(v, bool):
                    row_vals.append("True" if v else "False")
                else:
                    row_vals.append(str(v))
            fh.write("\t".join(row_vals) + "\n")
    print(f"\nwrote {out_md}\nwrote {out_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
