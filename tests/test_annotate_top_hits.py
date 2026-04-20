"""Regression tests for scripts/annotate_top_hits.py.

Covers the four bugs flagged in the 2026-04-20 reporting-surface review:

- `top_k_by` tie-break must match `evaluate_ranking` (candidate_id ascending).
- `annotate_gene` must treat UCSC transcript ends as half-open (pos == tx_end
  is outside the transcript, not gene_body).
- `annotate_cgi` must treat CGI ends as half-open (pos == cgi_end is outside
  the island, not inside).
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

# Load the script as a module so we can import its helpers cleanly.
_SCRIPT_PATH = Path(__file__).resolve().parent.parent / "scripts" / "annotate_top_hits.py"
_spec = importlib.util.spec_from_file_location("annotate_top_hits", _SCRIPT_PATH)
assert _spec and _spec.loader
annotate = importlib.util.module_from_spec(_spec)
sys.modules["annotate_top_hits"] = annotate
_spec.loader.exec_module(annotate)


# ---------- top_k_by tie-break ----------


def test_top_k_by_tie_break_is_candidate_id_ascending(tmp_path: Path):
    """Three records at tied score; K=2 must return candidate_ids 'a', 'b'
    (ascending), matching `evaluate_ranking`. Earlier implementation returned
    'c', 'b' via the bounded-heap's default tuple ordering."""

    def rec(cid: str, score: float) -> dict:
        return {
            "candidate": {
                "candidate_id": cid, "chrom": "chr1", "critical_c_pos": 10,
                "strand": "+", "pam": "ACGTCGA", "pam_family": "NNNNCGA",
                "is_cpg_pam": True, "local_seq_100bp": "",
            },
            "observation": {"candidate_id": cid, "cohort_name": "T",
                            "evidence_class": "exact"},
            "final_score": score,
            "probabilistic": None,
        }
    jsonl = tmp_path / "tied.jsonl"
    jsonl.write_text("\n".join(json.dumps(rec(c, 1.0)) for c in ("a", "b", "c")) + "\n")

    out = annotate.top_k_by("final_score", jsonl, 2)
    got = [r["candidate"]["candidate_id"] for r in out]
    assert got == ["a", "b"], (
        f"top_k_by returned {got!r}; expected ['a','b'] under candidate_id "
        f"ascending tie-break"
    )


# ---------- annotate_gene half-open ----------


def test_annotate_gene_transcript_end_is_exclusive():
    """Regression: pos == tx_end should NOT be inside the transcript under
    UCSC half-open semantics. Earlier implementation marked it as gene_body."""
    tx_list = [(1000, 2100, "+", "GENE")]
    # pos == tx_end (2100) is OUTSIDE the half-open interval [1000, 2100).
    gene, tss_d, feat = annotate.annotate_gene("chr1", 2100, tx_list)
    assert gene == "GENE"
    # The only transcript ends at 2100 — so pos is just past. TSS is at
    # tx_start = 1000 (strand +), distance = 1000 - 2100 = -1100 (downstream).
    # Feature class must NOT be "gene_body"; either "promoter" (if within
    # +/- 1 kb of TSS — here |−1100| > 1000, so no) or "intergenic".
    assert feat == "intergenic", f"expected intergenic for pos==tx_end, got {feat}"
    assert tss_d == -1100


def test_annotate_gene_promoter_within_1kb():
    """Sanity: pos within +/- 1 kb of TSS (and NOT overlapping any tx) is
    classified as promoter."""
    # Plus-strand tx [2000, 3000); TSS at 2000. pos 1500 → 500 bp upstream
    # of TSS, no overlap → promoter.
    tx_list = [(2000, 3000, "+", "FOO")]
    gene, tss_d, feat = annotate.annotate_gene("chr1", 1500, tx_list)
    assert gene == "FOO"
    assert feat == "promoter"
    assert tss_d == 500  # upstream positive


def test_annotate_gene_inside_is_gene_body():
    """Sanity: pos inside [tx_start, tx_end) far from TSS is gene_body."""
    tx_list = [(1000, 5000, "+", "FOO")]
    # pos 3000 is inside; 3000 - 1000 = 2000 bp from TSS → beyond 1 kb
    # promoter window → gene_body.
    _, _, feat = annotate.annotate_gene("chr1", 3000, tx_list)
    assert feat == "gene_body"


# ---------- annotate_cgi half-open ----------


def test_annotate_cgi_island_end_is_exclusive():
    """Regression: pos == cgi_end should classify as outside the island
    (distance 1 bp) rather than 'island'."""
    cgi_list = [(100, 200)]
    # pos 199 is last inside base of [100, 200) → island.
    ctx, d = annotate.annotate_cgi("chr1", 199, cgi_list)
    assert ctx == "island" and d == 0
    # pos 200 is first base outside → shore with distance 1.
    ctx, d = annotate.annotate_cgi("chr1", 200, cgi_list)
    assert ctx == "shore"
    assert d == 1


def test_annotate_cgi_shore_shelf_open_sea_bands():
    """Standard Illumina convention: shore +/- 2 kb of CGI edge, shelf
    +/- 2-4 kb, open_sea beyond."""
    cgi_list = [(10_000, 12_000)]
    # 1 kb after end → shore
    assert annotate.annotate_cgi("chr1", 13_000, cgi_list)[0] == "shore"
    # 3 kb after end → shelf
    assert annotate.annotate_cgi("chr1", 15_000, cgi_list)[0] == "shelf"
    # 10 kb before start → open_sea
    assert annotate.annotate_cgi("chr1", 0, cgi_list)[0] == "open_sea"
