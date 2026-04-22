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


def test_top_k_by_tie_break_handles_prefix_candidate_ids(tmp_path: Path):
    """Regression: byte-negation of UTF-8 bytes does not invert Python's
    native string ordering for prefix-related IDs (`a` < `ab` < `ac`).

    With all three records tied and K=2, the shortlist must keep
    `['a', 'ab']`, not `['ab', 'ac']`.
    """

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

    jsonl = tmp_path / "tied_prefix.jsonl"
    jsonl.write_text("\n".join(json.dumps(rec(c, 1.0)) for c in ("a", "ab", "ac")) + "\n")

    out = annotate.top_k_by("final_score", jsonl, 2)
    got = [r["candidate"]["candidate_id"] for r in out]
    assert got == ["a", "ab"], got


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


def test_annotate_gene_checks_long_earlier_transcript_overlap():
    """Regression: tx_list is sorted by tx_start, not tx_end. A short
    later-starting transcript can end before `pos` while an earlier-starting
    long transcript still overlaps. The scan must not terminate early on the
    short transcript's end coordinate.
    """

    tx_list = [
        (100, 5000, "+", "LONG"),
        (1000, 1100, "+", "SHORT_A"),
        (1500, 1600, "+", "SHORT_B"),
    ]
    gene, _, feat = annotate.annotate_gene("chr1", 3000, tx_list)
    assert gene == "LONG"
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


# ---------- annotate_repeat ----------


def test_annotate_repeat_inside_and_outside():
    rep_list = [
        (100, 200, "ALU", "SINE", "Alu"),
        (500, 700, "L1HS", "LINE", "L1"),
    ]
    # 150 is inside [100,200)
    in_rep, cls, fam, name = annotate.annotate_repeat("chr1", 150, rep_list)
    assert (in_rep, cls, fam, name) == (True, "SINE", "Alu", "ALU")
    # 200 is half-open exclusive end → outside
    in_rep, cls, fam, name = annotate.annotate_repeat("chr1", 200, rep_list)
    assert in_rep is False
    assert (cls, fam, name) == ("-", "-", "-")
    # between repeats → outside
    assert annotate.annotate_repeat("chr1", 300, rep_list)[0] is False
    # inside second
    assert annotate.annotate_repeat("chr1", 600, rep_list)[1:] == ("LINE", "L1", "L1HS")


def test_annotate_repeat_empty_list():
    assert annotate.annotate_repeat("chr1", 42, []) == (False, "-", "-", "-")


def test_annotate_repeat_finds_outer_nested_repeat():
    """Regression: pos may be covered by an outer repeat whose start lies
    several entries earlier than `idx` because of intervening shorter
    non-overlapping repeats. The earlier one-step-lookback implementation
    returned (False, ...) here, mislabeling the candidate as not in any
    repeat when it was actually inside an L1.
    """
    rep_list = [
        (100, 1000, "L1HS", "LINE", "L1"),    # outer; covers pos = 400
        (200, 250, "FLAM_C", "SINE", "Alu"),  # earlier; ends before pos
        (300, 350, "AluY", "SINE", "Alu"),    # nearest by start; ends before pos
    ]
    in_rep, cls, fam, name = annotate.annotate_repeat("chr1", 400, rep_list)
    assert (in_rep, cls, fam, name) == (True, "LINE", "L1", "L1HS")


def test_annotate_repeat_prefers_innermost_when_nested():
    """When pos is inside both an outer and an inner repeat, return the
    innermost (largest-start). The backward scan from `idx` encounters it
    first, so first-match-wins gives the right behavior.
    """
    rep_list = [
        (100, 1000, "L1HS", "LINE", "L1"),    # outer
        (400, 500, "AluY", "SINE", "Alu"),    # inner; covers pos = 450
    ]
    in_rep, cls, fam, name = annotate.annotate_repeat("chr1", 450, rep_list)
    assert (in_rep, cls, fam, name) == (True, "SINE", "Alu", "AluY")


# ---------- annotate_regulatory (DNase-HS clusters) ----------


def test_annotate_regulatory_inside_and_outside():
    reg_list = [
        (100, 200, 50),   # broad HS: 50 cell types
        (500, 600, 2),    # narrow HS: 2 cell types
    ]
    in_reg, count = annotate.annotate_regulatory("chr1", 150, reg_list)
    assert in_reg is True and count == 50
    # 200 is exclusive end → outside
    in_reg, count = annotate.annotate_regulatory("chr1", 200, reg_list)
    assert in_reg is False and count == 0
    # narrow HS returns its source count
    assert annotate.annotate_regulatory("chr1", 550, reg_list) == (True, 2)
    # far from any cluster
    assert annotate.annotate_regulatory("chr1", 10_000, reg_list) == (False, 0)


def test_annotate_regulatory_empty_list():
    assert annotate.annotate_regulatory("chr1", 42, []) == (False, 0)


# ---------- compute_flags (shortlist triage rules) ----------


def _card(**overrides) -> dict:
    base = {
        "feature_class": "intergenic",
        "cpg_island_context": "open_sea",
        "in_repeat": False,
        "repeat_family": "-",
        "in_dnase_cluster": False,
        "dnase_source_count": 0,
        "p_trust": 0.95,
        "delta_beta": 0.85,
    }
    base.update(overrides)
    return base


def test_compute_flags_island_localized_promoter_is_strong():
    flags = annotate.compute_flags(_card(
        feature_class="promoter", cpg_island_context="island",
    ))
    assert any(f.startswith("STRONG: island-localized promoter") for f in flags)


def test_compute_flags_active_promoter_requires_dnase_breadth():
    # Promoter + DNase in 15 cell types → STRONG
    flags = annotate.compute_flags(_card(
        feature_class="promoter", in_dnase_cluster=True, dnase_source_count=15,
    ))
    assert any("active promoter" in f and "15" in f for f in flags)
    # Promoter + DNase in only 3 cell types (below threshold) → no STRONG flag
    flags = annotate.compute_flags(_card(
        feature_class="promoter", in_dnase_cluster=True, dnase_source_count=3,
    ))
    assert not any("active promoter" in f for f in flags)


def test_compute_flags_in_repeat_is_caution():
    flags = annotate.compute_flags(_card(
        in_repeat=True, repeat_family="L1",
    ))
    assert any("CAUTION: overlaps L1 repeat" in f for f in flags)


def test_compute_flags_sparse_evidence_caution():
    flags = annotate.compute_flags(_card(p_trust=0.05))
    assert any("sparse evidence" in f for f in flags)
    # Threshold is strict <: 0.10 itself is not sparse.
    assert not any("sparse evidence" in f
                   for f in annotate.compute_flags(_card(p_trust=0.10)))


def test_compute_flags_small_differential_caution():
    flags = annotate.compute_flags(_card(delta_beta=0.15))
    assert any("small differential" in f for f in flags)
    # 0.30 is exactly the threshold and must NOT trigger (strict <).
    assert not any("small differential" in f
                   for f in annotate.compute_flags(_card(delta_beta=0.30)))


def test_compute_flags_gene_body_without_dnase_is_note():
    flags = annotate.compute_flags(_card(
        feature_class="gene_body", in_dnase_cluster=False,
    ))
    assert any("gene body with no DNase support" in f for f in flags)


def test_compute_flags_none_values_do_not_trigger():
    # When annotation layers are unavailable (None), the rule that depends
    # on them must not fire — tri-valued semantics matter.
    flags = annotate.compute_flags(_card(
        in_repeat=None, in_dnase_cluster=None,
        feature_class="gene_body",  # would be NOTE if in_dnase_cluster were False
    ))
    assert not any("overlaps" in f for f in flags)
    assert not any("gene body with no DNase" in f for f in flags)
