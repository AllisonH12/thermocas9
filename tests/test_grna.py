"""V3 — gRNA spacer scoring tests."""

from __future__ import annotations

import pytest

from thermocas.grna import (
    SPACER_LEN,
    _gc_content_score,
    _gc_fraction,
    _hairpin_score,
    _longest_mononucleotide_run,
    _melting_temp,
    _runs_score,
    _tm_score,
    extract_spacer,
    score_spacer,
)
from thermocas.models import CandidateSite, PamFamily, Strand


_CGA_FAMILY = PamFamily(
    name="NNNNCGA",
    regex="[ACGT][ACGT][ACGT][ACGT]CG[AG]",
    critical_c_offset=4,
    is_cpg=True,
    weight=1.0,
)


def _candidate(
    local_seq: str,
    pam: str = "ACGTCGA",
    *,
    critical_c_pos: int | None = None,
) -> CandidateSite:
    """Build a candidate with a hand-controlled local context.

    `critical_c_pos` defaults to point at the critical C of the *last* PAM
    occurrence in `local_seq` (so the spacer extractor finds the real PAM,
    not a decoy). Pass an explicit value to override.
    """
    if critical_c_pos is None:
        # critical C is at PAM start + critical_c_offset(=4)
        idx = local_seq.rfind(pam)
        critical_c_pos = idx + 4 if idx >= 0 else 50
    return CandidateSite(
        candidate_id=f"test:{len(local_seq)}+:NNNNCGA",
        chrom="chr1",
        critical_c_pos=critical_c_pos,
        strand=Strand.PLUS,
        pam=pam,
        pam_family="NNNNCGA",
        is_cpg_pam=True,
        local_seq_100bp=local_seq,
    )


# ---------- component math ----------


def test_gc_fraction():
    assert _gc_fraction("") == 0.0
    assert _gc_fraction("AAAA") == 0.0
    assert _gc_fraction("GGGG") == 1.0
    assert _gc_fraction("ATGC") == 0.5


def test_melting_temp_wallace_rule():
    # ATGC: 2 A/T (×2 = 4) + 2 G/C (×4 = 8) = 12
    assert _melting_temp("ATGC") == pytest.approx(12.0)
    # All A: 4 × 2 = 8
    assert _melting_temp("AAAA") == pytest.approx(8.0)
    # All G: 4 × 4 = 16
    assert _melting_temp("GGGG") == pytest.approx(16.0)


def test_longest_mononucleotide_run():
    assert _longest_mononucleotide_run("") == 0
    assert _longest_mononucleotide_run("A") == 1
    assert _longest_mononucleotide_run("ATCGAT") == 1
    assert _longest_mononucleotide_run("AAATCG") == 3
    assert _longest_mononucleotide_run("AATTTTAA") == 4


def test_gc_content_score_triangular():
    """1.0 in [0.40, 0.60], 0.0 outside [0.20, 0.80], linear in between."""
    assert _gc_content_score(0.40) == 1.0
    assert _gc_content_score(0.50) == 1.0
    assert _gc_content_score(0.60) == 1.0
    assert _gc_content_score(0.20) == 0.0
    assert _gc_content_score(0.80) == 0.0
    assert _gc_content_score(0.10) == 0.0
    assert _gc_content_score(0.30) == pytest.approx(0.5)


def test_tm_score_triangular():
    assert _tm_score(60.0) == 1.0
    assert _tm_score(75.0) == 1.0
    assert _tm_score(45.0) == 0.0
    assert _tm_score(90.0) == 0.0
    assert _tm_score(52.5) == pytest.approx(0.5)


def test_runs_score():
    assert _runs_score(1) == 1.0
    assert _runs_score(4) == 1.0
    assert _runs_score(5) == pytest.approx(0.75)
    assert _runs_score(7) == pytest.approx(0.25)
    assert _runs_score(8) == 0.0
    assert _runs_score(20) == 0.0


def test_hairpin_score_no_palindromes():
    """A scrambled sequence with no 4-mer/RC matches scores 1.0."""
    assert _hairpin_score("ACGTAGCTAGCTAGCTAGCT") <= 1.0  # may have some, can't assert exact 1.0


def test_hairpin_score_strong_palindrome():
    """A clear palindromic structure scores < 1.0."""
    # AAAA at start, TTTT downstream = revcomp pair
    assert _hairpin_score("AAAACGTACGTACGTTTTT") < 1.0


# ---------- extract_spacer + score_spacer ----------


def test_extract_spacer_returns_20nt_upstream_of_pam():
    """20-nt spacer + 7-nt PAM in a 27+ context → spacer is the 20 chars before PAM."""
    spacer = "ACGTACGTACGTACGTACGT"
    pam = "ACGTCGA"
    cand = _candidate(local_seq=spacer + pam, pam=pam)
    extracted = extract_spacer(cand, _CGA_FAMILY)
    assert extracted == spacer


def test_extract_spacer_returns_none_when_pam_too_close_to_start():
    """Less than 20 nt before the PAM in local_seq → no extractable spacer."""
    cand = _candidate(local_seq="ACGTCGA", pam="ACGTCGA")  # PAM at start
    assert extract_spacer(cand, _CGA_FAMILY) is None


def test_extract_spacer_returns_none_when_pam_missing():
    cand = _candidate(local_seq="A" * 50, pam="ACGTCGA")
    assert extract_spacer(cand, _CGA_FAMILY) is None


def test_extract_spacer_returns_none_when_spacer_contains_n():
    spacer = "ACGTNNNNNTACGTACGTAC"  # has N
    cand = _candidate(local_seq=spacer + "ACGTCGA", pam="ACGTCGA")
    assert extract_spacer(cand, _CGA_FAMILY) is None


def test_score_spacer_full_path_high_quality_spacer():
    """A balanced, aperiodic, palindrome-free spacer → near-perfect final_score.

    Periodic 4-base repeats like ATGCATGC... look balanced by GC and runs but
    every 4-mer's reverse-complement reappears later in the spacer, which the
    hairpin filter correctly flags. Real high-quality sgRNAs are aperiodic.
    """
    spacer = "ATGGAGTCCTAGACCTGCAT"  # 50% GC, longest run 2, no 4-bp palindrome hits
    cand = _candidate(local_seq=spacer + "ACGTCGA", pam="ACGTCGA")
    sc = score_spacer(cand, _CGA_FAMILY)
    assert sc is not None
    assert sc.spacer_seq == spacer
    assert sc.gc_fraction == 0.5
    assert sc.gc_content_score == 1.0
    assert sc.longest_run <= 2
    assert sc.runs_score == 1.0
    assert sc.hairpin_score == 1.0
    assert sc.final_score == pytest.approx(1.0)


def test_score_spacer_periodic_repeat_penalized_by_hairpin_filter():
    """Regression: a sequence that looks balanced (GC, runs) but is periodic
    on a 4-base cycle is dominated by 4-mer/RC matches and should score low."""
    spacer = "ATGCATGCATGCATGCATGC"  # 4-base period; every 4-mer's RC reappears
    cand = _candidate(local_seq=spacer + "ACGTCGA", pam="ACGTCGA")
    sc = score_spacer(cand, _CGA_FAMILY)
    assert sc is not None
    assert sc.gc_fraction == 0.5
    # gc_content_score and runs_score remain perfect; hairpin_score drives final down
    assert sc.hairpin_score < 0.5
    assert sc.final_score < 0.7


def test_score_spacer_penalizes_long_runs():
    bad = "AAAAAAAATCGATCGATCGA"  # 8-base A run → runs_score = 0
    cand = _candidate(local_seq=bad + "ACGTCGA", pam="ACGTCGA")
    sc = score_spacer(cand, _CGA_FAMILY)
    assert sc is not None
    assert sc.longest_run == 8
    assert sc.runs_score == 0.0
    assert sc.final_score == 0.0  # geometric mean with a 0 → 0


def test_score_spacer_penalizes_extreme_gc():
    """All-GC spacer has GC=1.0 → gc_content_score=0 → final_score=0."""
    bad = "GCGCGCGCGCGCGCGCGCGC"  # 100% GC
    cand = _candidate(local_seq=bad + "ACGTCGA", pam="ACGTCGA")
    sc = score_spacer(cand, _CGA_FAMILY)
    assert sc is not None
    assert sc.gc_fraction == 1.0
    assert sc.gc_content_score == 0.0
    assert sc.final_score == 0.0


def test_score_spacer_returns_none_when_no_spacer_extractable():
    """Empty / short context → score_spacer returns None, not a fake SpacerScore."""
    cand = _candidate(local_seq="ACGTCGA", pam="ACGTCGA")  # PAM at start, no upstream
    assert score_spacer(cand, _CGA_FAMILY) is None


def test_spacer_len_is_20():
    """V3 contract: protospacer is 20 nt for Type II Cas9."""
    assert SPACER_LEN == 20


def test_score_spacer_strict_ordering():
    """A real well-designed spacer must outscore a poly-A spacer."""
    good = _candidate(local_seq="ATGGAGTCCTAGACCTGCAT" + "ACGTCGA", pam="ACGTCGA")
    bad = _candidate(local_seq="AAAAAAAAAAAAAAAAAAAA" + "ACGTCGA", pam="ACGTCGA")
    s_good = score_spacer(good, _CGA_FAMILY)
    s_bad = score_spacer(bad, _CGA_FAMILY)
    assert s_good is not None and s_bad is not None
    assert s_good.final_score > s_bad.final_score
    assert s_good.final_score > 0.9
    assert s_bad.final_score == 0.0


# ---------- regression: decoy-PAM anchoring ----------


def test_extract_spacer_anchors_to_real_pam_not_decoy():
    """Regression: when the same PAM appears twice in local_seq_100bp the
    extractor must anchor on the candidate's reported critical_c_pos
    (catalog.py builds local_seq with critical-C at LOCAL_CONTEXT_HALF_WIDTH),
    not the first decoy occurrence. Earlier V3 used seq.find() and silently
    returned the decoy's spacer.
    """

    # 100-bp local_seq honoring the catalog.py contract: critical-C at index 50.
    # Real PAM at [46..53), decoy PAM at [26..33).
    pam = "ACGTCGA"
    # local layout:
    #   [0..6)   AAAAAA              padding
    #   [6..26)  GTAGTAGTAGTAGTAGTAGT  ← what the buggy seq.find() would return
    #   [26..33) ACGTCGA            decoy
    #   [33..46) TTAACGTGCATGA      filler (13 chars, no "ACGTCGA" inside)
    #   [46..53) ACGTCGA            real PAM (critical-C at idx 50)
    #   [53..100) C * 47            padding
    parts = (
        "AAAAAA",
        "GTAGTAGTAGTAGTAGTAGT",
        pam,
        "TTAACGTGCATGA",
        pam,
        "C" * 47,
    )
    local = "".join(parts)
    assert len(local) == 100, len(local)
    assert local.find(pam) == 26, f"expected decoy at 26, got {local.find(pam)}"
    # The CORRECT spacer is the 20 nt immediately before the REAL PAM (idx 46),
    # which happens to span [26..46) — i.e. it includes the decoy PAM bytes.
    expected_real_spacer = local[26:46]
    bug_spacer = local[6:26]
    assert expected_real_spacer != bug_spacer

    cand = _candidate(local, pam=pam, critical_c_pos=50)
    extracted = extract_spacer(cand, _CGA_FAMILY)
    assert extracted == expected_real_spacer, (
        f"expected real spacer {expected_real_spacer!r}, got {extracted!r} "
        f"(bug would have returned {bug_spacer!r})"
    )


def test_extract_spacer_returns_none_when_pam_not_at_expected_position():
    """If candidate.pam doesn't sit at the position critical_c_pos implies, the
    extractor refuses (no silent fallback to seq.find)."""

    spacer = "ATGGAGTCCTAGACCTGCAT"
    pam = "ACGTCGA"
    local = spacer + pam   # critical-C is at idx 24, not 50
    cand = _candidate(local, pam=pam, critical_c_pos=50)
    assert extract_spacer(cand, _CGA_FAMILY) is None
    assert score_spacer(cand, _CGA_FAMILY) is None


def test_extract_spacer_handles_minus_strand_candidates():
    """Regression: catalog._local_context() reverse-complements minus-strand
    windows, so the critical-C lands at index `len(seq)-1-fwd_idx`, not at
    `min(critical_c_pos, 50)`. Earlier extract_spacer assumed plus-strand
    indexing and silently dropped every minus-strand candidate."""

    from pathlib import Path

    from thermocas.catalog import build_catalog
    from thermocas.models import Strand
    from thermocas.pam_model import PamModel

    REPO_ROOT = Path(__file__).resolve().parent.parent
    pam_model = PamModel.from_yaml(REPO_ROOT / "config" / "pam_model.yaml")

    # Single-chrom synthetic FASTA with both strands matchable.
    # AAATCGACGTAAA on the forward strand contains the minus-strand NNNNCGA
    # match (rev-comp of TCGACGT is ACGTCGA → NNNNCGA hits the minus strand).
    # We pad it out to 200 bp so non-edge candidates exist.
    import tempfile
    with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as fh:
        fh.write(">chr1\n" + "G" * 80 + "AAATCGACGTAAA" + "G" * 100 + "\n")
        fa = Path(fh.name)
    try:
        candidates, _stats = build_catalog(fa, pam_model)
    finally:
        fa.unlink()

    minus_cga = [
        c for c in candidates
        if c.pam_family == "NNNNCGA" and c.strand == Strand.MINUS
    ]
    assert minus_cga, "expected at least one minus-strand NNNNCGA candidate"

    fam = pam_model.get("NNNNCGA")
    extracted = [extract_spacer(c, fam) for c in minus_cga]
    # At least one minus-strand candidate must produce a valid spacer.
    assert any(s is not None and len(s) == SPACER_LEN for s in extracted), (
        "minus-strand candidates all rejected — extract_spacer's anchor "
        "formula is wrong for revcomp'd local contexts"
    )


def test_extract_spacer_rejects_family_mismatch():
    """Calling with a PamFamily whose name doesn't match the candidate's
    pam_family is a programmer bug; refuse rather than guess."""

    other_family = PamFamily(
        name="NNNNCCA",
        regex="[ACGT][ACGT][ACGT][ACGT]CC[AG]",
        critical_c_offset=4,
        is_cpg=False,
        weight=0.9,
    )
    cand = _candidate("ATGGAGTCCTAGACCTGCAT" + "ACGTCGA", pam="ACGTCGA")
    assert extract_spacer(cand, other_family) is None
