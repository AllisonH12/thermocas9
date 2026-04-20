"""PAM model loader and matcher tests."""

from __future__ import annotations

from pathlib import Path

from thermocas.models import Strand
from thermocas.pam_model import PamModel, find_pam_matches, reverse_complement

REPO_ROOT = Path(__file__).resolve().parent.parent
PAM_YAML = REPO_ROOT / "config" / "pam_model.yaml"


def test_pam_model_loads_from_yaml():
    model = PamModel.from_yaml(PAM_YAML)
    names = {f.name for f in model.pam_families}
    assert {"NNNNCGA", "NNNNCCA"} <= names
    cga = model.get("NNNNCGA")
    assert cga.is_cpg is True
    assert cga.critical_c_offset == 4
    assert 0.0 < cga.weight <= 1.0


def test_reverse_complement_basic():
    assert reverse_complement("ACGTN") == "NACGT"
    # Case is preserved so callers can keep soft-masked repeat info.
    assert reverse_complement("aaaa") == "tttt"
    assert reverse_complement("AAAA") == "TTTT"


def test_find_pam_matches_forward_strand():
    model = PamModel.from_yaml(PAM_YAML)
    seq = "TTTTACGTCGAGGGG"   # NNNNCGA at index 4
    matches = find_pam_matches(seq, model)
    forward = [m for m in matches if m.strand == Strand.PLUS]
    assert any(m.family == "NNNNCGA" and m.critical_c_pos == 8 for m in forward)


def test_find_pam_matches_reverse_strand_coords_translate():
    model = PamModel.from_yaml(PAM_YAML)
    # Build a sequence so the forward strand has NO PAM but the reverse strand does.
    # Reverse-complement of NNNNCGA family member ACGTCGA is TCGACGT.
    # Embed TCGACGT on the forward strand → reverse strand reads ACGTCGA → matches NNNNCGA.
    seq = "AAATCGACGTAAA"  # forward has TCGACGT starting at index 3
    matches = find_pam_matches(seq, model)
    # critical C on the minus strand corresponds to a forward coord we can compute.
    rc = reverse_complement(seq)
    assert "ACGTCGA" in rc
    minus_matches = [m for m in matches if m.strand == Strand.MINUS]
    assert minus_matches, "expected at least one minus-strand PAM match"
    for m in minus_matches:
        # critical_c_pos must be within sequence bounds and on a C in the +strand sequence
        # (because reverse-complement of a G is C; the minus-strand PAM cytosine sits on a +strand G)
        assert 0 <= m.critical_c_pos < len(seq)
        assert seq[m.critical_c_pos] in "CG"  # forward position of critical C on minus strand


def test_find_pam_matches_no_false_positives():
    model = PamModel.from_yaml(PAM_YAML)
    # All-A sequence should match nothing
    matches = find_pam_matches("A" * 30, model)
    assert matches == []


def test_find_pam_matches_returns_overlapping_pams():
    """Regression: re.finditer on a fixed-width pattern silently drops overlapping
    matches. AAAACCACCA contains NNNNCCA at offsets 0 and 3 — both must be returned."""

    model = PamModel.from_yaml(PAM_YAML)
    seq = "AAAACCACCA"
    forward_ccA = [
        m for m in find_pam_matches(seq, model)
        if m.family == "NNNNCCA" and m.strand == Strand.PLUS
    ]
    starts = sorted(m.start for m in forward_ccA)
    assert starts == [0, 3], (
        f"expected NNNNCCA at offsets 0 and 3 on the forward strand, got {starts}"
    )
