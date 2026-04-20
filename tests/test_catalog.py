"""Catalog builder tests."""

from __future__ import annotations

from pathlib import Path

from thermocas.catalog import build_catalog, stream_catalog
from thermocas.models import Strand
from thermocas.pam_model import PamModel

REPO_ROOT = Path(__file__).resolve().parent.parent
PAM_YAML = REPO_ROOT / "config" / "pam_model.yaml"


def _write_fasta(tmp_path: Path, records: dict[str, str]) -> Path:
    p = tmp_path / "ref.fa"
    with open(p, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")
    return p


def test_build_catalog_finds_known_pam(tmp_path: Path):
    # Forward NNNNCGA at position 10 (C5* at index 14)
    seq = "AAAAAAAAAAACGTCGAGGGGG"
    fa = _write_fasta(tmp_path, {"chrSyn": seq})
    pam_model = PamModel.from_yaml(PAM_YAML)

    candidates, stats = build_catalog(fa, pam_model)

    assert stats.n_chromosomes == 1
    assert stats.n_filtered_out == 0
    forward_cga = [c for c in candidates if c.pam_family == "NNNNCGA" and c.strand == Strand.PLUS]
    assert any(c.critical_c_pos == 14 for c in forward_cga), (
        f"expected NNNNCGA C5 at pos 14; got {[c.critical_c_pos for c in forward_cga]}"
    )


def test_build_catalog_unique_ids(tmp_path: Path):
    fa = _write_fasta(tmp_path, {"chrA": "ACGTCGAACGTCGA", "chrB": "ACGTCGA"})
    pam_model = PamModel.from_yaml(PAM_YAML)
    candidates, _stats = build_catalog(fa, pam_model)
    ids = [c.candidate_id for c in candidates]
    assert len(ids) == len(set(ids)), "candidate_ids must be unique"


def test_region_filter_drops_candidates(tmp_path: Path):
    fa = _write_fasta(tmp_path, {"chr1": "ACGTCGAACGTCGAACGTCGA"})
    pam_model = PamModel.from_yaml(PAM_YAML)

    # keep only candidates with critical_c_pos < 10
    candidates, stats = build_catalog(
        fa, pam_model, region_filter=lambda chrom, pos: pos < 10
    )
    assert stats.n_filtered_out > 0
    assert all(c.critical_c_pos < 10 for c in candidates)


def test_local_seq_attached(tmp_path: Path):
    fa = _write_fasta(tmp_path, {"chr1": "G" * 200})
    fa = _write_fasta(tmp_path, {"chr1": "G" * 100 + "ACGTCGA" + "G" * 100})
    pam_model = PamModel.from_yaml(PAM_YAML)
    candidates, _ = build_catalog(fa, pam_model)
    forward_cga = [c for c in candidates if c.pam_family == "NNNNCGA" and c.strand == Strand.PLUS]
    assert forward_cga, "expected at least one NNNNCGA hit"
    site = forward_cga[0]
    # should be ~100 bp wide (50 each side, clipped at chrom ends)
    assert 90 <= len(site.local_seq_100bp) <= 100


def test_minus_strand_local_seq_contains_pam(tmp_path: Path):
    """Regression: local_seq_100bp must be on the *indicated strand*. For a
    minus-strand-only PAM hit the PAM string must appear in local_seq, not in
    its reverse complement."""

    # AAATCGACGTAAA: forward strand has no NNNNCGA; reverse complement is
    # TTTACGTCGATTT which contains ACGTCGA → minus-strand NNNNCGA.
    fa = _write_fasta(tmp_path, {"chr1": "AAATCGACGTAAA"})
    pam_model = PamModel.from_yaml(PAM_YAML)
    candidates, _ = build_catalog(fa, pam_model)
    minus_cga = [c for c in candidates if c.pam_family == "NNNNCGA" and c.strand == Strand.MINUS]
    assert minus_cga, "expected at least one minus-strand NNNNCGA"
    for c in minus_cga:
        assert c.pam in c.local_seq_100bp, (
            f"PAM {c.pam!r} must appear in strand-oriented local_seq {c.local_seq_100bp!r}"
        )


def test_stream_catalog_yields_same_as_build(tmp_path: Path):
    fa = _write_fasta(tmp_path, {"chr1": "ACGTCGAACGTCGA", "chr2": "TTTTCCATCGAA"})
    pam_model = PamModel.from_yaml(PAM_YAML)
    built, _ = build_catalog(fa, pam_model)
    streamed = list(stream_catalog(fa, pam_model))
    assert len(built) == len(streamed)
    assert {c.candidate_id for c in built} == {c.candidate_id for c in streamed}
