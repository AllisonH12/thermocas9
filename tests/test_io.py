"""IO helper tests — JSONL roundtrip, FASTA streaming, beta-matrix parsing."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from thermocas.io import (
    iter_fasta,
    read_beta_matrix,
    read_jsonl,
    read_tsv,
    write_jsonl,
    write_jsonl_atomic,
)
from thermocas.models import CandidateSite, Strand


def _candidate(cid: str, pos: int) -> CandidateSite:
    return CandidateSite(
        candidate_id=cid,
        chrom="chr1",
        critical_c_pos=pos,
        strand=Strand.PLUS,
        pam="AAAACGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def test_jsonl_roundtrip(tmp_path: Path):
    out = tmp_path / "candidates.jsonl"
    n = write_jsonl(out, [_candidate("a", 1), _candidate("b", 2)])
    assert n == 2
    rows = list(read_jsonl(out, CandidateSite))
    assert [r.candidate_id for r in rows] == ["a", "b"]


def test_jsonl_atomic_writes_then_replaces(tmp_path: Path):
    out = tmp_path / "candidates.jsonl"
    write_jsonl_atomic(out, [_candidate("x", 5)])
    assert out.exists()
    assert not (tmp_path / "candidates.jsonl.tmp").exists()


def test_iter_fasta_handles_multiline_records(tmp_path: Path):
    fa = tmp_path / "tiny.fa"
    fa.write_text(">chr1 description here\nACGT\nACGT\n>chr2\nTTTT\n")
    chroms = list(iter_fasta(fa))
    assert chroms == [("chr1", "ACGTACGT"), ("chr2", "TTTT")]


def test_iter_fasta_handles_gzip(tmp_path: Path):
    fa = tmp_path / "tiny.fa.gz"
    with gzip.open(fa, "wt") as f:
        f.write(">chrZ\nGGGGCCCC\n")
    assert list(iter_fasta(fa)) == [("chrZ", "GGGGCCCC")]


def test_read_tsv(tmp_path: Path):
    p = tmp_path / "probes.tsv"
    p.write_text("probe_id\tchrom\tpos\ncg001\tchr1\t100\ncg002\tchr2\t200\n")
    rows = list(read_tsv(p))
    assert rows == [
        {"probe_id": "cg001", "chrom": "chr1", "pos": "100"},
        {"probe_id": "cg002", "chrom": "chr2", "pos": "200"},
    ]


def test_read_beta_matrix_parses_and_handles_missing(tmp_path: Path):
    p = tmp_path / "betas.tsv"
    p.write_text(
        "probe_id\ts1\ts2\ts3\n"
        "cg001\t0.10\t0.20\t0.30\n"
        "cg002\t0.50\tNA\t0.70\n"
    )
    samples, betas = read_beta_matrix(p)
    assert samples == ["s1", "s2", "s3"]
    assert betas["cg001"] == [0.10, 0.20, 0.30]
    assert betas["cg002"] == [0.50, None, 0.70]


def test_read_beta_matrix_rejects_out_of_range(tmp_path: Path):
    p = tmp_path / "bad.tsv"
    p.write_text("probe_id\ts1\ncg001\t1.5\n")
    with pytest.raises(ValueError, match="not in"):
        read_beta_matrix(p)


def test_read_beta_matrix_rejects_row_length_mismatch(tmp_path: Path):
    p = tmp_path / "bad.tsv"
    p.write_text("probe_id\ts1\ts2\ncg001\t0.5\n")
    with pytest.raises(ValueError, match="expected 2"):
        read_beta_matrix(p)
