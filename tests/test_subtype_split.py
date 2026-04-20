"""V2 — LocalArrayBackend.split_by_subtype + io.read_sample_subtypes."""

from __future__ import annotations

from pathlib import Path

import pytest

from thermocas.io import read_sample_subtypes
from thermocas.methylation_backend import LocalArrayBackend


def _write(p: Path, content: str) -> Path:
    p.write_text(content)
    return p


def test_read_sample_subtypes(tmp_path: Path):
    p = _write(
        tmp_path / "subtypes.tsv",
        "sample_id\tsubtype\nT1\tLumA\nT2\tLumB\nT3\tBasal\nT4\tLumA\n",
    )
    m = read_sample_subtypes(p)
    assert m == {"T1": "LumA", "T2": "LumB", "T3": "Basal", "T4": "LumA"}


def test_read_sample_subtypes_rejects_missing_columns(tmp_path: Path):
    p = _write(tmp_path / "bad.tsv", "sample_id\tx\nT1\tLumA\n")
    with pytest.raises(ValueError, match="sample_id, subtype"):
        read_sample_subtypes(p)


def test_split_by_subtype_partitions_columns(tmp_path: Path):
    """Each subtype's backend summarizes only its own samples."""
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t100\n")
    # 4 tumor samples: T1+T4 LumA (low beta), T2+T3 LumB (high beta)
    tumor = _write(
        tmp_path / "tumor.tsv",
        "probe_id\tT1\tT2\tT3\tT4\n"
        "cg001\t0.05\t0.85\t0.90\t0.10\n",
    )
    normal = _write(
        tmp_path / "normal.tsv",
        "probe_id\tN1\tN2\ncg001\t0.95\t0.92\n",
    )
    subtypes = _write(
        tmp_path / "subtypes.tsv",
        "sample_id\tsubtype\nT1\tLumA\nT2\tLumB\nT3\tLumB\nT4\tLumA\n",
    )

    backends = LocalArrayBackend.split_by_subtype(probes, tumor, normal, subtypes)
    assert set(backends) == {"LumA", "LumB"}

    luma = backends["LumA"].tumor_summary("cg001")
    lumb = backends["LumB"].tumor_summary("cg001")
    assert luma is not None and lumb is not None
    assert luma.n_samples == 2 and lumb.n_samples == 2
    assert luma.mean is not None and luma.mean < 0.20  # T1 + T4
    assert lumb.mean is not None and lumb.mean > 0.80  # T2 + T3

    # Normal side is shared across subtype backends.
    n_luma = backends["LumA"].normal_summary("cg001")
    n_lumb = backends["LumB"].normal_summary("cg001")
    assert n_luma is not None and n_lumb is not None
    assert n_luma.mean == n_lumb.mean


def test_split_by_subtype_rejects_missing_sample_in_map(tmp_path: Path):
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t100\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\tT1\tT_orphan\ncg001\t0.5\t0.5\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tN1\ncg001\t0.5\n")
    subtypes = _write(tmp_path / "subtypes.tsv", "sample_id\tsubtype\nT1\tLumA\n")
    with pytest.raises(ValueError, match="missing from"):
        LocalArrayBackend.split_by_subtype(probes, tumor, normal, subtypes)


def test_split_by_subtype_empty_map_rejected(tmp_path: Path):
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t100\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\tT1\ncg001\t0.5\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tN1\ncg001\t0.5\n")
    empty = _write(tmp_path / "subtypes.tsv", "sample_id\tsubtype\n")
    with pytest.raises(ValueError, match="no sample"):
        LocalArrayBackend.split_by_subtype(probes, tumor, normal, empty)
