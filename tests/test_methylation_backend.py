"""Methylation backend tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from thermocas.methylation_backend import (
    GDCBackend,
    LocalArrayBackend,
    _summarize,
)


def _write(p: Path, content: str) -> Path:
    p.write_text(content)
    return p


def test_summarize_handles_missing_values():
    s = _summarize("cg001", [0.10, None, 0.30, None, 0.50])
    assert s.n_samples == 3
    assert s.mean == pytest.approx((0.10 + 0.30 + 0.50) / 3)
    # quantiles undefined but must satisfy q25 <= mean <= q75
    assert s.q25 is not None and s.q75 is not None
    assert s.q25 <= s.mean <= s.q75


def test_summarize_no_samples():
    s = _summarize("cg001", [None, None])
    assert s.n_samples == 0
    assert s.mean is None and s.q25 is None and s.q75 is None


def test_summarize_single_sample_collapses():
    s = _summarize("cg001", [0.42])
    assert s.n_samples == 1
    assert s.mean == 0.42
    assert s.q25 == 0.42 and s.q75 == 0.42


def test_local_array_backend_loads_and_summarizes(tmp_path: Path):
    probes = _write(
        tmp_path / "probes.tsv",
        "probe_id\tchrom\tpos\ncg001\tchr1\t100\ncg002\tchr1\t200\n",
    )
    tumor = _write(
        tmp_path / "tumor.tsv",
        "probe_id\ts1\ts2\ts3\ts4\n"
        "cg001\t0.05\t0.07\t0.10\t0.08\n"
        "cg002\t0.20\t0.30\t0.40\t0.35\n",
    )
    normal = _write(
        tmp_path / "normal.tsv",
        "probe_id\tn1\tn2\n"
        "cg001\t0.85\t0.90\n"
        "cg002\t0.50\t0.60\n",
    )

    backend = LocalArrayBackend(probes, tumor, normal)
    probe_list = list(backend.probes())
    assert {p.probe_id for p in probe_list} == {"cg001", "cg002"}

    t1 = backend.tumor_summary("cg001")
    n1 = backend.normal_summary("cg001")
    assert t1 is not None and n1 is not None
    assert t1.n_samples == 4
    assert t1.mean is not None and t1.mean < 0.15
    assert n1.n_samples == 2
    assert n1.mean is not None and n1.mean > 0.80


def test_local_array_backend_returns_none_for_unknown_probe(tmp_path: Path):
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t1\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\ts1\ncg001\t0.5\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\ncg001\t0.5\n")
    backend = LocalArrayBackend(probes, tumor, normal)
    assert backend.tumor_summary("missing") is None
    assert backend.normal_summary("missing") is None


def test_local_array_backend_rejects_bad_annotation(tmp_path: Path):
    bad = _write(tmp_path / "probes.tsv", "probe_id\tchrom\nMISSING_POS\tchr1\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\ts1\nMISSING_POS\t0.1\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\nMISSING_POS\t0.1\n")
    with pytest.raises(ValueError, match="probe_id, chrom, pos"):
        LocalArrayBackend(bad, tumor, normal)


def test_gdc_backend_uses_injected_opener_and_caches(tmp_path: Path):
    """V2 — exercise the real GDCBackend HTTP path with a mocked opener.

    We hand a fake `url_opener(url, body, headers) -> bytes` in. The first call
    should hit `/files`; subsequent `download_file` calls should hit `/data/...`
    once and then read from disk on a repeat.
    """
    import json

    calls: list[tuple[str, bytes | None]] = []

    files_payload = json.dumps({
        "data": {"hits": [
            {"file_id": "FA", "cases": [{"samples": [{"sample_type": "Primary Tumor"}]}]},
            {"file_id": "FB", "cases": [{"samples": [{"sample_type": "Primary Tumor"}]}]},
        ]}
    }).encode()
    file_a = b"Composite Element REF\tBeta_value\ncg001\t0.05\ncg002\t0.40\n"
    file_b = b"Composite Element REF\tBeta_value\ncg001\t0.10\ncg002\t0.50\n"

    def fake_opener(url: str, body: bytes | None, headers: dict[str, str]) -> bytes:
        calls.append((url, body))
        if url.endswith("/files"):
            return files_payload
        if url.endswith("/data/FA"):
            return file_a
        if url.endswith("/data/FB"):
            return file_b
        raise AssertionError(f"unexpected URL {url}")

    backend = GDCBackend(
        project_id="TCGA-BRCA",
        cache_dir=tmp_path,
        sample_type="Primary Tumor",
        url_opener=fake_opener,
    )
    summaries = backend.build_summaries()

    assert "cg001" in summaries and "cg002" in summaries
    assert summaries["cg001"].n_samples == 2
    assert summaries["cg001"].mean is not None
    # cached files written to disk
    assert (tmp_path / "FA.txt").exists()
    assert (tmp_path / "FB.txt").exists()

    # tumor_summary returns from build_summaries cache; normal_summary returns None
    assert backend.tumor_summary("cg001") is not None
    assert backend.normal_summary("cg001") is None

    # second call to download_file uses on-disk cache, no extra HTTP calls
    n_downloads_before = sum(1 for u, _ in calls if "/data/" in u)
    backend.download_file("FA")
    n_downloads_after = sum(1 for u, _ in calls if "/data/" in u)
    assert n_downloads_after == n_downloads_before


def test_gdc_backend_export_summaries(tmp_path: Path):
    """`export_summaries` writes a probe_id, n, mean, q25, q75 TSV."""
    import json

    files_payload = json.dumps({
        "data": {"hits": [{"file_id": "F1"}]}
    }).encode()
    file_payload = (
        b"Composite Element REF\tBeta_value\n"
        b"cg001\t0.20\ncg002\tNA\n"
    )

    def fake_opener(url: str, body: bytes | None, headers: dict[str, str]) -> bytes:
        return files_payload if url.endswith("/files") else file_payload

    backend = GDCBackend(
        project_id="TCGA-BRCA",
        cache_dir=tmp_path,
        sample_type="Primary Tumor",
        url_opener=fake_opener,
    )
    out = tmp_path / "summary.tsv"
    n = backend.export_summaries(out)
    assert n >= 1
    text = out.read_text()
    assert "probe_id\tn\tmean\tq25\tq75" in text
    assert "cg001" in text


def test_gdc_backend_normal_summary_when_sample_type_normal(tmp_path: Path):
    import json
    files_payload = json.dumps({"data": {"hits": [{"file_id": "F1"}]}}).encode()
    file_payload = b"Composite Element REF\tBeta_value\ncg001\t0.85\n"

    def fake_opener(url: str, body, headers):
        return files_payload if url.endswith("/files") else file_payload

    backend = GDCBackend(
        project_id="TCGA-BRCA",
        cache_dir=tmp_path,
        sample_type="Solid Tissue Normal",
        url_opener=fake_opener,
    )
    assert backend.normal_summary("cg001") is not None
    assert backend.tumor_summary("cg001") is None  # wrong side
