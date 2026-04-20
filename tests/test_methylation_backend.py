"""Methylation backend tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from thermocas.methylation_backend import (
    GDCBackend,
    LocalArrayBackend,
    LocalSummaryBackend,
    MethylationBackend,
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


def test_summarize_quartiles_stay_in_observed_range_for_small_n():
    """Regression: V3 method='exclusive' extrapolated outside [0.85, 0.90]
    (q25=0.8375, q75=0.9125) for two-sample input — invented support that
    never existed. The framework now uses 'inclusive' + a hard clamp to the
    empirical sample range."""

    s = _summarize("cg001", [0.85, 0.90])
    assert s.q25 is not None and s.q75 is not None
    assert 0.85 <= s.q25 <= 0.90
    assert 0.85 <= s.q75 <= 0.90

    s2 = _summarize("cg002", [0.05, 0.10])
    assert s2.q25 is not None and s2.q75 is not None
    assert 0.05 <= s2.q25 <= 0.10
    assert 0.05 <= s2.q75 <= 0.10


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


def test_gdc_backend_is_not_a_methylation_backend():
    """V3 — GDCBackend is intentionally not a MethylationBackend.

    The class can't satisfy the backend contract (probes() returns empty,
    one instance is single-sided). The supported live-data path is
    gdc-fetch → LocalSummaryBackend.
    """
    assert not issubclass(GDCBackend, MethylationBackend)


def test_gdc_backend_max_files_caps_downloads(tmp_path: Path):
    """V3 — max_files lets tests / quota-conscious users limit per-cohort downloads."""
    import json

    files_payload = json.dumps({
        "data": {"hits": [
            {"file_id": f"F{i}"} for i in range(1, 11)
        ]}
    }).encode()
    file_payload = b"Composite Element REF\tBeta_value\ncg001\t0.30\n"
    downloaded: list[str] = []

    def fake_opener(url: str, body: bytes | None, headers: dict[str, str]) -> bytes:
        if url.endswith("/files"):
            return files_payload
        downloaded.append(url)
        return file_payload

    backend = GDCBackend(
        project_id="TCGA-X", cache_dir=tmp_path,
        sample_type="Primary Tumor", max_files=3,
        url_opener=fake_opener,
    )
    backend.build_summaries()
    assert len(downloaded) == 3, f"max_files=3 should cap downloads to 3, got {len(downloaded)}"


# ---------- LocalSummaryBackend ----------


def test_local_summary_backend_round_trips_export(tmp_path: Path):
    """V3 — gdc-fetch's export_summaries output must be readable by LocalSummaryBackend."""
    import json

    files_payload = json.dumps({"data": {"hits": [{"file_id": "F1"}, {"file_id": "F2"}]}}).encode()
    file_a = b"Composite Element REF\tBeta_value\ncg001\t0.05\ncg002\t0.40\n"
    file_b = b"Composite Element REF\tBeta_value\ncg001\t0.10\ncg002\tNA\n"

    def fake_opener(url, body, headers):
        if url.endswith("/files"):
            return files_payload
        return file_a if url.endswith("/F1") else file_b

    # Step 1: gdc-fetch-style export
    tumor_summary = tmp_path / "tumor_summary.tsv"
    GDCBackend(
        project_id="TCGA-X", cache_dir=tmp_path,
        sample_type="Primary Tumor", url_opener=fake_opener,
    ).export_summaries(tumor_summary)

    # Reset cache by using a different cache dir for normal (real callers do per-side dirs;
    # here we just write a hand-crafted normal summary)
    normal_summary = _write(
        tmp_path / "normal_summary.tsv",
        "probe_id\tn\tmean\tq25\tq75\ncg001\t10\t0.85\t0.78\t0.92\ncg002\t10\t0.50\t0.45\t0.55\n",
    )
    probes = _write(
        tmp_path / "probes.tsv",
        "probe_id\tchrom\tpos\ncg001\tchr1\t100\ncg002\tchr1\t200\n",
    )

    # Step 2: LocalSummaryBackend reads the gdc-fetch output
    backend = LocalSummaryBackend(probes, tumor_summary, normal_summary)
    assert isinstance(backend, MethylationBackend)
    probe_list = list(backend.probes())
    assert {p.probe_id for p in probe_list} == {"cg001", "cg002"}
    t = backend.tumor_summary("cg001")
    n = backend.normal_summary("cg001")
    assert t is not None and n is not None
    assert t.n_samples == 2
    assert n.n_samples == 10
    assert n.mean is not None and n.mean > 0.80


def test_local_summary_backend_handles_na_summary_cells(tmp_path: Path):
    """NA cells in the summary TSV become None (not zeros)."""
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t1\n")
    s = _write(
        tmp_path / "summary.tsv",
        "probe_id\tn\tmean\tq25\tq75\ncg001\t0\tNA\tNA\tNA\n",
    )
    backend = LocalSummaryBackend(probes, s, s)
    res = backend.tumor_summary("cg001")
    assert res is not None
    assert res.mean is None and res.q25 is None and res.q75 is None


def test_probe_annotation_rejects_duplicate_probe_ids(tmp_path: Path):
    """Regression: LocalArrayBackend used to alias one probe_id across two
    genomic positions, fabricating EXACT evidence at the wrong locus."""
    probes = _write(
        tmp_path / "probes.tsv",
        "probe_id\tchrom\tpos\ncg001\tchr1\t10\ncg001\tchr1\t999\n",
    )
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\ts1\ncg001\t0.1\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\ncg001\t0.9\n")
    with pytest.raises(ValueError, match="duplicate probe_id"):
        LocalArrayBackend(probes, tumor, normal)


def test_local_summary_backend_rejects_duplicate_probe_rows(tmp_path: Path):
    """Regression: duplicate probe_id in summary TSV used to silently overwrite."""
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t1\n")
    s = _write(
        tmp_path / "summary.tsv",
        "probe_id\tn\tmean\tq25\tq75\n"
        "cg001\t10\t0.05\t0.02\t0.10\n"
        "cg001\t20\t0.85\t0.78\t0.92\n",
    )
    with pytest.raises(ValueError, match="duplicate probe_id"):
        LocalSummaryBackend(probes, s, s)


def test_local_summary_backend_rejects_missing_required_columns(tmp_path: Path):
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t1\n")
    bad = _write(tmp_path / "summary.tsv", "probe_id\tmean\ncg001\t0.5\n")
    with pytest.raises(ValueError, match="probe_id, n, mean"):
        LocalSummaryBackend(probes, bad, bad)


def test_gdc_backend_summarizes_normal_side(tmp_path: Path):
    """V3 — GDCBackend.build_summaries works for either side via sample_type."""
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
    summaries = backend.build_summaries()
    assert "cg001" in summaries
    assert summaries["cg001"].mean == pytest.approx(0.85)
