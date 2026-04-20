"""End-to-end CLI test: build-catalog → score-cohort → aggregate."""

from __future__ import annotations

from pathlib import Path

from thermocas.cli import main
from thermocas.io import read_jsonl
from thermocas.models import CandidateSite, PanCancerAggregate, ScoredCandidate

REPO_ROOT = Path(__file__).resolve().parent.parent
PAM_YAML = REPO_ROOT / "config" / "pam_model.yaml"


def _write(p: Path, content: str) -> Path:
    p.write_text(content)
    return p


def test_cli_full_pipeline(tmp_path: Path):
    # 1. tiny FASTA with a known NNNNCGA at C5* = 14
    fa = _write(tmp_path / "ref.fa", ">chr1\nAAAAAAAAAAACGTCGAGGGGGGGGG\n")

    # 2. cohort YAML
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n"
        "  name: TEST\n"
        "  tumor_dataset: TEST-T\n"
        "  normal_dataset: TEST-N\n"
        "  platform: HM450\n"
        "  min_samples_tumor: 1\n"
        "  min_samples_normal: 1\n",
    )

    # 3. probe + beta files (one probe at the candidate position)
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t14\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\ts1\ts2\ncg001\t0.05\t0.10\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\tn2\ncg001\t0.85\t0.90\n")

    catalog = tmp_path / "catalog.jsonl"
    rc = main([
        "build-catalog",
        "--reference", str(fa),
        "--pam-model", str(PAM_YAML),
        "--output", str(catalog),
    ])
    assert rc == 0
    catalog_records = list(read_jsonl(catalog, CandidateSite))
    assert len(catalog_records) >= 1
    assert any(c.critical_c_pos == 14 for c in catalog_records)

    scored = tmp_path / "scored.jsonl"
    rc = main([
        "score-cohort",
        "--catalog", str(catalog),
        "--cohort", str(cohort),
        "--pam-model", str(PAM_YAML),
        "--backend", "local",
        "--probe-annotation", str(probes),
        "--tumor-beta", str(tumor),
        "--normal-beta", str(normal),
        "--output", str(scored),
    ])
    assert rc == 0
    scored_records = list(read_jsonl(scored, ScoredCandidate))
    assert any(s.final_score > 0 for s in scored_records)

    pan = tmp_path / "pan.jsonl"
    rc = main([
        "aggregate",
        "--scored", f"TEST={scored}",
        "--high-score-threshold", "0.05",
        "--output", str(pan),
    ])
    assert rc == 0
    pan_records = list(read_jsonl(pan, PanCancerAggregate))
    assert len(pan_records) >= 1
    assert any(p.n_cohorts_high_score >= 1 for p in pan_records)


def test_cli_score_cohort_requires_local_inputs(tmp_path: Path, capsys):
    # missing --probe-annotation should fail with a clear error
    catalog = _write(tmp_path / "catalog.jsonl", "")
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n  name: T\n  tumor_dataset: a\n  normal_dataset: b\n  platform: HM450\n",
    )
    rc = main([
        "score-cohort",
        "--catalog", str(catalog),
        "--cohort", str(cohort),
        "--pam-model", str(PAM_YAML),
        "--backend", "local",
        "--output", str(tmp_path / "out.jsonl"),
    ])
    assert rc == 1
    err = capsys.readouterr().err
    assert "--probe-annotation" in err


def test_cli_score_cohort_with_probabilistic(tmp_path: Path):
    """V2 — --probabilistic flag attaches ProbabilisticScore to each ScoredCandidate."""

    fa = _write(tmp_path / "ref.fa", ">chr1\nAAAAAAAAAAACGTCGAGGGGGGGGGGGGGGGGGGGG\n")
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n  name: T\n  tumor_dataset: a\n  normal_dataset: b\n  platform: HM450\n"
        "  min_samples_tumor: 1\n  min_samples_normal: 1\n",
    )
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t14\n")
    tumor = _write(tmp_path / "tumor.tsv", "probe_id\ts1\ts2\ncg001\t0.05\t0.10\n")
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\tn2\ncg001\t0.85\t0.90\n")

    catalog = tmp_path / "catalog.jsonl"
    main([
        "build-catalog", "--reference", str(fa),
        "--pam-model", str(PAM_YAML), "--output", str(catalog),
    ])

    scored = tmp_path / "scored.jsonl"
    rc = main([
        "score-cohort",
        "--catalog", str(catalog),
        "--cohort", str(cohort),
        "--pam-model", str(PAM_YAML),
        "--backend", "local",
        "--probe-annotation", str(probes),
        "--tumor-beta", str(tumor),
        "--normal-beta", str(normal),
        "--probabilistic",
        "--output", str(scored),
    ])
    assert rc == 0
    records = list(read_jsonl(scored, ScoredCandidate))
    assert any(r.probabilistic is not None for r in records)
    # At least one record should have a non-zero probabilistic factor.
    nonzero = [r for r in records if r.probabilistic and r.probabilistic.p_targetable_tumor > 0]
    assert nonzero


def test_cli_score_cohort_by_subtype(tmp_path: Path):
    """V2 — --sample-subtypes fans out into per-subtype JSONLs with cohort_name suffixed."""

    fa = _write(tmp_path / "ref.fa", ">chr1\nAAAAAAAAAAACGTCGAGGGGGGGGGGGGGGGGGGGG\n")
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n  name: BRCA\n  tumor_dataset: a\n  normal_dataset: b\n  platform: HM450\n"
        "  min_samples_tumor: 1\n  min_samples_normal: 1\n",
    )
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t14\n")
    tumor = _write(
        tmp_path / "tumor.tsv",
        "probe_id\tT1\tT2\tT3\tT4\ncg001\t0.05\t0.10\t0.85\t0.90\n",
    )
    normal = _write(tmp_path / "normal.tsv", "probe_id\tn1\tn2\ncg001\t0.85\t0.90\n")
    subtypes = _write(
        tmp_path / "subtypes.tsv",
        "sample_id\tsubtype\nT1\tLumA\nT2\tLumA\nT3\tLumB\nT4\tLumB\n",
    )

    catalog = tmp_path / "catalog.jsonl"
    main([
        "build-catalog", "--reference", str(fa),
        "--pam-model", str(PAM_YAML), "--output", str(catalog),
    ])

    scored_prefix = tmp_path / "scored.jsonl"
    rc = main([
        "score-cohort",
        "--catalog", str(catalog),
        "--cohort", str(cohort),
        "--pam-model", str(PAM_YAML),
        "--backend", "local",
        "--probe-annotation", str(probes),
        "--tumor-beta", str(tumor),
        "--normal-beta", str(normal),
        "--sample-subtypes", str(subtypes),
        "--output", str(scored_prefix),
    ])
    assert rc == 0

    luma_path = tmp_path / "scored.LumA.jsonl"
    lumb_path = tmp_path / "scored.LumB.jsonl"
    assert luma_path.exists()
    assert lumb_path.exists()

    luma_records = list(read_jsonl(luma_path, ScoredCandidate))
    lumb_records = list(read_jsonl(lumb_path, ScoredCandidate))
    # cohort_name must be suffixed with the subtype
    assert all(r.observation.cohort_name == "BRCA::LumA" for r in luma_records)
    assert all(r.observation.cohort_name == "BRCA::LumB" for r in lumb_records)


def test_cli_gdc_backend_rejected_in_score_cohort(tmp_path: Path, capsys):
    """V2 — score-cohort still rejects --backend gdc; users should run gdc-fetch first."""
    catalog = _write(tmp_path / "catalog.jsonl", "")
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n  name: T\n  tumor_dataset: a\n  normal_dataset: b\n  platform: HM450\n",
    )
    rc = main([
        "score-cohort",
        "--catalog", str(catalog),
        "--cohort", str(cohort),
        "--pam-model", str(PAM_YAML),
        "--backend", "gdc",
        "--output", str(tmp_path / "out.jsonl"),
    ])
    assert rc == 1
    err = capsys.readouterr().err
    assert "gdc" in err.lower()
