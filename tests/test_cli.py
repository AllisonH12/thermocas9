"""End-to-end CLI test: build-catalog → score-cohort → aggregate."""

from __future__ import annotations

from pathlib import Path

import pytest

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


def test_cli_score_cohort_summary_backend_end_to_end(tmp_path: Path):
    """Regression: gdc-fetch writes per-probe summary TSVs, score-cohort
    --backend summary reads them. Earlier P1 bug: the formats didn't connect."""

    fa = _write(tmp_path / "ref.fa", ">chr1\nAAAAAAAAAAACGTCGAGGGGGGGGGGGGGGGGGGGG\n")
    cohort = _write(
        tmp_path / "cohort.yaml",
        "cohort:\n  name: T\n  tumor_dataset: a\n  normal_dataset: b\n  platform: HM450\n"
        "  min_samples_tumor: 1\n  min_samples_normal: 1\n",
    )
    probes = _write(tmp_path / "probes.tsv", "probe_id\tchrom\tpos\ncg001\tchr1\t14\n")
    # gdc-fetch-style summary TSV format: probe_id, n, mean, q25, q75 (NA-tolerant).
    tumor_summary = _write(
        tmp_path / "tumor_summary.tsv",
        "probe_id\tn\tmean\tq25\tq75\ncg001\t5\t0.05\t0.02\t0.10\n",
    )
    normal_summary = _write(
        tmp_path / "normal_summary.tsv",
        "probe_id\tn\tmean\tq25\tq75\ncg001\t5\t0.85\t0.78\t0.92\n",
    )

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
        "--backend", "summary",
        "--probe-annotation", str(probes),
        "--tumor-summary", str(tumor_summary),
        "--normal-summary", str(normal_summary),
        "--output", str(scored),
    ])
    assert rc == 0
    records = list(read_jsonl(scored, ScoredCandidate))
    # At least one candidate scored above zero — the summary backend correctly
    # produced a methylation observation, not UNOBSERVED.
    assert any(r.final_score > 0 for r in records)


def test_cli_score_cohort_summary_backend_requires_summary_inputs(tmp_path: Path, capsys):
    """Switching to --backend summary without --tumor-summary / --normal-summary fails clearly."""
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
        "--backend", "summary",
        "--output", str(tmp_path / "out.jsonl"),
    ])
    assert rc == 1
    assert "summary backend requires" in capsys.readouterr().err


def test_cli_inspect_top_zero_is_rejected_by_argparse(tmp_path: Path, capsys):
    """Regression: --top 0 used to crash with `IndexError: list index out of range`
    on an empty heap. argparse now rejects ≤0 with a clean error (SystemExit 2)."""
    out = _write(tmp_path / "tiny.jsonl", '{"x": 1}\n')
    with pytest.raises(SystemExit) as exc_info:
        main(["inspect", str(out), "--top", "0"])
    assert exc_info.value.code == 2
    assert "expected integer >= 1" in capsys.readouterr().err


def test_cli_inspect_reads_gzipped_jsonl(tmp_path: Path, capsys):
    """Regression: `_cmd_inspect` opened the file with Path.open() directly,
    failing on .jsonl.gz with a UTF-8 decode error on the gzip header byte."""
    import gzip
    from thermocas.io import write_jsonl
    from thermocas.models import CandidateSite, Strand

    cand = CandidateSite(
        candidate_id="x", chrom="chr1", critical_c_pos=10,
        strand=Strand.PLUS, pam="ACGTCGA", pam_family="NNNNCGA", is_cpg_pam=True,
    )
    out = tmp_path / "tiny.jsonl.gz"
    # write directly via gzip to be sure the file is *actually* gzipped
    with gzip.open(out, "wt") as f:
        f.write(cand.model_dump_json() + "\n")

    rc = main(["inspect", str(out)])
    assert rc == 0
    captured = capsys.readouterr().out
    assert "1 records (CandidateSite)" in captured


def test_cli_inspect_tie_break_matches_evaluate_ranking(tmp_path: Path, capsys):
    """Regression [P2]: inspect's top-K on ScoredCandidate JSONL used a
    bounded heap whose tuple comparison kept lex-larger candidate_ids on
    tied scores — the opposite of `evaluate_ranking`, which is
    `candidate_id` ascending.

    With three tied candidates a, b, c at the same final_score and
    --top 2, inspect must print ['a', 'b'], not ['c', 'b'].
    """
    import json

    # Three minimal ScoredCandidate records, all at final_score 1.0.
    # (We build them as plain dicts matching the model shape.)
    def rec(cid: str, score: float) -> str:
        return json.dumps({
            "candidate": {
                "candidate_id": cid, "chrom": "chr1", "critical_c_pos": 10,
                "strand": "+", "pam": "ACGTCGA", "pam_family": "NNNNCGA",
                "is_cpg_pam": True, "local_seq_100bp": "",
                "nearest_gene": None, "regulatory_context": None,
            },
            "observation": {
                "candidate_id": cid, "cohort_name": "T",
                "evidence_class": "exact", "evidence_distance_bp": 0,
                "probe_id": "cg001",
                "beta_tumor_mean": 0.05, "beta_tumor_q25": 0.02, "beta_tumor_q75": 0.10,
                "n_samples_tumor": 400,
                "beta_normal_mean": 0.85, "beta_normal_q25": 0.78, "beta_normal_q75": 0.92,
                "n_samples_normal": 80,
            },
            "components": {
                "sequence_score": 1.0, "selectivity_score": 1.0,
                "confidence_score": 1.0, "heterogeneity_penalty": 0.0,
                "low_coverage_penalty": 0.0,
            },
            "final_score": score,
            "probabilistic": None, "spacer": None,
        })

    out = tmp_path / "tied.jsonl"
    out.write_text("\n".join(rec(c, 1.0) for c in ("a", "b", "c")) + "\n")

    rc = main(["inspect", str(out), "--top", "2"])
    assert rc == 0
    printed = capsys.readouterr().out
    # `a` and `b` should appear in the output; `c` should NOT.
    # (Using a loose check — exact layout can change with formatting.)
    lines = [ln for ln in printed.splitlines() if ln.startswith("     1  ") or ln.startswith("     2  ")]
    assert len(lines) == 2
    assert " a " in lines[0] or lines[0].split()[1] == "a", lines
    assert " b " in lines[1] or lines[1].split()[1] == "b", lines
    # Explicit anti-check: `c` must not be present in the emitted top-2 block.
    top_block = "\n".join(lines)
    assert " c " not in top_block and "  c  " not in top_block, top_block


def test_cli_inspect_tie_break_handles_prefix_candidate_ids(tmp_path: Path, capsys):
    """Regression: prefix-related candidate IDs must still use native
    candidate_id ascending order inside tied scores.

    With tied `a`, `ab`, `ac` and --top 2, inspect must print `a`, `ab`.
    """
    import json

    def rec(cid: str, score: float) -> str:
        return json.dumps({
            "candidate": {
                "candidate_id": cid, "chrom": "chr1", "critical_c_pos": 10,
                "strand": "+", "pam": "ACGTCGA", "pam_family": "NNNNCGA",
                "is_cpg_pam": True, "local_seq_100bp": "",
                "nearest_gene": None, "regulatory_context": None,
            },
            "observation": {
                "candidate_id": cid, "cohort_name": "T",
                "evidence_class": "exact", "evidence_distance_bp": 0,
                "probe_id": "cg001",
                "beta_tumor_mean": 0.05, "beta_tumor_q25": 0.02, "beta_tumor_q75": 0.10,
                "n_samples_tumor": 400,
                "beta_normal_mean": 0.85, "beta_normal_q25": 0.78, "beta_normal_q75": 0.92,
                "n_samples_normal": 80,
            },
            "components": {
                "sequence_score": 1.0, "selectivity_score": 1.0,
                "confidence_score": 1.0, "heterogeneity_penalty": 0.0,
                "low_coverage_penalty": 0.0,
            },
            "final_score": score,
            "probabilistic": None, "spacer": None,
        })

    out = tmp_path / "tied_prefix.jsonl"
    out.write_text("\n".join(rec(c, 1.0) for c in ("a", "ab", "ac")) + "\n")

    rc = main(["inspect", str(out), "--top", "2"])
    assert rc == 0
    printed = capsys.readouterr().out
    lines = [ln for ln in printed.splitlines() if ln.startswith("     1  ") or ln.startswith("     2  ")]
    assert len(lines) == 2
    assert lines[0].split()[1] == "a", lines
    assert lines[1].split()[1] == "ab", lines


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
