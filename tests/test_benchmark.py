"""V3 — cross-validation benchmark tests."""

from __future__ import annotations

import pytest

from thermocas.benchmark import evaluate_ranking, split_by_chrom
from thermocas.models import (
    CandidateSite,
    EvidenceClass,
    MethylationObservation,
    ScoreComponents,
    ScoredCandidate,
    Strand,
)


def _candidate(cid: str, chrom: str = "chr1", pos: int = 100) -> CandidateSite:
    return CandidateSite(
        candidate_id=cid,
        chrom=chrom,
        critical_c_pos=pos,
        strand=Strand.PLUS,
        pam="ACGTCGA",
        pam_family="NNNNCGA",
        is_cpg_pam=True,
    )


def _scored(cid: str, score: float) -> ScoredCandidate:
    cand = _candidate(cid)
    obs = MethylationObservation(
        candidate_id=cid, cohort_name="TEST",
        evidence_class=EvidenceClass.EXACT,
        evidence_distance_bp=0, probe_id="cg001",
        beta_tumor_mean=0.05, beta_tumor_q25=0.02, beta_tumor_q75=0.10,
        beta_normal_mean=0.85, beta_normal_q25=0.78, beta_normal_q75=0.92,
        n_samples_tumor=400, n_samples_normal=80,
    )
    components = ScoreComponents(
        sequence_score=1.0, selectivity_score=score, confidence_score=1.0,
    )
    return ScoredCandidate(
        candidate=cand, observation=obs, components=components, final_score=score,
    )


# ---------- split_by_chrom ----------


def test_split_by_chrom_basic():
    cands = [
        _candidate("a", chrom="chr1"),
        _candidate("b", chrom="chr2"),
        _candidate("c", chrom="chr1"),
        _candidate("d", chrom="chrX"),
    ]
    train, test = split_by_chrom(cands, holdout_chroms={"chr2", "chrX"})
    assert [c.candidate_id for c in train] == ["a", "c"]
    assert [c.candidate_id for c in test] == ["b", "d"]


def test_split_by_chrom_empty_holdout():
    cands = [_candidate("a"), _candidate("b")]
    train, test = split_by_chrom(cands, holdout_chroms=set())
    assert len(train) == 2 and len(test) == 0


def test_split_by_chrom_preserves_order():
    cands = [_candidate(f"c{i}", chrom="chr1") for i in range(5)]
    train, _ = split_by_chrom(cands, holdout_chroms=set())
    assert [c.candidate_id for c in train] == ["c0", "c1", "c2", "c3", "c4"]


# ---------- evaluate_ranking ----------


def test_evaluate_perfect_ranking():
    """Positives at the top → P@K = R@K = 1.0, AUC = 1.0."""
    pos_ids = {"p1", "p2", "p3"}
    scored = [
        _scored("p1", 0.9), _scored("p2", 0.8), _scored("p3", 0.7),
        _scored("n1", 0.3), _scored("n2", 0.2), _scored("n3", 0.1),
    ]
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=3)
    assert r.n_total == 6
    assert r.n_positives == 3 and r.n_negatives == 3
    assert r.precision_at_k == pytest.approx(1.0)
    assert r.recall_at_k == pytest.approx(1.0)
    assert r.roc_auc == pytest.approx(1.0)


def test_evaluate_inverted_ranking():
    """Positives at the bottom → P@K = R@K = 0, AUC = 0."""
    pos_ids = {"p1", "p2", "p3"}
    scored = [
        _scored("n1", 0.9), _scored("n2", 0.8), _scored("n3", 0.7),
        _scored("p1", 0.3), _scored("p2", 0.2), _scored("p3", 0.1),
    ]
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=3)
    assert r.precision_at_k == pytest.approx(0.0)
    assert r.recall_at_k == pytest.approx(0.0)
    assert r.roc_auc == pytest.approx(0.0)


def test_evaluate_random_ranking_auc_is_half():
    """Interleaved positives/negatives → AUC ≈ 0.5."""
    pos_ids = {"p1", "p2"}
    scored = [
        _scored("p1", 0.9), _scored("n1", 0.8), _scored("p2", 0.7), _scored("n2", 0.6),
    ]
    # P@2 = 1/2 (one of top-2 is positive), R@2 = 1/2 (1 of 2 positives in top-2)
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=2)
    assert r.precision_at_k == pytest.approx(0.5)
    assert r.recall_at_k == pytest.approx(0.5)
    # P(p>n): p1>n1, p1>n2, p2<n1, p2>n2 → 3/4 = 0.75
    assert r.roc_auc == pytest.approx(0.75)


def test_evaluate_handles_ties_via_half_credit():
    """Tied scores between positive and negative contribute 0.5 to AUC."""
    pos_ids = {"p1"}
    scored = [_scored("p1", 0.5), _scored("n1", 0.5)]
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=1)
    assert r.roc_auc == pytest.approx(0.5)


def test_evaluate_returns_none_metrics_when_no_positives_in_set():
    """Edge case: positives list has IDs not present in scored → metrics None."""
    pos_ids = {"missing"}
    scored = [_scored("n1", 0.5), _scored("n2", 0.3)]
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=1)
    assert r.n_positives == 0
    assert r.precision_at_k is None
    assert r.recall_at_k is None
    assert r.roc_auc is None


def test_evaluate_top_k_clamped_to_n_total():
    pos_ids = {"p1"}
    scored = [_scored("p1", 0.5), _scored("n1", 0.3)]
    r = evaluate_ranking(scored, pos_ids, cohort_name="T", top_k=100)
    # P@K computed against k = min(100, 2) = 2 → 1/2
    assert r.precision_at_k == pytest.approx(0.5)


def test_evaluate_records_held_out_chromosomes():
    pos_ids = {"p1"}
    scored = [_scored("p1", 0.5)]
    r = evaluate_ranking(
        scored, pos_ids, cohort_name="T", top_k=1,
        held_out_chromosomes=["chr1", "chr2"],
    )
    assert r.held_out_chromosomes == ["chr1", "chr2"]


def _scored_on_chrom(cid: str, score: float, chrom: str) -> ScoredCandidate:
    """Variant of _scored that places the candidate on a specific chromosome."""
    sc = _scored(cid, score)
    cand = sc.candidate.model_copy(update={"chrom": chrom})
    return sc.model_copy(update={"candidate": cand})


def test_evaluate_held_out_filters_to_only_those_chromosomes():
    """Regression: held_out_chromosomes is the cross-validation control,
    not metadata. With a positive on chr2 and a negative on chr1, holding
    out chr2 must evaluate ONLY the chr2 candidate (n_total=1)."""

    scored = [
        _scored_on_chrom("p_chr2", 0.9, "chr2"),
        _scored_on_chrom("n_chr1", 0.5, "chr1"),
    ]
    r = evaluate_ranking(
        scored, {"p_chr2"}, cohort_name="T", top_k=1,
        held_out_chromosomes=["chr2"],
    )
    assert r.n_total == 1
    assert r.n_positives == 1
    assert r.n_negatives == 0


def test_evaluate_no_enforce_holdout_evaluates_full_set():
    """`enforce_holdout=False` keeps the metadata but evaluates everything."""

    scored = [
        _scored_on_chrom("p_chr2", 0.9, "chr2"),
        _scored_on_chrom("n_chr1", 0.5, "chr1"),
    ]
    r = evaluate_ranking(
        scored, {"p_chr2"}, cohort_name="T", top_k=1,
        held_out_chromosomes=["chr2"],
        enforce_holdout=False,
    )
    assert r.n_total == 2
    assert r.held_out_chromosomes == ["chr2"]


def test_evaluate_empty_holdout_evaluates_full_set():
    """No held-out chromosomes → no filter applied (intuitive default)."""

    scored = [
        _scored_on_chrom("p_chr2", 0.9, "chr2"),
        _scored_on_chrom("n_chr1", 0.5, "chr1"),
    ]
    r = evaluate_ranking(scored, {"p_chr2"}, cohort_name="T", top_k=1)
    assert r.n_total == 2


def test_evaluate_rejects_zero_top_k():
    with pytest.raises(ValueError, match="top_k must be"):
        evaluate_ranking([], set(), cohort_name="T", top_k=0)


def test_evaluate_score_field_dispatch_unknown():
    with pytest.raises(ValueError, match="unknown score_field"):
        evaluate_ranking(
            [_scored("p1", 0.5)],
            {"p1"},
            cohort_name="T", top_k=1,
            score_field="bogus",
        )


# ---------- missing-score policies ----------


def test_evaluate_missing_score_rank_last_keeps_candidates_in_total():
    """Regression: V3 score fields like spacer_final_score are missing on many
    candidates. rank_last (the default) keeps them in n_total/n_positives, never
    above a candidate that has the score. The pre-fix `drop` behavior silently
    shrank the evaluation set."""

    p1 = _scored("p1", 0.9)  # has final_score; spacer is None
    n1 = _scored("n1", 0.5)  # has final_score; spacer is None
    r = evaluate_ranking(
        [p1, n1], {"p1"}, cohort_name="T", top_k=1,
        score_field="spacer_final_score",
    )
    # Both candidates should count even though both lack spacer scores.
    assert r.n_total == 2
    assert r.n_positives == 1
    assert r.n_negatives == 1
    # AUC undefined here because both ranks tie at -inf → 0.5
    assert r.roc_auc == pytest.approx(0.5)


def test_evaluate_missing_score_drop_preserves_v30_behavior():
    """Opt-in to the V3.0 behavior for backward-compat numerics."""
    p1 = _scored("p1", 0.9)
    r = evaluate_ranking(
        [p1], {"p1"}, cohort_name="T", top_k=1,
        score_field="spacer_final_score",
        missing_score_policy="drop",
    )
    assert r.n_total == 0
    assert r.n_positives == 0
    assert r.precision_at_k is None


def test_evaluate_missing_score_error_raises():
    p1 = _scored("p1", 0.9)
    with pytest.raises(ValueError, match="has no 'spacer_final_score'"):
        evaluate_ranking(
            [p1], {"p1"}, cohort_name="T", top_k=1,
            score_field="spacer_final_score",
            missing_score_policy="error",
        )


def test_evaluate_naive_selectivity_score():
    """V3.1 — naive ablation: rank by (β_normal - β_tumor) only."""
    # Build two scored candidates with different methylation patterns.
    def _scored_with_methylation(cid, fs, bt_mean, bn_mean):
        sc = _scored(cid, fs)
        obs = sc.observation.model_copy(update={
            "beta_tumor_mean": bt_mean,
            "beta_normal_mean": bn_mean,
        })
        return sc.model_copy(update={"observation": obs})

    p1 = _scored_with_methylation("p1", 0.5, 0.05, 0.85)  # naive = 0.80
    n1 = _scored_with_methylation("n1", 0.5, 0.50, 0.55)  # naive = 0.05
    r = evaluate_ranking(
        [p1, n1], {"p1"}, cohort_name="T", top_k=1,
        score_field="naive_selectivity",
    )
    assert r.precision_at_k == pytest.approx(1.0)
    assert r.roc_auc == pytest.approx(1.0)


def test_evaluate_naive_selectivity_returns_none_for_unobserved():
    obs = MethylationObservation(
        candidate_id="x", cohort_name="T", evidence_class=EvidenceClass.UNOBSERVED,
    )
    sc = _scored("x", 0.0).model_copy(update={"observation": obs})
    r = evaluate_ranking(
        [sc], {"x"}, cohort_name="T", top_k=1,
        score_field="naive_selectivity",
        missing_score_policy="drop",
    )
    assert r.n_total == 0


def test_evaluate_missing_score_policy_validated():
    with pytest.raises(ValueError, match="missing_score_policy"):
        evaluate_ranking(
            [_scored("p1", 0.5)], {"p1"}, cohort_name="T", top_k=1,
            missing_score_policy="bogus",
        )


# ---------- P1 regression: deterministic tie-break + tie_band_size ----------


def test_evaluate_pak_does_not_depend_on_input_order_when_scores_tie():
    """P1 regression — from the code-review comment:

        with three candidates all scored 1.0, I reproduced P@1=1.0 when
        the positive is first in the input and P@1=0.0 when the same
        positive is last, while AUC correctly stays 0.5.

    After the fix, the secondary sort by candidate_id must make P@1 the
    same regardless of stream order."""
    p = _scored("apos", 1.0)
    n1 = _scored("bneg", 1.0)
    n2 = _scored("cneg", 1.0)

    r_pos_first = evaluate_ranking(
        [p, n1, n2], {"apos"}, cohort_name="T", top_k=1,
    )
    r_pos_last = evaluate_ranking(
        [n1, n2, p], {"apos"}, cohort_name="T", top_k=1,
    )
    r_pos_middle = evaluate_ranking(
        [n1, p, n2], {"apos"}, cohort_name="T", top_k=1,
    )

    assert r_pos_first.precision_at_k == r_pos_last.precision_at_k
    assert r_pos_last.precision_at_k == r_pos_middle.precision_at_k
    # And with the cid tie-break "apos" < "bneg" < "cneg", the positive
    # wins the tied top-1 deterministically.
    assert r_pos_first.precision_at_k == pytest.approx(1.0)


def test_evaluate_reports_tie_band_size_at_k():
    """BenchmarkResult.tie_band_size_at_k must report how many records
    share the K-th position's score. > 1 means P@K is tie-sensitive even
    under a deterministic secondary key."""
    # 5 candidates, all scored 1.0: whole set is one tied band.
    scored = [_scored(cid, 1.0) for cid in ("a", "b", "c", "d", "e")]
    r = evaluate_ranking(scored, {"a"}, cohort_name="T", top_k=3)
    assert r.tie_band_size_at_k == 5
    assert r.tie_break_policy == "candidate_id_asc"

    # All distinct scores → no ties at the cutoff.
    scored = [_scored(cid, float(i)) for i, cid in enumerate("abcde")]
    r = evaluate_ranking(scored, {"a"}, cohort_name="T", top_k=3)
    assert r.tie_band_size_at_k == 1


def test_evaluate_pak_min_max_bounds_tie_band_uncertainty():
    """P@K min/max must bracket the observed P@K and collapse to it when
    there is no tie band at the cutoff. Demonstrates the reporting
    contract: the [min, max] interval is the tie-band-aware uncertainty."""
    # 5 candidates all tied; 2 positives. K=3.
    # band_pos=2, k_band=3, drop_slots=2.
    # min: top-K positives = 0 + max(0, 2-2) = 0 → P@3 = 0.000
    # max: top-K positives = 0 + min(2, 3) = 2 → P@3 = 0.667
    scored = [_scored(cid, 1.0) for cid in ("a", "b", "c", "d", "e")]
    r = evaluate_ranking(scored, {"a", "b"}, cohort_name="T", top_k=3)
    assert r.tie_band_size_at_k == 5
    assert r.precision_at_k_min == pytest.approx(0.0)
    assert r.precision_at_k_max == pytest.approx(2.0 / 3.0)
    # Observed P@K must lie in [min, max].
    assert r.precision_at_k_min <= r.precision_at_k <= r.precision_at_k_max
    # Recall has denominator = n_pos = 2.
    assert r.recall_at_k_min == pytest.approx(0.0)
    assert r.recall_at_k_max == pytest.approx(1.0)


def test_evaluate_pak_min_max_collapse_on_no_tie():
    """With distinct scores, the tie band at the cutoff is size 1,
    so P@K_min == P@K == P@K_max."""
    scored = [_scored(cid, float(i)) for i, cid in enumerate("abcde")]
    r = evaluate_ranking(scored, {"e", "d"}, cohort_name="T", top_k=3)
    assert r.tie_band_size_at_k == 1
    assert r.precision_at_k_min == r.precision_at_k == r.precision_at_k_max
    assert r.recall_at_k_min == r.recall_at_k == r.recall_at_k_max
