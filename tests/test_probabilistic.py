"""V2 probabilistic scoring tests — CDF math + composition."""

from __future__ import annotations

import pytest

from thermocas.models import EvidenceClass, MethylationObservation
from thermocas.probabilistic import (
    _cdf,
    _fit_beta_method_of_moments,
    _piecewise_linear_cdf,
    p_observation_trustworthy,
    p_protected_normal,
    p_targetable_tumor,
    probabilistic_score,
    regularized_incomplete_beta,
)


def _obs(
    *,
    evidence: EvidenceClass = EvidenceClass.EXACT,
    bt_mean: float | None = 0.05,
    bt_q25: float | None = 0.02,
    bt_q75: float | None = 0.10,
    bn_mean: float | None = 0.85,
    bn_q25: float | None = 0.78,
    bn_q75: float | None = 0.92,
    n_t: int = 400,
    n_n: int = 80,
) -> MethylationObservation:
    if evidence == EvidenceClass.UNOBSERVED:
        return MethylationObservation(
            candidate_id="x", cohort_name="c", evidence_class=evidence,
        )
    return MethylationObservation(
        candidate_id="x", cohort_name="c",
        evidence_class=evidence, evidence_distance_bp=0, probe_id="cg001",
        beta_tumor_mean=bt_mean, beta_tumor_q25=bt_q25, beta_tumor_q75=bt_q75,
        beta_normal_mean=bn_mean, beta_normal_q25=bn_q25, beta_normal_q75=bn_q75,
        n_samples_tumor=n_t, n_samples_normal=n_n,
    )


def test_piecewise_cdf_anchors():
    """CDF must equal 0 at x=0 and 1 at x=1."""
    assert _piecewise_linear_cdf(0.0, 0.1, 0.5, 0.9) == 0.0
    assert _piecewise_linear_cdf(1.0, 0.1, 0.5, 0.9) == 1.0


def test_piecewise_cdf_hits_quantile_knots():
    """F(q25) == 0.25, F(mean) == 0.5, F(q75) == 0.75."""
    cdf = _piecewise_linear_cdf
    assert cdf(0.20, 0.20, 0.50, 0.80) == pytest.approx(0.25)
    assert cdf(0.50, 0.20, 0.50, 0.80) == pytest.approx(0.50)
    assert cdf(0.80, 0.20, 0.50, 0.80) == pytest.approx(0.75)


def test_piecewise_cdf_monotonic_random_grid():
    """For any monotonic knot set, F is non-decreasing."""
    knots = (0.05, 0.30, 0.70)
    xs = [i / 100 for i in range(0, 101)]
    fs = [_piecewise_linear_cdf(x, *knots) for x in xs]
    for a, b in zip(fs, fs[1:], strict=False):
        assert a <= b + 1e-12


def test_piecewise_cdf_mean_only_is_a_sharp_step():
    """Regression for P3: mean-only fallback must be a sharp step at mean,
    not a smoothly interpolated CDF (which would invent distributional shape
    that the data does not support)."""
    cdf = _piecewise_linear_cdf
    assert cdf(0.10, None, 0.50, None) == 0.0
    assert cdf(0.49, None, 0.50, None) == 0.0
    assert cdf(0.50, None, 0.50, None) == pytest.approx(0.5)
    assert cdf(0.51, None, 0.50, None) == 1.0
    assert cdf(0.90, None, 0.50, None) == 1.0


def test_piecewise_cdf_one_sided_quantile_uses_anchors():
    """When only q25 OR q75 is missing, the other knot still pins the curve."""
    # q75 known, q25 missing — F(q75) should still be ~0.75
    assert _piecewise_linear_cdf(0.80, None, 0.50, 0.80) == pytest.approx(0.75)
    # q25 known, q75 missing — F(q25) should still be ~0.25
    assert _piecewise_linear_cdf(0.20, 0.20, 0.50, None) == pytest.approx(0.25)


def test_p_targetable_tumor_high_when_tumor_hypomethylated():
    obs = _obs(bt_mean=0.05, bt_q25=0.02, bt_q75=0.10)
    p = p_targetable_tumor(obs, unmethylated_threshold=0.30)
    assert p > 0.75  # most of the tumor distribution is below 0.30


def test_p_targetable_tumor_low_when_tumor_methylated():
    obs = _obs(bt_mean=0.95, bt_q25=0.92, bt_q75=0.98)
    p = p_targetable_tumor(obs, unmethylated_threshold=0.30)
    # Piecewise-linear with (0,0) anchor over-counts tail mass relative to a
    # true Beta fit, but should still place clearly less than even odds.
    assert p < 0.15


def test_p_targetable_zero_when_unobserved():
    obs = _obs(evidence=EvidenceClass.UNOBSERVED)
    assert p_targetable_tumor(obs) == 0.0


def test_p_protected_normal_high_when_normal_methylated():
    obs = _obs(bn_mean=0.90, bn_q25=0.85, bn_q75=0.95)
    p = p_protected_normal(obs, methylated_threshold=0.50)
    assert p > 0.80


def test_p_protected_normal_low_when_normal_hypomethylated():
    obs = _obs(bn_mean=0.10, bn_q25=0.05, bn_q75=0.15)
    p = p_protected_normal(obs, methylated_threshold=0.50)
    assert p < 0.20


def test_p_protected_normal_strictly_greater_for_more_methylated():
    """Order matters: a more-methylated normal must always score as more protected."""
    high = p_protected_normal(_obs(bn_mean=0.90, bn_q25=0.85, bn_q75=0.95))
    low = p_protected_normal(_obs(bn_mean=0.10, bn_q25=0.05, bn_q75=0.15))
    assert high > low + 0.5  # comfortable margin


def test_p_targetable_strictly_greater_for_less_methylated():
    """Order matters: a less-methylated tumor must always score as more targetable."""
    low = p_targetable_tumor(_obs(bt_mean=0.95, bt_q25=0.92, bt_q75=0.98))
    high = p_targetable_tumor(_obs(bt_mean=0.05, bt_q25=0.02, bt_q75=0.10))
    assert high > low + 0.5


# ---------- V3: Beta-distribution CDF ----------


def test_regularized_incomplete_beta_boundaries():
    """Standard boundary identities for I_x(a, b)."""
    assert regularized_incomplete_beta(0.0, 2.0, 5.0) == 0.0
    assert regularized_incomplete_beta(1.0, 2.0, 5.0) == 1.0
    # Symmetry: Beta(a, a) is symmetric around 0.5.
    assert regularized_incomplete_beta(0.5, 3.0, 3.0) == pytest.approx(0.5, abs=1e-9)
    # Beta(1, 1) is uniform → I_x(1, 1) = x.
    for x in (0.1, 0.3, 0.7, 0.9):
        assert regularized_incomplete_beta(x, 1.0, 1.0) == pytest.approx(x, abs=1e-12)


def test_regularized_incomplete_beta_known_values():
    """Spot-check against known Beta CDF values."""
    # I_{0.5}(2, 2) = 0.5 (mean of Beta(2,2) is 0.5)
    assert regularized_incomplete_beta(0.5, 2.0, 2.0) == pytest.approx(0.5, abs=1e-9)
    # I_{0.5}(2, 1) = 0.25 (B(0.5; 2,1)/B(2,1) = (0.5²/2) / (1/2) = 0.25)
    assert regularized_incomplete_beta(0.5, 2.0, 1.0) == pytest.approx(0.25, abs=1e-9)
    # I_{0.5}(1, 2) = 0.75 by symmetry
    assert regularized_incomplete_beta(0.5, 1.0, 2.0) == pytest.approx(0.75, abs=1e-9)


def test_regularized_incomplete_beta_rejects_bad_params():
    with pytest.raises(ValueError, match="must be positive"):
        regularized_incomplete_beta(0.5, 0.0, 1.0)
    with pytest.raises(ValueError, match="must be positive"):
        regularized_incomplete_beta(0.5, 1.0, -1.0)


def test_fit_beta_method_of_moments_recovers_parameters():
    """Round-trip: pick (α, β) → derive μ, σ → re-fit → get back (α, β)."""
    alpha_true, beta_true = 5.0, 20.0
    mean = alpha_true / (alpha_true + beta_true)
    var = (alpha_true * beta_true) / (
        (alpha_true + beta_true) ** 2 * (alpha_true + beta_true + 1)
    )
    sigma = var ** 0.5
    fit = _fit_beta_method_of_moments(mean, sigma)
    assert fit is not None
    alpha_fit, beta_fit = fit
    assert alpha_fit == pytest.approx(alpha_true, rel=1e-9)
    assert beta_fit == pytest.approx(beta_true, rel=1e-9)


def test_fit_beta_returns_none_when_variance_too_large():
    """If σ² ≥ μ(1−μ) no Beta exists with those moments → None (PWL fallback)."""
    assert _fit_beta_method_of_moments(0.5, 0.6) is None  # σ² = 0.36 ≥ 0.25
    assert _fit_beta_method_of_moments(0.0, 0.1) is None  # mean at boundary
    assert _fit_beta_method_of_moments(1.0, 0.1) is None
    assert _fit_beta_method_of_moments(0.5, 0.0) is None  # zero stdev


def test_cdf_dispatch_uses_beta_when_quantiles_given():
    """The V3 dispatch picks Beta over PWL when quantiles let it fit."""
    # Heavily methylated: q25=0.92, μ=0.95, q75=0.98 → Beta gives ~1e-12, PWL gives 0.082.
    pwl = _piecewise_linear_cdf(0.30, 0.92, 0.95, 0.98)
    beta = _cdf(0.30, 0.92, 0.95, 0.98)
    assert pwl > 0.05  # historical PWL behavior
    assert beta < 1e-6  # Beta correctly puts ~0 mass below 0.30


def test_cdf_dispatch_falls_back_to_pwl_when_one_sided():
    """One-sided quantiles can't yield an IQR → no Beta fit → PWL is used."""
    only_q25 = _cdf(0.20, 0.20, 0.50, None)
    pwl_only = _piecewise_linear_cdf(0.20, 0.20, 0.50, None)
    assert only_q25 == pytest.approx(pwl_only)


def test_cdf_dispatch_falls_back_to_pwl_when_quantiles_collapse():
    """q25 == q75 means IQR = 0, Beta variance undefined → PWL."""
    collapsed = _cdf(0.40, 0.50, 0.50, 0.50)
    pwl = _piecewise_linear_cdf(0.40, 0.50, 0.50, 0.50)
    assert collapsed == pytest.approx(pwl)


def test_p_targetable_under_beta_is_essentially_zero_for_methylated_tumor():
    """Regression: V3 must place ≈0 mass on `targetable` when tumor is heavily
    methylated. V2 piecewise-linear gave ≈0.082 here, which was the issue."""
    obs = _obs(bt_mean=0.95, bt_q25=0.92, bt_q75=0.98)
    p = p_targetable_tumor(obs, unmethylated_threshold=0.30)
    assert p < 1e-6


def test_p_protected_under_beta_is_essentially_zero_for_hypomethylated_normal():
    """Symmetric V3 regression on the normal side."""
    obs = _obs(bn_mean=0.10, bn_q25=0.05, bn_q75=0.15)
    p = p_protected_normal(obs, methylated_threshold=0.50)
    assert p < 1e-3


def test_beta_cdf_monotonic_in_x():
    """For any fixed (α, β), I_x(a, b) is non-decreasing in x."""
    a, b = 4.0, 6.0
    xs = [i / 100 for i in range(0, 101)]
    fs = [regularized_incomplete_beta(x, a, b) for x in xs]
    for u, v in zip(fs, fs[1:], strict=False):
        assert u <= v + 1e-12


def test_p_observation_trustworthy_zero_for_unobserved():
    obs = _obs(evidence=EvidenceClass.UNOBSERVED)
    assert p_observation_trustworthy(obs) == 0.0


def test_p_observation_trustworthy_scales_with_samples():
    """Trust below the ramp threshold scales linearly with min(n_t, n_n)."""
    low = _obs(evidence=EvidenceClass.EXACT, n_t=5, n_n=5)
    high = _obs(evidence=EvidenceClass.EXACT, n_t=400, n_n=80)
    assert p_observation_trustworthy(low, ramp_n=30) < p_observation_trustworthy(high, ramp_n=30)


def test_p_observation_trustworthy_strict_ordering_by_evidence():
    """EXACT > PROXIMAL_CLOSE > PROXIMAL > REGIONAL > UNOBSERVED at saturation."""
    common = dict(n_t=400, n_n=80)
    weights = [
        p_observation_trustworthy(_obs(evidence=EvidenceClass.EXACT, **common)),
        p_observation_trustworthy(_obs(evidence=EvidenceClass.PROXIMAL_CLOSE, **common)),
        p_observation_trustworthy(_obs(evidence=EvidenceClass.PROXIMAL, **common)),
        p_observation_trustworthy(_obs(evidence=EvidenceClass.REGIONAL, **common)),
        p_observation_trustworthy(_obs(evidence=EvidenceClass.UNOBSERVED)),
    ]
    assert weights == sorted(weights, reverse=True)


def test_probabilistic_score_default_mode_is_tumor_only():
    """V2.4 — default mode is `tumor_only`; p_protected_normal is emitted for
    auditability but does NOT participate in p_therapeutic_selectivity.
    p_sel = p_targ × p_trust."""
    obs = _obs()
    ps = probabilistic_score(obs)
    assert ps.mode == "tumor_only"
    assert 0.0 < ps.p_targetable_tumor <= 1.0
    assert 0.0 < ps.p_protected_normal <= 1.0   # still emitted
    assert 0.0 < ps.p_observation_trustworthy <= 1.0
    expected = ps.p_targetable_tumor * ps.p_observation_trustworthy
    assert ps.p_therapeutic_selectivity == pytest.approx(expected)


def test_probabilistic_score_tumor_plus_normal_protection_mode():
    """Opt-in: p_sel = p_targ × p_prot × p_trust."""
    obs = _obs()
    ps = probabilistic_score(obs, mode="tumor_plus_normal_protection")
    assert ps.mode == "tumor_plus_normal_protection"
    expected = (
        ps.p_targetable_tumor
        * ps.p_protected_normal
        * ps.p_observation_trustworthy
    )
    assert ps.p_therapeutic_selectivity == pytest.approx(expected)


def test_probabilistic_score_mode_differs():
    """The two modes must produce different composites when p_prot ≠ 1.0."""
    obs = _obs(bn_mean=0.55, bn_q25=0.50, bn_q75=0.60)  # p_prot ~ moderate
    t = probabilistic_score(obs, mode="tumor_only")
    tp = probabilistic_score(obs, mode="tumor_plus_normal_protection")
    # tumor_plus scales down by p_prot, so it must be ≤ tumor_only
    assert tp.p_therapeutic_selectivity <= t.p_therapeutic_selectivity


def test_probabilistic_score_unknown_mode_rejected():
    from thermocas.probabilistic import probabilistic_score as ps_fn
    with pytest.raises(ValueError, match="unknown probabilistic_mode"):
        ps_fn(_obs(), mode="bogus")  # type: ignore[arg-type]


def test_probabilistic_score_zero_for_unobserved():
    ps = probabilistic_score(_obs(evidence=EvidenceClass.UNOBSERVED))
    assert ps.p_therapeutic_selectivity == 0.0


# ---------- V2.5 — differential-protection mode contract ----------


def _selective_obs() -> MethylationObservation:
    """Observation with a clear β_n − β_t gap — normal-mean 0.85, tumor-mean 0.05."""
    return _obs(
        bt_mean=0.05, bt_q25=0.02, bt_q75=0.10,
        bn_mean=0.85, bn_q25=0.78, bn_q75=0.92,
    )


def test_p_differential_protection_increases_with_gap():
    """Observation with a larger normal-vs-tumor gap must have higher p_diff."""
    from thermocas.probabilistic import p_differential_protection

    small_gap = _obs(bt_mean=0.40, bt_q25=0.35, bt_q75=0.45,
                     bn_mean=0.50, bn_q25=0.45, bn_q75=0.55)
    big_gap = _obs(bt_mean=0.05, bt_q25=0.02, bt_q75=0.08,
                   bn_mean=0.95, bn_q25=0.92, bn_q75=0.98)

    assert p_differential_protection(small_gap, delta=0.2) < p_differential_protection(big_gap, delta=0.2)


def test_p_differential_protection_half_at_breakpoint():
    """Under the normal approximation, when the observed gap equals δ exactly,
    `P(Δ > δ)` must equal 0.5 — the defining "no-information" point of the
    factor. Verifies the math, not a benchmark number."""
    from thermocas.probabilistic import p_differential_protection

    # symmetric spreads, gap = 0.20 = δ → z = 0 → p_diff = 0.5
    obs = _obs(bt_mean=0.30, bt_q25=0.25, bt_q75=0.35,
               bn_mean=0.50, bn_q25=0.45, bn_q75=0.55)
    assert p_differential_protection(obs, delta=0.2) == pytest.approx(0.5, abs=1e-12)


def test_p_differential_protection_unobserved_is_zero():
    """UNOBSERVED observations carry no β values — p_diff must be exactly 0,
    not some artifact of the normal approximation."""
    from thermocas.probabilistic import p_differential_protection
    obs = _obs(evidence=EvidenceClass.UNOBSERVED)
    assert p_differential_protection(obs, delta=0.2) == 0.0


def test_p_differential_protection_zero_when_either_side_missing():
    """No fabrication from one-sided data — return 0 if tumor or normal mean is None."""
    from thermocas.probabilistic import p_differential_protection
    obs_no_normal = MethylationObservation(
        candidate_id="x", cohort_name="c",
        evidence_class=EvidenceClass.EXACT, evidence_distance_bp=0, probe_id="cg001",
        beta_tumor_mean=0.1, beta_tumor_q25=0.08, beta_tumor_q75=0.12,
        n_samples_tumor=100, n_samples_normal=0,
    )
    assert p_differential_protection(obs_no_normal, delta=0.2) == 0.0


def test_probabilistic_score_differential_mode_differs_from_tumor_only():
    """V2.5 mode must produce a different composite than tumor_only when the
    normal-tumor gap is near the δ margin (so p_diff ≠ 1) — confirms the new
    factor is actually on the score path, not dropped."""
    # gap ≈ 0.10, δ = 0.2  → p_diff ≪ 1, composite must shrink vs tumor_only
    obs = _obs(bt_mean=0.40, bt_q25=0.35, bt_q75=0.45,
               bn_mean=0.50, bn_q25=0.45, bn_q75=0.55)
    t = probabilistic_score(obs, mode="tumor_only")
    d = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                            differential_delta=0.2)
    assert d.mode == "tumor_plus_differential_protection"
    # differential composite must be strictly less than tumor_only here
    assert d.p_therapeutic_selectivity < t.p_therapeutic_selectivity
    assert d.p_differential_protection is not None
    assert d.differential_delta is not None
    assert 0.0 < d.p_differential_protection < 1.0


def test_probabilistic_score_differential_records_delta():
    """The configured δ must appear on the emitted record so reviewers can
    tell which margin was in effect without consulting cohort config."""
    obs = _selective_obs()
    d = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                            differential_delta=0.35)
    assert d.differential_delta == pytest.approx(0.35)


def test_probabilistic_score_differential_higher_delta_is_stricter():
    """Higher δ = stricter criterion = smaller p_diff = smaller composite
    (p_targ and p_trust unchanged across δ)."""
    obs = _selective_obs()
    d_small = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                                  differential_delta=0.1)
    d_big = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                                differential_delta=0.6)
    assert d_big.p_differential_protection < d_small.p_differential_protection
    assert d_big.p_therapeutic_selectivity < d_small.p_therapeutic_selectivity


def test_probabilistic_score_differential_composite_formula():
    """V2.5 composite must be exactly p_targ × p_diff × p_trust — no silent
    extra factors."""
    obs = _selective_obs()
    d = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                            differential_delta=0.2)
    expected = (
        d.p_targetable_tumor
        * d.p_differential_protection
        * d.p_observation_trustworthy
    )
    assert d.p_therapeutic_selectivity == pytest.approx(expected)


def test_probabilistic_score_non_differential_leaves_diff_fields_none():
    """tumor_only / tumor_plus_normal_protection must NOT populate the
    differential audit fields — keeps record interpretation unambiguous."""
    obs = _selective_obs()
    for mode in ("tumor_only", "tumor_plus_normal_protection"):
        ps = probabilistic_score(obs, mode=mode)  # type: ignore[arg-type]
        assert ps.p_differential_protection is None
        assert ps.differential_delta is None


def test_probabilistic_score_rejects_composite_mismatch():
    """P2 regression: a ProbabilisticScore whose stored `p_therapeutic_selectivity`
    does not equal the factor product implied by `mode` must fail validation.
    Downstream ranking reads the stored scalar, so a silent mismatch would
    change rankings while still passing Pydantic.

    Reproduces the concrete failure case from the P2 review comment."""
    from pydantic import ValidationError

    from thermocas.models import ProbabilisticScore

    # mode=differential, factors imply 0.8*0.9*0.5 = 0.36 but stored = 0.01
    with pytest.raises(ValidationError, match="does not match the composite"):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_plus_differential_protection",
            p_targetable_tumor=0.8, p_protected_normal=0.5,
            p_observation_trustworthy=0.5,
            p_differential_protection=0.9, differential_delta=0.2,
            p_therapeutic_selectivity=0.01,
        )
    # Same for tumor_only
    with pytest.raises(ValidationError, match="does not match the composite"):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_only",
            p_targetable_tumor=0.4, p_protected_normal=0.0,
            p_observation_trustworthy=0.5,
            p_therapeutic_selectivity=0.99,  # expected 0.2
        )
    # And for tumor_plus_normal_protection
    with pytest.raises(ValidationError, match="does not match the composite"):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_plus_normal_protection",
            p_targetable_tumor=0.6, p_protected_normal=0.6,
            p_observation_trustworthy=0.5,
            p_therapeutic_selectivity=0.10,  # expected 0.18
        )


def test_probabilistic_score_accepts_consistent_composite_under_all_modes():
    """Positive-path test for the composite validator: the emitted
    ProbabilisticScore from `probabilistic_score` must always validate, under
    every mode. Uses the same observation so the only thing that changes is
    the composite arithmetic implied by mode."""
    obs = _obs()
    for mode in ("tumor_only", "tumor_plus_normal_protection",
                 "tumor_plus_differential_protection"):
        ps = probabilistic_score(obs, mode=mode)  # type: ignore[arg-type]
        # Round-trip through JSON to exercise deserialization too.
        from thermocas.models import ProbabilisticScore
        ProbabilisticScore.model_validate_json(ps.model_dump_json())


def test_probabilistic_score_rejects_inconsistent_diff_fields_at_validation():
    """A record claiming mode=differential without the diff audit fields, or
    vice versa, must fail validation — prevents malformed records from
    round-tripping through JSONL."""
    from pydantic import ValidationError

    from thermocas.models import ProbabilisticScore

    with pytest.raises(ValidationError):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_plus_differential_protection",
            p_targetable_tumor=0.5, p_protected_normal=0.5,
            p_observation_trustworthy=0.5,
            p_therapeutic_selectivity=0.125,
            # missing p_differential_protection and differential_delta
        )

    with pytest.raises(ValidationError):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_only",
            p_targetable_tumor=0.5, p_protected_normal=0.5,
            p_observation_trustworthy=0.5,
            p_therapeutic_selectivity=0.25,
            p_differential_protection=0.8,  # not allowed in tumor_only
            differential_delta=0.2,
        )


# ---------- tumor_plus_gap_sigmoid (shipping next tissue mode) ----------


def test_p_gap_sigmoid_monotone_in_gap():
    """Larger (β_normal − β_tumor) at fixed δ, σ_fixed → larger p_gap_sigmoid."""
    from thermocas.probabilistic import p_gap_sigmoid

    def obs_with_gap(gap):
        return _obs(bt_mean=0.1, bt_q25=0.08, bt_q75=0.12,
                    bn_mean=0.1 + gap, bn_q25=0.1 + gap - 0.02, bn_q75=0.1 + gap + 0.02)

    p_small = p_gap_sigmoid(obs_with_gap(0.1), delta=0.2, sigma_fixed=0.07)
    p_at = p_gap_sigmoid(obs_with_gap(0.2), delta=0.2, sigma_fixed=0.07)
    p_big = p_gap_sigmoid(obs_with_gap(0.5), delta=0.2, sigma_fixed=0.07)
    assert p_small < p_at < p_big
    # At gap == δ the sigmoid is at 0.5 exactly.
    assert p_at == pytest.approx(0.5)


def test_p_gap_sigmoid_missing_beta_returns_zero():
    """One-sided β → 0.0, same as p_differential_protection. No fabrication."""
    from thermocas.probabilistic import p_gap_sigmoid

    obs_no_normal = _obs(bt_mean=0.1, bt_q25=0.08, bt_q75=0.12,
                         bn_mean=None, bn_q25=None, bn_q75=None)
    assert p_gap_sigmoid(obs_no_normal, delta=0.2, sigma_fixed=0.07) == 0.0


def test_probabilistic_score_gap_sigmoid_mode_populates_fields():
    """tumor_plus_gap_sigmoid must populate p_gap_sigmoid + sigma_fixed +
    differential_delta; must NOT populate p_differential_protection."""
    obs = _selective_obs()
    r = probabilistic_score(obs, mode="tumor_plus_gap_sigmoid",
                            differential_delta=0.1, sigma_fixed=0.07)
    assert r.mode == "tumor_plus_gap_sigmoid"
    assert r.p_gap_sigmoid is not None
    assert 0.0 < r.p_gap_sigmoid <= 1.0
    assert r.sigma_fixed == pytest.approx(0.07)
    assert r.differential_delta == pytest.approx(0.1)
    assert r.p_differential_protection is None


def test_probabilistic_score_gap_sigmoid_composite_formula():
    """p_therapeutic_selectivity must equal p_targ × p_gap_sigmoid × p_trust."""
    obs = _selective_obs()
    r = probabilistic_score(obs, mode="tumor_plus_gap_sigmoid",
                            differential_delta=0.2, sigma_fixed=0.07)
    expected = r.p_targetable_tumor * r.p_gap_sigmoid * r.p_observation_trustworthy
    assert r.p_therapeutic_selectivity == pytest.approx(expected)


def test_probabilistic_score_gap_sigmoid_differs_from_differential():
    """gap_sigmoid with σ_fixed near σ_floor should produce a different
    composite than V2.5 differential on records where per-site σ_Δ differs
    from σ_fixed — confirms the new path is actually on the score, not a
    silent alias."""
    # Wide IQRs on both sides → V2.5 per-site σ_Δ ≈ 0.11 (well above σ_floor);
    # gap_sigmoid at σ_fixed = 0.05 sees a tighter bandwidth.
    obs = _obs(bt_mean=0.1, bt_q25=0.05, bt_q75=0.20,
               bn_mean=0.35, bn_q25=0.30, bn_q75=0.45)
    d = probabilistic_score(obs, mode="tumor_plus_differential_protection",
                            differential_delta=0.2)
    g = probabilistic_score(obs, mode="tumor_plus_gap_sigmoid",
                            differential_delta=0.2, sigma_fixed=0.05)
    assert d.p_therapeutic_selectivity != g.p_therapeutic_selectivity


def test_gap_sigmoid_score_rejects_differential_field():
    """A ProbabilisticScore record with mode=gap_sigmoid must not carry
    p_differential_protection — iff validator on the model."""
    from pydantic import ValidationError
    from thermocas.models import ProbabilisticScore

    with pytest.raises(ValidationError):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_plus_gap_sigmoid",
            p_targetable_tumor=0.8, p_protected_normal=0.0,
            p_observation_trustworthy=0.5,
            p_gap_sigmoid=0.7, sigma_fixed=0.07, differential_delta=0.2,
            p_differential_protection=0.6,  # not allowed in gap_sigmoid mode
            p_therapeutic_selectivity=0.8 * 0.7 * 0.5,
        )


def test_differential_score_rejects_gap_sigmoid_field():
    """And the reverse: a mode=differential record must not carry p_gap_sigmoid."""
    from pydantic import ValidationError
    from thermocas.models import ProbabilisticScore

    with pytest.raises(ValidationError):
        ProbabilisticScore(
            candidate_id="x", cohort_name="c",
            mode="tumor_plus_differential_protection",
            p_targetable_tumor=0.8, p_protected_normal=0.0,
            p_observation_trustworthy=0.5,
            p_differential_protection=0.9, differential_delta=0.2,
            p_gap_sigmoid=0.7,  # not allowed in differential mode
            p_therapeutic_selectivity=0.8 * 0.9 * 0.5,
        )


def test_p_gap_sigmoid_rejects_zero_sigma():
    """sigma_fixed = 0 would divide by zero — public function must reject it."""
    from thermocas.probabilistic import p_gap_sigmoid
    obs = _selective_obs()
    with pytest.raises(ValueError, match="sigma_fixed must be strictly positive"):
        p_gap_sigmoid(obs, delta=0.2, sigma_fixed=0.0)
    with pytest.raises(ValueError, match="sigma_fixed must be strictly positive"):
        p_gap_sigmoid(obs, delta=0.2, sigma_fixed=-0.1)


def test_cohort_config_rejects_zero_sigma_fixed():
    """CohortConfig.sigma_fixed has gt=0.0 — sigma_fixed: 0 in YAML must error."""
    from pydantic import ValidationError
    from thermocas.models import CohortConfig

    base = dict(
        name="x", tumor_dataset="t", normal_dataset="n", platform="HM450",
        min_samples_tumor=2, min_samples_normal=2,
        probabilistic_mode="tumor_plus_gap_sigmoid",
    )
    with pytest.raises(ValidationError):
        CohortConfig(**base, sigma_fixed=0.0)
    cfg = CohortConfig(**base, sigma_fixed=0.05)
    assert cfg.sigma_fixed == 0.05
