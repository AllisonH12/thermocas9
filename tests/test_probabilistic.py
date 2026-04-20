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
