"""Probabilistic targetability scoring.

Decomposes site-level targetability into three independent factors so reviewers
can audit the inputs to a clinical decision:

    P(therapeutic_selectivity)
        = P(targetable_in_tumor)
        × P(protected_in_normal)
        × P(observation_trustworthy)

Each factor is estimated from the per-probe `MethylationObservation` summary
(mean + q25 + q75 + n_samples + evidence_class).

CDF strategy (V3):

  1. **Method-of-moments Beta(α, β) fit** from the observed mean and an
     IQR-derived stdev (σ ≈ IQR / 1.349). Beta is the natural distribution
     for β-values bounded in [0, 1] and concentrates mass appropriately:
     a probe with q25=0.92, mean=0.95, q75=0.98 yields P(β<0.30) ≈ 1e-4 under
     Beta, vs ≈ 0.082 under the V2 piecewise-linear estimator.
     `I_x(α, β)` (regularized incomplete beta) is computed via Lentz's
     continued fraction — pure stdlib, no scipy.

  2. **Piecewise-linear fallback** for degenerate cases: missing quantiles,
     IQR == 0 (Beta variance undefined), or σ² ≥ μ(1−μ) (no Beta exists with
     that mean and variance). The V2 PWL estimator is preserved as the
     fallback path so behavior is well-defined for every observation.

  3. **Mean-only sharp step** when neither quantile is present — refuses to
     invent distributional shape from a single point estimate.
"""

from __future__ import annotations

import math

from thermocas.models import (
    EvidenceClass,
    MethylationObservation,
    ProbabilisticMode,
    ProbabilisticScore,
)

#: IQR → stdev factor for an approximately normal distribution. Used by the
#: method-of-moments Beta fit. Real β-value distributions are rarely normal,
#: but for the narrow IQRs typical of methylation arrays the approximation is
#: close enough to give Beta a reasonable second moment.
_IQR_TO_STDEV = 1.349

#: Default thresholds — tunable per call.
DEFAULT_UNMETHYLATED_THRESHOLD = 0.30
DEFAULT_METHYLATED_THRESHOLD = 0.50

#: V2.5 default differential margin (P(β_n − β_t > δ)). Matches the sweet
#: spot from the offline sweep on the MCF-7/MCF-10A surrogate.
DEFAULT_DIFFERENTIAL_DELTA = 0.2

#: V2.5 floor on per-side σ when deriving σ_Δ from IQRs. Prevents σ_Δ from
#: collapsing to zero at boundary β-values (where IQR==0 is common) and keeps
#: the normal-approximation well-defined for every observation.
DEFAULT_DIFFERENTIAL_SIGMA_FLOOR = 0.05

#: `tumor_plus_gap_sigmoid` default bandwidth. Chosen as √2 × σ_floor = 0.0707,
#: matching the σ_Δ V2.5 sees when its σ_floor binds on both sides — the modal
#: case on low-`n` matched cell-line cohorts per PAPER.md §3.5 binding-rate
#: diagnostic. Bandwidth-robust on tissue across {0.05, 0.0707, 0.10}
#: (PAPER.md §5.2.1 sweep).
import math as _math
DEFAULT_GAP_SIGMOID_SIGMA_FIXED = _math.sqrt(2) * DEFAULT_DIFFERENTIAL_SIGMA_FLOOR

#: Default per-EvidenceClass trust ceiling. With enough samples,
#: P(trustworthy) saturates at this value.
_TRUST_BASE: dict[EvidenceClass, float] = {
    EvidenceClass.EXACT: 0.95,
    EvidenceClass.PROXIMAL_CLOSE: 0.75,
    EvidenceClass.PROXIMAL: 0.45,
    EvidenceClass.REGIONAL: 0.15,
    EvidenceClass.UNOBSERVED: 0.0,
}

#: Sample count above which trust no longer scales (saturation point).
DEFAULT_TRUST_RAMP_N = 30


# ---------- public API ----------


def p_targetable_tumor(
    obs: MethylationObservation,
    unmethylated_threshold: float = DEFAULT_UNMETHYLATED_THRESHOLD,
) -> float:
    """P(β_tumor < unmethylated_threshold) — the cancer cell's PAM cytosine
    is unmethylated, so ThermoCas9 can bind."""

    if obs.beta_tumor_mean is None:
        return 0.0
    return _cdf(
        unmethylated_threshold,
        obs.beta_tumor_q25,
        obs.beta_tumor_mean,
        obs.beta_tumor_q75,
    )


def p_protected_normal(
    obs: MethylationObservation,
    methylated_threshold: float = DEFAULT_METHYLATED_THRESHOLD,
) -> float:
    """P(β_normal ≥ methylated_threshold) — the normal cell's PAM cytosine
    is methylated, so ThermoCas9 cannot bind. This is the "protection" factor."""

    if obs.beta_normal_mean is None:
        return 0.0
    return 1.0 - _cdf(
        methylated_threshold,
        obs.beta_normal_q25,
        obs.beta_normal_mean,
        obs.beta_normal_q75,
    )


def p_differential_protection(
    obs: MethylationObservation,
    delta: float = DEFAULT_DIFFERENTIAL_DELTA,
    sigma_floor: float = DEFAULT_DIFFERENTIAL_SIGMA_FLOOR,
) -> float:
    """V2.5 — `P(β_normal − β_tumor > δ)` under an independent-normal approximation.

    Replaces the threshold-based `p_protected_normal` (which encodes the
    *biologically-specialized* assumption that normal tissue is methylated
    above 0.5) with a differential factor that makes no static-threshold
    claim about the normal side — it asks only whether the normal-vs-tumor
    gap exceeds δ.

    Model:
        σ_k ≈ IQR_k / 1.349 (normal approximation to the β-value spread)
        σ_Δ = sqrt(max(σ_t, floor)² + max(σ_n, floor)²)
        P(Δ > δ) = 1 − Φ((δ − (μ_n − μ_t)) / σ_Δ)

    Returns 0 when either side lacks a mean — no fabrication of signal from
    one-sided data.
    """

    mu_t = obs.beta_tumor_mean
    mu_n = obs.beta_normal_mean
    if mu_t is None or mu_n is None:
        return 0.0

    q25_t = obs.beta_tumor_q25
    q75_t = obs.beta_tumor_q75
    q25_n = obs.beta_normal_q25
    q75_n = obs.beta_normal_q75

    sigma_t = (q75_t - q25_t) / _IQR_TO_STDEV if (q25_t is not None and q75_t is not None) else 0.0
    sigma_n = (q75_n - q25_n) / _IQR_TO_STDEV if (q25_n is not None and q75_n is not None) else 0.0
    floor = max(0.0, sigma_floor)
    sigma_sq = max(sigma_t, floor) ** 2 + max(sigma_n, floor) ** 2
    z = (delta - (mu_n - mu_t)) / math.sqrt(sigma_sq)
    return 0.5 * (1.0 - math.erf(z / math.sqrt(2.0)))


def p_gap_sigmoid(
    obs: MethylationObservation,
    delta: float = DEFAULT_DIFFERENTIAL_DELTA,
    sigma_fixed: float = DEFAULT_GAP_SIGMOID_SIGMA_FIXED,
) -> float:
    """Fixed-bandwidth sigmoid gap factor: `sigmoid((Δβ − δ) / σ_fixed)`.

    Alternative to `p_differential_protection` whose per-site σ_Δ collapses
    to σ_floor on essentially every record in low-`n` matched cell-line
    cohorts (§3.5 binding-rate diagnostic). On tissue cohorts, where
    per-site σ_Δ actually varies, `p_diff`'s per-site σ adaptation
    underperforms a fixed-bandwidth response (PAPER.md §5.2.1 ablation +
    §5.2.2 genome-wide gating).

    Returns 0 when either side lacks a β mean — no fabrication from
    one-sided data.
    """

    if sigma_fixed <= 0.0:
        raise ValueError(
            f"p_gap_sigmoid: sigma_fixed must be strictly positive, got {sigma_fixed!r}. "
            "Use a small positive value (e.g. 1e-6) for an effectively-hard threshold."
        )
    mu_t = obs.beta_tumor_mean
    mu_n = obs.beta_normal_mean
    if mu_t is None or mu_n is None:
        return 0.0
    x = ((mu_n - mu_t) - delta) / sigma_fixed
    # Numerically-stable logistic.
    if x >= 0:
        e = math.exp(-x)
        return 1.0 / (1.0 + e)
    e = math.exp(x)
    return e / (1.0 + e)


def p_observation_trustworthy(
    obs: MethylationObservation,
    ramp_n: int = DEFAULT_TRUST_RAMP_N,
) -> float:
    """Trust scales with `EvidenceClass` and saturates with sample count.

    The sample-count factor is the smaller of the two cohorts' n divided by the
    ramp threshold (clamped to 1.0). UNOBSERVED → 0.0 regardless of n.
    """

    base = _TRUST_BASE[obs.evidence_class]
    if base == 0.0:
        return 0.0
    sample_floor = min(obs.n_samples_tumor, obs.n_samples_normal)
    if ramp_n <= 0:
        return base
    return base * min(1.0, sample_floor / ramp_n)


def probabilistic_score(
    obs: MethylationObservation,
    *,
    mode: ProbabilisticMode = "tumor_only",
    unmethylated_threshold: float = DEFAULT_UNMETHYLATED_THRESHOLD,
    methylated_threshold: float = DEFAULT_METHYLATED_THRESHOLD,
    trust_ramp_n: int = DEFAULT_TRUST_RAMP_N,
    differential_delta: float = DEFAULT_DIFFERENTIAL_DELTA,
    sigma_fixed: float = DEFAULT_GAP_SIGMOID_SIGMA_FIXED,
) -> ProbabilisticScore:
    """Compose the three (or four) factors into a `ProbabilisticScore`.

    Pure function — no IO, deterministic given inputs. Suitable for direct
    unit testing and vectorized application across genome-wide candidates.

    Args:
        mode: which factors participate in the composite. See
            `ProbabilisticScore.mode` for supported policies. Default
            `tumor_only` — safer because p_protected_normal inverts on
            cohorts where the normal comparator doesn't methylate target
            promoters.
        differential_delta: margin for `tumor_plus_differential_protection`
            and `tumor_plus_gap_sigmoid`. Ignored by the other modes.
        sigma_fixed: fixed-bandwidth σ for `tumor_plus_gap_sigmoid`.
            Ignored by the other modes.
    """

    p_t = p_targetable_tumor(obs, unmethylated_threshold)
    p_p = p_protected_normal(obs, methylated_threshold)
    p_r = p_observation_trustworthy(obs, trust_ramp_n)
    p_diff: float | None = None
    p_gap: float | None = None
    delta_recorded: float | None = None
    sigma_recorded: float | None = None

    if mode == "tumor_only":
        p_sel = p_t * p_r
    elif mode == "tumor_plus_normal_protection":
        p_sel = p_t * p_p * p_r
    elif mode == "tumor_plus_differential_protection":
        p_diff = p_differential_protection(obs, differential_delta)
        delta_recorded = differential_delta
        p_sel = p_t * p_diff * p_r
    elif mode == "tumor_plus_gap_sigmoid":
        p_gap = p_gap_sigmoid(obs, delta=differential_delta, sigma_fixed=sigma_fixed)
        delta_recorded = differential_delta
        sigma_recorded = sigma_fixed
        p_sel = p_t * p_gap * p_r
    else:
        raise ValueError(f"unknown probabilistic_mode: {mode!r}")

    return ProbabilisticScore(
        candidate_id=obs.candidate_id,
        cohort_name=obs.cohort_name,
        mode=mode,
        p_targetable_tumor=p_t,
        p_protected_normal=p_p,
        p_observation_trustworthy=p_r,
        p_differential_protection=p_diff,
        differential_delta=delta_recorded,
        p_gap_sigmoid=p_gap,
        sigma_fixed=sigma_recorded,
        p_therapeutic_selectivity=p_sel,
    )


# ---------- CDF dispatch ----------


def _cdf(
    x: float,
    q25: float | None,
    mean: float | None,
    q75: float | None,
) -> float:
    """V3 CDF dispatch: try Beta(α, β) first, fall back to PWL on degenerate fits.

    Strategy:
      * mean missing → 0
      * mean only → sharp step (delegated to `_piecewise_linear_cdf`)
      * mean + ≥1 quantile + IQR > 0 → method-of-moments Beta fit; if Beta is
        well-defined (σ² < μ(1−μ)), use the regularized incomplete beta
      * otherwise → V2 piecewise-linear estimator (also handles q25==q75)
    """

    if mean is None:
        return 0.0
    if q25 is None and q75 is None:
        return _piecewise_linear_cdf(x, q25, mean, q75)

    # Need both quantiles to derive an IQR-based stdev.
    if q25 is not None and q75 is not None and q75 > q25:
        iqr = q75 - q25
        sigma = iqr / _IQR_TO_STDEV
        params = _fit_beta_method_of_moments(mean, sigma)
        if params is not None:
            alpha, beta = params
            return regularized_incomplete_beta(x, alpha, beta)

    # One-sided quantile, zero IQR, or Beta moments out of range → PWL.
    return _piecewise_linear_cdf(x, q25, mean, q75)


# ---------- Beta-distribution CDF ----------


def _fit_beta_method_of_moments(mean: float, sigma: float) -> tuple[float, float] | None:
    """Method-of-moments Beta(α, β) from mean and stdev.

    The Beta distribution has mean μ = α / (α+β) and variance
    σ² = αβ / ((α+β)² (α+β+1)).  Solving for (α, β):

        α + β = μ(1−μ)/σ² − 1
        α     = μ * (α + β)
        β     = (1−μ) * (α + β)

    Returns None when no Beta exists with the requested moments
    (i.e. σ² ≥ μ(1−μ)) or when the mean is at the {0, 1} boundary.
    """

    if not (0.0 < mean < 1.0):
        return None
    if sigma <= 0.0:
        return None
    var = sigma * sigma
    max_var = mean * (1.0 - mean)
    if var >= max_var:
        return None
    nu = max_var / var - 1.0
    if nu <= 0.0:
        return None
    alpha = mean * nu
    beta = (1.0 - mean) * nu
    if alpha <= 0.0 or beta <= 0.0:
        return None
    return alpha, beta


def regularized_incomplete_beta(x: float, a: float, b: float) -> float:
    """Compute I_x(a, b), the regularized incomplete beta function.

    Definition:  I_x(a, b) = B(x; a, b) / B(a, b)

    Implementation: Lentz's continued fraction expansion (Numerical Recipes,
    section 6.4). Switches branches at x < (a+1)/(a+b+2) for fast convergence.

    Boundary handling:
      * x ≤ 0  → 0.0
      * x ≥ 1  → 1.0
      * a, b > 0 required; raises ValueError otherwise.
    """

    if a <= 0.0 or b <= 0.0:
        raise ValueError(f"Beta parameters must be positive: a={a}, b={b}")
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    # log B(a, b) prefactor
    ln_pref = (
        math.lgamma(a + b)
        - math.lgamma(a)
        - math.lgamma(b)
        + a * math.log(x)
        + b * math.log(1.0 - x)
    )
    pref = math.exp(ln_pref)

    cf_threshold = (a + 1.0) / (a + b + 2.0)
    if x < cf_threshold:
        return pref * _beta_continued_fraction(x, a, b) / a
    # Use the symmetry I_x(a, b) = 1 − I_{1−x}(b, a) for fast convergence on
    # the upper branch.
    return 1.0 - pref * _beta_continued_fraction(1.0 - x, b, a) / b


def _beta_continued_fraction(x: float, a: float, b: float) -> float:
    """Lentz-modified continued fraction for the incomplete beta.

    Returns the value of the CF that, multiplied by the prefactor in
    `regularized_incomplete_beta`, gives the (un-regularized) incomplete beta.
    """

    MAX_ITER = 200
    EPS = 3.0e-16
    FPMIN = 1.0e-300

    qab = a + b
    qap = a + 1.0
    qam = a - 1.0

    c = 1.0
    d = 1.0 - qab * x / qap
    if abs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    h = d

    for m in range(1, MAX_ITER + 1):
        m2 = 2 * m
        # even step
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c
        # odd step
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if abs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < EPS:
            return h

    # Shouldn't happen for sensible inputs; return last estimate.
    return h


# ---------- piecewise-linear CDF (V2, kept as fallback) ----------


def _piecewise_linear_cdf(
    x: float,
    q25: float | None,
    mean: float | None,
    q75: float | None,
) -> float:
    """Estimate F(x) for β ∈ [0, 1] using anchored quantile knots.

    With both quantiles present, knots are (0, 0), (q25, 0.25), (mean, 0.5),
    (q75, 0.75), (1, 1) and the result is piecewise-linear between them.

    With only `mean` known, the estimator collapses to a **sharp step** at
    `mean` — 0 below, 0.5 at, 1.0 above. Linear interpolation between (0, 0),
    (mean, 0.5), (1, 1) would fabricate distributional shape that the data
    does not support, so we refuse to do that.

    With one quantile known but not the other, the knots present are used,
    plus the (0, 0) and (1, 1) anchors. With nothing known, returns 0.
    """

    if mean is None:
        return 0.0

    have_quantiles = q25 is not None or q75 is not None
    if not have_quantiles:
        # Mean-only fallback: sharp step at mean.
        if x < mean:
            return 0.0
        if x > mean:
            return 1.0
        return 0.5

    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0

    knots: list[tuple[float, float]] = [(0.0, 0.0)]
    if q25 is not None:
        knots.append((q25, 0.25))
    knots.append((mean, 0.5))
    if q75 is not None:
        knots.append((q75, 0.75))
    knots.append((1.0, 1.0))

    # Sort by β; collapse exact-tie β values to their highest CDF value so the
    # estimator stays a function (not a relation) at boundaries.
    knots.sort(key=lambda t: t[0])
    deduped: list[tuple[float, float]] = []
    for beta, cdf in knots:
        if deduped and deduped[-1][0] == beta:
            deduped[-1] = (beta, max(deduped[-1][1], cdf))
        else:
            deduped.append((beta, cdf))

    # Find the segment containing x and linearly interpolate.
    for (b0, f0), (b1, f1) in zip(deduped, deduped[1:], strict=False):
        if b0 <= x <= b1:
            if b1 == b0:
                return max(f0, f1)
            t = (x - b0) / (b1 - b0)
            return f0 + t * (f1 - f0)

    # Shouldn't reach here given the (0,0) and (1,1) anchors, but fall through.
    return 0.0
