#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["numpy>=2.0", "scipy>=1.13"]
# ///
"""limma-eBayes moderated t-statistic per probe, pure-python.

Implements the Smyth (2004) empirical-Bayes moderated t-statistic. For
each probe i we observe samples from two groups, compute the ordinary
residual variance s²_i and mean difference Δ_i, and then *moderate*
s²_i by borrowing strength across probes via a scaled-inverse-χ²
prior fit by method-of-moments on log(s²_i).

Output per probe: moderated t = Δ_i / sqrt(s̃²_i × (1/n1 + 1/n2)),
where s̃²_i = (d0·s0² + d·s²_i) / (d0 + d) and the total degrees of
freedom are d0 + d.

Also implements a small CLI for applying the transform to a β matrix
TSV (probe_id + per-sample β columns) given sample-group labels.

Reference: Smyth, G.K. (2004). Linear models and empirical Bayes
methods for assessing differential expression in microarray
experiments. Stat. Appl. Genet. Mol. Biol., 3(1), Article 3. Core
formulas: §6 (moderated t) and Appendix A (variance prior fit by
fitFDist / moment matching on log-scale).

Caveat: at tiny per-probe d (e.g. 2 samples per side → d = 2), the
log-scale variance of s²_i approaches the irreducible noise floor
trigamma(d/2), the prior-df estimator diverges to ∞, and the
moderated-t collapses toward a signed effect-size / prior-sd
statistic. This script handles that limit deterministically and
emits d0 = inf when it happens.
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path

import numpy as np
from scipy.special import digamma, polygamma
from scipy.stats import t as tdist


def fit_variance_prior(s_sq: np.ndarray, d: float, tol: float = 1e-8, max_iter: int = 500) -> tuple[float, float]:
    """Smyth (2004) MoM fit of a scaled-inverse-χ² prior to s²_i at fixed df d.

    Returns (s0_sq, d0). When the log-variance is at/below the
    irreducible trigamma(d/2) floor (over-dispersion not present),
    returns (exp(mean_log_s_sq − digamma(d/2) + log(d/2)), +inf) —
    indicating the prior dominates and the moderated-t becomes an
    effect-size statistic.
    """
    s_sq = np.asarray(s_sq, dtype=np.float64)
    s_sq = s_sq[np.isfinite(s_sq) & (s_sq > 0)]
    if s_sq.size < 10:
        raise ValueError(f"need at least 10 non-zero probes to fit a variance prior; got {s_sq.size}")

    log_s_sq = np.log(s_sq)
    # E[log(s²_i)] = log(s0²) − log(d/2) + digamma(d/2) + log(d0/2) − digamma(d0/2)
    #             = log(s0²) + (digamma(d/2) − log(d/2)) − (digamma(d0/2) − log(d0/2))
    # Var[log(s²_i)] = trigamma(d/2) + trigamma(d0/2)
    mean_log = float(np.mean(log_s_sq))
    var_log = float(np.var(log_s_sq, ddof=1))

    target_trigamma = var_log - polygamma(1, d / 2)
    if target_trigamma <= 0:
        # Prior dominates; no over-dispersion. d0 → ∞, s0² = exp(mean_log − digamma(d/2) + log(d/2)).
        s0_sq = float(np.exp(mean_log - digamma(d / 2) + np.log(d / 2)))
        return s0_sq, float("inf")

    # Solve trigamma(d0_half) = target_trigamma for d0_half = d0/2.
    # trigamma is strictly decreasing on (0, ∞); Newton's method is stable.
    d0_half = 1.0
    for _ in range(max_iter):
        f = polygamma(1, d0_half) - target_trigamma
        fp = polygamma(2, d0_half)  # derivative (negative)
        if fp == 0.0:
            break
        step = f / fp
        new_d0_half = d0_half - step
        if new_d0_half <= 0:
            new_d0_half = d0_half / 2.0
        if abs(new_d0_half - d0_half) < tol:
            d0_half = new_d0_half
            break
        d0_half = new_d0_half
    d0 = 2.0 * d0_half
    # Log-scale identity: log(s0²) = mean_log + (digamma(d0_half) − log(d0_half)) − (digamma(d/2) − log(d/2))
    s0_sq = float(np.exp(mean_log - (digamma(d / 2) - np.log(d / 2)) + (digamma(d0_half) - np.log(d0_half))))
    return s0_sq, float(d0)


def limma_ebayes(
    betas: np.ndarray,
    group: np.ndarray,
) -> dict[str, np.ndarray]:
    """Compute moderated-t for every row of `betas`.

    Args:
        betas: (n_probes, n_samples) array of β values. Missing → np.nan.
        group: (n_samples,) array of 0/1 group labels. 0 = reference
            (e.g. normal), 1 = test (e.g. tumor). `Δ = mean(1) − mean(0)`.

    Returns:
        dict with keys:
          - "delta":   (n_probes,) mean β difference, group1 − group0
          - "s_sq":    (n_probes,) ordinary pooled variance
          - "s_tilde_sq": (n_probes,) moderated variance
          - "t_mod":   (n_probes,) moderated t-statistic
          - "df_tot":  (n_probes,) total df = d + d0
          - "p_value": (n_probes,) two-sided p from moderated-t null
          - "s0_sq":   scalar prior variance
          - "d0":      scalar prior df (or +inf)
          - "d":       scalar per-probe df
    """
    betas = np.asarray(betas, dtype=np.float64)
    group = np.asarray(group)
    if betas.ndim != 2:
        raise ValueError(f"betas must be 2D; got shape {betas.shape}")
    if group.shape != (betas.shape[1],):
        raise ValueError(f"group must have len = n_samples = {betas.shape[1]}; got {group.shape}")

    mask0 = group == 0
    mask1 = group == 1
    n0 = int(mask0.sum()); n1 = int(mask1.sum())
    if n0 < 2 or n1 < 2:
        raise ValueError(f"need at least 2 samples per group; got n0={n0}, n1={n1}")

    g0 = betas[:, mask0]; g1 = betas[:, mask1]

    # Per-probe means and residual-variance across NaN-aware samples.
    # For the sample-size we report the actual non-NaN count per side per
    # probe; probes with < 2 valid samples either side get s²_i = NaN and
    # are excluded from the prior fit.
    n0_vec = np.sum(~np.isnan(g0), axis=1)
    n1_vec = np.sum(~np.isnan(g1), axis=1)
    mean0 = np.nanmean(g0, axis=1)
    mean1 = np.nanmean(g1, axis=1)
    delta = mean1 - mean0

    # Pooled variance at the complete-case df count.
    # var0 = sum((x−mean0)²) / (n0 − 1); same for var1.
    ss0 = np.nansum((g0 - mean0[:, None]) ** 2, axis=1)
    ss1 = np.nansum((g1 - mean1[:, None]) ** 2, axis=1)
    # Treat any row with n0_vec < 2 or n1_vec < 2 as degenerate (NaN).
    degenerate = (n0_vec < 2) | (n1_vec < 2)
    d_vec = (n0_vec + n1_vec - 2).astype(np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        s_sq = (ss0 + ss1) / d_vec
    s_sq[degenerate] = np.nan

    # Majority d (use the mode df; we'll fit the prior on probes that
    # share that df for the cleanest Smyth closed-form). In practice all
    # probes have the same d if all samples have no missing values at
    # that probe.
    d_values, d_counts = np.unique(d_vec[~degenerate], return_counts=True)
    d_mode = float(d_values[np.argmax(d_counts)])

    # Fit prior on probes at d_mode.
    mask_mode = (d_vec == d_mode) & ~degenerate & np.isfinite(s_sq) & (s_sq > 0)
    s0_sq, d0 = fit_variance_prior(s_sq[mask_mode], d=d_mode)

    # Moderated variance per probe. When d0 = inf, s̃² → s0² uniformly.
    if np.isinf(d0):
        s_tilde_sq = np.full_like(s_sq, s0_sq)
    else:
        with np.errstate(divide="ignore", invalid="ignore"):
            s_tilde_sq = (d0 * s0_sq + d_vec * s_sq) / (d0 + d_vec)
        s_tilde_sq[degenerate] = np.nan

    # Moderated t-statistic.
    # Standard error: sqrt(s̃² · (1/n0 + 1/n1)) — computed per probe with valid counts.
    with np.errstate(divide="ignore", invalid="ignore"):
        se = np.sqrt(s_tilde_sq * (1.0 / n0_vec + 1.0 / n1_vec))
    t_mod = delta / se

    # Degrees of freedom of the moderated t.
    df_tot = d_vec + (0 if np.isinf(d0) else d0)  # d + d0; if d0 inf, df_tot → inf → normal limit
    # Two-sided p-value. When d0 = inf, use the standard normal limit.
    if np.isinf(d0):
        from scipy.stats import norm as _norm
        p = 2.0 * (1.0 - _norm.cdf(np.abs(t_mod)))
    else:
        p = 2.0 * (1.0 - tdist.cdf(np.abs(t_mod), df=df_tot))

    return {
        "delta": delta,
        "s_sq": s_sq,
        "s_tilde_sq": s_tilde_sq,
        "t_mod": t_mod,
        "df_tot": df_tot,
        "p_value": p,
        "s0_sq": s0_sq,
        "d0": d0,
        "d": d_mode,
        "n0": n0,
        "n1": n1,
    }


# ---------- CLI ----------


def _open_maybe_gz(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _read_beta_matrix(path: Path) -> tuple[list[str], list[str], np.ndarray]:
    """Read a TSV with header `probe_id<TAB>sample1<TAB>sample2<TAB>…`
    and body `probeid<TAB>β<TAB>β<TAB>…` (NA-tolerant)."""
    with _open_maybe_gz(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        samples = header[1:]
        probes: list[str] = []
        rows: list[list[float]] = []
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            probes.append(parts[0])
            vals = []
            for v in parts[1:]:
                if v == "" or v == "NA" or v.lower() == "nan":
                    vals.append(np.nan)
                else:
                    try:
                        vals.append(float(v))
                    except ValueError:
                        vals.append(np.nan)
            rows.append(vals)
    return probes, samples, np.asarray(rows, dtype=np.float64)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--beta-matrix", type=Path,
                    help="TSV (optionally .gz): probe_id + per-sample β columns")
    ap.add_argument("--group", type=Path,
                    help="TSV: sample_id<TAB>group (0 = normal / ref, 1 = tumor / test)")
    ap.add_argument("--output", type=Path,
                    help="Output TSV: probe_id, delta, s_sq, s_tilde_sq, t_mod, df_tot, p_value")
    ap.add_argument("--self-test", action="store_true",
                    help="Run a self-test on synthetic data and exit.")
    args = ap.parse_args()

    if args.self_test:
        return _self_test()
    if not (args.beta_matrix and args.group and args.output):
        ap.error("--beta-matrix, --group, --output are required (or pass --self-test)")

    probes, samples, betas = _read_beta_matrix(args.beta_matrix)
    print(f"loaded {len(probes):,} probes × {len(samples)} samples from {args.beta_matrix.name}", flush=True)

    # Load groups.
    group_map: dict[str, int] = {}
    with args.group.open() as fh:
        for i, line in enumerate(fh):
            parts = line.rstrip("\n").split("\t")
            if i == 0 and (parts[0].lower() == "sample_id" or parts[0].lower() == "sample"):
                continue
            if len(parts) < 2:
                continue
            group_map[parts[0]] = int(parts[1])
    group = np.array([group_map[s] for s in samples if s in group_map])
    keep_cols = np.array([s in group_map for s in samples])
    betas_used = betas[:, keep_cols]
    samples_used = [s for s in samples if s in group_map]
    n0 = int((group == 0).sum()); n1 = int((group == 1).sum())
    print(f"  using {len(samples_used)} samples: {n0} group=0, {n1} group=1", flush=True)

    r = limma_ebayes(betas_used, group)
    print(f"  d = {r['d']}; prior: s0² = {r['s0_sq']:.4g}, d0 = {r['d0']!r}", flush=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as fh:
        fh.write("probe_id\tdelta_beta\ts_sq\ts_tilde_sq\tt_mod\tdf_tot\tp_value\n")
        for i, pid in enumerate(probes):
            fh.write(
                f"{pid}\t"
                f"{r['delta'][i]:.6g}\t"
                f"{r['s_sq'][i]:.6g}\t"
                f"{r['s_tilde_sq'][i]:.6g}\t"
                f"{r['t_mod'][i]:.6g}\t"
                f"{r['df_tot'][i]:.6g}\t"
                f"{r['p_value'][i]:.6g}\n"
            )
    print(f"wrote {args.output}", flush=True)

    # Emit a small companion JSON with global fit parameters for provenance.
    meta_path = args.output.with_suffix(".meta.json")
    meta_path.write_text(json.dumps({
        "n_probes": int(len(probes)),
        "n_samples_used": int(len(samples_used)),
        "n0": n0, "n1": n1,
        "d": r["d"],
        "s0_sq": r["s0_sq"],
        "d0": None if np.isinf(r["d0"]) else r["d0"],
        "d0_is_infinite": bool(np.isinf(r["d0"])),
        "notes": ("d0 = inf means the log-variance across probes did not "
                  "exceed the irreducible trigamma(d/2) noise floor — the "
                  "moderated-t collapses to a signed-effect-size-over-prior-sd "
                  "statistic. This is the expected regime at n0 = n1 = 2 per "
                  "the Smyth (2004) prior-fit limit."),
    }, indent=2))
    print(f"wrote {meta_path}", flush=True)
    return 0


# ---------- self-test ----------


def _self_test() -> int:
    """Validate limma-eBayes against two independent checks:

    1. **d0 → ∞ limit agrees with an unmoderated ordinary-t formula.** When
       the prior dominates (all probes share a single variance estimate s0²),
       the moderated-t reduces to `Δ / sqrt(s0² · (1/n0 + 1/n1))`. Check
       closed-form.
    2. **Null p-values are uniform under the null.** Pure Beta noise with
       no spikes: the empirical p-value CDF should be ≈ identity (KS < 0.1
       at n_probes = 5000).
    """
    rng = np.random.default_rng(42)
    n_probes = 5000
    n0, n1 = 5, 5

    # Pure-null β matrix; same distribution in both groups.
    base = rng.beta(2.0, 2.0, size=(n_probes, n0 + n1))
    group = np.concatenate([np.zeros(n0), np.ones(n1)]).astype(int)
    r = limma_ebayes(base, group)

    print(f"self-test: n_probes={n_probes}, n0=n1={n0}")
    print(f"  d = {r['d']}, d0 = {r['d0']!r}, s0² = {r['s0_sq']:.4f}")

    t = r["t_mod"]
    delta = r["delta"]

    # (1) Closed-form agreement when d0 = inf.
    if np.isinf(r["d0"]):
        expected_se = np.sqrt(r["s0_sq"] * (1.0 / n0 + 1.0 / n1))
        expected_t = delta / expected_se
        max_abs_err = float(np.nanmax(np.abs(t - expected_t)))
        print(f"  d0 = ∞; |t − Δ/√(s0²(1/n0+1/n1))| max = {max_abs_err:.2e}")
        assert max_abs_err < 1e-9, f"closed-form check failed: max_abs_err = {max_abs_err:.4g}"
    else:
        # Finite d0: moderated variance s̃² = (d0·s0² + d·s²)/(d0+d); check one probe.
        i = 0
        expected = (r["d0"] * r["s0_sq"] + r["d"] * r["s_sq"][i]) / (r["d0"] + r["d"])
        max_abs_err = float(abs(r["s_tilde_sq"][i] - expected))
        print(f"  d0 = {r['d0']:.3f}; per-probe s̃² closed-form check: |err| = {max_abs_err:.2e}")
        assert max_abs_err < 1e-9

    # (2) Null p-values are uniform.
    p = r["p_value"]
    valid = np.isfinite(p)
    p_sorted = np.sort(p[valid])
    uniform = np.linspace(0, 1, len(p_sorted))
    ks = float(np.max(np.abs(p_sorted - uniform)))
    print(f"  null p-value KS distance to uniform: {ks:.3f}")
    assert ks < 0.10, f"null p-values deviate too far from uniform (KS={ks:.3f})"

    # Bonus sanity: Δ should be mean-zero under the null (no spike).
    mean_delta = float(np.mean(delta[valid]))
    print(f"  mean(Δ) under null: {mean_delta:+.4f} (expected ≈ 0)")
    assert abs(mean_delta) < 0.01, f"null Δ mean drift {mean_delta:+.4f}"

    print("self-test PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
