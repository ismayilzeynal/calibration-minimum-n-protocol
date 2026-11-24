"""Shared helpers for the Monte Carlo precision analyses (r2_E_*) and the
pre-specified-grid crossings (r2_F_grid.py).

Crossing rule: identical to reliability.run().min_n_cross (first grid point
at or above the level, log-linear interpolation from the previous point),
applied to any grid. A curve that never reaches the level within the grid
returns inf ("> 200").
"""
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
PAPER = os.path.dirname(REV)
RESULTS = os.path.join(REV, "results")
PARTS = os.path.join(RESULTS, "r2_E_parts")
SEED = 42
SIDS = {"t12": 1, "t34": 2}
LAB = {"t12": "S1", "t34": "S2"}
FULL_GRID = [5, 10, 12, 15, 20, 25, 30, 40, 50, 100, 200]
FROZEN_GRID = [5, 10, 20, 50, 100, 200]
THRS = (0.70, 0.80, 0.90)
BIG = 1.0e6
S2_MEDIAN_VOLUME = 31  # replication criterion (b), analysis plan section 7


def cross(rvec, thr, grid=FULL_GRID):
    lg = np.log(grid)
    rv = np.asarray(rvec, float)
    if np.all(np.isnan(rv)):
        return np.nan
    for i, N in enumerate(grid):
        if not np.isnan(rv[i]) and rv[i] >= thr:
            if i == 0:
                return float(grid[0])
            r0, r1 = rv[i - 1], rv[i]
            if np.isnan(r0) or r1 == r0:
                return float(N)
            f = (thr - r0) / (r1 - r0)
            return float(np.exp(lg[i - 1] + f * (lg[i] - lg[i - 1])))
    return np.inf


def on_frozen(rvec_full):
    """Values of a full-grid curve at the six pre-specified grid points."""
    rv = np.asarray(rvec_full, float)
    return rv[[FULL_GRID.index(n) for n in FROZEN_GRID]]


def interp15(rvec_full):
    """Frozen H2 quantity: log-linear interpolation between N = 10 and 20."""
    i10, i20 = FULL_GRID.index(10), FULL_GRID.index(20)
    w = (np.log(15) - np.log(10)) / (np.log(20) - np.log(10))
    return rvec_full[i10] + w * (rvec_full[i20] - rvec_full[i10])


def fnum(v):
    if v is None:
        return None
    v = float(v)
    if np.isnan(v):
        return None
    if np.isinf(v) or v > 200:
        return "inf"
    return v


def delta_se_cross(rvec, se, thr, grid=FULL_GRID):
    """Delta-method SE of the crossing from independent per-N SEs.
    None when the crossing is not interior (first grid point or no crossing)."""
    lg = np.log(grid)
    rv = np.asarray(rvec, float)
    for i, N in enumerate(grid):
        if not np.isnan(rv[i]) and rv[i] >= thr:
            if i == 0 or np.isnan(rv[i - 1]) or rv[i] == rv[i - 1]:
                return None
            r0, r1 = rv[i - 1], rv[i]
            f = (thr - r0) / (r1 - r0)
            c = np.exp(lg[i - 1] + f * (lg[i] - lg[i - 1]))
            dl = lg[i] - lg[i - 1]
            d0 = c * dl * (thr - r1) / (r1 - r0) ** 2
            d1 = -c * dl * (thr - r0) / (r1 - r0) ** 2
            return float(np.sqrt(d0 ** 2 * se[i - 1] ** 2 + d1 ** 2 * se[i] ** 2))
    return None


def draws_summary(d, extra_gt=None):
    """Summary of crossing values over Monte Carlo draws; a draw that never
    reaches the level within 200 is ranked above 200 (as ci_all_draws)."""
    d = np.asarray(d, float)
    fin = np.isfinite(d)
    dd = np.where(fin, d, BIG)
    lo, med, hi = np.percentile(dd, [2.5, 50, 97.5])
    out = {"n_draws": int(len(d)), "share_reaching_within_200": float(fin.mean()),
           "median": fnum(med), "p2_5": fnum(lo), "p97_5": fnum(hi),
           "min": fnum(dd.min()), "max": fnum(dd.max())}
    if fin.all():
        out["mean"] = float(d.mean())
        out["sd"] = float(d.std(ddof=1)) if len(d) > 1 else None
    else:
        out["mean"] = None
        out["sd"] = None
        out["sd_of_draws_reaching"] = float(d[fin].std(ddof=1)) if fin.sum() > 1 else None
    if extra_gt is not None:
        out[f"share_above_{extra_gt}"] = float(np.mean(dd > extra_gt))
    return out


def resample_split_crossings(per_split, thrs, rng, B, grid=FULL_GRID):
    """Resample the per-split values with replacement within each N and
    recompute the mean curve and its crossings (split Monte Carlo error only).
    per_split: list over grid of 1-D arrays (NaN allowed)."""
    G = len(grid)
    curves = np.empty((B, G))
    for gi in range(G):
        v = np.asarray(per_split[gi], float)
        n = len(v)
        idx = rng.integers(0, n, (B, n))
        with np.errstate(invalid="ignore"):
            curves[:, gi] = np.nanmean(v[idx], axis=1)
    out = {}
    for thr in thrs:
        cr = np.array([cross(curves[b], thr, grid) for b in range(B)])
        fin = np.isfinite(cr)
        out[f"{int(round(thr * 100))}"] = {
            "resample_sd": float(cr[fin].std(ddof=1)) if fin.sum() > 1 else None,
            "resample_share_reaching_within_200": float(fin.mean()),
            "resample_p2_5": fnum(np.percentile(np.where(fin, cr, BIG), 2.5)),
            "resample_p97_5": fnum(np.percentile(np.where(fin, cr, BIG), 97.5))}
    return out, curves


def mc_se(v):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    return float(np.std(v, ddof=1) / np.sqrt(len(v))) if len(v) > 1 else None
