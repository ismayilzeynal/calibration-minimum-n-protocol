"""Revision, gate fix round 1 (post-review, exploratory additions).

Design A is rerun with the submitted seeds so that the same subsamples,
random splits, temporal splits and joint bootstrap draws are regenerated:

  subsample + 50 random splits : default_rng([42, 5, N, sid])  (reliability.py)
  joint student bootstrap      : default_rng([42, 77, sid]), integers(0, n, n)

The submitted cache keeps the bootstrap draws for random halves only. This
script also keeps them for temporal halves, and adds:

  (1) crossing CIs as percentiles over ALL 1,000 draws, a draw that never
      crosses within 200 being ranked above 200 (the submitted code took the
      2.5th percentile of the crossing draws only and printed "> 200" for the
      upper bound whenever any draw failed to cross); both rules are written;
  (2) CIs for the temporal curve at every N, including N = 200;
  (3) the Monte Carlo SE of the mean over the 50 random splits at every N;
  (4) bias and Brier under the guess-corrected coding c' = .25 + .75 c
      (same subsamples, splits and draws);
  (5) a pool sensitivity analysis without students whose rated answers (all
      of them) have accuracy >= .99 or carry one constant rating; their rows
      are removed AFTER subsampling and splitting, so every other student
      keeps the submitted subsample; the bootstrap over the kept students
      uses the new tag [42, 401, sid];
  (6) AUROC under both rules: tie convention (submitted code) and plan rule
      (undefined wherever gamma is undefined).

Regression gates: random-half curves of every key against the submitted
curves (baseline and plan-rule reruns of stream A), and the regenerated
random-half bootstrap draws against both cached draw arrays.

Usage: python -B r2_A_rerun.py        (both samples; about 10 minutes)
Writes: revision-2026-11/results/r2_A.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
PAPER = os.path.dirname(REV)
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import indices as IX  # noqa: E402
import reliability as REL  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")
GRID = REL.GRID
R = REL.R_SPLITS
B = REL.B_BOOT
SEED = 42
SIDS = {"t12": 1, "t34": 2}
LEVELS_G = 0.25 + 0.75 * IX.LEVELS
KEYS = ["bias", "brier", "gamma", "auroc", "auroc_pr", "conf", "acc", "bias_g", "brier_g"]
BOOTK = ["bias", "brier", "gamma", "auroc", "auroc_pr", "bias_g", "brier_g"]
EXCK = ["bias", "brier"]
BIG = 1.0e6
THRS = (0.70, 0.80, 0.90)


def idx_all(cnt):
    """Submitted indices plus the plan-rule AUROC and the guess-coded bias and Brier."""
    o = IX.indices_from_counts(cnt)
    ntot = cnt.sum(axis=(1, 2))
    with np.errstate(divide="ignore", invalid="ignore"):
        nk = cnt.sum(axis=2)
        conf_g = (nk * LEVELS_G).sum(axis=1) / ntot
        brier_g = ((cnt[:, :, 0] * LEVELS_G ** 2).sum(axis=1)
                   + (cnt[:, :, 1] * (LEVELS_G - 1.0) ** 2).sum(axis=1)) / ntot
    out = {k: o[k] for k in ("bias", "brier", "gamma", "auroc", "conf", "acc")}
    out["auroc_pr"] = np.where(np.isnan(o["gamma"]), np.nan, o["auroc"])
    out["bias_g"] = conf_g - o["acc"]
    out["brier_g"] = brier_g
    return out


def halves(C, X, rng):
    """RNG call order identical to reliability.half_indices."""
    n, N = C.shape
    hA = N // 2
    A = {k: np.empty((n, R)) for k in KEYS}
    Bh = {k: np.empty((n, R)) for k in KEYS}
    for s in range(R):
        keys = rng.random((n, N))
        rank = np.argsort(np.argsort(keys, axis=1), axis=1)
        mA = rank < hA
        iA = idx_all(IX.count_tables(C, X, mA))
        iB = idx_all(IX.count_tables(C, X, ~mA))
        for k in KEYS:
            A[k][:, s] = iA[k]
            Bh[k][:, s] = iB[k]
    mT = np.zeros((n, N), dtype=bool)
    mT[:, :hA] = True
    tA = idx_all(IX.count_tables(C, X, mT))
    tB = idx_all(IX.count_tables(C, X, ~mT))
    full = idx_all(IX.count_tables(C, X, np.ones((n, N), bool)))
    TA = {k: tA[k][:, None] for k in KEYS}
    TB = {k: tB[k][:, None] for k in KEYS}
    return A, Bh, TA, TB, full


def cross(rvec, thr):
    """First crossing, log-linear between grid points (reliability.min_n_cross)."""
    logN = np.log(GRID)
    rv = np.asarray(rvec, float)
    if np.all(np.isnan(rv)):
        return np.nan
    for i, N in enumerate(GRID):
        if not np.isnan(rv[i]) and rv[i] >= thr:
            if i == 0:
                return float(GRID[0])
            r0, r1 = rv[i - 1], rv[i]
            if np.isnan(r0) or r1 == r0:
                return float(N)
            f = (thr - r0) / (r1 - r0)
            return float(np.exp(logN[i - 1] + f * (logN[i] - logN[i - 1])))
    return np.inf


def interp15(rvec):
    i10, i20 = GRID.index(10), GRID.index(20)
    w = (np.log(15) - np.log(10)) / (np.log(20) - np.log(10))
    return rvec[i10] + w * (rvec[i20] - rvec[i10])


def fnum(v):
    if v is None:
        return None
    v = float(v)
    if np.isnan(v):
        return None
    if v > 200 or np.isinf(v):
        return "inf"
    return v


def ci_all(draws):
    """Percentiles over all draws; a non-crossing draw ranks above 200."""
    d = np.where(np.isfinite(draws), draws, BIG)
    lo, hi = np.percentile(d, [2.5, 97.5])
    return [fnum(lo), fnum(hi)]


def ci_finite_only(draws):
    """The submitted rule (reliability.ci_or_partial), kept for comparison."""
    fin = draws[np.isfinite(draws)]
    if len(fin) == len(draws):
        return [float(v) for v in np.percentile(fin, [2.5, 97.5])]
    if len(fin):
        return [float(np.percentile(fin, 2.5)), "inf"]
    return ["inf", "inf"]


def crossing_block(point_curve, draws_mat):
    out = {}
    for thr in THRS:
        t = f"{int(round(thr * 100))}"
        d = np.array([cross(draws_mat[b], thr) for b in range(draws_mat.shape[0])])
        out[t] = {"point": fnum(cross(point_curve, thr)),
                  "ci_all_draws": ci_all(d),
                  "ci_finite_draws_only_submitted_rule": ci_finite_only(d),
                  "share_draws_crossing": float(np.mean(np.isfinite(d))),
                  "n_draws_crossing": int(np.isfinite(d).sum())}
    return out


def run(sample):
    t0 = time.time()
    sid = SIDS[sample]
    df = common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    poolA = counts >= max(GRID)
    pool = np.where(poolA)[0]
    nA = len(pool)

    # pool description on all rated answers of each pool student
    acc_full = np.array([x[offsets[i]:offsets[i + 1]].mean() for i in pool])
    const_full = np.array([codes[offsets[i]:offsets[i + 1]].min() == codes[offsets[i]:offsets[i + 1]].max()
                           for i in pool])
    excl = (acc_full >= 0.99) | const_full
    keep = ~excl
    pool_desc = {"n_pool": int(nA),
                 "all_rated_correct": int((acc_full == 1).sum()),
                 "rated_accuracy_ge_099": int((acc_full >= 0.99).sum()),
                 "constant_rating_all_rated": int(const_full.sum()),
                 "excluded_union": int(excl.sum()),
                 "excluded_union_share": float(excl.mean()),
                 "n_kept": int(keep.sum()),
                 "median_rated_accuracy": float(np.median(acc_full))}

    base = pd.read_csv(os.path.join(REV, "results", "baseline", "results",
                                    f"reliability_curves_{sample}.csv"))
    pre = pd.read_csv(os.path.join(REV, "results", "r1_A", "prereg", "results",
                                   f"reliability_curves_{sample}.csv"))
    bnpz = np.load(os.path.join(REV, "results", "baseline", "data", "cache", f"designA_{sample}.npz"))
    pnpz = np.load(os.path.join(REV, "results", "r1_A", "prereg", "data", "cache", f"designA_{sample}.npz"))

    G = len(GRID)
    pt = {k: {"random": np.full(G, np.nan), "temporal": np.full(G, np.nan)} for k in KEYS}
    mcse = {k: np.full(G, np.nan) for k in KEYS}
    split_sd = {k: np.full(G, np.nan) for k in KEYS}
    und_full = {k: np.full(G, np.nan) for k in ("gamma", "auroc", "auroc_pr")}
    exc_pt = {k: {"random": np.full(G, np.nan), "temporal": np.full(G, np.nan)} for k in EXCK}
    H = {}
    regress = {}
    for gi, N in enumerate(GRID):
        rngN = np.random.default_rng([SEED, 5, N, sid])
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolA, N, rngN)
        A, Bh, TA, TB, full = halves(C, X, rngN)
        H[N] = (A, Bh, TA, TB)
        for k in KEYS:
            cols = REL.sb_corr_cols(A[k], Bh[k])
            pt[k]["random"][gi] = float(np.nanmean(cols))
            pt[k]["temporal"][gi] = REL.sb_mean(TA[k], TB[k])
            v = cols[~np.isnan(cols)]
            split_sd[k][gi] = float(np.std(v, ddof=1))
            mcse[k][gi] = float(np.std(v, ddof=1) / np.sqrt(len(v)))
        for k in ("gamma", "auroc", "auroc_pr"):
            und_full[k][gi] = float(np.mean(np.isnan(full[k])))
        for k in EXCK:
            exc_pt[k]["random"][gi] = REL.sb_mean(A[k][keep], Bh[k][keep])
            exc_pt[k]["temporal"][gi] = REL.sb_mean(TA[k][keep], TB[k][keep])
        for k, src, kk in [("bias", base, "bias"), ("brier", base, "brier"), ("gamma", base, "gamma"),
                           ("auroc", base, "auroc"), ("conf", base, "conf"), ("acc", base, "acc"),
                           ("auroc_pr", pre, "auroc")]:
            for sp in ("random", "temporal"):
                row = src[(src.design == "A") & (src.split == sp) & (src["index"] == kk) & (src.N == N)]
                d = abs(float(row.r.iloc[0]) - pt[k][sp][gi])
                regress[f"curve_{k}_{sp}"] = max(regress.get(f"curve_{k}_{sp}", 0.0), d)
        print(f"[{sample}] N={N} halves done {time.time() - t0:.0f}s", flush=True)

    # joint bootstrap, submitted draws
    brng = np.random.default_rng([SEED, 77, sid])
    bs = {k: {"random": np.full((B, G), np.nan), "temporal": np.full((B, G), np.nan)} for k in BOOTK}
    for b in range(B):
        idx = brng.integers(0, nA, nA)
        for gi, N in enumerate(GRID):
            A, Bh, TA, TB = H[N]
            for k in BOOTK:
                bs[k]["random"][b, gi] = REL.sb_mean(A[k][idx], Bh[k][idx])
                bs[k]["temporal"][b, gi] = REL.sb_mean(TA[k][idx], TB[k][idx])
    print(f"[{sample}] bootstrap done {time.time() - t0:.0f}s", flush=True)
    for k in ("bias", "brier", "gamma", "auroc"):
        regress[f"boot_random_{k}_vs_submitted_cache"] = float(np.nanmax(np.abs(bs[k]["random"] - bnpz[f"boot_{k}"])))
    regress["boot_random_auroc_pr_vs_planrule_cache"] = float(
        np.nanmax(np.abs(bs["auroc_pr"]["random"] - pnpz["boot_auroc"])))

    # bootstrap over the kept students (pool sensitivity), new tag
    keep_idx = np.where(keep)[0]
    nK = len(keep_idx)
    erng = np.random.default_rng([SEED, 401, sid])
    ebs = {k: {"random": np.full((B, G), np.nan), "temporal": np.full((B, G), np.nan)} for k in EXCK}
    for b in range(B):
        idx = keep_idx[erng.integers(0, nK, nK)]
        for gi, N in enumerate(GRID):
            A, Bh, TA, TB = H[N]
            for k in EXCK:
                ebs[k]["random"][b, gi] = REL.sb_mean(A[k][idx], Bh[k][idx])
                ebs[k]["temporal"][b, gi] = REL.sb_mean(TA[k][idx], TB[k][idx])
    del H

    out ={"sample": sample, "n_poolA": int(nA), "grid": GRID, "R_splits": R, "B_boot": B,
           "pool_description": pool_desc, "indices": {}, "pool_sensitivity": {}}
    for k in KEYS:
        blk = {"curve_random_point": {str(N): float(pt[k]["random"][gi]) for gi, N in enumerate(GRID)},
               "curve_temporal_point": {str(N): float(pt[k]["temporal"][gi]) for gi, N in enumerate(GRID)},
               "mc_se_over_splits": {str(N): float(mcse[k][gi]) for gi, N in enumerate(GRID)},
               "sd_over_splits": {str(N): float(split_sd[k][gi]) for gi, N in enumerate(GRID)}}
        if k in BOOTK:
            for sp in ("random", "temporal"):
                D = bs[k][sp]
                lo = np.nanpercentile(D, 2.5, axis=0)
                hi = np.nanpercentile(D, 97.5, axis=0)
                blk[f"curve_{sp}_ci"] = {str(N): [float(lo[gi]), float(hi[gi])] for gi, N in enumerate(GRID)}
                blk[f"crossings_{sp}"] = crossing_block(pt[k][sp], D)
            r15i = np.array([interp15(bs[k]["random"][b]) for b in range(B)])
            r15d = bs[k]["random"][:, GRID.index(15)]
            blk["r15"] = {"interpolated_10_20_point": float(interp15(pt[k]["random"])),
                          "direct_grid_point": float(pt[k]["random"][GRID.index(15)]),
                          "share_draws_interpolated_ge_070": float(np.mean(r15i >= 0.70)),
                          "share_draws_direct_ge_070": float(np.mean(r15d >= 0.70))}
        if k in und_full:
            blk["undefined_full_set"] = {str(N): float(und_full[k][gi]) for gi, N in enumerate(GRID)}
        out["indices"][k] = blk
    for k in EXCK:
        blk = {"curve_random_point": {str(N): float(exc_pt[k]["random"][gi]) for gi, N in enumerate(GRID)},
               "curve_temporal_point": {str(N): float(exc_pt[k]["temporal"][gi]) for gi, N in enumerate(GRID)}}
        for sp in ("random", "temporal"):
            D = ebs[k][sp]
            lo = np.nanpercentile(D, 2.5, axis=0)
            hi = np.nanpercentile(D, 97.5, axis=0)
            blk[f"curve_{sp}_ci"] = {str(N): [float(lo[gi]), float(hi[gi])] for gi, N in enumerate(GRID)}
            blk[f"crossings_{sp}"] = crossing_block(exc_pt[k][sp], D)
        out["pool_sensitivity"][k] = blk
    out["regression_max_abs_diff"] = regress
    out["regression_pass"] = bool(max(regress.values()) <= 1e-12)
    out["seeds"] = {"subsample_and_splits": "[42, 5, N, sid] (as submitted)",
                    "joint_bootstrap": "[42, 77, sid] (as submitted)",
                    "pool_sensitivity_bootstrap": "[42, 401, sid] (new)", "sid": sid}
    out["runtime_seconds"] = round(time.time() - t0, 1)
    print(sample, "done", out["runtime_seconds"], "regression pass", out["regression_pass"],
          json.dumps(regress), flush=True)
    return out


if __name__ == "__main__":
    res = {"_definitions": {
        "ci_all_draws": "2.5th and 97.5th percentiles over all 1,000 joint bootstrap draws; a draw whose curve never crosses within 200 is ranked above 200; 'inf' = bound above 200",
        "ci_finite_draws_only_submitted_rule": "submitted rule: percentiles of the crossing draws only; upper bound 'inf' whenever any draw fails to cross",
        "mc_se_over_splits": "SD of the 50 per-split Spearman-Brown values divided by sqrt(50): Monte Carlo SE of the reported mean",
        "bias_g, brier_g": "guess-corrected coding c' = .25 + .75 c (rating 0 read as the 25% guessing rate of a four-option item)",
        "auroc": "tie convention (submitted code); auroc_pr = plan rule, undefined wherever gamma is undefined",
        "pool_sensitivity": "Design A without students whose rated answers have accuracy >= .99 or one constant rating (computed on all of their rated answers); rows removed after subsampling and splitting",
        "r15": "interpolated_10_20 = frozen H2 quantity (log-linear between N = 10 and 20); direct = grid point N = 15 (addendum 1)"}}
    for s, lab in (("t12", "S1"), ("t34", "S2")):
        res[lab] = run(s)
    with open(os.path.join(REV, "results", "r2_A.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("written results/r2_A.json")
