"""Split-half reliability machinery (Designs A/B, random + temporal splits),
bootstrap CIs, and minimum-N estimation. Writes results/ CSV/JSON + cache npz.

Run from the paper root:  python scripts/reliability.py t12
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import indices as IX

GRID = [5, 10, 12, 15, 20, 25, 30, 40, 50, 100, 200]  # Design A (densified, addendum 1)
GRID_B = [5, 10, 20, 50, 100, 200]                     # Design B (original frozen grid)
FLAG_N = [10, 15, 20, 30]                              # flag-consistency grid (addendum 5)
R_SPLITS = 50
B_BOOT = 1000
MASTER_SEED = 42
MIN_POOL = 300
STORE_N = {20, 100}  # Design A half/full caches kept for the confound analysis

IDX_KEYS = IX.CONFIRMATORY + IX.EXPLORATORY + ["acc", "conf"]


def student_arrays(df):
    """Tagged responses grouped per student, date-sorted. Returns
    (user_ids, offsets, codes, x) with codes int8 level 0-4."""
    tg = common.tagged(df).sort_values(["UserId", "DateAnswered"], kind="stable")
    codes = (tg["Confidence"].to_numpy() / 25.0).round().astype(np.int8)
    x = tg["IsCorrect"].to_numpy().astype(np.int8)
    uid = tg["UserId"].to_numpy()
    users, starts = np.unique(uid, return_index=True)
    offsets = np.append(starts, len(uid))
    return users, offsets, codes, x


def subsample_matrix(users, offsets, codes, x, pool_mask, N, rng):
    """(n_pool, N) code/x matrices; columns date-ordered within student."""
    pool = np.where(pool_mask)[0]
    C = np.empty((len(pool), N), dtype=np.int8)
    X = np.empty((len(pool), N), dtype=np.int8)
    for r, i in enumerate(pool):
        s, e = offsets[i], offsets[i + 1]
        sel = np.sort(rng.choice(e - s, size=N, replace=False)) + s
        C[r] = codes[sel]
        X[r] = x[sel]
    return C, X, users[pool]


def half_indices(C, X, rng):
    """Random-split halves (R_SPLITS) + one temporal split.
    Returns dict key -> (A, B) arrays of shape (n, R_SPLITS) / (n, 1)."""
    n, N = C.shape
    hA = N // 2
    outA = {k: np.empty((n, R_SPLITS)) for k in IDX_KEYS}
    outB = {k: np.empty((n, R_SPLITS)) for k in IDX_KEYS}
    for s in range(R_SPLITS):
        keys = rng.random((n, N))
        rank = np.argsort(np.argsort(keys, axis=1), axis=1)
        mA = rank < hA
        iA = IX.indices_from_counts(IX.count_tables(C, X, mA))
        iB = IX.indices_from_counts(IX.count_tables(C, X, ~mA))
        for k in IDX_KEYS:
            outA[k][:, s] = iA[k]
            outB[k][:, s] = iB[k]
    # temporal: first hA columns (date-ordered) vs rest
    mT = np.zeros((n, N), dtype=bool)
    mT[:, :hA] = True
    tA = IX.indices_from_counts(IX.count_tables(C, X, mT))
    tB = IX.indices_from_counts(IX.count_tables(C, X, ~mT))
    outTA = {k: tA[k][:, None] for k in IDX_KEYS}
    outTB = {k: tB[k][:, None] for k in IDX_KEYS}
    full = IX.indices_from_counts(IX.count_tables(C, X, np.ones((n, N), bool)))
    return outA, outB, outTA, outTB, full


def sb_corr_cols(A, B, min_pairs=10):
    """Column-wise nan-aware Pearson + Spearman-Brown. A,B (n, S) -> (S,)"""
    M = ~(np.isnan(A) | np.isnan(B))
    Az = np.where(M, A, 0.0)
    Bz = np.where(M, B, 0.0)
    n = M.sum(axis=0).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        sA, sB = Az.sum(0), Bz.sum(0)
        cov = (Az * Bz).sum(0) - sA * sB / n
        vA = (Az * Az).sum(0) - sA ** 2 / n
        vB = (Bz * Bz).sum(0) - sB ** 2 / n
        r = cov / np.sqrt(vA * vB)
        sb = 2.0 * r / (1.0 + r)
        sb = np.where(n >= min_pairs, sb, np.nan)
    return sb


def sb_mean(A, B):
    return float(np.nanmean(sb_corr_cols(A, B)))


def run(sample):
    rng = np.random.default_rng([MASTER_SEED, {"t12": 1, "t34": 2}[sample]])
    df = common.build_responses(sample)
    users, offsets, codes, x = student_arrays(df)
    counts = np.diff(offsets)
    os.makedirs("results", exist_ok=True)
    os.makedirs(os.path.join("data", "cache"), exist_ok=True)

    rows = []
    bootA = {k: np.full((B_BOOT, len(GRID)), np.nan) for k in IX.CONFIRMATORY}
    bootA_T = {k: np.full((B_BOOT, len(GRID)), np.nan) for k in IX.CONFIRMATORY}
    cacheA = {}

    # ---------- Design A: population-constant pool (>= max grid) ----------
    poolA = counts >= max(GRID)
    nA = int(poolA.sum())
    hA_data = {}
    flags = []
    sid = {"t12": 1, "t34": 2}[sample]
    for N in GRID:
        rngN = np.random.default_rng([MASTER_SEED, 5, N, sid])
        C, X, uidA = subsample_matrix(users, offsets, codes, x, poolA, N, rngN)
        A, B, TA, TB, full = half_indices(C, X, rngN)
        hA_data[N] = (A, B, TA, TB)
        if N in STORE_N:
            cacheA[N] = (A, B, full, uidA)
        for k in IDX_KEYS:
            und = float(np.mean(np.isnan(A[k])))
            rows.append(dict(sample=sample, design="A", split="random", index=k,
                             N=N, n_pool=nA, r=sb_mean(A[k], B[k]),
                             pct_undefined=und,
                             pct_undefined_full=float(np.mean(np.isnan(full[k])))))
            rows.append(dict(sample=sample, design="A", split="temporal", index=k,
                             N=N, n_pool=nA, r=sb_mean(TA[k], TB[k]),
                             pct_undefined=float(np.mean(np.isnan(TA[k]))),
                             pct_undefined_full=float(np.mean(np.isnan(full[k])))))
        if N in FLAG_N:  # addendum 5: raw-bias flag consistency across halves
            k10 = max(int(round(nA * 0.10)), 1)
            jac, agree = [], []
            for s in range(R_SPLITS):
                a, b = A["bias"][:, s], B["bias"][:, s]
                tA_ = set(np.argsort(-a)[:k10])
                tB_ = set(np.argsort(-b)[:k10])
                jac.append(len(tA_ & tB_) / len(tA_ | tB_))
                agree.append(len(tA_ & tB_) / k10)
            flags.append(dict(N=N, jaccard=float(np.mean(jac)),
                              pct_flags_retained=float(np.mean(agree))))
    # joint bootstrap across the grid (same student resample per draw)
    brng = np.random.default_rng([MASTER_SEED, 77, {"t12": 1, "t34": 2}[sample]])
    for b in range(B_BOOT):
        idx = brng.integers(0, nA, nA)
        for gi, N in enumerate(GRID):
            A, B, TA, TB = hA_data[N]
            for k in IX.CONFIRMATORY:
                bootA[k][b, gi] = sb_mean(A[k][idx], B[k][idx])
                bootA_T[k][b, gi] = sb_mean(TA[k][idx], TB[k][idx])
    del hA_data

    # ---------- Design B: all-eligible per N ----------
    for N in GRID_B:
        poolB = counts >= N
        nB = int(poolB.sum())
        if nB < MIN_POOL:
            continue
        rngN = np.random.default_rng([MASTER_SEED, 6, N, sid])
        C, X, _ = subsample_matrix(users, offsets, codes, x, poolB, N, rngN)
        A, B, TA, TB, _ = half_indices(C, X, rngN)
        brngB = np.random.default_rng([MASTER_SEED, 78, N,
                                       {"t12": 1, "t34": 2}[sample]])
        for k in IDX_KEYS:
            r_pt = sb_mean(A[k], B[k])
            if k in IX.CONFIRMATORY:
                draws = np.empty(B_BOOT)
                for b in range(B_BOOT):
                    idx = brngB.integers(0, nB, nB)
                    draws[b] = sb_mean(A[k][idx], B[k][idx])
                lo, hi = np.nanpercentile(draws, [2.5, 97.5])
            else:
                lo = hi = np.nan
            rows.append(dict(sample=sample, design="B", split="random", index=k,
                             N=N, n_pool=nB, r=r_pt, ci_lo=lo, ci_hi=hi,
                             pct_undefined=float(np.mean(np.isnan(A[k])))))
            rows.append(dict(sample=sample, design="B", split="temporal", index=k,
                             N=N, n_pool=nB, r=sb_mean(TA[k], TB[k]),
                             pct_undefined=float(np.mean(np.isnan(TA[k])))))

    # ---------- assemble curves + CIs for Design A ----------
    curves = pd.DataFrame(rows)
    for k in IX.CONFIRMATORY:
        lo = np.nanpercentile(bootA[k], 2.5, axis=0)
        hi = np.nanpercentile(bootA[k], 97.5, axis=0)
        loT = np.nanpercentile(bootA_T[k], 2.5, axis=0)
        hiT = np.nanpercentile(bootA_T[k], 97.5, axis=0)
        for gi, N in enumerate(GRID):
            m = ((curves["sample"] == sample) & (curves.design == "A")
                 & (curves["index"] == k) & (curves.N == N))
            curves.loc[m & (curves.split == "random"), "ci_lo"] = lo[gi]
            curves.loc[m & (curves.split == "random"), "ci_hi"] = hi[gi]
            curves.loc[m & (curves.split == "temporal"), "ci_lo"] = loT[gi]
            curves.loc[m & (curves.split == "temporal"), "ci_hi"] = hiT[gi]
    curves.to_csv(os.path.join("results", f"reliability_curves_{sample}.csv"),
                  index=False)

    # ---------- min-N + r(15) from Design A ----------
    logN = np.log(GRID)

    def min_n_cross(rvec, thr):
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

    def ci_or_partial(draws):
        fin = draws[np.isfinite(draws)]
        if len(fin) == B_BOOT:
            return [float(v) for v in np.percentile(fin, [2.5, 97.5])]
        if len(fin):
            return [float(np.percentile(fin, 2.5)), None]
        return [None, None]

    minn = {}
    pt = {k: curves[(curves.design == "A") & (curves.split == "random")
                    & (curves["index"] == k)].sort_values("N").r.to_numpy()
          for k in IX.CONFIRMATORY}
    ptT = {k: curves[(curves.design == "A") & (curves.split == "temporal")
                     & (curves["index"] == k)].sort_values("N").r.to_numpy()
           for k in IX.CONFIRMATORY}
    ptB = {k: curves[(curves.design == "B") & (curves.split == "random")
                     & (curves["index"] == k)].sort_values("N").r.to_numpy()
           for k in IX.CONFIRMATORY}

    def cross_on(rvec, thr, grid):
        lg = np.log(grid)
        rv = np.asarray(rvec, float)
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

    for k in IX.CONFIRMATORY:
        d70 = np.array([min_n_cross(bootA[k][b], 0.70) for b in range(B_BOOT)])
        d80 = np.array([min_n_cross(bootA[k][b], 0.80) for b in range(B_BOOT)])
        d70T = np.array([min_n_cross(bootA_T[k][b], 0.70) for b in range(B_BOOT)])
        d80T = np.array([min_n_cross(bootA_T[k][b], 0.80) for b in range(B_BOOT)])
        r15 = np.array([interp15(bootA[k][b]) for b in range(B_BOOT)])
        r20 = bootA[k][:, GRID.index(20)]
        r20pt = pt[k][GRID.index(20)]
        minn[k] = {
            "minN70_point": min_n_cross(pt[k], 0.70),
            "minN80_point": min_n_cross(pt[k], 0.80),
            "minN70_ci": ci_or_partial(d70),
            "minN80_ci": ci_or_partial(d80),
            "share_minN70_gt200": float(np.mean(~np.isfinite(d70))),
            "share_minN80_gt200": float(np.mean(~np.isfinite(d80))),
            "minN70_temporal_point": min_n_cross(ptT[k], 0.70),
            "minN80_temporal_point": min_n_cross(ptT[k], 0.80),
            "minN70_temporal_ci": ci_or_partial(d70T),
            "minN80_temporal_ci": ci_or_partial(d80T),
            "minN70_designB_point": cross_on(ptB[k], 0.70, GRID_B),
            "minN80_designB_point": cross_on(ptB[k], 0.80, GRID_B),
            "prophecy_minN70_from20": float(
                20 * (0.70 / 0.30) * ((1 - r20pt) / r20pt)) if r20pt > 0 else None,
            "r15_point": float(interp15(pt[k])),
            "r15_direct_point": float(pt[k][GRID.index(15)]),
            "r15_ci": [float(v) for v in np.nanpercentile(r15, [2.5, 97.5])],
            "p_r15_ge_070": float(np.mean(r15 >= 0.70)),
            "r20_point": float(r20pt),
            "r20_ci": [float(v) for v in np.nanpercentile(r20, [2.5, 97.5])],
            "r20_upper_1sided95": float(np.nanpercentile(r20, 95)),
        }
    minn["flag_consistency_designA_bias"] = flags
    minn["_meta"] = dict(sample=sample, n_poolA=nA, grid=GRID, grid_B=GRID_B,
                         R=R_SPLITS, B=B_BOOT, seed=MASTER_SEED)
    with open(os.path.join("results", f"minN_{sample}.json"), "w") as fh:
        json.dump(minn, fh, indent=1)

    np.savez_compressed(
        os.path.join("data", "cache", f"designA_{sample}.npz"),
        **{f"boot_{k}": bootA[k] for k in IX.CONFIRMATORY},
        **{f"A20_{k}": cacheA[20][0][k] for k in IDX_KEYS},
        **{f"B20_{k}": cacheA[20][1][k] for k in IDX_KEYS},
        **{f"A100_{k}": cacheA[100][0][k] for k in IDX_KEYS},
        **{f"B100_{k}": cacheA[100][1][k] for k in IDX_KEYS},
        **{f"full20_{k}": cacheA[20][2][k] for k in IDX_KEYS},
        **{f"full100_{k}": cacheA[100][2][k] for k in IDX_KEYS},
        uid20=cacheA[20][3], uid100=cacheA[100][3])
    print(f"[{sample}] done: curves + minN + cache written; poolA n={nA}")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
