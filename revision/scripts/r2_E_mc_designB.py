"""Monte Carlo precision of the Design B estimates (revision round 2).

Design B includes at each N of the pre-specified grid {5, 10, 20, 50, 100,
200} every student with at least N rated answers (reliability.py). Its
crossings are point estimates (no joint bootstrap). This script:

  (1) reruns draw 0 with the submitted seeds default_rng([42, 6, N, sid]),
      keeps the 50 per-split values at every N and reports the split
      Monte Carlo SE of every value and crossing (delta method, checked by
      resampling the per-split values, B = 2,000, default_rng([42, 504, sid]));
  (2) reruns the same estimator on K = 50 further subsample draws,
      draw k using default_rng([42, 502, N, sid, k]); the SD over the K draws
      is the total Monte Carlo SD of a single-draw value.

Regression gate: draw 0 curves against the submitted Design B curves
(results/baseline; plan-rule AUROC against results/r1_A/prereg), <= 1e-12.

Usage: python -B r2_E_mc_designB.py t12   (then t34)
Writes: revision-2026-11/results/r2_E_parts/designB_<sample>.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_E_common as EC  # noqa: E402
import r2_A_rerun as RA  # noqa: E402

REL = RA.REL
GRID_B = list(REL.GRID_B)
assert GRID_B == EC.FROZEN_GRID
K_DRAWS = 50
B_RS = 2000
KEYS = ["bias", "brier", "gamma", "auroc", "auroc_pr", "conf", "acc"]
CROSS_KEYS = ["bias", "brier", "gamma", "auroc", "auroc_pr"]


def one_draw(arrs, counts, seedfn, keep_splits):
    users, offsets, codes, x = arrs
    G = len(GRID_B)
    rnd = {k: np.full(G, np.nan) for k in KEYS}
    splits = {k: [] for k in KEYS}
    npool = []
    for gi, N in enumerate(GRID_B):
        poolB = counts >= N
        npool.append(int(poolB.sum()))
        if poolB.sum() < REL.MIN_POOL:
            for k in KEYS:
                splits[k].append(np.full(REL.R_SPLITS, np.nan))
            continue
        rng = np.random.default_rng(seedfn(N))
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolB, N, rng)
        A, Bh, TA, TB, full = RA.halves(C, X, rng)
        for k in KEYS:
            cols = REL.sb_corr_cols(A[k], Bh[k])
            rnd[k][gi] = float(np.nanmean(cols))
            if keep_splits:
                splits[k].append(cols)
    return rnd, (splits if keep_splits else None), npool


def cd(v):
    return {str(N): (None if np.isnan(v[i]) else float(v[i])) for i, N in enumerate(GRID_B)}


def run(sample):
    t0 = time.time()
    sid = EC.SIDS[sample]
    lab = EC.LAB[sample]
    df = RA.common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    arrs = (users, offsets, codes, x)
    base = pd.read_csv(os.path.join(EC.RESULTS, "baseline", "results", f"reliability_curves_{sample}.csv"))
    pre = pd.read_csv(os.path.join(EC.RESULTS, "r1_A", "prereg", "results", f"reliability_curves_{sample}.csv"))

    r0, s0, npool = one_draw(arrs, counts, lambda N: [EC.SEED, 6, N, sid], True)
    reg = {}
    for k in KEYS:
        src, kk = (pre, "auroc") if k == "auroc_pr" else (base, k)
        ref = []
        for N in GRID_B:
            row = src[(src.design == "B") & (src.split == "random") & (src["index"] == kk) & (src.N == N)]
            ref.append(float(row.r.iloc[0]) if len(row) else np.nan)
        reg[f"curve_{k}"] = float(np.nanmax(np.abs(r0[k] - np.array(ref))))
    reg_pass = bool(max(reg.values()) <= 1e-12)
    print(f"[{sample}] B draw 0 {time.time() - t0:.0f}s regression pass {reg_pass}", flush=True)
    if not reg_pass:
        print(json.dumps(reg, indent=1))
        raise SystemExit("draw 0 does not reproduce the submitted Design B curves")

    rs_rng = np.random.default_rng([EC.SEED, 504, sid])
    split_part = {}
    for k in KEYS:
        se = np.array([EC.mc_se(v) if np.any(~np.isnan(v)) else np.nan for v in s0[k]], dtype=float)
        blk = {"per_N_split_se": cd(se), "max_split_se_all_N": float(np.nanmax(se)),
               "argmax_N": int(GRID_B[int(np.nanargmax(se))])}
        if k in CROSS_KEYS:
            rs, _ = EC.resample_split_crossings(s0[k], EC.THRS, rs_rng, B_RS, GRID_B)
            cr = {}
            for thr in EC.THRS:
                t = f"{int(round(thr * 100))}"
                cr[t] = {"point": EC.fnum(EC.cross(r0[k], thr, GRID_B)),
                         "split_se_delta": EC.delta_se_cross(r0[k], se, thr, GRID_B)}
                cr[t].update(rs[t])
            blk["crossings"] = cr
        split_part[k] = blk

    draws = []
    for kk in range(1, K_DRAWS + 1):
        rk, _, _ = one_draw(arrs, counts, lambda N, q=kk: [EC.SEED, 502, N, sid, q], False)
        draws.append(rk)
        if kk % 10 == 0:
            print(f"[{sample}] B draw {kk}/{K_DRAWS} {time.time() - t0:.0f}s", flush=True)
    curves = {}
    crossings = {}
    for k in KEYS:
        M = np.array([d[k] for d in draws])
        curves[k] = {"draw0": cd(r0[k]), "mean_over_draws": cd(M.mean(axis=0)),
                     "sd_over_draws": cd(M.std(axis=0, ddof=1))}
        if k in CROSS_KEYS:
            for thr in EC.THRS:
                t = f"{int(round(thr * 100))}"
                dv = np.array([EC.cross(row, thr, GRID_B) for row in M])
                blk = EC.draws_summary(dv, extra_gt=EC.S2_MEDIAN_VOLUME if (sample == "t34" and thr == 0.70) else None)
                blk["draw0"] = EC.fnum(EC.cross(r0[k], thr, GRID_B))
                blk["crossing_of_mean_curve"] = EC.fnum(EC.cross(M.mean(axis=0), thr, GRID_B))
                crossings[f"{k}_{t}"] = blk
    out = {"_definitions": {
               "split_error": "per index, printed draw (draw 0): per-N Monte Carlo SE of the mean over the 50 random splits, its maximum, and per crossing the point value, the delta-method split SE and the SD, share reaching and 2.5th / 97.5th percentiles over 2,000 resamples of the per-split values within each N",
               "draw_curves": "per index: draw 0 value, mean and SD over the K further subsample draws (random halves)",
               "draw_crossings": "<index>_<level>: over the K further draws, median and 2.5th / 97.5th percentiles (a draw that never reaches the level within 200 ranked above 200), mean and SD only when every draw reaches the level, share reaching within 200, draw 0 value, crossing of the mean curve, share above 31 (S2, .70)",
               "n_pool_by_N": "Design B pool at each N: students with at least N rated answers",
               "regression_max_abs_diff": "draw 0 against the submitted Design B curves (results/baseline; plan-rule AUROC against results/r1_A/prereg)",
               "index names": "bias, brier, gamma, auroc (tie convention), auroc_pr (plan rule), conf, acc"},
           "sample": lab, "grid": GRID_B, "n_pool_by_N": dict(zip(map(str, GRID_B), npool)),
           "R_splits": REL.R_SPLITS, "K_further_draws": K_DRAWS, "B_split_resamples": B_RS,
           "seeds": {"draw0_subsample_and_splits": "[42, 6, N, sid] (as submitted)",
                     "further_draw_k_subsample_and_splits": "[42, 502, N, sid, k], k = 1..K",
                     "split_resampling": "[42, 504, sid]", "sid": sid},
           "regression_max_abs_diff": reg, "regression_pass": reg_pass,
           "split_error": split_part, "draw_curves": curves, "draw_crossings": crossings,
           "runtime_seconds": round(time.time() - t0, 1)}
    with open(os.path.join(EC.PARTS, f"designB_{sample}.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"[{sample}] B written {out['runtime_seconds']}s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
