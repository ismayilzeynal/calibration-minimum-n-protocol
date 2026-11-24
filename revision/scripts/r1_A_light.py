"""Stream A (E5) light pass, no bootstrap:

 1. Full-set incomputability for Design A and Design B at every grid N, under
    both AUROC rules. The subsamples are re-drawn with the exact seeds of
    reliability.py ([42, 5, N, sid] for A, [42, 6, N, sid] for B); the
    subsample is the first use of each generator, so the full N-response sets
    are identical to the submitted run. Design A values are checked against
    the baseline curves CSV (pct_undefined_full).
 2. Persistence of the constant-confidence state across split halves at
    N = 20 and N = 100 (Design A, baseline cache): among student x split
    cells where AUROC is defined in both halves, how often a half is
    "forced" (constant confidence, AUROC = .5), and the phi correlation
    between forced-in-half-A and forced-in-half-B.
 3. Tagged-response count distribution of reporters (share with >= k tagged
    answers, k = 1..400) for threshold coverage.

Writes revision-2026-11/results/r1_A/light_{t12,t34}.json.
Usage: python -B r1_A_light.py t12|t34
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUTDIR = os.path.join(REV, "results", "r1_A")
BASE = os.path.join(REV, "results", "baseline")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import indices as IX  # noqa: E402
import reliability as RL  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")


def shares(full):
    g_nan = np.isnan(full["gamma"])
    a_nan = np.isnan(full["auroc"])
    forced = g_nan & ~a_nan
    return {"n_students": int(len(g_nan)),
            "gamma_undefined": float(g_nan.mean()),
            "auroc_undefined_submitted_rule": float(a_nan.mean()),
            "auroc_undefined_prereg_rule": float((g_nan | a_nan).mean()),
            "auroc_forced_05_share": float(forced.mean()),
            "check_gamma_nan_equals_auroc_nan_or_forced": bool(
                np.array_equal(g_nan, a_nan | forced)),
            "check_no_auroc_nan_where_gamma_defined": bool(
                not np.any(a_nan & ~g_nan))}


def persistence(cache, N):
    gA, gB = cache[f"A{N}_gamma"], cache[f"B{N}_gamma"]
    aA, aB = cache[f"A{N}_auroc"], cache[f"B{N}_auroc"]
    fA = np.isnan(gA) & ~np.isnan(aA)
    fB = np.isnan(gB) & ~np.isnan(aB)
    both = ~np.isnan(aA) & ~np.isnan(aB)
    phis, pa, pb_given_a = [], [], []
    for s in range(aA.shape[1]):
        m = both[:, s]
        u, v = fA[m, s].astype(float), fB[m, s].astype(float)
        if u.std() > 0 and v.std() > 0:
            phis.append(float(np.corrcoef(u, v)[0, 1]))
        pa.append(float(u.mean()))
        if u.sum() > 0:
            pb_given_a.append(float(v[u == 1].mean()))
    nonforced = ~np.isnan(aA) & ~fA
    vals = aA[nonforced]
    return {"N": N,
            "cells_auroc_defined_both_halves": int(both.sum()),
            "share_cells_any_half_forced": float((both & (fA | fB)).sum() / both.sum()),
            "p_forced_halfA_among_defined_both": float(np.mean(pa)),
            "p_forced_halfB_given_forced_halfA": float(np.mean(pb_given_a)),
            "phi_forcedA_forcedB_mean_over_splits": float(np.mean(phis)),
            "halfA_auroc_nonforced_mean": float(vals.mean()),
            "halfA_auroc_nonforced_sd": float(vals.std(ddof=1)),
            "halfA_auroc_nonforced_share_exactly_05": float(np.mean(vals == 0.5))}


def main(sample):
    t0 = time.time()
    sid = {"t12": 1, "t34": 2}[sample]
    df = common.build_responses(sample)
    users, offsets, codes, x = RL.student_arrays(df)
    counts = np.diff(offsets)
    out = {"sample": sample}

    curves = pd.read_csv(os.path.join(BASE, "results",
                                      f"reliability_curves_{sample}.csv"))
    # Design A full sets
    poolA = counts >= max(RL.GRID)
    blkA, maxdiff = {}, 0.0
    for N in RL.GRID:
        rngN = np.random.default_rng([RL.MASTER_SEED, 5, N, sid])
        C, X, _ = RL.subsample_matrix(users, offsets, codes, x, poolA, N, rngN)
        full = IX.indices_from_counts(IX.count_tables(C, X, np.ones(C.shape, bool)))
        blkA[str(N)] = shares(full)
        for k, key in (("gamma", "gamma_undefined"),
                       ("auroc", "auroc_undefined_submitted_rule")):
            ref = curves[(curves.design == "A") & (curves.split == "random")
                         & (curves["index"] == k) & (curves.N == N)]
            maxdiff = max(maxdiff, abs(float(ref.pct_undefined_full.iloc[0])
                                       - blkA[str(N)][key]))
    out["designA_full"] = blkA
    out["designA_full_check_vs_baseline_csv_max_abs_diff"] = maxdiff

    # Design B full sets
    blkB = {}
    for N in RL.GRID_B:
        poolB = counts >= N
        if int(poolB.sum()) < RL.MIN_POOL:
            continue
        rngN = np.random.default_rng([RL.MASTER_SEED, 6, N, sid])
        C, X, _ = RL.subsample_matrix(users, offsets, codes, x, poolB, N, rngN)
        full = IX.indices_from_counts(IX.count_tables(C, X, np.ones(C.shape, bool)))
        blkB[str(N)] = shares(full)
    out["designB_full"] = blkB

    # persistence of the constant-confidence state (baseline cache)
    cache = np.load(os.path.join(BASE, "data", "cache", f"designA_{sample}.npz"))
    out["forced_persistence_designA"] = {str(N): persistence(cache, N)
                                         for N in (20, 100)}

    # reporter coverage
    out["reporters"] = int(len(counts))
    out["share_reporters_ge_k"] = {str(k): float((counts >= k).mean())
                                   for k in range(1, 401)}
    out["wall_seconds"] = time.time() - t0
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, f"light_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("share_reporters_ge_k",)}, indent=1))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "t12")
