"""H3: empirical-Bayes shrinkage vs raw bias for low-N person-ranking stability.

Ecological setting per the frozen plan: all students with >= 10 tagged
responses, variable n_i, random half-split per student, EB fit within each
half. Metric: Spearman correlation of person estimates across halves
(ordinal ranks; ties broken by stable order under the fixed seed).

Run:  python scripts/h3_shrinkage.py t12
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from reliability import student_arrays

R_SPLITS = 50
B_BOOT = 1000
MASTER_SEED = 42
MIN_TAGGED = 10


def eb_shrink(bias, s2, n):
    """Normal-normal EB toward the grand mean; tau2 by method of moments."""
    samp_var = s2 / n
    tau2 = max(float(np.var(bias, ddof=1) - np.mean(samp_var)), 0.0)
    w = tau2 / (tau2 + samp_var) if tau2 > 0 else np.zeros_like(samp_var)
    gm = float(np.mean(bias))
    return w * bias + (1.0 - w) * gm


def half_bias(d_list, rng):
    """Per-student raw + EB bias for both halves of one random split."""
    n = len(d_list)
    bA = np.empty(n); bB = np.empty(n)
    sA = np.empty(n); sB = np.empty(n)
    nA = np.empty(n); nB = np.empty(n)
    for i, d in enumerate(d_list):
        perm = rng.permutation(len(d))
        h = len(d) // 2
        a, b = d[perm[:h]], d[perm[h:]]
        bA[i], bB[i] = a.mean(), b.mean()
        sA[i], sB[i] = a.var(ddof=1), b.var(ddof=1)
        nA[i], nB[i] = len(a), len(b)
    ebA = eb_shrink(bA, sA, nA)
    ebB = eb_shrink(bB, sB, nB)
    return bA, bB, ebA, ebB


def spearman_cols(A, B):
    """Column-wise Spearman via ordinal ranks (A,B same shape (n,S))."""
    rA = np.argsort(np.argsort(A, axis=0), axis=0).astype(float)
    rB = np.argsort(np.argsort(B, axis=0), axis=0).astype(float)
    rA -= rA.mean(axis=0); rB -= rB.mean(axis=0)
    num = (rA * rB).sum(axis=0)
    den = np.sqrt((rA ** 2).sum(axis=0) * (rB ** 2).sum(axis=0))
    return num / den


def run(sample):
    rng = np.random.default_rng([MASTER_SEED, 3, {"t12": 1, "t34": 2}[sample]])
    df = common.build_responses(sample)
    users, offsets, codes, x = student_arrays(df)
    counts = np.diff(offsets)
    keep = np.where(counts >= MIN_TAGGED)[0]
    d_all = codes.astype(np.float64) * 0.25 - x  # c - x per response
    d_list = [d_all[offsets[i]:offsets[i + 1]] for i in keep]
    n = len(d_list)

    rawA = np.empty((n, R_SPLITS)); rawB = np.empty((n, R_SPLITS))
    ebA = np.empty((n, R_SPLITS)); ebB = np.empty((n, R_SPLITS))
    jac_raw = []; jac_eb = []
    k10 = max(int(round(n * 0.10)), 1)
    for s in range(R_SPLITS):
        a, b, ea, eb = half_bias(d_list, rng)
        rawA[:, s], rawB[:, s], ebA[:, s], ebB[:, s] = a, b, ea, eb
        tA = set(np.argsort(-a)[:k10]); tB = set(np.argsort(-b)[:k10])
        jac_raw.append(len(tA & tB) / len(tA | tB))
        tA = set(np.argsort(-ea)[:k10]); tB = set(np.argsort(-eb)[:k10])
        jac_eb.append(len(tA & tB) / len(tA | tB))

    rho_raw = spearman_cols(rawA, rawB)
    rho_eb = spearman_cols(ebA, ebB)
    point_raw = float(rho_raw.mean()); point_eb = float(rho_eb.mean())
    # artifact-immune validation (addendum 4b): EB estimate from one half must
    # predict the OTHER half's raw score; symmetrized over the two directions
    rho_eb_val = float(0.5 * (spearman_cols(ebA, rawB)
                              + spearman_cols(ebB, rawA)).mean())

    # flag validation vs raw other-half criterion
    k10 = max(int(round(n * 0.10)), 1)
    jv_raw, jv_eb = [], []
    for s in range(R_SPLITS):
        critB = set(np.argsort(-rawB[:, s])[:k10])
        tRaw = set(np.argsort(-rawA[:, s])[:k10])
        tEb = set(np.argsort(-ebA[:, s])[:k10])
        jv_raw.append(len(tRaw & critB) / len(tRaw | critB))
        jv_eb.append(len(tEb & critB) / len(tEb | critB))

    brng = np.random.default_rng([MASTER_SEED, 79, {"t12": 1, "t34": 2}[sample]])
    diffs = np.empty(B_BOOT)
    diffs_val = np.empty(B_BOOT)
    for bdraw in range(B_BOOT):
        idx = brng.integers(0, n, n)
        dr = spearman_cols(rawA[idx], rawB[idx]).mean()
        de = spearman_cols(ebA[idx], ebB[idx]).mean()
        dv = 0.5 * (spearman_cols(ebA[idx], rawB[idx])
                    + spearman_cols(ebB[idx], rawA[idx])).mean()
        diffs[bdraw] = de - dr
        diffs_val[bdraw] = dv - dr
    out = {
        "sample": sample, "n_students": n, "min_tagged": MIN_TAGGED,
        "rho_raw": point_raw, "rho_eb": point_eb,
        "diff_point": point_eb - point_raw,
        "diff_ci95": [float(v) for v in np.percentile(diffs, [2.5, 97.5])],
        "diff_lower_1sided95": float(np.percentile(diffs, 5)),
        "p_onesided_diff_le_0": float(np.mean(diffs <= 0)),
        "rho_eb_validation": rho_eb_val,
        "diff_validation_point": rho_eb_val - point_raw,
        "diff_validation_ci95": [float(v) for v in
                                 np.percentile(diffs_val, [2.5, 97.5])],
        "diff_validation_lower_1sided95": float(np.percentile(diffs_val, 5)),
        "p_onesided_validation_le_0": float(np.mean(diffs_val <= 0)),
        "jaccard_top_decile_raw": float(np.mean(jac_raw)),
        "jaccard_top_decile_eb": float(np.mean(jac_eb)),
        "jaccard_validation_raw": float(np.mean(jv_raw)),
        "jaccard_validation_eb": float(np.mean(jv_eb)),
        "R": R_SPLITS, "B": B_BOOT, "seed": MASTER_SEED,
    }
    os.makedirs("results", exist_ok=True)
    with open(os.path.join("results", f"h3_shrinkage_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(sample, "rho_raw", round(point_raw, 4), "rho_eb", round(point_eb, 4),
          "diff", round(point_eb - point_raw, 4),
          "p<=0", out["p_onesided_diff_le_0"])


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
