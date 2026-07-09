"""Addenda 3 and 6: temporal-gap descriptive + bias-variance decomposition.

Run AFTER the reliability rerun (cache must include conf):
    python scripts/rev_extras.py t12
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from reliability import GRID, MASTER_SEED, sb_corr_cols

R_SPLITS = 50


def temporal_gap(sample, rng):
    """Median elapsed days between median dates of temporal halves, per N."""
    df = common.build_responses(sample)
    tg = common.tagged(df).sort_values(["UserId", "DateAnswered"], kind="stable")
    day = (tg["DateAnswered"] - tg["DateAnswered"].min()).dt.days.to_numpy()
    uid = tg["UserId"].to_numpy()
    users, starts = np.unique(uid, return_index=True)
    offsets = np.append(starts, len(uid))
    counts = np.diff(offsets)
    pool = np.where(counts >= max(GRID))[0]
    out = {}
    for N in GRID:
        gaps = []
        h = N // 2
        for i in pool:
            s, e = offsets[i], offsets[i + 1]
            sel = np.sort(rng.choice(e - s, size=N, replace=False)) + s
            d = day[sel]
            gaps.append(np.median(d[h:]) - np.median(d[:h]))
        out[str(N)] = {"median_gap_days": float(np.median(gaps)),
                       "iqr": [float(np.percentile(gaps, 25)),
                               float(np.percentile(gaps, 75))]}
    return out


def decomposition(sample):
    """What is the reliable bias variance made of? Cross-half correlations of
    bias with conf and acc, and the partialed-bias vs conf correlation."""
    cache = np.load(os.path.join("data", "cache", f"designA_{sample}.npz"))
    out = {}
    for N in (20, 100):
        bA, bB = cache[f"A{N}_bias"], cache[f"B{N}_bias"]
        cA, cB = cache[f"A{N}_conf"], cache[f"B{N}_conf"]
        aA, aB = cache[f"A{N}_acc"], cache[f"B{N}_acc"]
        blk = {
            "sb_bias": float(np.nanmean(sb_corr_cols(bA, bB))),
            "sb_conf": float(np.nanmean(sb_corr_cols(cA, cB))),
            "sb_acc": float(np.nanmean(sb_corr_cols(aA, aB))),
            # raw cross-half correlations (no SB): what predicts bias in the
            # other half
            "r_biasA_biasB": float(np.nanmean(_xcorr(bA, bB))),
            "r_biasA_confB": float(np.nanmean(_xcorr(bA, cB))),
            "r_biasA_accB": float(np.nanmean(_xcorr(bA, aB))),
        }
        # accuracy-partialed bias half-scores vs mean confidence (same half)
        rs = []
        for s in range(bA.shape[1]):
            b, a, c = bA[:, s], aA[:, s], cA[:, s]
            m = ~(np.isnan(b) | np.isnan(a) | np.isnan(c))
            beta = np.cov(a[m], b[m])[0, 1] / np.var(a[m])
            resid = b[m] - beta * (a[m] - a[m].mean())
            rs.append(np.corrcoef(resid, c[m])[0, 1])
        blk["r_partialedBias_conf_samehalf"] = float(np.mean(rs))
        out[f"N{N}"] = blk
    return out


def _xcorr(A, B):
    M = ~(np.isnan(A) | np.isnan(B))
    Az = np.where(M, A, 0.0); Bz = np.where(M, B, 0.0)
    n = M.sum(axis=0).astype(float)
    sA, sB = Az.sum(0), Bz.sum(0)
    cov = (Az * Bz).sum(0) - sA * sB / n
    vA = (Az * Az).sum(0) - sA ** 2 / n
    vB = (Bz * Bz).sum(0) - sB ** 2 / n
    return cov / np.sqrt(vA * vB)


def run(sample):
    rng = np.random.default_rng([MASTER_SEED, 91, {"t12": 1, "t34": 2}[sample]])
    out = {"sample": sample,
           "temporal_gap_days": temporal_gap(sample, rng),
           "bias_decomposition": decomposition(sample)}
    with open(os.path.join("results", f"rev_extras_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(sample, "rev_extras written")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
