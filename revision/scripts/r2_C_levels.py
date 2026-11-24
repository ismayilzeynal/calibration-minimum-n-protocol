"""Revision, gate fix round 1 (post-review, exploratory addition): the LEVEL
of signed bias by item difficulty, content strand and approximate school term,
and how often an absolute +.15 flag changes between halves formed by
difficulty or by term.

(1) Answer-level means over rated answers (mean confidence, accuracy, bias =
    their difference) by difficulty quintile (global quintiles of the
    leave-one-out item proportion correct, as in stream C), by level-1 strand
    and by approximate school term; for all rated answers of the sample and
    for the Design A pool.
(2) Design A pool, N = 20, the subsamples and halves of stream C
    (r1_C_splits.py; subsample and random halves [42, 5, N, sid], term halves
    [42, 332, N, sid], harder-versus-easier halves [42, 333, N, sid]):
    per student, the bias of each 10-answer half, the share of students whose
    +.15 flag differs between the two halves (mean over the 50 repetitions),
    and the mean harder-minus-easier bias difference.

Usage: python -B r2_C_levels.py
Writes: revision-2026-11/results/r2_C.json
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
sys.path.insert(0, HERE)
import r1_C_common as RC  # noqa: E402
import r1_C_splits as SP  # noqa: E402
import reliability as RL  # noqa: E402

N = 20
MAIN_STRANDS = ["Number", "Algebra", "Geometry and Measure", "Data and Statistics"]


def level_table(t, key):
    g = t.groupby(key, observed=True)
    out = {}
    for lab, sub in g:
        out[str(lab)] = {"n_answers": int(len(sub)),
                         "share_answers": float(len(sub) / len(t)),
                         "mean_confidence": float(sub.c.mean()),
                         "accuracy": float(sub.IsCorrect.mean()),
                         "mean_bias": float(sub.d.mean())}
    return out


def levels(tg):
    t = tg.copy()
    t["c"] = (t["Confidence"].to_numpy() / 25.0).round() * 0.25
    t["d"] = t["c"] - t["IsCorrect"].astype(float)
    t["strand_main"] = np.where(t["strand"].astype(str).isin(MAIN_STRANDS),
                                t["strand"].astype(str), "Combined or other")
    t["dq"] = np.where(t["dstrat"] < RC.N_DIFF_STRATA, t["dstrat"], -1)
    t = t[t.dq >= 0]
    return {"overall": {"n_answers": int(len(t)), "mean_bias": float(t.d.mean())},
            "difficulty_quintile": level_table(t, "dq"),
            "strand": level_table(t, "strand_main"),
            "term": level_table(t, "term")}


def run(sample):
    t0 = time.time()
    sid = RC.SID[sample]
    df = RC.load(sample)
    tg = df[df.tag].sort_values(["UserId", "DateAnswered"], kind="stable")
    codes = (tg["Confidence"].to_numpy() / 25.0).round().astype(np.int8)
    x = tg["IsCorrect"].to_numpy().astype(np.int8)
    uid = tg["UserId"].to_numpy()
    users, starts = np.unique(uid, return_index=True)
    offsets = np.append(starts, len(uid))
    counts = np.diff(offsets)
    pool = np.where(counts >= max(SP.GRID))[0]
    nA = len(pool)
    pool_users = set(users[pool].tolist())

    lev = {"all_rated": levels(tg),
           "designA_pool": levels(tg[tg.UserId.isin(pool_users)]),
           "difficulty_quintile_edges_p_loo": df.attrs["dstrat_edges"]}

    term_c, _ = pd.factorize(tg["term"].astype(str))
    p_t = tg["p_loo"].to_numpy().astype(float)
    p_t = np.where(np.isnan(p_t), np.nanmean(df["p_loo"].to_numpy()), p_t)
    del df

    rngN = np.random.default_rng([RC.MASTER_SEED, 5, N, sid])
    S = SP.subsample_idx(offsets, pool, N, rngN)
    C, X = codes[S], x[S]
    mR = SP.random_masks(nA, N, rngN)
    mK = SP.keyed_masks(term_c[S], np.random.default_rng([RC.MASTER_SEED, 332, N, sid]))
    mD = SP.ordered_masks(p_t[S], np.random.default_rng([RC.MASTER_SEED, 333, N, sid]))

    ref = json.load(open(os.path.join(RC.OUT_DIR, f"r1_C_splits_{sample}.json")))["curves"]
    i20 = SP.GRID.index(N)
    C64, X64 = C.astype(np.int64), X.astype(np.int64)

    def half_stats(masks):
        flag_diff, dbias, adbias = [], [], []
        A_b = np.empty((nA, len(masks)))
        B_b = np.empty((nA, len(masks)))
        for j, m in enumerate(masks):
            nh_a, nh_b = m.sum(1), (~m).sum(1)
            sc_a, sc_b = (C64 * m).sum(1), (C64 * ~m).sum(1)
            sx_a, sx_b = (X64 * m).sum(1), (X64 * ~m).sum(1)
            ba = (0.25 * sc_a - sx_a) / nh_a
            bb = (0.25 * sc_b - sx_b) / nh_b
            fa = 25 * sc_a - 100 * sx_a > 15 * nh_a
            fb = 25 * sc_b - 100 * sx_b > 15 * nh_b
            flag_diff.append((fa != fb).mean())
            dbias.append(ba - bb)
            adbias.append(np.abs(ba - bb))
            A_b[:, j], B_b[:, j] = ba, bb
        dmean = np.mean(dbias, axis=0)
        return {"share_students_flag15_differs": float(np.mean(flag_diff)),
                "mean_bias_halfA_minus_halfB": float(dmean.mean()),
                "sd_over_students_bias_halfA_minus_halfB": float(dmean.std(ddof=1)),
                "mean_abs_bias_difference": float(np.mean(adbias)),
                "share_flagged_halfA": float(np.mean([(25 * (C64 * m).sum(1) - 100 * (X64 * m).sum(1) > 15 * m.sum(1)).mean() for m in masks])),
                "share_flagged_halfB": float(np.mean([(25 * (C64 * ~m).sum(1) - 100 * (X64 * ~m).sum(1) > 15 * (~m).sum(1)).mean() for m in masks])),
                "bias_r_SB": RL.sb_mean(A_b, B_b)}

    within = {"random": half_stats(mR), "term": half_stats(mK), "difficulty_harder_A": half_stats(mD)}
    gate = max(abs(within["random"]["bias_r_SB"] - ref["random"]["bias"][i20]),
               abs(within["term"]["bias_r_SB"] - ref["term"]["bias"][i20]),
               abs(within["difficulty_harder_A"]["bias_r_SB"] - ref["difficulty"]["bias"][i20]))
    out = {"sample": sample, "n_poolA": int(nA), "N": N, "levels": lev, "within_student_halves_N20": within,
           "regression_bias_r20_vs_stream_C_max_abs_diff": float(gate),
           "runtime_seconds": round(time.time() - t0, 1)}
    print(sample, "done", out["runtime_seconds"], "gate", gate, flush=True)
    return out


if __name__ == "__main__":
    res = {"_definitions": {
        "levels": "answer-level means over rated answers; bias = mean confidence (0-1) minus accuracy; difficulty quintile 0 = hardest (lowest leave-one-out item proportion correct)",
        "within_student_halves_N20": "Design A pool, 20-answer subsample, two 10-answer halves; share of students whose +.15 flag (bias above +.15) differs between halves, mean over 50 repetitions; difficulty_harder_A: half A = harder half",
        "term": "approximate English school term from the calendar month (Autumn Sep-Dec, Spring Jan-Mar, Summer Apr-Aug)"}}
    for s, lab in (("t12", "S1"), ("t34", "S2")):
        res[lab] = run(s)
    with open(os.path.join(REV, "results", "r2_C.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("written results/r2_C.json")
