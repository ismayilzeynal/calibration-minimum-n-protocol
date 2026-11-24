"""Revision, gate fix round 1 (post-review, exploratory additions): flag
stability at N = 20 on the disjoint-set design of stream D.

Same pool (Design A, >= 200 rated answers), same disjoint N-answer sets
(seed [42, 92, N, sid], RNG call order of rev_round2.disjoint_flags) and the
same joint student bootstrap ([42, 77, sid]) as r1_D_construct_flags.py.

  (1) matched prevalence: each absolute rule (+.15, +.20, 95% interval above
      zero) against a rank rule that flags the same share of students (top k
      by bias, k = round(share flagged by the absolute rule x n)); paired
      bootstrap of the kappa difference;
  (2) guess-corrected coding c' = .25 + .75 c: share flagged, retention and
      kappa of the +.15 rule (exact integer comparison);
  (3) pool sensitivity: the +.15 rule and the top decile without students
      whose rated answers have accuracy >= .99 or one constant rating
      (columns removed after the sets are drawn).

Usage: python -B r2_D_flags.py
Writes: revision-2026-11/results/r2_D.json
"""
import json
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import r1_D_construct_flags as D  # noqa: E402

N = 20
D.NGRID = [N]
B = D.B
SIDS = {"t12": 1, "t34": 2}


def run(sample):
    t0 = time.time()
    sid = SIDS[sample]
    df = D.common.build_responses(sample)
    users, offsets, codes, x = D.REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    pool = np.where(counts >= 200)[0]
    n = len(pool)
    acc_full = np.array([x[offsets[i]:offsets[i + 1]].mean() for i in pool])
    const_full = np.array([codes[offsets[i]:offsets[i + 1]].min() == codes[offsets[i]:offsets[i + 1]].max()
                           for i in pool])
    keep = ~((acc_full >= 0.99) | const_full)

    ds = D.disjoint_sets(offsets, codes, x, pool, sid)[N]
    cA, aA, fA, _ = D.derive(N, ds["sA"])
    cB, aB, fB, _ = D.derive(N, ds["sB"])
    k10 = max(int(round(n * 0.10)), 1)
    fA["decile"] = D.top_mask(ds["bA"], k10)
    fB["decile"] = D.top_mask(ds["bB"], k10)

    # regression against stream D
    ref = json.load(open(os.path.join(REV, "tmp", f"r1_D_{sample}.json")))["per_N"][str(N)]["flags"]
    gate = 0.0
    for rl in ("decile", "thr15", "thr20", "lcb_z"):
        m = D.flag_metrics(fA[rl], fB[rl])
        for kk in ("share_flagged", "retention", "kappa"):
            gate = max(gate, abs(m[kk] - ref[rl][kk]))

    # (1) matched prevalence
    matched = {}
    brng = np.random.default_rng([D.SEED, 77, sid])
    IDX = [brng.integers(0, n, n) for _ in range(B)]
    for rl in ("thr15", "thr20", "lcb_z"):
        share = D.flag_metrics(fA[rl], fB[rl])["share_flagged"]
        k = max(int(round(share * n)), 1)
        tA = D.top_mask(ds["bA"], k)
        tB = D.top_mask(ds["bB"], k)
        ma = D.flag_metrics(fA[rl], fB[rl])
        mt = D.flag_metrics(tA, tB)
        dk = np.empty(B)
        dr = np.empty(B)
        for b, idx in enumerate(IDX):
            m1 = D.flag_metrics(fA[rl][:, idx], fB[rl][:, idx])
            m2 = D.flag_metrics(D.top_mask(ds["bA"][:, idx], k), D.top_mask(ds["bB"][:, idx], k))
            dk[b] = m1["kappa"] - m2["kappa"]
            dr[b] = m1["retention"] - m2["retention"]
        matched[rl] = {"share_flagged_absolute": ma["share_flagged"], "k_rank": k,
                       "share_flagged_rank": mt["share_flagged"],
                       "absolute": {kk: ma[kk] for kk in ("retention", "kappa", "chance_retention")},
                       "rank_matched": {kk: mt[kk] for kk in ("retention", "kappa", "chance_retention")},
                       "kappa_absolute_minus_rank": {"point": ma["kappa"] - mt["kappa"],
                                                     "ci": D.ci(dk)},
                       "retention_absolute_minus_rank": {"point": ma["retention"] - mt["retention"],
                                                         "ci": D.ci(dr)}}

    # (2) guess-corrected coding: bias' = .25 + .75 conf - acc; +.15 rule
    #     1600 N bias' = 400 N + 300 sc - 1600 sx  ->  bias' > .15  <=>  15 sc - 80 sx > -8 N
    def g15(s):
        return 15 * s[..., 0] - 80 * s[..., 1] > -8 * N
    gA, gB = g15(ds["sA"]), g15(ds["sB"])
    mg = D.flag_metrics(gA, gB)
    ms = D.flag_metrics(fA["thr15"], fB["thr15"])
    bias_g_mean = float(np.mean(((0.25 + 0.75 * cA - aA).mean() + (0.25 + 0.75 * cB - aB).mean()) / 2))
    guess = {"thr15_linear": {kk: ms[kk] for kk in ("share_flagged", "retention", "kappa", "chance_retention")},
             "thr15_guess": {kk: mg[kk] for kk in ("share_flagged", "retention", "kappa", "chance_retention")},
             "share_students_flag_differs_between_codings": float(
                 np.mean(np.concatenate([(gA != fA["thr15"]).mean(1), (gB != fB["thr15"]).mean(1)]))),
             "mean_bias_linear": float(np.mean(((cA - aA).mean() + (cB - aB).mean()) / 2)),
             "mean_bias_guess": bias_g_mean}

    # (3) pool sensitivity
    kk_ = np.where(keep)[0]
    nk = len(kk_)
    k10k = max(int(round(nk * 0.10)), 1)
    exc = {"n_kept": int(nk)}
    m15 = D.flag_metrics(fA["thr15"][:, kk_], fB["thr15"][:, kk_])
    mdec = D.flag_metrics(D.top_mask(ds["bA"][:, kk_], k10k), D.top_mask(ds["bB"][:, kk_], k10k))
    exc["thr15"] = {kk: m15[kk] for kk in ("share_flagged", "retention", "kappa", "chance_retention")}
    exc["decile"] = {kk: mdec[kk] for kk in ("share_flagged", "retention", "kappa", "chance_retention")}
    exc["full_pool_thr15"] = {kk: ms[kk] for kk in ("share_flagged", "retention", "kappa")}
    out = {"sample": sample, "N": N, "n_pool": int(n), "regression_vs_stream_D_max_abs_diff": gate,
           "matched_prevalence": matched, "guess_coding": guess, "pool_sensitivity": exc,
           "runtime_seconds": round(time.time() - t0, 1)}
    print(sample, "done", out["runtime_seconds"], "gate", gate, flush=True)
    return out


if __name__ == "__main__":
    res = {"_definitions": {
        "retention": "Dice overlap of the two flag sets (share of flagged students flagged again), as in stream D",
        "matched_prevalence": "rank rule flags the top k students by bias, k = round(share flagged by the absolute rule x n); same sets and bootstrap draws",
        "guess_coding": "c' = .25 + .75 c; +.15 rule evaluated by exact integer comparison",
        "pool_sensitivity": "students with rated accuracy >= .99 or one constant rating on all their rated answers removed after the sets were drawn"}}
    for s, lab in (("t12", "S1"), ("t34", "S2")):
        res[lab] = run(s)
    with open(os.path.join(REV, "results", "r2_D.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print("written results/r2_D.json")
