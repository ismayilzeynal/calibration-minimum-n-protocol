"""Stream A (E5): establish exactly how the submitted indices.py treats
constant-confidence and constant-correctness response sets, and verify that
its AUROC is the standard Mann-Whitney estimator with 0.5 tie credit and its
gamma is Goodman-Kruskal gamma.

Checks:
 1. hand-built count tables for the edge cases;
 2. 5,000 random response vectors (skewed level use, many ties, many
    constant-confidence cases) against brute-force pair counting and
    scipy.stats.mannwhitneyu (U / (n1 * n0));
 3. every S2 reporter's full tagged set (real data) against the same.

Writes revision-2026-11/results/r1_A/indices_check.json.
Usage: python -B r1_A_indices_check.py
"""
import json
import os
import sys
import time

import numpy as np
from scipy.stats import mannwhitneyu

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUT = os.path.join(REV, "results", "r1_A", "indices_check.json")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import indices as IX  # noqa: E402
import reliability  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")


def code_indices(codes, x):
    cnt = IX.count_tables(codes[None, :], x[None, :].astype(np.int8),
                          np.ones((1, len(codes)), bool))
    v = IX.indices_from_counts(cnt)
    return float(v["auroc"][0]), float(v["gamma"][0])


def brute(codes, x):
    c1 = codes[x == 1].astype(float)
    c0 = codes[x == 0].astype(float)
    n1, n0 = len(c1), len(c0)
    if n1 == 0 or n0 == 0:
        return np.nan, np.nan, np.nan
    gt = (c1[:, None] > c0[None, :]).sum()
    lt = (c1[:, None] < c0[None, :]).sum()
    eq = (c1[:, None] == c0[None, :]).sum()
    auroc = (gt + 0.5 * eq) / (n1 * n0)
    gamma = (gt - lt) / (gt + lt) if gt + lt > 0 else np.nan
    u = mannwhitneyu(c1, c0, alternative="two-sided").statistic / (n1 * n0)
    return float(auroc), float(gamma), float(u)


def same(a, b, tol=1e-12):
    if np.isnan(a) or np.isnan(b):
        return bool(np.isnan(a) and np.isnan(b))
    return abs(a - b) <= tol


def main():
    t0 = time.time()
    out = {"source": "scripts/indices.py lines 37-46 (unchanged)"}

    # 1. hand-built edge cases
    cases = {
        "constant_conf_mixed_correct": (np.array([3] * 6), np.array([1, 1, 1, 0, 0, 1])),
        "constant_conf_all_correct": (np.array([4] * 5), np.array([1] * 5)),
        "varied_conf_all_correct": (np.array([0, 2, 4, 3]), np.array([1, 1, 1, 1])),
        "varied_conf_all_incorrect": (np.array([0, 2, 4, 3]), np.array([0, 0, 0, 0])),
        "perfect_resolution": (np.array([4, 4, 1, 0]), np.array([1, 1, 0, 0])),
        "inverse_resolution": (np.array([0, 1, 4, 4]), np.array([1, 1, 0, 0])),
        "partial_ties": (np.array([4, 3, 3, 2, 3]), np.array([1, 1, 0, 0, 1])),
    }
    hand = {}
    for name, (c, xx) in cases.items():
        a, g = code_indices(c, xx)
        ba, bg, bu = brute(c, xx)
        hand[name] = {"code_auroc": None if np.isnan(a) else a,
                      "code_gamma": None if np.isnan(g) else g,
                      "brute_auroc": None if np.isnan(ba) else ba,
                      "brute_gamma": None if np.isnan(bg) else bg,
                      "scipy_U_over_n1n0": None if np.isnan(bu) else bu,
                      "auroc_match": same(a, ba), "gamma_match": same(g, bg)}
    out["hand_cases"] = hand

    # 2. random vectors
    rng = np.random.default_rng([42, 300])
    n_rand, n_bad_a, n_bad_g, n_bad_u, n_const = 5000, 0, 0, 0, 0
    for _ in range(n_rand):
        n = int(rng.integers(2, 41))
        p = rng.dirichlet(np.full(5, rng.choice([0.1, 0.5, 2.0])))
        c = rng.choice(5, size=n, p=p)
        xx = (rng.random(n) < rng.uniform(0.05, 0.95)).astype(int)
        a, g = code_indices(c, xx)
        ba, bg, bu = brute(c, xx)
        n_const += int(len(np.unique(c)) == 1 and 0 < xx.sum() < n)
        n_bad_a += int(not same(a, ba))
        n_bad_g += int(not same(g, bg))
        n_bad_u += int(not same(a, bu))
    out["random_vectors"] = {"n": n_rand, "n_constant_conf_mixed_correct": n_const,
                             "auroc_mismatch_vs_bruteforce": n_bad_a,
                             "gamma_mismatch_vs_bruteforce": n_bad_g,
                             "auroc_mismatch_vs_scipy_mannwhitneyu": n_bad_u}

    # 3. real data: every S2 reporter's full tagged set
    df = common.build_responses("t34")
    users, offsets, codes, x = reliability.student_arrays(df)
    bad_a = bad_g = bad_u = 0
    n_cc = n_cx = 0
    cc_values = set()
    for i in range(len(users)):
        s, e = offsets[i], offsets[i + 1]
        c, xx = codes[s:e], x[s:e]
        a, g = code_indices(c, xx)
        ba, bg, bu = brute(c, xx)
        bad_a += int(not same(a, ba))
        bad_g += int(not same(g, bg))
        bad_u += int(not same(a, bu))
        if len(np.unique(c)) == 1 and 0 < xx.sum() < len(xx):
            n_cc += 1
            cc_values.add(a)
            n_cx += int(np.isnan(g))
    out["real_S2_full_sets"] = {
        "n_students": int(len(users)),
        "auroc_mismatch_vs_bruteforce": bad_a,
        "gamma_mismatch_vs_bruteforce": bad_g,
        "auroc_mismatch_vs_scipy_mannwhitneyu": bad_u,
        "n_constant_conf_mixed_correct": n_cc,
        "auroc_values_in_those_cases": sorted(cc_values),
        "gamma_nan_in_those_cases": n_cx}
    out["conclusion"] = (
        "Submitted code: AUROC = Mann-Whitney U/(n1*n0) with 0.5 tie credit; "
        "defined whenever both correct and incorrect responses exist; equals "
        "exactly .5 when confidence is constant and correctness varies. Gamma "
        "is NaN iff C + D = 0, i.e. constant confidence or constant correctness.")
    out["wall_seconds"] = time.time() - t0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
