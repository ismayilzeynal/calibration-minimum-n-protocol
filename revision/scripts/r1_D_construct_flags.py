"""Stream D, revision round 1 (post-review, exploratory).

(1) What the signed bias score is made of: split-half reliability of mean
    confidence and of accuracy at the same N as bias, and a cross-set
    (true-score) decomposition of the between-student variance of bias into
    a confidence part, an accuracy part and their covariance.
(2) Absolute-criterion flag stability: "overconfident" = bias above a fixed
    threshold (+.15, +.20) or 95% lower confidence bound above 0, compared
    across two disjoint N-response sets per student, next to the existing
    top-decile flag.
(3) Inputs for the decision table (assembled by r1_D_report.py).

Pool: Design A (students with >= 200 tagged responses). Two designs are used:
  - split-half (paper machinery, reliability.subsample_matrix + half_indices,
    seed [42, 5, N, sid]): reliability of an N-response score, SB-corrected;
  - disjoint sets (rev_round2.disjoint_flags design, seed [42, 92, N, sid]):
    two disjoint sets of N responses per student, 50 repeats; flags and the
    cross-set decomposition use these N-response scores directly.
Student bootstrap: the joint draws of reliability.py ([42, 77, sid]),
B = 1000, applied to every statistic below.

Usage: python -B r1_D_construct_flags.py t12   (or t34)
Writes revision-2026-11/tmp/r1_D_<sample>.json (merged by r1_D_report.py).
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy import stats

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import reliability as REL  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")

NGRID = [10, 15, 20, 30, 40, 50, 60, 100]
R = 50
B = 1000
SEED = 42
Z975 = float(stats.norm.ppf(0.975))
THR_PTS = [15, 20]            # fixed thresholds in confidence-scale points
REL_KEYS = ["bias", "conf", "acc", "brier"]
BOOT_REL_KEYS = ["bias", "conf", "acc"]
RULES = ["decile", "thr15", "thr20", "lcb_z", "thr15_ge", "thr20_ge", "lcb_t"]
BOOT_RULES = ["decile", "thr15", "thr20", "lcb_z"]


# ---------------------------------------------------------------- part 1a
def splithalf(users, offsets, codes, x, poolA, sid):
    point, halves = {}, {}
    for N in NGRID:
        rngN = np.random.default_rng([SEED, 5, N, sid])
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolA, N, rngN)
        A, Bh, TA, TB, _ = REL.half_indices(C, X, rngN)
        point[N] = {k: {"random": REL.sb_mean(A[k], Bh[k]),
                        "temporal": REL.sb_mean(TA[k], TB[k])}
                    for k in REL_KEYS}
        halves[N] = {k: (A[k], Bh[k], TA[k], TB[k]) for k in BOOT_REL_KEYS}
    return point, halves


# ---------------------------------------------------------------- part 1b/2
def disjoint_sets(offsets, codes, x, pool, sid):
    """Two disjoint N-sets per student per repeat, RNG call order identical
    to rev_round2.disjoint_flags. Returns per N: float bias (as rev_round2
    computes it) and integer sums of codes, correctness and (4d)^2."""
    d_all = codes.astype(np.float64) * 0.25 - x
    c64 = codes.astype(np.int64)
    x64 = x.astype(np.int64)
    d4 = c64 - 4 * x64
    tri = np.stack([c64, x64, d4 * d4], axis=1)
    n = len(pool)
    out = {}
    for N in NGRID:
        rngN = np.random.default_rng([SEED, 92, N, sid])
        bA = np.empty((R, n)); bB = np.empty((R, n))
        sA = np.empty((R, n, 3), np.int64); sB = np.empty((R, n, 3), np.int64)
        for rep in range(R):
            for r, i in enumerate(pool):
                s, e = offsets[i], offsets[i + 1]
                sel = rngN.choice(e - s, size=2 * N, replace=False) + s
                a_, b_ = sel[:N], sel[N:]
                bA[rep, r] = d_all[a_].mean()
                bB[rep, r] = d_all[b_].mean()
                sA[rep, r] = tri[a_].sum(0)
                sB[rep, r] = tri[b_].sum(0)
        out[N] = dict(bA=bA, bB=bB, sA=sA, sB=sB)
    return out


def derive(N, s):
    """Per-set quantities from integer sums s (R, n, 3)."""
    sc, sx, sq = s[..., 0], s[..., 1], s[..., 2]
    conf = sc * 0.25 / N
    acc = sx / N
    bint = 25 * sc - 100 * sx              # = 100 * N * bias, exact integer
    sd4 = sc - 4 * sx                      # sum of 4d
    q = N * sq - sd4 * sd4                 # = N (N - 1) var(4d), exact integer
    root = np.sqrt(q / (N - 1.0))          # lcb > 0  <=>  sd4 > z * root
    tq = float(stats.t.ppf(0.975, N - 1))
    flags = {
        "thr15": bint > 15 * N, "thr20": bint > 20 * N,
        "thr15_ge": bint >= 15 * N, "thr20_ge": bint >= 20 * N,
        "lcb_z": sd4 > Z975 * root, "lcb_t": sd4 > tq * root,
        "ucb_z": sd4 < -Z975 * root,
    }
    extra = {"tie15": bint == 15 * N, "tie20": bint == 20 * N,
             "sd0": q == 0, "sd0_pos": (q == 0) & (sd4 > 0)}
    return conf, acc, flags, extra


def top_mask(b, k):
    m = np.zeros(b.shape, bool)
    order = np.argsort(-b, axis=1)[:, :k]
    np.put_along_axis(m, order, True, axis=1)
    return m


def flag_metrics(fA, fB):
    n = fA.shape[1]
    nA = fA.sum(1).astype(float); nB = fB.sum(1).astype(float)
    nAB = (fA & fB).sum(1).astype(float); nU = (fA | fB).sum(1).astype(float)
    pA, pB = nA / n, nB / n
    po = (fA == fB).mean(1)
    pe = pA * pB + (1 - pA) * (1 - pB)
    with np.errstate(divide="ignore", invalid="ignore"):
        dice = np.where(nA + nB > 0, 2 * nAB / (nA + nB), np.nan)
        rAB = np.where(nA > 0, nAB / nA, np.nan)
        rBA = np.where(nB > 0, nAB / nB, np.nan)
        jac = np.where(nU > 0, nAB / nU, np.nan)
        kap = np.where(pe < 1, (po - pe) / (1 - pe), np.nan)
        chance = np.where(nA + nB > 0, 2 * nA * nB / (n * (nA + nB)), np.nan)
    return {
        "share_flagged": float(np.mean((pA + pB) / 2)),
        "retention": float(np.nanmean(dice)),
        "retention_directional": float(np.nanmean(np.concatenate([rAB, rBA]))),
        "jaccard": float(np.nanmean(jac)),
        "kappa": float(np.nanmean(kap)),
        "percent_agreement": float(np.mean(po)),
        "chance_retention": float(np.nanmean(chance)),
        "repeats_with_no_flag": int(np.sum(nA + nB == 0)),
    }


def three_way(oA, uA, oB, uB):
    cA = oA.astype(int) - uA.astype(int)
    cB = oB.astype(int) - uB.astype(int)
    po = (cA == cB).mean(1)
    pe = sum((cA == v).mean(1) * (cB == v).mean(1) for v in (-1, 0, 1))
    with np.errstate(divide="ignore", invalid="ignore"):
        kap = np.where(pe < 1, (po - pe) / (1 - pe), np.nan)
    return {
        "percent_agreement": float(np.mean(po)),
        "kappa": float(np.nanmean(kap)),
        "share_over": float(np.mean(((cA == 1).mean(1) + (cB == 1).mean(1)) / 2)),
        "share_under": float(np.mean(((cA == -1).mean(1) + (cB == -1).mean(1)) / 2)),
        "share_unclassified": float(np.mean(((cA == 0).mean(1) + (cB == 0).mean(1)) / 2)),
        "share_opposite_class": float(np.mean((cA * cB == -1).mean(1))),
    }


def cov_rows(X, Y):
    return ((X - X.mean(1, keepdims=True)) * (Y - Y.mean(1, keepdims=True))).mean(1)


def corr_rows(X, Y):
    return cov_rows(X, Y) / np.sqrt(cov_rows(X, X) * cov_rows(Y, Y))


def decomposition(cA, aA, cB, aB):
    """Cross-set covariances estimate true-score (co)variances at volume N
    because the two sets share no response. Covariances are averaged over the
    50 repeats before ratios are taken."""
    bA, bB = cA - aA, cB - aB
    m = lambda v: float(np.mean(v))  # noqa: E731
    Vb = m((cov_rows(bA, bA) + cov_rows(bB, bB)) / 2)
    Vc = m((cov_rows(cA, cA) + cov_rows(cB, cB)) / 2)
    Va = m((cov_rows(aA, aA) + cov_rows(aB, aB)) / 2)
    Vca = m((cov_rows(cA, aA) + cov_rows(cB, aB)) / 2)
    Tb = m(cov_rows(bA, bB))
    Tc = m(cov_rows(cA, cB))
    Ta = m(cov_rows(aA, aB))
    Tca = m((cov_rows(cA, aB) + cov_rows(aA, cB)) / 2)
    Eb, Ec, Ea, Eca = Vb - Tb, Vc - Tc, Va - Ta, Vca - Tca
    return {
        "var_obs_bias": Vb, "var_obs_conf": Vc, "var_obs_acc": Va,
        "cov_obs_conf_acc": Vca,
        "obs_share_conf": Vc / Vb, "obs_share_acc": Va / Vb,
        "obs_share_cross": -2 * Vca / Vb,
        "var_true_bias": Tb, "var_true_conf": Tc, "var_true_acc": Ta,
        "cov_true_conf_acc": Tca,
        "identity_residual_true": Tb - (Tc + Ta - 2 * Tca),
        "true_share_conf": Tc / Tb, "true_share_acc": Ta / Tb,
        "true_share_cross": -2 * Tca / Tb,
        "true_attrib_conf": (Tc - Tca) / Tb, "true_attrib_acc": (Ta - Tca) / Tb,
        "obs_attrib_conf": (Vc - Vca) / Vb, "obs_attrib_acc": (Va - Vca) / Vb,
        "var_err_bias": Eb, "var_err_conf": Ec, "var_err_acc": Ea,
        "cov_err_conf_acc": Eca,
        "err_share_conf": Ec / Eb, "err_share_acc": Ea / Eb,
        "err_share_cross": -2 * Eca / Eb,
        "share_obs_var_bias_true": Tb / Vb,
        "rho_true_conf_acc": Tca / np.sqrt(Tc * Ta),
        "rho_true_bias_conf": (Tc - Tca) / np.sqrt(Tb * Tc),
        "rho_true_bias_acc": (Tca - Ta) / np.sqrt(Tb * Ta),
        "r_obs_bias_conf_sameset": m((corr_rows(bA, cA) + corr_rows(bB, cB)) / 2),
        "r_obs_bias_acc_sameset": m((corr_rows(bA, aA) + corr_rows(bB, aB)) / 2),
        "r_obs_conf_acc_sameset": m((corr_rows(cA, aA) + corr_rows(cB, aB)) / 2),
        "r_cross_bias_conf": m((corr_rows(bA, cB) + corr_rows(bB, cA)) / 2),
        "r_cross_bias_acc": m((corr_rows(bA, aB) + corr_rows(bB, aA)) / 2),
        "r_parallel_bias": m(corr_rows(bA, bB)),
        "r_parallel_conf": m(corr_rows(cA, cB)),
        "r_parallel_acc": m(corr_rows(aA, aB)),
    }


DEC_BOOT = ["true_attrib_conf", "true_attrib_acc",
            "true_share_conf", "true_share_acc", "true_share_cross",
            "rho_true_bias_conf", "rho_true_bias_acc", "rho_true_conf_acc",
            "r_parallel_bias", "r_parallel_conf", "r_parallel_acc",
            "share_obs_var_bias_true"]


def ci(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return [None, None]
    return [float(x) for x in np.percentile(v, [2.5, 97.5])]


def run(sample):
    t0 = time.time()
    sid = {"t12": 1, "t34": 2}[sample]
    df = common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    poolA = counts >= max(REL.GRID)
    pool = np.where(counts >= 200)[0]
    assert np.array_equal(np.where(poolA)[0], pool)
    n = len(pool)
    k10 = max(int(round(n * 0.10)), 1)
    t_load = time.time() - t0

    t1 = time.time()
    sh_point, halves = splithalf(users, offsets, codes, x, poolA, sid)
    t_split = time.time() - t1

    t1 = time.time()
    ds = disjoint_sets(offsets, codes, x, pool, sid)
    t_disj = time.time() - t1

    # ---------------- point estimates on disjoint sets
    t1 = time.time()
    per_n = {}
    cache = {}
    for N in NGRID:
        d = ds[N]
        cA, aA, fA, eA = derive(N, d["sA"])
        cB, aB, fB, eB = derive(N, d["sB"])
        # rev_round2-identical decile flags (float bias, argsort tie order)
        fA["decile"] = top_mask(d["bA"], k10)
        fB["decile"] = top_mask(d["bB"], k10)
        cache[N] = (cA, aA, cB, aB, fA, fB)
        blk = {"flags": {r: flag_metrics(fA[r], fB[r]) for r in RULES}}
        bAi, bBi = cA - aA, cB - aB
        for r in RULES:
            nf = fA[r].sum() + fB[r].sum()
            rev = ((fA[r] & (bBi <= 0)).sum() + (fB[r] & (bAi <= 0)).sum())
            blk["flags"][r]["share_flagged_with_other_set_bias_le_0"] = (
                float(rev / nf) if nf else None)
        blk["ci_three_way_z"] = three_way(fA["lcb_z"], fA["ucb_z"],
                                          fB["lcb_z"], fB["ucb_z"])
        blk["set_properties"] = {
            "share_bias_exactly_015": float((eA["tie15"].mean() + eB["tie15"].mean()) / 2),
            "share_bias_exactly_020": float((eA["tie20"].mean() + eB["tie20"].mean()) / 2),
            "share_sd_zero": float((eA["sd0"].mean() + eB["sd0"].mean()) / 2),
            "share_sd_zero_and_mean_pos": float((eA["sd0_pos"].mean() + eB["sd0_pos"].mean()) / 2),
            "mean_bias": float(((cA - aA).mean() + (cB - aB).mean()) / 2),
            "sd_bias_between_students": float(np.mean(
                (np.std(cA - aA, axis=1) + np.std(cB - aB, axis=1)) / 2)),
            "t_crit_975": float(stats.t.ppf(0.975, N - 1)),
            "max_abs_float_vs_int_bias": float(max(
                np.abs(d["bA"] - (cA - aA)).max(), np.abs(d["bB"] - (cB - aB)).max())),
        }
        blk["decomposition"] = decomposition(cA, aA, cB, aB)
        per_n[N] = blk
    t_point = time.time() - t1

    # ---------------- joint student bootstrap (reliability.py draws)
    t1 = time.time()
    brng = np.random.default_rng([SEED, 77, sid])
    G = len(NGRID)
    bs_rel = {k: {sp: np.full((B, G), np.nan) for sp in ("random", "temporal")}
              for k in BOOT_REL_KEYS}
    bs_flag = {r: {m: np.full((B, G), np.nan) for m in ("retention", "kappa", "share_flagged")}
               for r in BOOT_RULES}
    bs_dec = {k: np.full((B, G), np.nan) for k in DEC_BOOT}
    for b in range(B):
        idx = brng.integers(0, n, n)
        for gi, N in enumerate(NGRID):
            for k in BOOT_REL_KEYS:
                A_, B_, TA_, TB_ = halves[N][k]
                bs_rel[k]["random"][b, gi] = REL.sb_mean(A_[idx], B_[idx])
                bs_rel[k]["temporal"][b, gi] = REL.sb_mean(TA_[idx], TB_[idx])
            cA, aA, cB, aB, fA, fB = cache[N]
            for r in BOOT_RULES:
                if r == "decile":
                    d = ds[N]
                    mA = top_mask(d["bA"][:, idx], k10)
                    mB = top_mask(d["bB"][:, idx], k10)
                else:
                    mA, mB = fA[r][:, idx], fB[r][:, idx]
                fm = flag_metrics(mA, mB)
                for m_ in ("retention", "kappa", "share_flagged"):
                    bs_flag[r][m_][b, gi] = fm[m_]
            dec = decomposition(cA[:, idx], aA[:, idx], cB[:, idx], aB[:, idx])
            for k in DEC_BOOT:
                bs_dec[k][b, gi] = dec[k]
    t_boot = time.time() - t1

    # ---------------- assemble
    for gi, N in enumerate(NGRID):
        blk = per_n[N]
        blk["splithalf_SB"] = {
            k: {sp: {"r": sh_point[N][k][sp],
                     "ci": ci(bs_rel[k][sp][:, gi]) if k in BOOT_REL_KEYS else None}
                for sp in ("random", "temporal")}
            for k in REL_KEYS}
        for r in BOOT_RULES:
            for m_ in ("retention", "kappa", "share_flagged"):
                blk["flags"][r][m_ + "_ci"] = ci(bs_flag[r][m_][:, gi])
        for k in DEC_BOOT:
            blk["decomposition"][k + "_ci"] = ci(bs_dec[k][:, gi])

    # ---------------- regression gates against canonical files
    base = pd.read_csv(os.path.join(REV, "results", "baseline", "results",
                                    f"reliability_curves_{sample}.csv"))
    base = base[base.design == "A"]
    gate_r, gate_ci = 0.0, 0.0
    checked = []
    for gi, N in enumerate(NGRID):
        if N not in REL.GRID:
            continue
        for k in REL_KEYS:
            for sp in ("random", "temporal"):
                row = base[(base["index"] == k) & (base.N == N) & (base.split == sp)].iloc[0]
                gate_r = max(gate_r, abs(row.r - sh_point[N][k][sp]))
                if k == "bias":
                    lo, hi = ci(bs_rel[k][sp][:, gi])
                    gate_ci = max(gate_ci, abs(row.ci_lo - lo), abs(row.ci_hi - hi))
        checked.append(N)
    with open(os.path.join(PAPER, "results", f"rev_round2_{sample}.json")) as fh:
        r2 = json.load(fh)
    gate_flag = 0.0
    flag_checked = []
    for f in r2["disjoint_flag_consistency"]:
        N = f["N"]
        if N in per_n:
            mine = per_n[N]["flags"]["decile"]
            gate_flag = max(gate_flag, abs(mine["jaccard"] - f["jaccard"]),
                            abs(mine["retention"] - f["pct_flags_retained"]))
            flag_checked.append(N)

    out = {
        "sample": sample,
        "n_pool": n,
        "k_top_decile": k10,
        "grid_N": NGRID,
        "per_N": {str(N): per_n[N] for N in NGRID},
        "regression_gates": {
            "splithalf_r_vs_baseline_curves_max_abs_diff": gate_r,
            "splithalf_N_checked": checked,
            "bias_bootstrap_ci_vs_baseline_curves_max_abs_diff": gate_ci,
            "decile_flags_vs_rev_round2_max_abs_diff": gate_flag,
            "decile_N_checked": flag_checked,
        },
        "seeds": {
            "splithalf": "[42, 5, N, sid] (reliability.py Design A; N = 60 is a new grid point under the same rule)",
            "disjoint_sets": "[42, 92, N, sid] (rev_round2.disjoint_flags; N = 40, 60 are new grid points under the same rule)",
            "bootstrap": "[42, 77, sid] (reliability.py joint draws, integers(0, n, n) per draw)",
            "sid": sid,
        },
        "runtime_seconds": {
            "load": round(t_load, 1), "splithalf": round(t_split, 1),
            "disjoint_draws": round(t_disj, 1), "point_estimates": round(t_point, 1),
            "bootstrap": round(t_boot, 1), "total": round(time.time() - t0, 1),
        },
    }
    os.makedirs(os.path.join(REV, "tmp"), exist_ok=True)
    with open(os.path.join(REV, "tmp", f"r1_D_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(sample, "done", out["runtime_seconds"], out["regression_gates"])


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
