"""Pre-specified-grid crossings, alpha-type checks without Monte Carlo error,
and AUROC/gamma diagnostics (revision round 2).

Part G (grid). The analysis plan (section 3) fixed N in {5, 10, 20, 50, 100,
200}; addendum 1 (2026-07-09, before submission) added 12, 15, 25, 30 and 40
to Design A as an interpolation-sensitivity check. The printed Design A
crossings use the full 11-point grid. This part recomputes every Design A
crossing (random and temporal halves, every index, coding and AUROC rule,
both samples) on the 6 pre-specified points only, with the same crossing
rule. CIs come from the submitted joint student bootstrap
(default_rng([42, 77, sid]), 1,000 draws), regenerated here for random and
temporal halves and checked against the cached random-half draws
(results/baseline and results/r1_A/prereg caches, exact) and against the
temporal intervals in results/r2_A.json (<= 1e-12). Design B already uses
the pre-specified grid.

Part H (alpha). On the crossed S2 subset (14 quizzes with >= 300 students
who rated every item), alpha minus ICC(1,k) and Phi minus ICC(1,k) carry no
Monte Carlo error. ICC(1,k), alpha and Phi are recomputed from the data with
the r1_B_alpha.py routines and checked against results/r1_B_parts (1e-7).
The printed differences against the per-student split-half estimator are
medians over quizzes of 50-split means; their Monte Carlo error is
estimated by rerunning the 50 splits per quiz with K = 200 fresh seeds
default_rng([42, 505, 2, QuizId, k]); draw 0 ([42, 103, 2, QuizId]) is
checked against the stored split-half values.

Part I (AUROC). AUROC minus gamma reliability at N <= 50 under both AUROC
rules (results/r2_A.json); the share of Design A students that enter the
gamma and AUROC curves (both halves defined), and a diagnostic of the
non-monotone AUROC curve between N = 15 and 20, from
results/r2_E_parts/designA_<sample>.json (run r2_E_mc_designA.py first).

Usage: python -B r2_F_grid.py
Writes: revision-2026-11/results/r2_F_grid.json
"""
import json
import os
import sys
import time

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_E_common as EC  # noqa: E402
import r2_A_rerun as RA  # noqa: E402
import r1_B_alpha as BA  # noqa: E402  (helpers only)

REL = RA.REL
GRID = EC.FULL_GRID
BOOTK = list(RA.BOOTK)  # bias brier gamma auroc auroc_pr bias_g brier_g
B = REL.B_BOOT
K_ALPHA = 200


def eq(a, b):
    if isinstance(a, float) and isinstance(b, float):
        return abs(a - b) <= 1e-9
    return a == b


def ci_all(d):
    dd = np.where(np.isfinite(d), d, EC.BIG)
    lo, hi = np.percentile(dd, [2.5, 97.5])
    return [EC.fnum(lo), EC.fnum(hi)]


def part_grid(sample, r2a):
    t0 = time.time()
    sid = EC.SIDS[sample]
    df = RA.common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    poolA = counts >= max(GRID)
    nA = int(poolA.sum())
    H = {}
    for N in GRID:
        rng = np.random.default_rng([EC.SEED, 5, N, sid])
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolA, N, rng)
        A, Bh, TA, TB, _ = RA.halves(C, X, rng)
        H[N] = (A, Bh, TA, TB)
    G = len(GRID)
    pt = {k: {"random": np.array([REL.sb_mean(H[N][0][k], H[N][1][k]) for N in GRID]),
              "temporal": np.array([REL.sb_mean(H[N][2][k], H[N][3][k]) for N in GRID])} for k in BOOTK}
    bs = {k: {"random": np.full((B, G), np.nan), "temporal": np.full((B, G), np.nan)} for k in BOOTK}
    brng = np.random.default_rng([EC.SEED, 77, sid])
    for b in range(B):
        idx = brng.integers(0, nA, nA)
        for gi, N in enumerate(GRID):
            A, Bh, TA, TB = H[N]
            for k in BOOTK:
                bs[k]["random"][b, gi] = REL.sb_mean(A[k][idx], Bh[k][idx])
                bs[k]["temporal"][b, gi] = REL.sb_mean(TA[k][idx], TB[k][idx])
    del H
    bnpz = np.load(os.path.join(EC.RESULTS, "baseline", "data", "cache", f"designA_{sample}.npz"))
    pnpz = np.load(os.path.join(EC.RESULTS, "r1_A", "prereg", "data", "cache", f"designA_{sample}.npz"))
    reg = {}
    for k in ("bias", "brier", "gamma", "auroc"):
        reg[f"boot_random_{k}_vs_cache"] = float(np.nanmax(np.abs(bs[k]["random"] - bnpz[f"boot_{k}"])))
    reg["boot_random_auroc_pr_vs_cache"] = float(np.nanmax(np.abs(bs["auroc_pr"]["random"] - pnpz["boot_auroc"])))
    for k in BOOTK:
        for sp in ("random", "temporal"):
            ref_pt = np.array([r2a["indices"][k][f"curve_{sp}_point"][str(N)] for N in GRID])
            reg[f"curve_{sp}_{k}"] = float(np.nanmax(np.abs(pt[k][sp] - ref_pt)))
            D = bs[k][sp]
            lo = np.nanpercentile(D, 2.5, axis=0)
            hi = np.nanpercentile(D, 97.5, axis=0)
            ref = np.array([r2a["indices"][k][f"curve_{sp}_ci"][str(N)] for N in GRID])
            reg[f"ci_{sp}_{k}"] = float(np.nanmax(np.abs(np.c_[lo, hi] - ref)))
    reg_pass = bool(max(reg.values()) <= 1e-12)
    print(f"[{sample}] grid bootstrap regenerated {time.time() - t0:.0f}s regression pass {reg_pass}", flush=True)
    if not reg_pass:
        print(json.dumps(reg, indent=1))
        raise SystemExit("regenerated draws do not reproduce the cached draws / r2_A.json")

    out = {}
    match_full = True
    for k in BOOTK:
        for sp in ("random", "temporal"):
            for thr in EC.THRS:
                t = f"{int(round(thr * 100))}"
                D = bs[k][sp]
                d_full = np.array([EC.cross(D[b], thr) for b in range(B)])
                d_fz = np.array([EC.cross(EC.on_frozen(D[b]), thr, EC.FROZEN_GRID) for b in range(B)])
                full_pt = EC.cross(pt[k][sp], thr)
                fz_pt = EC.cross(EC.on_frozen(pt[k][sp]), thr, EC.FROZEN_GRID)
                ref = r2a["indices"][k][f"crossings_{sp}"][t]
                mine_full = {"point": EC.fnum(full_pt), "ci_all_draws": ci_all(d_full)}
                same = eq(mine_full["point"], ref["point"]) and \
                    all(eq(a, b) for a, b in zip(mine_full["ci_all_draws"], ref["ci_all_draws"]))
                match_full &= bool(same)
                out[f"{sp}_{k}_{t}"] = {
                    "full_grid": {"point": mine_full["point"], "ci_all_draws": mine_full["ci_all_draws"],
                                  "share_draws_reaching": float(np.mean(np.isfinite(d_full)))},
                    "frozen_grid": {"point": EC.fnum(fz_pt), "ci_all_draws": ci_all(d_fz),
                                    "share_draws_reaching": float(np.mean(np.isfinite(d_fz)))},
                    "frozen_minus_full": (float(fz_pt - full_pt) if np.isfinite(fz_pt) and np.isfinite(full_pt) else None)}
    # S2 replication criterion (b) and the S1 median comparison on both grids
    extra = {}
    if sample == "t34":
        for gname, sel, grid in (("full_grid", lambda v: v, GRID),
                                 ("frozen_grid", EC.on_frozen, EC.FROZEN_GRID)):
            for rule, akey in (("plan_rule", "auroc_pr"), ("tie_convention", "auroc")):
                vals = {k: EC.cross(sel(pt[k]["random"]), 0.70, grid) for k in ("bias", "brier", "gamma", akey)}
                hits = {k: bool((not np.isfinite(v)) or v > EC.S2_MEDIAN_VOLUME) for k, v in vals.items()}
                extra[f"criterion_b_{gname}_{rule}"] = {"minN70": {k: EC.fnum(v) for k, v in vals.items()},
                                                        "above_31": hits, "met": bool(sum(hits.values()) >= 3)}
    else:
        for gname, sel, grid in (("full_grid", lambda v: v, GRID),
                                 ("frozen_grid", EC.on_frozen, EC.FROZEN_GRID)):
            v = EC.cross(sel(pt["bias"]["random"]), 0.70, grid)
            extra[f"bias70_vs_S1_median_15_{gname}"] = {"minN70": EC.fnum(v), "above_15": bool(v > 15)}
    return {"n_poolA": nA, "regression_max_abs_diff": reg, "regression_pass": reg_pass,
            "full_grid_reproduces_r2_A": bool(match_full), "crossings": out, "comparisons": extra,
            "runtime_seconds": round(time.time() - t0, 1)}


def part_alpha():
    t0 = time.time()
    sample, sid = "t34", 2
    with open(os.path.join(EC.RESULTS, "r1_B_parts", "alpha_t34.json"), encoding="utf-8") as fh:
        ref = json.load(fh)
    refq = {r["QuizId"]: r for r in ref["quizzes"] if r["n_students"] >= BA.PRIMARY_STUD}
    df = BA.common.build_responses(sample)
    df = df[df.QuizId != -1]
    qpq = df.groupby("QuizId").QuestionId.nunique().rename("k")
    tg = BA.common.tagged(df)
    g = tg.groupby(["QuizId", "UserId"]).agg(nq=("QuestionId", "nunique"), n=("QuestionId", "size")).reset_index()
    g = g.merge(qpq.reset_index(), on="QuizId")
    comp = g[(g.nq == g.k) & (g.n == g.nq)]
    per_quiz = comp.groupby("QuizId").agg(nstud=("UserId", "size"), k=("k", "first"))
    prim = per_quiz[(per_quiz.nstud >= BA.PRIMARY_STUD) & (per_quiz.k >= BA.MIN_K)]
    assert sorted(prim.index.tolist()) == sorted(refq.keys())
    tgc = tg.merge(comp[["QuizId", "UserId"]], on=["QuizId", "UserId"])
    rows = []
    maxdiff = {"icc1_k": 0.0, "alpha": 0.0, "phi": 0.0, "sb_per_student": 0.0}
    for qid in sorted(prim.index.tolist()):
        sub = tgc[tgc.QuizId == qid]
        items = np.sort(sub.QuestionId.unique())
        uids = np.sort(sub.UserId.unique())
        C = (sub.pivot(index="UserId", columns="QuestionId", values="Confidence").loc[uids, items]
             .to_numpy() / 25.0).round().astype(np.int8)
        X = sub.pivot(index="UserId", columns="QuestionId", values="IsCorrect").loc[uids, items].to_numpy().astype(np.int8)
        c = REL.IX.LEVELS[C.astype(np.int64)]
        Ys = {"bias": c - X, "brier": (c - X) ** 2}
        rng0 = np.random.default_rng([EC.SEED, 103, sid, int(qid)])
        A0, B0, _, _, _ = REL.half_indices(C, X, rng0)
        reruns = []
        for kk in range(1, K_ALPHA + 1):
            rr = np.random.default_rng([EC.SEED, 505, sid, int(qid), kk])
            Ak, Bk, _, _, _ = REL.half_indices(C, X, rr)
            reruns.append({key: REL.sb_mean(Ak[key], Bk[key]) for key in ("bias", "brier")})
        row = {"QuizId": int(qid), "n_students": int(C.shape[0]), "k_items": int(C.shape[1])}
        for key in ("bias", "brier"):
            Y = Ys[key]
            av = BA.anova2(Y)
            est = {"alpha": float(BA.alpha_of(Y)), "phi": float(av["phi"]), "icc1_k": float(av["icc1_k"]),
                   "sb_per_student_draw0": REL.sb_mean(A0[key], B0[key]),
                   "sb_per_student_split_se_draw0": EC.mc_se(REL.sb_corr_cols(A0[key], B0[key])),
                   "sb_per_student_reruns": [r[key] for r in reruns]}
            for m, mm in (("icc1_k", "icc1_k"), ("alpha", "alpha"), ("phi", "phi"), ("sb_per_student_draw0", "sb_per_student")):
                maxdiff[mm] = max(maxdiff[mm], abs(est[m] - refq[qid][key][mm]))
            row[key] = est
        rows.append(row)
    ok = maxdiff["icc1_k"] <= 1e-7 and maxdiff["alpha"] <= 1e-7 and maxdiff["phi"] <= 1e-7 \
        and maxdiff["sb_per_student"] <= 1e-12
    print(f"[alpha] 14 quizzes recomputed {time.time() - t0:.0f}s checks pass {ok} {maxdiff}", flush=True)
    if not ok:
        raise SystemExit("alpha-subset recomputation does not match results/r1_B_parts")

    def med_rng(v):
        v = np.asarray(v, float)
        return {"median": float(np.median(v)), "min": float(v.min()), "max": float(v.max())}

    out = {"n_quizzes": len(rows), "n_quiz_student_rows": int(sum(r["n_students"] for r in rows)),
           "k_range": [int(min(r["k_items"] for r in rows)), int(max(r["k_items"] for r in rows))],
           "checks_max_abs_diff_vs_r1_B_parts": maxdiff, "checks_pass": bool(ok), "K_split_reruns": K_ALPHA}
    printed = {"bias": {"alpha_minus_sb": 0.022, "phi_minus_sb": 0.007},
               "brier": {"alpha_minus_sb": 0.039, "phi_minus_sb": 0.018}}
    for key in ("bias", "brier"):
        a = np.array([r[key]["alpha"] for r in rows])
        p = np.array([r[key]["phi"] for r in rows])
        i1 = np.array([r[key]["icc1_k"] for r in rows])
        s0 = np.array([r[key]["sb_per_student_draw0"] for r in rows])
        S = np.array([r[key]["sb_per_student_reruns"] for r in rows])  # quizzes x K
        blk = {"alpha_minus_icc1_k": med_rng(a - i1), "phi_minus_icc1_k": med_rng(p - i1),
               "alpha_median": float(np.median(a)), "phi_median": float(np.median(p)),
               "icc1_k_median": float(np.median(i1)),
               "n_quizzes_alpha_gt_icc1_k": int(np.sum(a > i1)), "n_quizzes_phi_gt_icc1_k": int(np.sum(p > i1)),
               "printed_draw0": {"alpha_minus_sb_median": float(np.median(a - s0)),
                                 "phi_minus_sb_median": float(np.median(p - s0)),
                                 "sb_minus_icc1_k_median": float(np.median(s0 - i1))},
               "max_split_se_per_quiz_draw0": float(max(r[key]["sb_per_student_split_se_draw0"] for r in rows))}
        for lab_, base in (("alpha_minus_sb", a), ("phi_minus_sb", p)):
            meds = np.median(base[:, None] - S, axis=0)
            pv = printed[key][lab_]
            blk[f"{lab_}_median_over_reruns"] = {
                "mean": float(meds.mean()), "sd": float(meds.std(ddof=1)),
                "p2_5": float(np.percentile(meds, 2.5)), "p97_5": float(np.percentile(meds, 97.5)),
                "printed_value": pv, "share_reruns_rounding_to_printed": float(np.mean(np.round(meds, 3) == pv)),
                "share_reruns_at_or_above_draw0": float(np.mean(meds >= blk["printed_draw0"][f"{lab_}_median"]))}
        sbm = np.median(S - i1[:, None], axis=0)
        blk["sb_minus_icc1_k_median_over_reruns"] = {"mean": float(sbm.mean()), "sd": float(sbm.std(ddof=1))}
        out[key] = blk
    out["per_quiz"] = [{"QuizId": r["QuizId"], "n_students": r["n_students"], "k_items": r["k_items"],
                        **{f"{key}_{m}": r[key][m] for key in ("bias", "brier")
                           for m in ("alpha", "phi", "icc1_k", "sb_per_student_draw0")}} for r in rows]
    out["runtime_seconds"] = round(time.time() - t0, 1)
    return out


def part_auroc(r2a_all):
    out = {}
    small = [N for N in GRID if N <= 50]
    for sample in ("t12", "t34"):
        lab = EC.LAB[sample]
        ix = r2a_all[lab]["indices"]
        g = np.array([ix["gamma"]["curve_random_point"][str(N)] for N in small])
        blk = {}
        for rule, k in (("plan_rule", "auroc_pr"), ("tie_convention", "auroc")):
            a = np.array([ix[k]["curve_random_point"][str(N)] for N in small])
            d = a - g
            blk[f"auroc_minus_gamma_{rule}"] = {"min": float(d.min()), "max": float(d.max()),
                                                "N_at_min": int(small[int(d.argmin())]),
                                                "N_at_max": int(small[int(d.argmax())]),
                                                "per_N": {str(N): float(d[i]) for i, N in enumerate(small)}}
        with open(os.path.join(EC.PARTS, f"designA_{sample}.json"), encoding="utf-8") as fh:
            pa = json.load(fh)
        dc0 = pa["draw_curves"]
        Mg = np.array(dc0["random_gamma"]["per_draw"], float)[:, :len(small)]
        for rule, k in (("plan_rule", "auroc_pr"), ("tie_convention", "auroc")):
            Ma = np.array(dc0[f"random_{k}"]["per_draw"], float)[:, :len(small)]
            D = Ma - Mg
            mD = D.mean(axis=0)
            blk[f"auroc_minus_gamma_{rule}"]["over_further_draws"] = {
                "mean_per_N": {str(N): float(mD[i]) for i, N in enumerate(small)},
                "sd_per_N": {str(N): float(D[:, i].std(ddof=1)) for i, N in enumerate(small)},
                "min_of_means": float(mD.min()), "max_of_means": float(mD.max()),
                "K": int(D.shape[0])}
        ks = pa["kept_shares"]
        blk["students_entering_resolution_curves"] = {
            "definition": "share of Design A students with both halves defined, mean over the 50 random splits; "
                          "sb_corr_cols drops a student when either half is undefined",
            "gamma_and_auroc_plan_rule_draw0": ks["gamma"]["random_both_halves_defined"]["draw0"],
            "auroc_plan_rule_equals_gamma": bool(ks["gamma"]["random_both_halves_defined"]["draw0"]
                                                 == ks["auroc_pr"]["random_both_halves_defined"]["draw0"]),
            "auroc_tie_convention_draw0": ks["auroc"]["random_both_halves_defined"]["draw0"],
            "gamma_one_half_undefined_draw0": ks["gamma"]["random_half_undefined"]["draw0"],
            "gamma_temporal_both_halves_defined_draw0": ks["gamma"]["temporal_both_halves_defined"]["draw0"],
            "gamma_and_auroc_plan_rule_mean_over_further_draws":
                ks["gamma"]["random_both_halves_defined"]["mean_over_draws"]}
        dc = pa["draw_curves"]
        diag = {}
        for k in ("auroc_pr", "auroc", "gamma"):
            c = dc[f"random_{k}"]
            diag[k] = {str(N): {"draw0": c["draw0"][str(N)], "mean_over_draws": c["mean_over_draws"][str(N)],
                                "sd_over_draws": c["sd_over_draws"][str(N)]} for N in (10, 12, 15, 20, 25, 30)}
        M = np.array(dc["random_auroc_pr"]["per_draw"], float)
        i15, i20 = GRID.index(15), GRID.index(20)
        diag["auroc_pr_share_draws_r15_gt_r20"] = float(np.mean(M[:, i15] > M[:, i20]))
        diag["auroc_pr_draw0_r15_minus_r20"] = float(dc["random_auroc_pr"]["draw0"]["15"] - dc["random_auroc_pr"]["draw0"]["20"])
        diag["auroc_pr_mean_r15_minus_r20"] = float(M[:, i15].mean() - M[:, i20].mean())
        diag["auroc_pr_sd_of_r15_minus_r20_over_draws"] = float((M[:, i15] - M[:, i20]).std(ddof=1))
        diag["half_sizes"] = {"15": [7, 8], "20": [10, 10]}
        diag["auroc_pr_mean_curve_monotone_N_le_50"] = bool(np.all(np.diff(
            [dc["random_auroc_pr"]["mean_over_draws"][str(N)] for N in small]) > 0))
        blk["auroc_dip_diagnostic"] = diag
        out[lab] = blk
    return out


if __name__ == "__main__":
    with open(os.path.join(EC.RESULTS, "r2_A.json"), encoding="utf-8") as fh:
        R2A = json.load(fh)
    res = {"_definitions": {
        "full_grid": "Design A grid N = 5, 10, 12, 15, 20, 25, 30, 40, 50, 100, 200 (pre-specified points plus the five points added by addendum 1); the printed crossings",
        "frozen_grid": "the pre-specified grid N = 5, 10, 20, 50, 100, 200 (analysis plan section 3); same curves and crossing rule, the five added points ignored",
        "ci_all_draws": "2.5th and 97.5th percentiles over all 1,000 joint student bootstrap draws (seed [42, 77, sid]); a draw that never reaches the level within 200 is ranked above 200; 'inf' = above 200",
        "share_draws_reaching": "share of the 1,000 bootstrap draws whose curve reaches the level within 200",
        "criterion_b": "S2 replication criterion (b): minimum N for .70 above the S2 median of 31 rated answers for at least three of bias, Brier, gamma and AUROC",
        "alpha_minus_icc1_k": "per-quiz coefficient alpha minus the one-way ICC(1,k) (items nested in students), median and range over the 14 quizzes; deterministic, no Monte Carlo error",
        "phi_minus_icc1_k": "per-quiz G-theory absolute coefficient Phi minus ICC(1,k); deterministic",
        "x_minus_sb_median_over_reruns": "median over the 14 quizzes of (alpha or Phi) minus the 50-split per-student split-half mean, recomputed on K = 200 fresh split seeds [42, 505, 2, QuizId, k]; mean, SD and percentiles over the K reruns",
        "students_entering_resolution_curves": "share of Design A students with both halves defined (mean over the 50 random splits), i.e. the students whose half-scores enter the gamma or AUROC correlation at that N",
        "auroc_minus_gamma": "Design A random-half reliability of AUROC minus that of gamma, at N <= 50"}}
    for s in ("t12", "t34"):
        res[f"grid_{EC.LAB[s]}"] = part_grid(s, R2A[EC.LAB[s]])
    res["alpha_subset_S2"] = part_alpha()
    res["auroc_gamma"] = part_auroc(R2A)
    with open(os.path.join(EC.RESULTS, "r2_F_grid.json"), "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print("written results/r2_F_grid.json")
