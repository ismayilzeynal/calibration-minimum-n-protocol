"""Monte Carlo precision of the Design A estimates (revision round 2).

Two Monte Carlo components sit in every printed Design A value:

  (1) split error: each reliability is the mean of 50 Spearman-Brown values,
      one per random split of ONE fixed N-answer subsample per student;
  (2) subsample-draw error: the printed curves rest on one seeded subsample
      of N rated answers per student at each N.

Part 1 (split error) reruns draw 0 with the submitted seeds
default_rng([42, 5, N, sid]) through the same routines as r2_A_rerun.py,
keeps the 50 per-split values at every N, and reports their Monte Carlo SE
for every index, coding and rule, plus the split SE of every crossing
(delta method, checked by resampling the 50 per-split values within each N,
B = 2,000, default_rng([42, 503, sid])). Temporal halves have no split
component (one deterministic split per subsample); neither does the
one-way variance-components (ICC) curve.

Part 2 (subsample draws) reruns the frozen Design A estimator (same pool,
same N grid, 50 random splits plus the temporal split) on K = 100 further
subsample draws, draw k using default_rng([42, 501, N, sid, k]) for both the
subsample and its splits, exactly as draw 0 uses [42, 5, N, sid]. The SD over
the K further draws is the total Monte Carlo SD of a single-draw value
(subsample plus split component); draw 0 is not among the K draws.
No bootstrap is run per draw.

Regression gates: draw 0 curves of every key against results/r2_A.json
(<= 1e-12); per-N split SE against r2_A.json mc_se_over_splits (<= 1e-12);
pool-sensitivity curves against r2_A.json; ICC(1,N) against
results/r1_B_parts/vc_<sample>.json (<= 1e-12); full-set undefined shares
against r2_A.json.

Usage: python -B r2_E_mc_designA.py t12   (then t34)
Writes: revision-2026-11/results/r2_E_parts/designA_<sample>.json
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
import r2_A_rerun as RA  # noqa: E402  (imports the paper's common / indices / reliability)

REL = RA.REL
IX = RA.IX
GRID = EC.FULL_GRID
assert list(REL.GRID) == GRID
K_DRAWS = 100
B_RS = 2000
KEYS = list(RA.KEYS)          # bias brier gamma auroc auroc_pr conf acc bias_g brier_g
CROSS_KEYS = ["bias", "brier", "gamma", "auroc", "auroc_pr", "bias_g", "brier_g"]
RES_KEYS = ["gamma", "auroc_pr", "auroc"]
EXCK = ["bias", "brier"]


def oneway_icc(C, X):
    c = IX.LEVELS[C.astype(np.int64)]
    d = c - X
    out = {}
    for k, Y in (("bias", d), ("brier", d ** 2)):
        n, N = Y.shape
        m = Y.mean(axis=1)
        ssw = ((Y - m[:, None]) ** 2).sum(axis=1)
        msb = N * np.var(m, ddof=1)
        msw = ssw.sum() / (n * (N - 1))
        out[k] = float((msb - msw) / msb)
    return out


def one_draw(arrs, poolA, keep, seedfn, keep_splits):
    users, offsets, codes, x = arrs
    G = len(GRID)
    rnd = {k: np.full(G, np.nan) for k in KEYS}
    tmp = {k: np.full(G, np.nan) for k in KEYS}
    icc = {k: np.full(G, np.nan) for k in EXCK}
    kept = {k: {"random_both_halves_defined": np.full(G, np.nan),
                "random_half_undefined": np.full(G, np.nan),
                "temporal_both_halves_defined": np.full(G, np.nan)} for k in RES_KEYS}
    und_full = {k: np.full(G, np.nan) for k in RES_KEYS}
    splits = {k: [] for k in KEYS} if keep_splits else None
    splits_keep = {k: [] for k in EXCK} if keep_splits else None
    keep_curves = {k: {"random": np.full(G, np.nan), "temporal": np.full(G, np.nan)} for k in EXCK}
    for gi, N in enumerate(GRID):
        rng = np.random.default_rng(seedfn(N))
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolA, N, rng)
        A, Bh, TA, TB, full = RA.halves(C, X, rng)
        for k in KEYS:
            cols = REL.sb_corr_cols(A[k], Bh[k])
            rnd[k][gi] = float(np.nanmean(cols))
            tmp[k][gi] = REL.sb_mean(TA[k], TB[k])
            if keep_splits:
                splits[k].append(cols)
        for k in EXCK:
            if keep_splits:
                splits_keep[k].append(REL.sb_corr_cols(A[k][keep], Bh[k][keep]))
            keep_curves[k]["random"][gi] = REL.sb_mean(A[k][keep], Bh[k][keep])
            keep_curves[k]["temporal"][gi] = REL.sb_mean(TA[k][keep], TB[k][keep])
        ic = oneway_icc(C, X)
        for k in EXCK:
            icc[k][gi] = ic[k]
        for k in RES_KEYS:
            both = ~(np.isnan(A[k]) | np.isnan(Bh[k]))
            kept[k]["random_both_halves_defined"][gi] = float(both.mean())
            kept[k]["random_half_undefined"][gi] = float(np.mean(np.isnan(A[k])))
            kept[k]["temporal_both_halves_defined"][gi] = float(np.mean(~(np.isnan(TA[k]) | np.isnan(TB[k]))))
            und_full[k][gi] = float(np.mean(np.isnan(full[k])))
    return dict(random=rnd, temporal=tmp, icc=icc, kept=kept, und_full=und_full,
                splits=splits, splits_keep=splits_keep, keep_curves=keep_curves)


def curve_dict(v):
    return {str(N): (None if np.isnan(v[i]) else float(v[i])) for i, N in enumerate(GRID)}


def run(sample):
    t0 = time.time()
    sid = EC.SIDS[sample]
    lab = EC.LAB[sample]
    df = RA.common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    poolA = counts >= max(GRID)
    pool = np.where(poolA)[0]
    acc_full = np.array([x[offsets[i]:offsets[i + 1]].mean() for i in pool])
    const_full = np.array([codes[offsets[i]:offsets[i + 1]].min() == codes[offsets[i]:offsets[i + 1]].max()
                           for i in pool])
    keep = ~((acc_full >= 0.99) | const_full)
    arrs = (users, offsets, codes, x)
    with open(os.path.join(EC.RESULTS, "r2_A.json"), encoding="utf-8") as fh:
        r2a = json.load(fh)[lab]
    with open(os.path.join(EC.RESULTS, "r1_B_parts", f"vc_{sample}.json"), encoding="utf-8") as fh:
        vc = json.load(fh)

    # ---------- draw 0: submitted seeds, per-split values kept ----------
    d0 = one_draw(arrs, poolA, keep, lambda N: [EC.SEED, 5, N, sid], True)
    reg = {}
    for k in KEYS:
        ref_r = np.array([r2a["indices"][k]["curve_random_point"][str(N)] for N in GRID])
        ref_t = np.array([r2a["indices"][k]["curve_temporal_point"][str(N)] for N in GRID])
        ref_se = np.array([r2a["indices"][k]["mc_se_over_splits"][str(N)] for N in GRID])
        se0 = np.array([EC.mc_se(v) for v in d0["splits"][k]])
        reg[f"curve_random_{k}"] = float(np.nanmax(np.abs(d0["random"][k] - ref_r)))
        reg[f"curve_temporal_{k}"] = float(np.nanmax(np.abs(d0["temporal"][k] - ref_t)))
        reg[f"split_se_{k}"] = float(np.nanmax(np.abs(se0 - ref_se)))
    for k in EXCK:
        for sp in ("random", "temporal"):
            ref = np.array([r2a["pool_sensitivity"][k][f"curve_{sp}_point"][str(N)] for N in GRID])
            reg[f"pool_sensitivity_{sp}_{k}"] = float(np.nanmax(np.abs(d0["keep_curves"][k][sp] - ref)))
        ref_icc = np.array([r["%s" % k]["icc1_N"] for r in vc["by_N"]])
        reg[f"icc_{k}"] = float(np.nanmax(np.abs(d0["icc"][k] - ref_icc)))
    for k in RES_KEYS:
        ref = np.array([r2a["indices"][k]["undefined_full_set"][str(N)] for N in GRID])
        reg[f"undefined_full_{k}"] = float(np.nanmax(np.abs(d0["und_full"][k] - ref)))
    reg_pass = bool(max(reg.values()) <= 1e-12)
    print(f"[{sample}] draw 0 done {time.time() - t0:.0f}s regression pass {reg_pass}", flush=True)
    if not reg_pass:
        print(json.dumps(reg, indent=1))
        raise SystemExit("draw 0 does not reproduce r2_A.json")

    # ---------- part 1: split error ----------
    rs_rng = np.random.default_rng([EC.SEED, 503, sid])
    split_part = {}
    fams = [(k, d0["splits"][k], d0["random"][k]) for k in KEYS] + \
           [(f"pool_sensitivity_{k}", d0["splits_keep"][k], d0["keep_curves"][k]["random"]) for k in EXCK]
    for name, sp, curve in fams:
        se = np.array([EC.mc_se(v) for v in sp])
        blk = {"per_N_split_se": curve_dict(se),
               "max_split_se_all_N": float(np.nanmax(se)),
               "argmax_N": int(GRID[int(np.nanargmax(se))]),
               "max_split_se_N_ge_20": float(np.nanmax(se[GRID.index(20):])),
               "max_split_se_N_ge_50": float(np.nanmax(se[GRID.index(50):])),
               "per_N_split_values_n_defined": {str(N): int(np.sum(~np.isnan(sp[i]))) for i, N in enumerate(GRID)},
               "per_split_values": {str(N): [None if np.isnan(v) else float(v) for v in sp[i]]
                                    for i, N in enumerate(GRID)}}
        base_key = name.replace("pool_sensitivity_", "")
        if base_key in CROSS_KEYS:
            rs, _ = EC.resample_split_crossings(sp, EC.THRS, rs_rng, B_RS)
            cr = {}
            for thr in EC.THRS:
                t = f"{int(round(thr * 100))}"
                pt = EC.cross(curve, thr)
                cr[t] = {"point": EC.fnum(pt), "split_se_delta": EC.delta_se_cross(curve, se, thr)}
                cr[t].update(rs[t])
                fz = EC.on_frozen(curve)
                cr[t]["frozen_grid_point"] = EC.fnum(EC.cross(fz, thr, EC.FROZEN_GRID))
                cr[t]["frozen_grid_split_se_delta"] = EC.delta_se_cross(fz, EC.on_frozen(se), thr, EC.FROZEN_GRID)
            blk["crossings"] = cr
            w = (np.log(15) - np.log(10)) / (np.log(20) - np.log(10))
            i10, i20 = GRID.index(10), GRID.index(20)
            blk["r15_interpolated"] = {"point": float(EC.interp15(curve)),
                                       "split_se": float(np.sqrt((1 - w) ** 2 * se[i10] ** 2 + w ** 2 * se[i20] ** 2))}
            r20 = curve[i20]
            blk["spearman_brown_projection_070_from_N20"] = {
                "point": float(20 * (0.7 / 0.3) * (1 - r20) / r20) if r20 > 0 else None,
                "split_se_delta": float(20 * (0.7 / 0.3) * se[i20] / r20 ** 2) if r20 > 0 else None}
        split_part[name] = blk
    print(f"[{sample}] part 1 done {time.time() - t0:.0f}s", flush=True)

    # ---------- part 2: K further subsample draws ----------
    draws = []
    for k in range(1, K_DRAWS + 1):
        dk = one_draw(arrs, poolA, keep, lambda N, kk=k: [EC.SEED, 501, N, sid, kk], False)
        draws.append(dk)
        if k % 10 == 0:
            print(f"[{sample}] draw {k}/{K_DRAWS} {time.time() - t0:.0f}s", flush=True)

    def stack(field, k, sub=None):
        if sub is None:
            return np.array([d[field][k] for d in draws])
        return np.array([d[field][k][sub] for d in draws])

    curves_part = {}
    for field in ("random", "temporal"):
        for k in KEYS:
            M = stack(field, k)
            curves_part[f"{field}_{k}"] = {
                "draw0": curve_dict(d0[field][k]),
                "mean_over_draws": curve_dict(M.mean(axis=0)),
                "sd_over_draws": curve_dict(M.std(axis=0, ddof=1)),
                "per_draw": [[None if np.isnan(v) else float(v) for v in row] for row in M]}
    for k in EXCK:
        M = stack("icc", k)
        curves_part[f"icc_{k}"] = {"draw0": curve_dict(d0["icc"][k]),
                                   "mean_over_draws": curve_dict(M.mean(axis=0)),
                                   "sd_over_draws": curve_dict(M.std(axis=0, ddof=1)),
                                   "per_draw": [[float(v) for v in row] for row in M]}
    kept_part = {}
    for k in RES_KEYS:
        kept_part[k] = {}
        for sub in ("random_both_halves_defined", "random_half_undefined", "temporal_both_halves_defined"):
            M = stack("kept", k, sub)
            kept_part[k][sub] = {"draw0": curve_dict(d0["kept"][k][sub]),
                                 "mean_over_draws": curve_dict(M.mean(axis=0)),
                                 "sd_over_draws": curve_dict(M.std(axis=0, ddof=1))}
        M = stack("und_full", k)
        kept_part[k]["full_set_undefined"] = {"draw0": curve_dict(d0["und_full"][k]),
                                              "mean_over_draws": curve_dict(M.mean(axis=0)),
                                              "sd_over_draws": curve_dict(M.std(axis=0, ddof=1))}

    cross_part = {}
    for field in ("random", "temporal", "icc"):
        keys = CROSS_KEYS if field != "icc" else EXCK
        for k in keys:
            M = stack(field, k)
            c0 = d0[field][k]
            for gname, grid in (("full_grid", GRID), ("frozen_grid", EC.FROZEN_GRID)):
                sel = (lambda v: v) if gname == "full_grid" else EC.on_frozen
                for thr in EC.THRS:
                    t = f"{int(round(thr * 100))}"
                    dv = np.array([EC.cross(sel(row), thr, grid) for row in M])
                    blk = EC.draws_summary(dv, extra_gt=EC.S2_MEDIAN_VOLUME if (sample == "t34" and thr == 0.70) else None)
                    blk["draw0"] = EC.fnum(EC.cross(sel(c0), thr, grid))
                    blk["crossing_of_mean_curve"] = EC.fnum(EC.cross(sel(M.mean(axis=0)), thr, grid))
                    v0 = EC.cross(sel(c0), thr, grid)
                    blk["share_draws_below_draw0"] = (float(np.mean(np.where(np.isfinite(dv), dv, EC.BIG) < v0))
                                                      if np.isfinite(v0) else None)
                    if gname == "full_grid":
                        blk["per_draw"] = [EC.fnum(v) for v in dv]
                    cross_part[f"{field}_{k}_{gname}_{t}"] = blk
    # pool sensitivity crossings per draw are not run (draw 0 only)

    # confirmatory and replication quantities per draw (informational)
    bias_r = stack("random", "bias")
    i20 = GRID.index(20)
    info = {}
    r15 = np.array([EC.interp15(row) for row in bias_r])
    info["bias_r20"] = {"draw0": float(d0["random"]["bias"][i20]), "mean": float(bias_r[:, i20].mean()),
                        "sd": float(bias_r[:, i20].std(ddof=1)),
                        "share_below_070": float(np.mean(bias_r[:, i20] < 0.70))}
    info["bias_r15_interpolated"] = {"draw0": float(EC.interp15(d0["random"]["bias"])), "mean": float(r15.mean()),
                                     "sd": float(r15.std(ddof=1)), "min": float(r15.min()), "max": float(r15.max()),
                                     "share_ge_070": float(np.mean(r15 >= 0.70))}
    if sample == "t34":
        crit = {}
        for gname, grid in (("full_grid", GRID), ("frozen_grid", EC.FROZEN_GRID)):
            sel = (lambda v: v) if gname == "full_grid" else EC.on_frozen
            for rule, akey in (("plan_rule", "auroc_pr"), ("tie_convention", "auroc")):
                hits_all = []
                for d in [d0] + draws:
                    hits = 0
                    for k in ("bias", "brier", "gamma", akey):
                        v = EC.cross(sel(d["random"][k]), 0.70, grid)
                        hits += int((not np.isfinite(v)) or v > EC.S2_MEDIAN_VOLUME)
                    hits_all.append(hits)
                hits_all = np.array(hits_all)
                crit[f"{gname}_{rule}"] = {"draw0_indices_above_31": int(hits_all[0]),
                                           "draw0_met": bool(hits_all[0] >= 3),
                                           "share_further_draws_met": float(np.mean(hits_all[1:] >= 3))}
            br = np.array([EC.cross(sel(d["random"]["brier"]), 0.70, grid) for d in draws])
            crit[f"{gname}_brier70"] = {"draw0": EC.fnum(EC.cross(sel(d0["random"]["brier"]), 0.70, grid)),
                                        "share_further_draws_above_31": float(np.mean(np.where(np.isfinite(br), br, EC.BIG) > 31)),
                                        "mean": float(br.mean()) if np.isfinite(br).all() else None,
                                        "median": float(np.median(br))}
        crit["criterion_a_bias_r20_below_070"] = {"draw0": float(d0["random"]["bias"][i20]),
                                                  "share_further_draws_below_070": float(np.mean(bias_r[:, i20] < 0.70))}
        info["replication_criteria_S2"] = crit

    out = {"_definitions": {
               "split_error": "per index or pool-sensitivity variant, printed draw (draw 0): per-N Monte Carlo SE of the mean over the 50 random splits (SD of the per-split Spearman-Brown values / sqrt(50)), its maximum, the per-split values, and per crossing the point value, the delta-method split SE (split_se_delta) and the SD, share reaching and 2.5th / 97.5th percentiles over 2,000 resamples of the per-split values within each N; frozen_grid_* = the same on the pre-specified grid N = 5, 10, 20, 50, 100, 200",
               "draw_curves": "per curve (random_<index>, temporal_<index>, icc_<index>): draw 0 value, mean and SD over the K further subsample draws, and the value in each further draw",
               "draw_crossings": "<halves>_<index>_<grid>_<level>: over the K further draws, median and 2.5th / 97.5th percentiles (a draw that never reaches the level within 200 ranked above 200), mean and SD only when every draw reaches the level, share reaching within 200, draw 0 value, crossing of the mean curve, share of draws below the draw 0 value, share above 31 (S2, .70)",
               "kept_shares": "per resolution index: share of students with both halves defined (random halves, mean over the 50 splits; temporal halves), share with half A undefined, share undefined on the full N-answer set; draw 0 and mean / SD over the further draws",
               "draw_quantities": "bias r(20), bias r(15) interpolated between N = 10 and 20, and (S2) the replication criteria, per further draw; informational",
               "regression_max_abs_diff": "draw 0 against results/r2_A.json (curves, split SE, pool-sensitivity curves, full-set undefined shares) and results/r1_B_parts (ICC curve)",
               "index names": "bias, brier, gamma, auroc (tie convention), auroc_pr (plan rule), conf, acc, bias_g and brier_g (guess-corrected coding c' = .25 + .75 c)"},
           "sample": lab, "n_poolA": int(poolA.sum()), "n_pool_sensitivity_kept": int(keep.sum()),
           "grid": GRID, "frozen_grid": EC.FROZEN_GRID, "R_splits": REL.R_SPLITS,
           "K_further_draws": K_DRAWS, "B_split_resamples": B_RS,
           "seeds": {"draw0_subsample_and_splits": "[42, 5, N, sid] (as submitted)",
                     "further_draw_k_subsample_and_splits": "[42, 501, N, sid, k], k = 1..K",
                     "split_resampling": "[42, 503, sid]", "sid": sid},
           "regression_max_abs_diff": reg, "regression_pass": reg_pass,
           "split_error": split_part, "draw_curves": curves_part, "draw_crossings": cross_part,
           "kept_shares": kept_part, "draw_quantities": info,
           "runtime_seconds": round(time.time() - t0, 1)}
    os.makedirs(EC.PARTS, exist_ok=True)
    with open(os.path.join(EC.PARTS, f"designA_{sample}.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"[{sample}] written, {out['runtime_seconds']}s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
