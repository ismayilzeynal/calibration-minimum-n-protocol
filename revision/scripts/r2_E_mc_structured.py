"""Split Monte Carlo error of the structured-halves and reweighting results
(the bias table of content, difficulty, time and selection; revision round 2).

Draw 0 of r1_C_splits.py is regenerated without the bootstrap, with the same
seeds and the same helper functions (imported, not copied):

  subsample + random splits : [42, 5, N, sid]
  difficulty-balanced subset: [42, 320, N, sid]
  strand / topic / term     : [42, 330 / 331 / 332, N, sid]
  harder versus easier      : [42, 333, N, sid]
  shuffled-weight nulls     : [42, 350, sid], [42, 351, sid]

Every variant is a mean over 50 splits (random splits, random category
orders, or random tie-breaks), so each has a split Monte Carlo SE; the
temporal split has none. The 50 per-split values are kept at every N and
the split SE of each crossing is computed by the delta method and checked by
resampling the per-split values (B = 2,000, default_rng([42, 506, sid])).

Regression gate: point curves of every variant and the target-subset curves
against results/C/r1_C_splits_<sample>.json (<= 1e-12).

Usage: python -B r2_E_mc_structured.py t12   (then t34)
Writes: revision-2026-11/results/r2_E_parts/structured_<sample>.json
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
import r1_C_splits as CS  # noqa: E402  (helpers only; run() is not called)

RC = CS.RC
RL = CS.RL
GRID = EC.FULL_GRID
assert list(CS.GRID) == GRID
B_RS = 2000
KEYS = ["bias", "brier"]
VARIANTS = ["random", "strand", "topic", "term", "difficulty", "w_diff_random", "w_null_diff_random",
            "w_diffstrand_random", "w_null_diffstrand_random", "balanced_random"]
TARGET_VARIANTS = ["random", "w_diff_random", "w_diffstrand_random", "w_null_diff_random",
                   "w_null_diffstrand_random", "balanced_random"]


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
    pool = np.where(counts >= max(GRID))[0]
    nA = len(pool)
    strand_c, strand_lab = pd.factorize(tg["strand"].astype(str))
    topic_c, _ = pd.factorize(tg["topic"].astype(str))
    term_c, _ = pd.factorize(tg["term"].astype(str))
    p_t = tg["p_loo"].to_numpy().astype(float)
    p_t = np.where(np.isnan(p_t), np.nanmean(df["p_loo"].to_numpy()), p_t)
    ds_all = df["dstrat"].to_numpy().astype(np.int64)
    ds_t = tg["dstrat"].to_numpy().astype(np.int64)
    sc_all = pd.Categorical(df["strand"].astype(str), categories=list(strand_lab)).codes.astype(np.int64)
    sc_all = np.where(sc_all < 0, 9999, sc_all)
    dss_all = ds_all * 10000 + sc_all
    dss_t = ds_t * 10000 + strand_c.astype(np.int64)
    Wd, _, ht_d, _, _, _ = CS.strata_weights(df, None, uid, ds_all, ds_t)
    Wds, _, ht_ds, _, _, _ = CS.strata_weights(df, None, uid, dss_all, dss_t)
    pool_users = users[pool]
    ht_pool = ht_d.reindex(pool_users).fillna(False).to_numpy().astype(bool)
    ht_pool2 = ht_pool & ht_ds.reindex(pool_users).fillna(False).to_numpy().astype(bool)
    pr1 = np.random.default_rng([RC.MASTER_SEED, 350, sid]).random(len(uid))
    Wd_null = Wd[np.lexsort((pr1, uid))]
    pr2 = np.random.default_rng([RC.MASTER_SEED, 351, sid]).random(len(uid))
    Wds_null = Wds[np.lexsort((pr2, uid))]
    del df
    print(f"[{sample}] prep {time.time() - t0:.0f}s", flush=True)

    with open(os.path.join(EC.RESULTS, "C", f"r1_C_splits_{sample}.json"), encoding="utf-8") as fh:
        ref = json.load(fh)

    per_split = {v: {k: [] for k in KEYS} for v in VARIANTS}
    per_split_t = {v: {k: [] for k in KEYS} for v in TARGET_VARIANTS}
    temporal = {k: [] for k in KEYS}
    for N in GRID:
        rngN = np.random.default_rng([RC.MASTER_SEED, 5, N, sid])
        S = CS.subsample_idx(offsets, pool, N, rngN)
        C, X = codes[S], x[S]
        hA = N // 2
        mR = CS.random_masks(nA, N, rngN)
        mT = np.zeros((nA, N), bool)
        mT[:, :hA] = True
        HV = {"random": CS.halves(C, X, mR)}
        HT = CS.halves(C, X, [mT])
        HV["w_diff_random"] = CS.halves(C, X, mR, Wd[S])
        HV["w_diffstrand_random"] = CS.halves(C, X, mR, Wds[S])
        HV["w_null_diff_random"] = CS.halves(C, X, mR, Wd_null[S])
        HV["w_null_diffstrand_random"] = CS.halves(C, X, mR, Wds_null[S])
        rngB = np.random.default_rng([RC.MASTER_SEED, 320, N, sid])
        Sb = np.empty_like(S)
        for r, i in enumerate(pool):
            s, e = offsets[i], offsets[i + 1]
            pr = np.maximum(Wd[s:e], CS.W_FLOOR)
            Sb[r] = np.sort(rngB.choice(e - s, size=N, replace=False, p=pr / pr.sum())) + s
        mB = CS.random_masks(nA, N, rngB)
        HV["balanced_random"] = CS.halves(codes[Sb], x[Sb], mB)
        Ks = {"strand": strand_c[S], "topic": topic_c[S], "term": term_c[S]}
        for j, key in enumerate(["strand", "topic", "term"]):
            rngK = np.random.default_rng([RC.MASTER_SEED, 330 + j, N, sid])
            HV[key] = CS.halves(C, X, CS.keyed_masks(Ks[key], rngK))
        rngD = np.random.default_rng([RC.MASTER_SEED, 333, N, sid])
        HV["difficulty"] = CS.halves(C, X, CS.ordered_masks(p_t[S], rngD))
        for v in VARIANTS:
            A, B = HV[v]
            for k in KEYS:
                per_split[v][k].append(RL.sb_corr_cols(A[k], B[k]))
        for v in TARGET_VARIANTS:
            A, B = HV[v]
            for k in KEYS:
                per_split_t[v][k].append(RL.sb_corr_cols(A[k][ht_pool2], B[k][ht_pool2]))
        for k in KEYS:
            temporal[k].append(CS.sbm(HT[0][k], HT[1][k]))
        print(f"[{sample}] N={N} {time.time() - t0:.0f}s", flush=True)

    reg = {}
    for v in VARIANTS:
        for k in KEYS:
            mine = np.array([np.nanmean(c) for c in per_split[v][k]])
            r = np.array([np.nan if z is None else z for z in ref["curves"][v][k]], float)
            reg[f"curve_{v}_{k}"] = float(np.nanmax(np.abs(mine - r)))
    for v in TARGET_VARIANTS:
        for k in KEYS:
            mine = np.array([np.nanmean(c) for c in per_split_t[v][k]])
            r = np.array([t[f"{v}_{k}"] for t in ref["target_subset_curves"]], float)
            reg[f"target_{v}_{k}"] = float(np.nanmax(np.abs(mine - r)))
    for k in KEYS:
        r = np.array([np.nan if z is None else z for z in ref["curves"]["temporal"][k]], float)
        reg[f"curve_temporal_{k}"] = float(np.nanmax(np.abs(np.array(temporal[k]) - r)))
    reg_pass = bool(max(reg.values()) <= 1e-12)
    print(f"[{sample}] regression pass {reg_pass}", flush=True)
    if not reg_pass:
        print(json.dumps(reg, indent=1))
        raise SystemExit("draw 0 does not reproduce results/C curves")

    rs_rng = np.random.default_rng([EC.SEED, 506, sid])
    i20 = GRID.index(20)

    def block(sp):
        curve = np.array([np.nanmean(c) for c in sp])
        se = np.array([EC.mc_se(c) for c in sp])
        rs, _ = EC.resample_split_crossings(sp, (0.70, 0.80), rs_rng, B_RS)
        cr = {}
        for thr in (0.70, 0.80):
            t = f"{int(round(thr * 100))}"
            cr[t] = {"point": EC.fnum(EC.cross(curve, thr)), "split_se_delta": EC.delta_se_cross(curve, se, thr)}
            cr[t].update(rs[t])
        return {"curve": {str(N): float(curve[i]) for i, N in enumerate(GRID)},
                "per_N_split_se": {str(N): float(se[i]) for i, N in enumerate(GRID)},
                "max_split_se_all_N": float(np.nanmax(se)), "r20": float(curve[i20]),
                "r20_split_se": float(se[i20]), "crossings": cr}

    out_v = {v: {k: block(per_split[v][k]) for k in KEYS} for v in VARIANTS}
    out_t = {v: {k: block(per_split_t[v][k]) for k in KEYS} for v in TARGET_VARIANTS}
    out = {"_definitions": {
               "variants": "per variant and index: curve (mean over the 50 splits, orders or tie-breaks), per-N split SE, its maximum, r(20) and its split SE, and per crossing (.70, .80) the point value, the delta-method split SE and the SD, share reaching and 2.5th / 97.5th percentiles over 2,000 resamples of the per-split values",
               "target_subset": "the same for the Design A students with a weighting target (at least 20 unrated answers and a matched difficulty-by-strand mix)",
               "temporal_no_split_error": "temporal-half curves (one deterministic split; no split error)",
               "variant names": "random = random halves; strand, topic, term = halves by whole category; difficulty = harder versus easier items; w_diff_random / w_diffstrand_random = weighted to the unrated difficulty (difficulty-by-strand) mix; w_null_* = the same weights shuffled within student; balanced_random = difficulty-balanced subset, equal weights",
               "regression_max_abs_diff": "draw 0 against results/C/r1_C_splits_<sample>.json"},
           "sample": EC.LAB[sample], "n_poolA": int(nA), "n_target_subset": int(ht_pool2.sum()),
           "grid": GRID, "R_splits": RL.R_SPLITS, "B_split_resamples": B_RS,
           "seeds": {"subsample_and_random_splits": "[42, 5, N, sid]", "balanced_subset": "[42, 320, N, sid]",
                     "strand_topic_term": "[42, 330 / 331 / 332, N, sid]", "difficulty": "[42, 333, N, sid]",
                     "null_weights": "[42, 350, sid], [42, 351, sid]", "split_resampling": "[42, 506, sid]",
                     "sid": sid},
           "regression_max_abs_diff": reg, "regression_pass": reg_pass,
           "variants": out_v, "target_subset": out_t,
           "temporal_no_split_error": {k: {str(N): float(temporal[k][i]) for i, N in enumerate(GRID)} for k in KEYS},
           "runtime_seconds": round(time.time() - t0, 1)}
    with open(os.path.join(EC.PARTS, f"structured_{sample}.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)
    print(f"[{sample}] structured written {out['runtime_seconds']}s", flush=True)


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t34")
