"""Stream C parts (2) and (3): difficulty reweighting sensitivity of the
signed-bias minimum N, and structured (content / difficulty / term) splits.

Everything runs on the Design A pool and the SAME per-N subsamples as the
submitted analysis (seed [42, 5, N, sid], original subsample routine), and
the random and temporal splits are regenerated from the same generator, so
the "random" and "temporal" curves reproduce the submitted ones exactly
(checked below). The joint student bootstrap reuses the submitted draws
(seed [42, 77, sid]).

Variants (all Spearman-Brown corrected split-half, mean over splits):
  random, temporal             - as submitted (regression check)
  w_diff_random / _temporal    - same splits, each tagged response weighted
                                 so that the student's tagged difficulty-
                                 stratum mix (global quintiles of the
                                 leave-one-out item p-value) matches the mix
                                 of that student's untagged answers
  w_diffstrand_random          - same, strata = difficulty quintile x strand
  balanced_random              - N responses drawn with probability
                                 proportional to the difficulty weight
                                 (difficulty-balanced subset), then 50
                                 random splits, unweighted indices
  strand, topic, term          - halves formed by whole categories: each
                                 student's categories are put in a random
                                 order and the first N//2 responses in that
                                 order form half A (only one category can
                                 straddle the cut); 50 random orders
  difficulty                   - harder half vs easier half by item p-value
                                 (ties broken at random; 50 repetitions)
  w_null_diff_random /         - null: the student's own w_diff (w_diffstrand)
  w_null_diffstrand_random       weights permuted across that student's rated
                                 answers (same Kish ESS, no difficulty link)

Usage: python r1_C_splits.py t34   (or t12)
Writes: revision-2026-11/results/C/r1_C_splits_{sample}.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_C_common as RC  # noqa: E402
import indices as IX  # noqa: E402
import reliability as RL  # noqa: E402

GRID = RL.GRID
R = RL.R_SPLITS
BB = RL.B_BOOT
IDXK = ["bias", "brier", "conf", "acc"]
BOOTK = ["bias", "brier"]
MIN_UNTAG_TARGET = 20
W_FLOOR = 1e-6
BASE = os.path.join(RC.REV, "results", "baseline")


def subsample_idx(offsets, pool, N, rng):
    S = np.empty((len(pool), N), dtype=np.int64)
    for r, i in enumerate(pool):
        s, e = offsets[i], offsets[i + 1]
        S[r] = np.sort(rng.choice(e - s, size=N, replace=False)) + s
    return S


def halves(C, X, masks, W=None):
    n = C.shape[0]
    A = {k: np.empty((n, len(masks))) for k in IDXK}
    B = {k: np.empty((n, len(masks))) for k in IDXK}
    mn = 1e-9 if W is not None else 1
    for j, m in enumerate(masks):
        wa = m if W is None else m * W
        wb = ~m if W is None else (~m) * W
        iA = IX.indices_from_counts(IX.count_tables(C, X, wa), min_n=mn)
        iB = IX.indices_from_counts(IX.count_tables(C, X, wb), min_n=mn)
        for k in IDXK:
            A[k][:, j] = iA[k]
            B[k][:, j] = iB[k]
    return A, B


def random_masks(n, N, rng):
    hA = N // 2
    out = []
    for _ in range(R):
        keys = rng.random((n, N))
        rank = np.argsort(np.argsort(keys, axis=1), axis=1)
        out.append(rank < hA)
    return out


def rank_from_order(order):
    n, N = order.shape
    rank = np.empty_like(order)
    np.put_along_axis(rank, order, np.broadcast_to(np.arange(N), (n, N)).copy(), axis=1)
    return rank


def keyed_masks(K, rng):
    n, N = K.shape
    hA = N // 2
    ncat = int(K.max()) + 1
    rows = np.arange(n)[:, None]
    out = []
    for _ in range(R):
        u = rng.random((n, ncat))
        val = u[rows, K]
        tie = rng.random((n, N))
        order = np.lexsort((tie, val), axis=1)
        out.append(rank_from_order(order) < hA)
    return out


def ordered_masks(P, rng):
    n, N = P.shape
    hA = N // 2
    out = []
    for _ in range(R):
        tie = rng.random((n, N))
        order = np.lexsort((tie, P), axis=1)
        out.append(rank_from_order(order) < hA)  # half A = harder (lower p)
    return out


def cat_counts(K, w, ncat):
    n, N = K.shape
    flat = (np.arange(n)[:, None] * ncat + K).ravel()
    return np.bincount(flat, weights=np.asarray(w, float).ravel(),
                       minlength=n * ncat).reshape(n, ncat)


def shared_share(K, m, ncat):
    """Per student: share of the N responses whose category occurs in both halves."""
    ca = cat_counts(K, m, ncat)
    cb = cat_counts(K, ~m, ncat)
    sh = (ca > 0) & (cb > 0)
    return ((ca + cb) * sh).sum(axis=1) / K.shape[1]


def sb_cols(A, B):
    return RL.sb_corr_cols(A, B)


def sbm(A, B):
    v = sb_cols(A, B)
    return float(np.nanmean(v)) if np.any(~np.isnan(v)) else float("nan")


def min_n_cross(rvec, thr):
    """Copy of the nested function in reliability.run (lines 208-221)."""
    logN = np.log(GRID)
    rv = np.asarray(rvec, float)
    if np.all(np.isnan(rv)):
        return np.nan
    for i, N in enumerate(GRID):
        if not np.isnan(rv[i]) and rv[i] >= thr:
            if i == 0:
                return float(GRID[0])
            r0, r1 = rv[i - 1], rv[i]
            if np.isnan(r0) or r1 == r0:
                return float(N)
            f = (thr - r0) / (r1 - r0)
            return float(np.exp(logN[i - 1] + f * (logN[i] - logN[i - 1])))
    return np.inf


def ci_or_partial(draws):
    fin = draws[np.isfinite(draws)]
    if len(fin) == len(draws):
        return [float(v) for v in np.percentile(fin, [2.5, 97.5])]
    if len(fin):
        return [float(np.percentile(fin, 2.5)), None]
    return [None, None]


def fnum(v):
    if v is None:
        return None
    v = float(v)
    if np.isinf(v):
        return "inf"
    if np.isnan(v):
        return None
    return v


def strata_weights(df, tg_index, uid_t, strat_all, strat_t):
    """Per tagged response weight P_u(s) / P_t(s) within student; students
    with < MIN_UNTAG_TARGET untagged answers keep weight 1. Weights are
    renormalised to mean 1 over the student's tagged responses."""
    d = pd.DataFrame({"u": df.UserId.to_numpy(), "t": df.tag.to_numpy(), "s": strat_all})
    ct = d.groupby(["u", "t", "s"]).size().unstack("s", fill_value=0)
    tt = ct.xs(True, level="t")
    uu = ct.xs(False, level="t").reindex(tt.index).fillna(0)
    nu = uu.sum(axis=1)
    Pt = tt.div(tt.sum(axis=1), axis=0)
    Pu = uu.div(nu.where(nu > 0, 1), axis=0)
    ratio = (Pu / Pt.where(Pt > 0)).fillna(0.0)
    unmatched = (Pu * (Pt == 0)).sum(axis=1)
    has_target = (nu >= MIN_UNTAG_TARGET) & (unmatched < 1 - 1e-9)
    norm = (1 - unmatched).where(has_target, 1)
    ratio = ratio.div(norm, axis=0)
    cols = {c: j for j, c in enumerate(ratio.columns)}
    rmat = ratio.to_numpy()
    rowpos = pd.Index(ratio.index).get_indexer(uid_t)
    colpos = np.array([cols[s] for s in strat_t]) if not np.issubdtype(np.asarray(strat_t).dtype, np.integer) \
        else pd.Index(ratio.columns).get_indexer(strat_t)
    w = rmat[rowpos, colpos]
    ht = has_target.reindex(ratio.index).to_numpy()[rowpos]
    w = np.where(ht, w, 1.0)
    info = {"n_reporters": int(len(tt)), "n_with_target": int(has_target.sum()),
            "unmatched_mass_mean_with_target": float(unmatched[has_target].mean()),
            "unmatched_mass_q90_with_target": float(unmatched[has_target].quantile(0.9))}
    Pu_mat = Pu.reindex(columns=ratio.columns).to_numpy()
    return w, info, has_target, Pu_mat, ratio.index, list(ratio.columns)


def run(sample):
    t0 = time.time()
    sid = RC.SID[sample]
    df = RC.load(sample)
    # tagged arrays exactly as reliability.student_arrays
    tg = df[df.tag].sort_values(["UserId", "DateAnswered"], kind="stable")
    codes = (tg["Confidence"].to_numpy() / 25.0).round().astype(np.int8)
    x = tg["IsCorrect"].to_numpy().astype(np.int8)
    uid = tg["UserId"].to_numpy()
    users, starts = np.unique(uid, return_index=True)
    offsets = np.append(starts, len(uid))
    u0, o0, c0, x0 = RL.student_arrays(df)
    assert np.array_equal(u0, users) and np.array_equal(o0, offsets)
    assert np.array_equal(c0, codes) and np.array_equal(x0, x)
    counts = np.diff(offsets)
    pool = np.where(counts >= max(GRID))[0]
    nA = len(pool)

    # keys per tagged response
    strand_c, strand_lab = pd.factorize(tg["strand"].astype(str))
    topic_c, topic_lab = pd.factorize(tg["topic"].astype(str))
    term_c, term_lab = pd.factorize(tg["term"].astype(str))
    p_t = tg["p_loo"].to_numpy().astype(float)
    p_t = np.where(np.isnan(p_t), np.nanmean(df["p_loo"].to_numpy()), p_t)
    ds_all = df["dstrat"].to_numpy().astype(np.int64)
    ds_t = tg["dstrat"].to_numpy().astype(np.int64)
    sc_all = pd.Categorical(df["strand"].astype(str), categories=list(strand_lab)).codes.astype(np.int64)
    # strand codes of untagged-only categories become -1 -> own code
    sc_all = np.where(sc_all < 0, 9999, sc_all)
    dss_all = ds_all * 10000 + sc_all
    dss_t = ds_t * 10000 + strand_c.astype(np.int64)
    Wd, info_d, ht_d, Pu_d, Pu_index, Pu_cols = strata_weights(df, None, uid, ds_all, ds_t)
    Wds, info_ds, ht_ds, Pu_ds, Pu_index2, Pu_cols2 = strata_weights(df, None, uid, dss_all, dss_t)
    tgt2 = Pu_ds[Pu_index2.get_indexer(users[np.where(counts >= max(GRID))[0]])]
    dsj2 = pd.Series(dss_t).map({c: j for j, c in enumerate(Pu_cols2)}).fillna(-1).to_numpy().astype(np.int64)
    nc2 = len(Pu_cols2)
    pool_users = users[pool]
    ht_pool = ht_d.reindex(pool_users).fillna(False).to_numpy().astype(bool)
    ht_pool2 = ht_pool & ht_ds.reindex(pool_users).fillna(False).to_numpy().astype(bool)
    prow = Pu_index.get_indexer(pool_users)
    tgt = Pu_d[prow]
    colpos = {c: j for j, c in enumerate(Pu_cols)}
    dsj = pd.Series(ds_t).map(colpos).fillna(-1).to_numpy().astype(np.int64)
    nc = len(Pu_cols)
    target_subset = []
    # null weights: each student's own weight vector permuted across that student's rated answers
    # (same Kish ESS, link to difficulty / strand broken); uid is sorted, so blocks are contiguous
    pr1 = np.random.default_rng([RC.MASTER_SEED, 350, sid]).random(len(uid))
    Wd_null = Wd[np.lexsort((pr1, uid))]
    pr2 = np.random.default_rng([RC.MASTER_SEED, 351, sid]).random(len(uid))
    Wds_null = Wds[np.lexsort((pr2, uid))]
    assert np.allclose(np.add.reduceat(Wd_null, offsets[:-1]), np.add.reduceat(Wd, offsets[:-1]))
    t_prep = time.time() - t0

    # baseline curves for the regression check
    base = pd.read_csv(os.path.join(BASE, "results", f"reliability_curves_{sample}.csv"))
    bnpz = np.load(os.path.join(BASE, "data", "cache", f"designA_{sample}.npz"))

    brng = np.random.default_rng([RC.MASTER_SEED, 77, sid])
    IDXS = [brng.integers(0, nA, nA) for _ in range(BB)]

    variants = ["random", "temporal", "w_diff_random", "w_diff_temporal", "w_diffstrand_random",
                "w_null_diff_random", "w_null_diffstrand_random",
                "balanced_random", "balanced_ds_random", "strand", "topic", "term", "difficulty"]
    boot_vars = [v for v in variants if v != "w_diff_temporal"]
    curves = {v: {k: [] for k in IDXK} for v in variants}
    undef = {v: [] for v in ["w_diff_random", "w_diffstrand_random"]}
    boot = {v: {k: np.full((BB, len(GRID)), np.nan) for k in BOOTK} for v in boot_vars}
    diag = {"strand": [], "topic": [], "term": [], "difficulty": [], "weights": [], "balanced": []}
    restricted = []
    regress = {}
    rows = np.arange(nA)[:, None]
    for gi, N in enumerate(GRID):
        tN = time.time()
        rngN = np.random.default_rng([RC.MASTER_SEED, 5, N, sid])
        S = subsample_idx(offsets, pool, N, rngN)
        C, X = codes[S], x[S]
        Cchk, Xchk, _ = RL.subsample_matrix(users, offsets, codes, x, counts >= max(GRID), N,
                                            np.random.default_rng([RC.MASTER_SEED, 5, N, sid]))
        assert np.array_equal(Cchk, C) and np.array_equal(Xchk, X)
        hA = N // 2
        mR = random_masks(nA, N, rngN)
        mT = np.zeros((nA, N), bool)
        mT[:, :hA] = True
        HV = {}
        HV["random"] = halves(C, X, mR)
        HV["temporal"] = halves(C, X, [mT])
        W1 = Wd[S]
        W2 = Wds[S]
        HV["w_diff_random"] = halves(C, X, mR, W1)
        HV["w_diff_temporal"] = halves(C, X, [mT], W1)
        HV["w_diffstrand_random"] = halves(C, X, mR, W2)
        HV["w_null_diff_random"] = halves(C, X, mR, Wd_null[S])
        HV["w_null_diffstrand_random"] = halves(C, X, mR, Wds_null[S])
        undef["w_diff_random"].append(float(np.mean(np.isnan(HV["w_diff_random"][0]["bias"]))))
        undef["w_diffstrand_random"].append(float(np.mean(np.isnan(HV["w_diffstrand_random"][0]["bias"]))))
        # difficulty-balanced subset
        rngB = np.random.default_rng([RC.MASTER_SEED, 320, N, sid])
        Sb = np.empty_like(S)
        n_floor = 0
        for r, i in enumerate(pool):
            s, e = offsets[i], offsets[i + 1]
            w = Wd[s:e]
            if (w > 0).sum() < N:
                n_floor += 1
            pr = np.maximum(w, W_FLOOR)
            Sb[r] = np.sort(rngB.choice(e - s, size=N, replace=False, p=pr / pr.sum())) + s
        mB = random_masks(nA, N, rngB)
        HV["balanced_random"] = halves(codes[Sb], x[Sb], mB)
        # difficulty x strand balanced subset
        rngB2 = np.random.default_rng([RC.MASTER_SEED, 321, N, sid])
        Sb2 = np.empty_like(S)
        n_floor2 = 0
        for r, i in enumerate(pool):
            s, e = offsets[i], offsets[i + 1]
            w = Wds[s:e]
            if (w > 0).sum() < N:
                n_floor2 += 1
            pr = np.maximum(w, W_FLOOR)
            Sb2[r] = np.sort(rngB2.choice(e - s, size=N, replace=False, p=pr / pr.sum())) + s
        mB2 = random_masks(nA, N, rngB2)
        HV["balanced_ds_random"] = halves(codes[Sb2], x[Sb2], mB2)

        # composition check vs target (students with a target only)
        def tvd_to_target(SS, dj, tg_, ncc):
            K = dj[SS]
            cc = cat_counts(np.where(K < 0, 0, K), np.ones_like(K, float), ncc) / N
            return 0.5 * np.abs(cc - tg_).sum(axis=1)
        ht_p2 = ht_ds.reindex(pool_users).fillna(False).to_numpy().astype(bool)
        diag["balanced"].append(dict(
            N=N, n_students_floor_needed=int(n_floor), n_students_floor_needed_ds=int(n_floor2),
            tvd_to_untagged_mix_plain_mean=float(tvd_to_target(S, dsj, tgt, nc)[ht_pool].mean()),
            tvd_to_untagged_mix_balanced_mean=float(tvd_to_target(Sb, dsj, tgt, nc)[ht_pool].mean()),
            ds_tvd_to_untagged_mix_plain_mean=float(tvd_to_target(S, dsj2, tgt2, nc2)[ht_p2].mean()),
            ds_tvd_to_untagged_mix_balanced_ds_mean=float(tvd_to_target(Sb2, dsj2, tgt2, nc2)[ht_p2].mean())))
        ess = (W1.sum(1) ** 2) / (W1 ** 2).sum(1) / N
        with np.errstate(divide="ignore", invalid="ignore"):
            ess2 = (W2.sum(1) ** 2) / (W2 ** 2).sum(1) / N
        diag["weights"].append(dict(N=N, kish_ess_ratio_mean=float(ess.mean()),
                                    kish_ess_ratio_median=float(np.median(ess)),
                                    kish_ess_ratio_mean_target_students=float(ess[ht_pool].mean()),
                                    weight_max=float(W1.max()), weight_q99=float(np.quantile(W1, .99)),
                                    diffstrand_kish_ess_ratio_mean=float(np.nanmean(ess2)),
                                    diffstrand_kish_ess_ratio_mean_target_students=float(np.nanmean(ess2[ht_pool2])),
                                    diffstrand_weight_max=float(W2.max())))
        # keyed splits
        Ks = {"strand": strand_c[S], "topic": topic_c[S], "term": term_c[S]}
        for j, key in enumerate(["strand", "topic", "term"]):
            rngK = np.random.default_rng([RC.MASTER_SEED, 330 + j, N, sid])
            K = Ks[key]
            mk = keyed_masks(K, rngK)
            HV[key] = halves(C, X, mk)
            ncat = int(K.max()) + 1
            cnt = cat_counts(K, np.ones_like(K, float), ncat)
            largest = cnt.max(axis=1) / N
            ncats = (cnt > 0).sum(axis=1)
            sh_k = np.mean([shared_share(K, m, ncat) for m in mk], axis=0)
            sh_r = np.mean([shared_share(K, m, ncat) for m in mR], axis=0)
            diag[key].append(dict(N=N, n_categories_mean=float(ncats.mean()),
                                  n_categories_median=float(np.median(ncats)),
                                  share_students_single_category=float((ncats == 1).mean()),
                                  largest_category_share_mean=float(largest.mean()),
                                  shared_category_response_share_keyed=float(sh_k.mean()),
                                  shared_category_response_share_random=float(sh_r.mean())))
            sub = largest <= 0.5
            rr = {"N": N, "key": key, "n_students": int(sub.sum())}
            if sub.sum() >= 10:
                for k in BOOTK:
                    rr[f"{k}_random"] = sbm(HV["random"][0][k][sub], HV["random"][1][k][sub])
                    rr[f"{k}_keyed"] = sbm(HV[key][0][k][sub], HV[key][1][k][sub])
                    rr[f"{k}_temporal"] = sbm(HV["temporal"][0][k][sub], HV["temporal"][1][k][sub])
            restricted.append(rr)
        rngD = np.random.default_rng([RC.MASTER_SEED, 333, N, sid])
        P = p_t[S]
        md = ordered_masks(P, rngD)
        HV["difficulty"] = halves(C, X, md)
        pa = np.mean([np.nanmean(np.where(m, P, np.nan), axis=1) for m in md], axis=0)
        pb = np.mean([np.nanmean(np.where(~m, P, np.nan), axis=1) for m in md], axis=0)
        pra = np.mean([np.abs(np.nanmean(np.where(m, P, np.nan), axis=1)
                              - np.nanmean(np.where(~m, P, np.nan), axis=1)) for m in mR], axis=0)
        diag["difficulty"].append(dict(N=N, mean_p_harder_half=float(pa.mean()),
                                       mean_p_easier_half=float(pb.mean()),
                                       mean_abs_p_gap_random_split=float(pra.mean())))
        # point curves
        for v in variants:
            A, B = HV[v]
            for k in IDXK:
                curves[v][k].append(sbm(A[k], B[k]))
        # same comparison restricted to pool students who have a weighting target
        ts = {"N": N, "n_students": int(ht_pool2.sum())}
        for v in ["random", "w_diff_random", "w_diffstrand_random", "w_null_diff_random", "w_null_diffstrand_random", "balanced_random", "balanced_ds_random"]:
            A, B = HV[v]
            for k in BOOTK:
                ts[f"{v}_{k}"] = sbm(A[k][ht_pool2], B[k][ht_pool2])
        target_subset.append(ts)
        # regression check vs baseline curves
        for v, split in [("random", "random"), ("temporal", "temporal")]:
            for k in IDXK:
                row = base[(base.design == "A") & (base.split == split) & (base["index"] == k) & (base.N == N)]
                d = abs(float(row.r.iloc[0]) - curves[v][k][-1])
                key = f"{split}_{k}"
                regress[key] = max(regress.get(key, 0.0), d)
        # bootstrap
        for v in boot_vars:
            A, B = HV[v]
            for k in BOOTK:
                a, b = A[k], B[k]
                bv = boot[v][k]
                for bi in range(BB):
                    idx = IDXS[bi]
                    bv[bi, gi] = sbm(a[idx], b[idx])
        print(f"[{sample}] N={N} done in {time.time() - tN:.1f}s", flush=True)
    regress["boot_random_bias_vs_baseline_npz"] = float(np.nanmax(np.abs(boot["random"]["bias"] - bnpz["boot_bias"])))
    regress["boot_random_brier_vs_baseline_npz"] = float(np.nanmax(np.abs(boot["random"]["brier"] - bnpz["boot_brier"])))

    # min-N and r(20) summaries
    summ = {}
    i20 = GRID.index(20)
    for v in variants:
        summ[v] = {}
        for k in IDXK:
            cv = np.array(curves[v][k])
            blk = {"minN70_point": fnum(min_n_cross(cv, 0.70)), "minN80_point": fnum(min_n_cross(cv, 0.80)),
                   "r20_point": fnum(cv[i20])}
            if v in boot_vars and k in BOOTK:
                bv = boot[v][k]
                d70 = np.array([min_n_cross(bv[b], 0.70) for b in range(BB)])
                d80 = np.array([min_n_cross(bv[b], 0.80) for b in range(BB)])
                blk["minN70_ci"] = ci_or_partial(d70)
                blk["minN80_ci"] = ci_or_partial(d80)
                blk["share_minN70_gt200"] = float(np.mean(~np.isfinite(d70)))
                blk["share_minN80_gt200"] = float(np.mean(~np.isfinite(d80)))
                blk["r20_ci"] = [float(q) for q in np.nanpercentile(bv[:, i20], [2.5, 97.5])]
                if v != "random":
                    rb = boot["random"][k]
                    r70 = np.array([min_n_cross(rb[b], 0.70) for b in range(BB)])
                    r80 = np.array([min_n_cross(rb[b], 0.80) for b in range(BB)])
                    for thr, dv, dr in [("70", d70, r70), ("80", d80, r80)]:
                        both = np.isfinite(dv) & np.isfinite(dr)
                        dd = dv[both] - dr[both]
                        pt_v = min_n_cross(cv, int(thr) / 100)
                        pt_r = min_n_cross(np.array(curves["random"][k]), int(thr) / 100)
                        blk[f"minN{thr}_minus_random"] = {
                            "point": fnum(pt_v - pt_r) if np.isfinite(pt_v) and np.isfinite(pt_r) else None,
                            "ci_finite_draws": [float(q) for q in np.percentile(dd, [2.5, 97.5])] if len(dd) else None,
                            "share_draws_both_finite": float(both.mean())}
                    dr20 = bv[:, i20] - rb[:, i20]
                    blk["r20_minus_random"] = {"point": fnum(cv[i20] - curves["random"][k][i20]),
                                               "ci": [float(q) for q in np.nanpercentile(dr20, [2.5, 97.5])]}
            summ[v][k] = blk

    # full-set bias level: raw vs difficulty-reweighted (all tagged responses)
    dres = codes.astype(float) * 0.25 - x
    st = offsets[:-1]
    sw = np.add.reduceat(Wd, st)
    raw = np.add.reduceat(dres, st) / counts
    wb_d = np.add.reduceat(dres * Wd, st) / sw
    wb_ds = np.add.reduceat(dres * Wds, st) / np.add.reduceat(Wds, st)
    wb_dn = np.add.reduceat(dres * Wd_null, st) / np.add.reduceat(Wd_null, st)
    wb_dsn = np.add.reduceat(dres * Wds_null, st) / np.add.reduceat(Wds_null, st)
    htu = ht_d.reindex(users).fillna(False).to_numpy().astype(bool)
    htu2 = ht_ds.reindex(users).fillna(False).to_numpy().astype(bool)
    full = {}
    for lab, m in [("designA_pool", counts >= 200), ("ge20_tagged", counts >= 20), ("ge5_tagged", counts >= 5)]:
        blk = {}
        for wl, wv, hm in [("diff", wb_d, htu), ("diffstrand", wb_ds, htu2),
                           ("null_diff", wb_dn, htu), ("null_diffstrand", wb_dsn, htu2)]:
            mm = m & hm
            a, b = raw[mm], wv[mm]
            k10 = max(int(round(0.10 * mm.sum())), 1)
            ta, tb = set(np.argsort(-a)[:k10]), set(np.argsort(-b)[:k10])
            blk[wl] = {"n_students": int(mm.sum()), "mean_raw_bias": float(a.mean()),
                       "mean_weighted_bias": float(b.mean()),
                       "weighted_minus_raw": RC.summ(b - a, seed=[RC.MASTER_SEED, 340, sid]),
                       "pearson_raw_weighted": float(np.corrcoef(a, b)[0, 1]),
                       "top_decile_jaccard_raw_vs_weighted": len(ta & tb) / len(ta | tb)}
        full[lab] = blk

    out = {"sample": sample, "n_poolA": int(nA), "grid": GRID, "R_splits": R, "B_boot": BB,
           "seeds": {"subsample_and_random_splits": "[42, 5, N, sid] (as submitted)",
                     "joint_bootstrap": "[42, 77, sid] (as submitted)",
                     "balanced_subset": "[42, 320, N, sid]", "balanced_ds_subset": "[42, 321, N, sid]", "strand": "[42, 330, N, sid]",
                     "topic": "[42, 331, N, sid]", "term": "[42, 332, N, sid]",
                     "difficulty": "[42, 333, N, sid]", "bias_shift_bootstrap": "[42, 340, sid]",
                     "null_weight_permutation": "[42, 350, sid] (difficulty), [42, 351, sid] (difficulty x strand)"},
           "regression_max_abs_diff_vs_baseline": regress,
           "curves": {v: {k: [fnum(r) for r in curves[v][k]] for k in IDXK} for v in variants},
           "summary": summ,
           "weighting": {"difficulty": info_d, "difficulty_x_strand": info_ds,
                         "pool_students_with_target": int(ht_pool.sum()),
                         "min_untagged_for_target": MIN_UNTAG_TARGET,
                         "share_undefined_half_bias_w_diff_random_by_N": undef["w_diff_random"],
                         "share_undefined_half_bias_w_diffstrand_random_by_N": undef["w_diffstrand_random"],
                         "difficulty_strata_edges_p_loo": df.attrs["dstrat_edges"]},
           "diagnostics": diag, "restricted_largest_category_le_half": restricted,
           "target_subset_curves": target_subset,
           "target_subset_minN": {f"{v}_{k}": {"minN70_point": fnum(min_n_cross([t[f"{v}_{k}"] for t in target_subset], 0.70)),
                                               "minN80_point": fnum(min_n_cross([t[f"{v}_{k}"] for t in target_subset], 0.80))}
                                  for v in ["random", "w_diff_random", "w_diffstrand_random", "w_null_diff_random", "w_null_diffstrand_random", "balanced_random", "balanced_ds_random"]
                                  for k in BOOTK},
           "fullset_bias_shift": full,
           "category_labels": {"strand": list(map(str, strand_lab)), "term": list(map(str, term_lab)),
                               "n_topic": int(len(topic_lab))},
           "runtime_seconds": {"prep": t_prep, "total": time.time() - t0}}
    os.makedirs(RC.OUT_DIR, exist_ok=True)
    with open(os.path.join(RC.OUT_DIR, f"r1_C_splits_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(sample, "splits done", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t34")
