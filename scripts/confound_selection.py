"""Pre-specified exploratory analyses: accuracy confound (plan Sec 5),
reporting selection (Sec 6), temporal test-retest (Sec 3).

Run after reliability.py:  python scripts/confound_selection.py t12
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
import indices as IX
from reliability import student_arrays, subsample_matrix, half_indices, \
    sb_corr_cols, sb_mean, GRID, MASTER_SEED

B_BOOT = 1000


def nan_pearson(a, b):
    m = ~(np.isnan(a) | np.isnan(b))
    if m.sum() < 10:
        return np.nan
    return float(np.corrcoef(a[m], b[m])[0, 1])


def nan_spearman(a, b):
    m = ~(np.isnan(a) | np.isnan(b))
    ra = pd.Series(a[m]).rank().to_numpy()
    rb = pd.Series(b[m]).rank().to_numpy()
    return float(np.corrcoef(ra, rb)[0, 1])


def residualize_cols(Y, Xc):
    """Residualize each column of Y on the same column of Xc (nan-aware OLS)."""
    out = np.full_like(Y, np.nan)
    for s in range(Y.shape[1]):
        y, xc = Y[:, s], Xc[:, s]
        m = ~(np.isnan(y) | np.isnan(xc))
        if m.sum() < 10:
            continue
        xm, ym = xc[m], y[m]
        b = np.cov(xm, ym)[0, 1] / np.var(xm)
        out[m, s] = ym - b * (xm - xm.mean())
    return out


def logit_irls(X, y, iters=25):
    beta = np.zeros(X.shape[1])
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-X @ beta))
        W = p * (1 - p)
        H = X.T @ (X * W[:, None])
        g = X.T @ (y - p)
        step = np.linalg.solve(H, g)
        beta += step
        if np.max(np.abs(step)) < 1e-10:
            break
    se = np.sqrt(np.diag(np.linalg.inv(H)))
    return beta, se, p


def run(sample):
    out = {"sample": sample}
    df = common.build_responses(sample)
    cache = np.load(os.path.join("data", "cache", f"designA_{sample}.npz"))

    # ---------- Sec 5: accuracy confound ----------
    conf = {}
    for N in (20, 100):
        blk = {}
        for k in IX.CONFIRMATORY + ["bias3"]:
            full_k = cache[f"full{N}_{k}"]
            full_acc = cache[f"full{N}_acc"]
            blk[k] = {"pearson_vs_acc": nan_pearson(full_k, full_acc),
                      "spearman_vs_acc": nan_spearman(full_k, full_acc)}
            A, B = cache[f"A{N}_{k}"], cache[f"B{N}_{k}"]
            accA, accB = cache[f"A{N}_acc"], cache[f"B{N}_acc"]
            rA = float(np.nanmean(sb_corr_cols(A, B)))
            Ares = residualize_cols(A, accA)
            Bres = residualize_cols(B, accB)
            blk[k]["reliability_raw"] = rA
            blk[k]["reliability_acc_partialed"] = float(
                np.nanmean(sb_corr_cols(Ares, Bres)))
        # top-decile membership change after partialing (bias, full set)
        fb, fa = cache[f"full{N}_bias"], cache[f"full{N}_acc"]
        m = ~(np.isnan(fb) | np.isnan(fa))
        b = np.cov(fa[m], fb[m])[0, 1] / np.var(fa[m])
        resid = fb[m] - b * (fa[m] - fa[m].mean())
        k10 = max(int(round(m.sum() * 0.10)), 1)
        traw = set(np.argsort(-fb[m])[:k10])
        tres = set(np.argsort(-resid)[:k10])
        blk["top_decile_jaccard_raw_vs_partialed_bias"] = \
            len(traw & tres) / len(traw | tres)
        conf[f"N{N}"] = blk
    out["accuracy_confound"] = conf

    # ---------- Sec 6: selection ----------
    sel = {}
    df["tag"] = df["Confidence"].notna() & (df["Confidence"] >= 0) & \
        (df["Confidence"] <= 100)
    g = df.groupby("UserId").agg(
        n=("IsCorrect", "size"), acc=("IsCorrect", "mean"),
        ntag=("tag", "sum"),
        acc_tag=("IsCorrect", lambda s: np.nan),  # placeholder
        first=("DateAnswered", "min"))
    tag_acc = df[df.tag].groupby("UserId")["IsCorrect"].mean()
    untag_acc = df[~df.tag].groupby("UserId")["IsCorrect"].mean()
    g["acc_tag"] = tag_acc.reindex(g.index)
    g["acc_untag"] = untag_acc.reindex(g.index)
    g["reporter"] = (g["ntag"] > 0).astype(int)

    # within-student tagged vs untagged accuracy (>=5 of each)
    w = g[(g.ntag >= 5) & ((g.n - g.ntag) >= 5)].copy()
    diffs = (w.acc_tag - w.acc_untag).to_numpy()
    rng = np.random.default_rng([MASTER_SEED, 80])
    bs = np.array([np.mean(diffs[rng.integers(0, len(diffs), len(diffs))])
                   for _ in range(B_BOOT)])
    sel["within_student"] = {
        "n_students": int(len(w)),
        "mean_acc_tagged_minus_untagged": float(np.mean(diffs)),
        "ci95": [float(v) for v in np.percentile(bs, [2.5, 97.5])]}

    # student-level standardized differences
    def cohens_d(a, b):
        sp = np.sqrt(((len(a) - 1) * np.var(a, ddof=1)
                      + (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
        return float((np.mean(a) - np.mean(b)) / sp)
    rep, non = g[g.reporter == 1], g[g.reporter == 0]
    sel["student_level"] = {
        "n_reporters": int(len(rep)), "n_nonreporters": int(len(non)),
        "acc_mean_rep": float(rep.acc.mean()), "acc_mean_non": float(non.acc.mean()),
        "d_accuracy": cohens_d(rep.acc.to_numpy(), non.acc.to_numpy()),
        "logvol_d": cohens_d(np.log10(rep.n.to_numpy()),
                             np.log10(non.n.to_numpy())),
        "vol_median_rep": float(rep.n.median()),
        "vol_median_non": float(non.n.median())}

    # logistic propensity: reporter ~ z(acc) + z(log vol) + z(first month)
    t0 = g["first"].min()
    month = ((g["first"] - t0).dt.days / 30.44).to_numpy()
    acc = g["acc"].to_numpy()
    lv = np.log10(g["n"].to_numpy())
    Z = np.column_stack([acc, lv, month])
    Zs = (Z - Z.mean(0)) / Z.std(0)
    X = np.column_stack([np.ones(len(g)), Zs])
    y = g["reporter"].to_numpy().astype(float)
    beta, se, p = logit_irls(X, y)
    zstat = beta / se
    sel["propensity_model"] = {
        "coef_std": {"intercept": beta[0], "acc": beta[1],
                     "log10_vol": beta[2], "first_month": beta[3]},
        "se": {"intercept": se[0], "acc": se[1],
               "log10_vol": se[2], "first_month": se[3]},
        "z": {"acc": zstat[1], "log10_vol": zstat[2], "first_month": zstat[3]},
        "note": "covariates standardized; coefficients are per-SD log-odds"}
    for kk in ("coef_std", "se", "z"):
        sel["propensity_model"][kk] = {a: float(b) for a, b in
                                       sel["propensity_model"][kk].items()}

    # IPW sensitivity: Design B N=20 reliability, weighted by 1/p among eligible
    users, offsets, codes, x = student_arrays(df)
    counts = np.diff(offsets)
    pool = counts >= 20
    rngW = np.random.default_rng([MASTER_SEED, 81, {"t12": 1, "t34": 2}[sample]])
    C, Xm, uidp = subsample_matrix(users, offsets, codes, x, pool, 20, rngW)
    A, B, _, _, _ = half_indices(C, Xm, rngW)
    phat = pd.Series(p, index=g.index).reindex(uidp).to_numpy()
    wgt = np.mean(phat) / phat
    wgt = np.minimum(wgt, np.percentile(wgt, 99))

    def weighted_sb(Amat, Bmat, wv):
        rs = []
        for s in range(Amat.shape[1]):
            a, b = Amat[:, s], Bmat[:, s]
            m = ~(np.isnan(a) | np.isnan(b))
            if m.sum() < 10:
                continue
            ww = wv[m] / wv[m].sum()
            ma, mb = (ww * a[m]).sum(), (ww * b[m]).sum()
            cov = (ww * (a[m] - ma) * (b[m] - mb)).sum()
            va = (ww * (a[m] - ma) ** 2).sum()
            vb = (ww * (b[m] - mb) ** 2).sum()
            r = cov / np.sqrt(va * vb)
            rs.append(2 * r / (1 + r))
        return float(np.mean(rs))
    ipw = {}
    for k in IX.CONFIRMATORY:
        ipw[k] = {"unweighted": sb_mean(A[k], B[k]),
                  "ipw": weighted_sb(A[k], B[k], wgt)}
    sel["ipw_designB_N20"] = ipw
    out["selection"] = sel

    # ---------- Sec 3: temporal test-retest ----------
    tg = common.tagged(df).sort_values(["UserId", "DateAnswered"], kind="stable")
    tg = tg.assign(code=(tg.Confidence / 25).round().astype(np.int8))
    first_tag = tg.groupby("UserId")["DateAnswered"].transform("min")
    days = (tg["DateAnswered"] - first_tag).dt.days
    winA = tg[days <= 28]
    winB = tg[days >= 57]

    def win_indices(win):
        res = {}
        for uid, grp in win.groupby("UserId"):
            if len(grp) < 20:
                continue
            cnt = IX.count_tables(grp["code"].to_numpy()[None, :],
                                  grp["IsCorrect"].to_numpy(np.int8)[None, :],
                                  np.ones((1, len(grp)), bool))
            vals = IX.indices_from_counts(cnt)
            res[uid] = {k: float(vals[k][0]) for k in IX.CONFIRMATORY}
            res[uid]["n"] = int(len(grp))
        return pd.DataFrame(res).T
    ia, ib = win_indices(winA), win_indices(winB)
    both = ia.join(ib, how="inner", lsuffix="_a", rsuffix="_b")
    curves = pd.read_csv(os.path.join("results",
                                      f"reliability_curves_{sample}.csv"))
    bcurve = curves[(curves.design == "B") & (curves.split == "random")]

    def interp_rel(k, N):
        sub = bcurve[bcurve["index"] == k].sort_values("N")
        return float(np.interp(np.log(N), np.log(sub.N), sub.r))
    retest = {"n_students": int(len(both)),
              "median_n_winA": float(both["n_a"].median()) if len(both) else None,
              "median_n_winB": float(both["n_b"].median()) if len(both) else None}
    for k in IX.CONFIRMATORY:
        if not len(both):
            break
        r_obs = nan_pearson(both[f"{k}_a"].to_numpy(float),
                            both[f"{k}_b"].to_numpy(float))
        relA = interp_rel(k, both["n_a"].median())
        relB = interp_rel(k, both["n_b"].median())
        retest[k] = {"r_obs": r_obs,
                     "rel_at_median_n": [relA, relB],
                     "r_disattenuated": (r_obs / np.sqrt(relA * relB)
                                         if relA > 0 and relB > 0 else None)}
    out["test_retest"] = retest

    with open(os.path.join("results", f"confound_selection_{sample}.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(sample, "confound/selection/retest written")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
