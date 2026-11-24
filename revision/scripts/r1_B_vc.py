"""Stream B, part 2: variance-components (one-way random-effects) reliability
for the additive indices (signed bias, Brier) on the SAME Design A subsamples
as the submitted split-half curves, plus SEM.

Item-level scores per tagged response (c = LEVELS[code], x = correctness):
  d = c - x       (signed bias is the mean of d)
  q = (c - x)^2   (Brier is the mean of q)
Responses are nested in students (each student answers different items), so
the design is r:p (one-way). For a balanced n x N matrix:
  MSB = N * var(row means, ddof 1), MSW = pooled within-student variance
  sigma_p^2 = (MSB - MSW) / N, sigma_e^2 = MSW
  ICC(1,N) = (MSB - MSW) / MSB = sigma_p^2 / (sigma_p^2 + sigma_e^2 / N)
Implied n for threshold t: n_t = t / (1 - t) * sigma_e^2 / sigma_p^2.

Gates (recorded in the output):
  G1 split-half r recomputed on the same generator == baseline curves (<= 1e-12)
  G2 row means of d and q == full-set bias and Brier from indices.py (<= 1e-12)
  G3 re-implemented min_n_cross reproduces baseline minN70/80_point

Usage:  python -B r1_B_vc.py t12   (or t34)
Writes: revision-2026-11/results/r1_B_parts/vc_<sample>.json
"""
import json
import os
import sys
import time

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUTDIR = os.path.join(REV, "results", "r1_B_parts")
BASE = os.path.join(REV, "results", "baseline", "results")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common  # noqa: E402
import indices as IX  # noqa: E402
import reliability as REL  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")

GRID = REL.GRID
B_BOOT = 1000
THRS = (0.70, 0.80)
VC_KEYS = ("bias", "brier")
CI_N = (15, 20, 30, 50)


def min_n_cross(rvec, thr, grid=GRID):
    """Exact re-implementation of reliability.run().min_n_cross."""
    lg = np.log(grid)
    rv = np.asarray(rvec, float)
    if np.all(np.isnan(rv)):
        return np.nan
    for i, N in enumerate(grid):
        if not np.isnan(rv[i]) and rv[i] >= thr:
            if i == 0:
                return float(grid[0])
            r0, r1 = rv[i - 1], rv[i]
            if np.isnan(r0) or r1 == r0:
                return float(N)
            f = (thr - r0) / (r1 - r0)
            return float(np.exp(lg[i - 1] + f * (lg[i] - lg[i - 1])))
    return np.inf


def ci_or_partial(draws):
    """Same rule as reliability.run().ci_or_partial (finite-draw percentiles)."""
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
    if np.isnan(v):
        return None
    if np.isinf(v):
        return "inf"
    return v


def oneway_from_stats(m, ssw, N):
    n = len(m)
    msb = N * np.var(m, ddof=1)
    msw = ssw.sum() / (n * (N - 1))
    sp2 = (msb - msw) / N
    icc = (msb - msw) / msb
    return msb, msw, sp2, icc


def run(sample):
    t0 = time.time()
    sid = {"t12": 1, "t34": 2}[sample]
    df = common.build_responses(sample)
    users, offsets, codes, x = REL.student_arrays(df)
    del df
    counts = np.diff(offsets)
    poolA = counts >= max(GRID)
    nA = int(poolA.sum())

    base = pd.read_csv(os.path.join(BASE, f"reliability_curves_{sample}.csv"))
    with open(os.path.join(BASE, f"minN_{sample}.json")) as fh:
        base_minn = json.load(fh)

    def base_curve(k, split="random"):
        s = base[(base.design == "A") & (base.split == split)
                 & (base["index"] == k)].sort_values("N")
        return s

    # G3: re-implemented crossing rule vs baseline
    g3 = {}
    for k in IX.CONFIRMATORY:
        pt = base_curve(k).r.to_numpy()
        for thr, key in ((0.70, "minN70_point"), (0.80, "minN80_point")):
            mine = min_n_cross(pt, thr)
            ref = base_minn[k][key]
            ref = np.inf if ref in (None, "Infinity") else float(ref)
            same = (np.isinf(mine) and np.isinf(ref)) or abs(mine - ref) <= 1e-12
            g3[f"{k}_{key}"] = dict(reimpl=fnum(mine), baseline=fnum(ref),
                                    match=bool(same))

    stats = {}      # N -> key -> (m, ssw)
    out_N = []
    g1_max = {k: 0.0 for k in VC_KEYS}
    g2_max = {k: 0.0 for k in VC_KEYS}
    for N in GRID:
        rngN = np.random.default_rng([REL.MASTER_SEED, 5, N, sid])
        C, X, _ = REL.subsample_matrix(users, offsets, codes, x, poolA, N, rngN)
        A, Bh, TA, TB, full = REL.half_indices(C, X, rngN)
        c = IX.LEVELS[C.astype(np.int64)]
        d = c - X
        q = d ** 2
        row = {"N": N, "n_pool": nA}
        for k, Y in (("bias", d), ("brier", q)):
            m = Y.mean(axis=1)
            ssw = ((Y - m[:, None]) ** 2).sum(axis=1)
            stats.setdefault(N, {})[k] = (m, ssw)
            g2_max[k] = max(g2_max[k], float(np.max(np.abs(m - full[k]))))
            r_split = REL.sb_mean(A[k], Bh[k])
            r_base = float(base_curve(k).set_index("N").loc[N, "r"])
            g1_max[k] = max(g1_max[k], abs(r_split - r_base))
            msb, msw, sp2, icc = oneway_from_stats(m, ssw, N)
            sd_obs = float(np.std(m, ddof=1))
            s_i = np.sqrt(ssw / (N - 1))
            vc = dict(
                split_half_sb=r_split,
                split_half_ci_baseline=[float(base_curve(k).set_index("N").loc[N, "ci_lo"]),
                                        float(base_curve(k).set_index("N").loc[N, "ci_hi"])],
                MSB=float(msb), MSW=float(msw),
                sigma_p2=float(sp2), sigma_e2=float(msw),
                icc1_single=float(sp2 / (sp2 + msw)) if sp2 > 0 else None,
                icc1_N=float(icc),
                gap_icc_minus_split=float(icc - r_split),
                implied_n70=fnum(0.70 / 0.30 * msw / sp2) if sp2 > 0 else None,
                implied_n80=fnum(0.80 / 0.20 * msw / sp2) if sp2 > 0 else None,
                sd_obs=sd_obs,
                sem_split=float(sd_obs * np.sqrt(max(1.0 - r_split, 0.0))),
                sem_vc=float(np.sqrt(msw / N)),
                median_ci_halfwidth=float(np.median(1.96 * s_i / np.sqrt(N))),
                share_sd_zero=float(np.mean(s_i == 0)),
            )
            row[k] = vc
        out_N.append(row)
        print(f"[{sample}] N={N} done ({time.time() - t0:.1f}s)", flush=True)
    t_curves = time.time() - t0

    # joint student bootstrap across the grid for ICC(1,N)
    brng = np.random.default_rng([REL.MASTER_SEED, 101, sid])
    boot = {k: np.full((B_BOOT, len(GRID)), np.nan) for k in VC_KEYS}
    boot_sp_nonpos = {k: 0 for k in VC_KEYS}
    boot_n70_20 = {k: np.full(B_BOOT, np.nan) for k in VC_KEYS}
    for b in range(B_BOOT):
        idx = brng.integers(0, nA, nA)
        for gi, N in enumerate(GRID):
            for k in VC_KEYS:
                m, ssw = stats[N][k]
                msb, msw, sp2, icc = oneway_from_stats(m[idx], ssw[idx], N)
                boot[k][b, gi] = icc
                if sp2 <= 0:
                    boot_sp_nonpos[k] += 1
                elif N == 20:
                    boot_n70_20[k][b] = 0.70 / 0.30 * msw / sp2
    t_boot = time.time() - t0 - t_curves

    curves = {}
    for k in VC_KEYS:
        icc_pt = np.array([r[k]["icc1_N"] for r in out_N])
        sh_pt = np.array([r[k]["split_half_sb"] for r in out_N])
        lo = np.percentile(boot[k], 2.5, axis=0)
        hi = np.percentile(boot[k], 97.5, axis=0)
        for gi, r in enumerate(out_N):
            r[k]["icc1_N_ci"] = [float(lo[gi]), float(hi[gi])]
        d70 = np.array([min_n_cross(boot[k][b], 0.70) for b in range(B_BOOT)])
        d80 = np.array([min_n_cross(boot[k][b], 0.80) for b in range(B_BOOT)])
        gaps = icc_pt - sh_pt
        odd = np.array([N % 2 == 1 for N in GRID])
        r20 = base_minn[k]["r20_point"]
        curves[k] = dict(
            icc_minN70_point=fnum(min_n_cross(icc_pt, 0.70)),
            icc_minN80_point=fnum(min_n_cross(icc_pt, 0.80)),
            icc_minN70_ci=ci_or_partial(d70),
            icc_minN80_ci=ci_or_partial(d80),
            icc_share_minN70_gt200=float(np.mean(~np.isfinite(d70))),
            split_minN70_point=fnum(base_minn[k]["minN70_point"]),
            split_minN80_point=fnum(base_minn[k]["minN80_point"]),
            split_minN70_ci=base_minn[k]["minN70_ci"],
            split_minN80_ci=base_minn[k]["minN80_ci"],
            vc_implied_n70_at20=out_N[GRID.index(20)][k]["implied_n70"],
            vc_implied_n80_at20=out_N[GRID.index(20)][k]["implied_n80"],
            vc_implied_n70_at20_ci=[fnum(v) for v in np.nanpercentile(boot_n70_20[k], [2.5, 97.5])],
            split_prophecy_n70_from20=fnum(base_minn[k]["prophecy_minN70_from20"]),
            split_prophecy_n80_from20_computed_here=fnum(
                20 * (0.80 / 0.20) * ((1 - r20) / r20)) if r20 > 0 else None,
            max_abs_gap_icc_vs_split=float(np.max(np.abs(gaps))),
            N_at_max_abs_gap=int(GRID[int(np.argmax(np.abs(gaps)))]),
            mean_gap_icc_minus_split=float(np.mean(gaps)),
            mean_abs_gap_odd_N=float(np.mean(np.abs(gaps[odd]))),
            mean_abs_gap_even_N=float(np.mean(np.abs(gaps[~odd]))),
            boot_draws_sigma_p2_nonpositive=int(boot_sp_nonpos[k]),
        )

    # all-record moment estimator (pool students, all their tagged responses)
    pool_idx = np.where(poolA)[0]
    rec = {}
    for k in VC_KEYS:
        mi = np.empty(nA)
        vi = np.empty(nA)
        ni = np.empty(nA)
        for j, i in enumerate(pool_idx):
            s, e = offsets[i], offsets[i + 1]
            y = IX.LEVELS[codes[s:e].astype(np.int64)] - x[s:e]
            if k == "brier":
                y = y ** 2
            mi[j] = y.mean()
            vi[j] = y.var(ddof=1)
            ni[j] = e - s
        rec[k] = (mi, vi, ni)
    rrng = np.random.default_rng([REL.MASTER_SEED, 102, sid])
    record = {}
    for k in VC_KEYS:
        mi, vi, ni = rec[k]

        def comp(ix):
            sp2 = np.var(mi[ix], ddof=1) - np.mean(vi[ix] / ni[ix])
            se2 = np.mean(vi[ix])
            return sp2, se2
        sp2, se2 = comp(np.arange(nA))
        draws70 = np.full(B_BOOT, np.nan)
        draws80 = np.full(B_BOOT, np.nan)
        nonpos = 0
        for b in range(B_BOOT):
            ix = rrng.integers(0, nA, nA)
            a, e2 = comp(ix)
            if a <= 0:
                nonpos += 1
                continue
            draws70[b] = 0.70 / 0.30 * e2 / a
            draws80[b] = 0.80 / 0.20 * e2 / a
        rel_at = {str(N): float(sp2 / (sp2 + se2 / N)) for N in GRID}
        record[k] = dict(
            n_students=nA,
            median_responses=float(np.median(ni)),
            sigma_p2=float(sp2), sigma_e2=float(se2),
            implied_reliability_at_N=rel_at,
            implied_n70=fnum(0.70 / 0.30 * se2 / sp2) if sp2 > 0 else None,
            implied_n80=fnum(0.80 / 0.20 * se2 / sp2) if sp2 > 0 else None,
            implied_n70_ci=[fnum(v) for v in np.nanpercentile(draws70, [2.5, 97.5])],
            implied_n80_ci=[fnum(v) for v in np.nanpercentile(draws80, [2.5, 97.5])],
            boot_draws_sigma_p2_nonpositive=nonpos,
        )

    res = dict(
        sample=sample, design="A", n_poolA=nA, grid=GRID, B=B_BOOT,
        seeds=dict(subsample_and_splits="default_rng([42, 5, N, sid]) (same as reliability.py)",
                   icc_joint_bootstrap=f"default_rng([42, 101, {sid}])",
                   record_bootstrap=f"default_rng([42, 102, {sid}])"),
        gates=dict(
            G1_split_half_vs_baseline_max_abs_diff={k: g1_max[k] for k in VC_KEYS},
            G1_pass=all(v <= 1e-12 for v in g1_max.values()),
            G2_rowmean_vs_indices_max_abs_diff={k: g2_max[k] for k in VC_KEYS},
            G2_pass=all(v <= 1e-12 for v in g2_max.values()),
            G3_min_n_cross_reimpl=g3,
            G3_pass=all(v["match"] for v in g3.values()),
        ),
        by_N=out_N,
        summary=curves,
        all_record_moment_estimator=record,
        runtime_seconds=dict(curves_and_gates=t_curves, bootstrap=t_boot,
                             total=time.time() - t0),
    )
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, f"vc_{sample}.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"[{sample}] gates G1={res['gates']['G1_pass']} G2={res['gates']['G2_pass']} "
          f"G3={res['gates']['G3_pass']}; total {res['runtime_seconds']['total']:.1f}s")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
