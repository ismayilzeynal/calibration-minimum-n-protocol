"""Round-2 additions (addenda 8, 9/4c, 10):
- disjoint full-N flag consistency (Design A pools)
- permutation-null certification of flag metrics (H3 pools)
- descriptives artifact (threshold shares + Table 1 values)

Run:  python scripts/rev_round2.py t12
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from reliability import student_arrays, MASTER_SEED
from h3_shrinkage import half_bias, MIN_TAGGED

R_SPLITS = 50
FLAG_N2 = [10, 15, 20, 30, 50, 100]


def disjoint_flags(users, offsets, codes, x, counts, sid, rng):
    d_all = codes.astype(np.float64) * 0.25 - x
    pool = np.where(counts >= 200)[0]
    n = len(pool)
    k10 = max(int(round(n * 0.10)), 1)
    out = []
    for N in FLAG_N2:
        rngN = np.random.default_rng([MASTER_SEED, 92, N, sid])
        jac, ret = [], []
        for rep in range(R_SPLITS):
            bA = np.empty(n); bB = np.empty(n)
            for r, i in enumerate(pool):
                s, e = offsets[i], offsets[i + 1]
                sel = rngN.choice(e - s, size=2 * N, replace=False) + s
                bA[r] = d_all[sel[:N]].mean()
                bB[r] = d_all[sel[N:]].mean()
            tA = set(np.argsort(-bA)[:k10])
            tB = set(np.argsort(-bB)[:k10])
            jac.append(len(tA & tB) / len(tA | tB))
            ret.append(len(tA & tB) / k10)
        out.append(dict(N=N, jaccard=float(np.mean(jac)),
                        pct_flags_retained=float(np.mean(ret))))
    # analytic chance for two independent decile sets: J = p/(2-p), p=.10
    return out, dict(chance_jaccard=0.1 / 1.9, chance_retention=0.10,
                     n_pool=n)


def null_flags(users, offsets, codes, x, counts, sid):
    rng = np.random.default_rng([MASTER_SEED, 90, sid])  # same as rev_h3null
    keep = np.where(counts >= MIN_TAGGED)[0]
    d_all = codes.astype(np.float64) * 0.25 - x
    blocks = [d_all[offsets[i]:offsets[i + 1]] for i in keep]
    pooled = np.concatenate(blocks)
    perm = rng.permutation(pooled)
    cuts = np.cumsum([len(b) for b in blocks])[:-1]
    d_list = np.split(perm, cuts)
    n = len(d_list)
    k10 = max(int(round(n * 0.10)), 1)
    j_raw, j_ebval = [], []
    for s in range(R_SPLITS):
        a, b, ea, eb = half_bias(d_list, rng)
        tA = set(np.argsort(-a)[:k10]); tB = set(np.argsort(-b)[:k10])
        j_raw.append(len(tA & tB) / len(tA | tB))
        tE = set(np.argsort(-ea)[:k10])
        j_ebval.append(len(tE & tB) / len(tE | tB))
    return dict(null_jaccard_raw_raw=float(np.mean(j_raw)),
                null_jaccard_eb_vs_rawcriterion=float(np.mean(j_ebval)),
                n_students=n)


def descriptives(df):
    tg = common.tagged(df)
    cnt = tg.groupby("UserId").size()
    return {
        "reporters": int(cnt.size),
        "median_tagged": float(cnt.median()),
        "mean_tagged": float(cnt.mean()),
        "share_ge": {str(k): float((cnt >= k).mean())
                     for k in (15, 18, 20, 26, 33, 50, 100, 200)},
        "mean_confidence": float(tg.Confidence.mean()),
        "accuracy_tagged": float(tg.IsCorrect.mean()),
        "accuracy_all": float(df.IsCorrect.mean()),
    }


def run(sample):
    sid = {"t12": 1, "t34": 2}[sample]
    df = common.build_responses(sample)
    users, offsets, codes, x = student_arrays(df)
    counts = np.diff(offsets)
    rng = np.random.default_rng([MASTER_SEED, 93, sid])
    flags, meta = disjoint_flags(users, offsets, codes, x, counts, sid, rng)
    out = {"sample": sample,
           "disjoint_flag_consistency": flags,
           "flag_meta": meta,
           "null_flags": null_flags(users, offsets, codes, x, counts, sid),
           "descriptives": descriptives(df)}
    with open(os.path.join("results", f"rev_round2_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(sample, "round-2 additions written")
    print(" flags:", [(f["N"], round(f["jaccard"], 3),
                       round(f["pct_flags_retained"], 3)) for f in flags])
    print(" null flags:", {k: round(v, 4) if isinstance(v, float) else v
                           for k, v in out["null_flags"].items()})


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
