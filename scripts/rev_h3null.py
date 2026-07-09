"""Addendum 4: H3 certification null. Permute response-level (c - x) values
across students (preserving each student's volume n_i), destroying person
effects; rerun the H3 shrinkage machinery unchanged. Expected EB gain ~ 0.

Run:  python scripts/rev_h3null.py t12
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common
from reliability import student_arrays
from h3_shrinkage import half_bias, spearman_cols, R_SPLITS, MASTER_SEED, \
    MIN_TAGGED


def run(sample):
    rng = np.random.default_rng([MASTER_SEED, 90, {"t12": 1, "t34": 2}[sample]])
    df = common.build_responses(sample)
    users, offsets, codes, x = student_arrays(df)
    counts = np.diff(offsets)
    keep = np.where(counts >= MIN_TAGGED)[0]
    d_all = codes.astype(np.float64) * 0.25 - x
    blocks = [d_all[offsets[i]:offsets[i + 1]] for i in keep]
    pooled = np.concatenate(blocks)
    perm = rng.permutation(pooled)
    sizes = [len(b) for b in blocks]
    cuts = np.cumsum(sizes)[:-1]
    d_list = np.split(perm, cuts)
    n = len(d_list)

    rawA = np.empty((n, R_SPLITS)); rawB = np.empty((n, R_SPLITS))
    ebA = np.empty((n, R_SPLITS)); ebB = np.empty((n, R_SPLITS))
    for s in range(R_SPLITS):
        a, b, ea, eb = half_bias(d_list, rng)
        rawA[:, s], rawB[:, s], ebA[:, s], ebB[:, s] = a, b, ea, eb
    rho_raw = float(spearman_cols(rawA, rawB).mean())
    rho_eb = float(spearman_cols(ebA, ebB).mean())
    rho_eb_val = float(0.5 * (spearman_cols(ebA, rawB)
                              + spearman_cols(ebB, rawA)).mean())
    out = {"sample": sample, "n_students": n,
           "null_rho_raw": rho_raw, "null_rho_eb": rho_eb,
           "null_gain": rho_eb - rho_raw,
           "null_rho_eb_validation": rho_eb_val,
           "note": "person effects destroyed by cross-student permutation of "
                   "(c - x); volumes n_i preserved; same machinery as H3; "
                   "validation metric (EB_A vs raw_B) is expected ~0 under "
                   "this null while the symmetric EB-vs-EB metric is inflated "
                   "by shared shrinkage weights (addendum 4b)"}
    with open(os.path.join("results", f"h3_null_{sample}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(sample, "NULL: rho_raw", round(rho_raw, 4), "rho_eb", round(rho_eb, 4),
          "gain", round(rho_eb - rho_raw, 4))


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t12")
