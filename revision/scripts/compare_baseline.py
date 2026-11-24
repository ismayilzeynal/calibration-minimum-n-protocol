"""Compare the baseline rerun (revision-2026-11/results/baseline/results/)
with the submitted results (paper results/). Prints and writes a JSON of
max absolute differences per file and field.

Usage: python compare_baseline.py
"""
import json
import math
import os

import numpy as np
import pandas as pd

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
ORIG = os.path.join(PAPER, "results")
NEW = os.path.join(REV, "results", "baseline", "results")
OUT = os.path.join(REV, "results", "baseline", "compare_baseline.json")


def flatten(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{prefix}.{k}" if prefix else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{prefix}[{i}]")
    else:
        yield prefix, obj


def num(v):
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    return v


def cmp_json(a_path, b_path):
    a = dict(flatten(json.load(open(a_path))))
    b = dict(flatten(json.load(open(b_path))))
    keys = sorted(set(a) | set(b))
    maxd, worst, mism = 0.0, None, []
    n_num = 0
    for k in keys:
        va, vb = num(a.get(k)), num(b.get(k))
        if isinstance(va, float) and isinstance(vb, float):
            if math.isinf(va) or math.isinf(vb) or math.isnan(va) or math.isnan(vb):
                same = (math.isinf(va) and math.isinf(vb) and va == vb) or \
                       (math.isnan(va) and math.isnan(vb))
                if not same:
                    mism.append([k, va, vb])
                continue
            n_num += 1
            d = abs(va - vb)
            if d > maxd:
                maxd, worst = d, k
        elif va != vb:
            mism.append([k, va, vb])
    return {"numeric_fields": n_num, "max_abs_diff": maxd, "argmax_field": worst,
            "nonnumeric_or_special_mismatches": mism[:50],
            "n_nonnumeric_or_special_mismatches": len(mism)}


def cmp_curves(a_path, b_path):
    key = ["sample", "design", "split", "index", "N"]
    a = pd.read_csv(a_path)
    b = pd.read_csv(b_path)
    m = a.merge(b, on=key, how="outer", suffixes=("_o", "_n"), indicator=True)
    out = {"rows_orig": int(len(a)), "rows_new": int(len(b)),
           "rows_unmatched": int((m["_merge"] != "both").sum()), "fields": {}}
    for col in ["n_pool", "r", "ci_lo", "ci_hi", "pct_undefined", "pct_undefined_full"]:
        if f"{col}_o" not in m:
            continue
        o, n = m[f"{col}_o"].to_numpy(float), m[f"{col}_n"].to_numpy(float)
        nan_mismatch = int((np.isnan(o) != np.isnan(n)).sum())
        ok = ~(np.isnan(o) | np.isnan(n))
        d = np.abs(o[ok] - n[ok])
        rec = {"compared": int(ok.sum()), "nan_pattern_mismatches": nan_mismatch,
               "max_abs_diff": float(d.max()) if d.size else None}
        if d.size and d.max() > 0:
            i = np.argmax(d)
            row = m[ok].iloc[i]
            rec["argmax_row"] = {k: (row[k].item() if hasattr(row[k], "item") else row[k])
                                 for k in key}
        out["fields"][col] = rec
    return out


def main():
    res = {}
    for s in ("t12", "t34"):
        mj_new = os.path.join(NEW, f"minN_{s}.json")
        cv_new = os.path.join(NEW, f"reliability_curves_{s}.csv")
        if os.path.exists(mj_new):
            res[f"minN_{s}.json"] = cmp_json(os.path.join(ORIG, f"minN_{s}.json"), mj_new)
        if os.path.exists(cv_new):
            res[f"reliability_curves_{s}.csv"] = cmp_curves(
                os.path.join(ORIG, f"reliability_curves_{s}.csv"), cv_new)
    with open(OUT, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
