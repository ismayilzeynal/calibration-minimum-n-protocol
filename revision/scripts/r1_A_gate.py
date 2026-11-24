"""Stream A (E5) regression gate and positive control.

Regression gate: the pre-registered-rule rerun (results/r1_A/prereg) must
reproduce every NON-AUROC output of the baseline rerun (results/baseline,
itself equal to the submitted results up to 3.0e-15 in one curve cell)
with max abs diff <= 1e-12 and an identical NaN pattern. For the confound
and confirmatory JSONs the reference is the submitted paper results/.

Positive control: under the pre-registered rule AUROC must be undefined in
exactly the cells where gamma is undefined (otherwise the patch did not apply
and "no change" would be meaningless).

Writes revision-2026-11/results/r1_A/gate.json.
Usage: python -B r1_A_gate.py
"""
import json
import math
import os

import numpy as np
import pandas as pd

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
BASE = os.path.join(REV, "results", "baseline")
PRE = os.path.join(REV, "results", "r1_A", "prereg")
OUT = os.path.join(REV, "results", "r1_A", "gate.json")
TOL = 1e-12


def flatten(obj, prefix=""):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{prefix}.{k}" if prefix else k)
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{prefix}[{i}]")
    else:
        yield prefix, obj


def cmp_values(a, b, skip):
    keys = sorted(set(a) | set(b))
    maxd, worst, mism, n, skipped = 0.0, None, [], 0, 0
    for k in keys:
        if skip(k):
            skipped += 1
            continue
        va, vb = a.get(k), b.get(k)
        if isinstance(va, bool) or isinstance(vb, bool) or not (
                isinstance(va, (int, float)) and isinstance(vb, (int, float))):
            if va != vb:
                mism.append([k, va, vb])
            continue
        va, vb = float(va), float(vb)
        if math.isnan(va) or math.isnan(vb) or math.isinf(va) or math.isinf(vb):
            if not ((math.isnan(va) and math.isnan(vb)) or va == vb):
                mism.append([k, va, vb])
            continue
        n += 1
        d = abs(va - vb)
        if d > maxd:
            maxd, worst = d, k
    return {"numeric_fields": n, "skipped_auroc_fields": skipped,
            "max_abs_diff": maxd, "argmax": worst,
            "mismatches": mism[:20], "n_mismatches": len(mism),
            "pass": bool(maxd <= TOL and not mism)}


def cmp_json(pa, pb, skip):
    a = dict(flatten(json.load(open(pa))))
    b = dict(flatten(json.load(open(pb))))
    return cmp_values(a, b, skip)


def is_auroc(k):
    return "auroc" in k.lower()


def cmp_curves(pa, pb):
    key = ["sample", "design", "split", "index", "N"]
    a, b = pd.read_csv(pa), pd.read_csv(pb)
    m = a.merge(b, on=key, how="outer", suffixes=("_o", "_n"), indicator=True)
    res = {"rows_unmatched": int((m["_merge"] != "both").sum())}
    na = m[m["index"] != "auroc"]
    worst, nanmis = 0.0, 0
    for col in ["n_pool", "r", "ci_lo", "ci_hi", "pct_undefined", "pct_undefined_full"]:
        o, n = na[f"{col}_o"].to_numpy(float), na[f"{col}_n"].to_numpy(float)
        nanmis += int((np.isnan(o) != np.isnan(n)).sum())
        ok = ~(np.isnan(o) | np.isnan(n))
        if ok.any():
            worst = max(worst, float(np.abs(o[ok] - n[ok]).max()))
    res["non_auroc_rows"] = int(len(na))
    res["non_auroc_max_abs_diff"] = worst
    res["non_auroc_nan_pattern_mismatches"] = nanmis
    res["pass"] = bool(worst <= TOL and nanmis == 0 and res["rows_unmatched"] == 0)
    au = m[m["index"] == "auroc"]
    d = (au["r_n"] - au["r_o"]).to_numpy(float)
    res["auroc_rows"] = int(len(au))
    res["auroc_r_changed_rows"] = int(np.sum(np.abs(np.nan_to_num(d)) > TOL))
    # positive control: prereg AUROC undefined share == gamma undefined share
    g = b[b["index"] == "gamma"].set_index(["design", "split", "N"])
    p = b[b["index"] == "auroc"].set_index(["design", "split", "N"])
    j = g.join(p, lsuffix="_g", rsuffix="_a")
    dd = np.abs(j.pct_undefined_g - j.pct_undefined_a).to_numpy(float)
    ddf = np.abs(j.pct_undefined_full_g - j.pct_undefined_full_a).to_numpy(float)
    res["positive_control_rows"] = int(len(j))
    res["positive_control_max_diff_pct_undefined"] = float(np.nanmax(dd))
    res["positive_control_max_diff_pct_undefined_full"] = float(np.nanmax(ddf)) \
        if np.isfinite(ddf).any() else None
    res["positive_control_pass"] = bool(np.nanmax(dd) == 0 and
                                        (not np.isfinite(ddf).any() or np.nanmax(ddf) == 0))
    return res


def cmp_npz(pa, pb):
    a, b = np.load(pa), np.load(pb)
    worst, nanmis, n = 0.0, 0, 0
    for k in a.files:
        if "auroc" in k:
            continue
        x, y = a[k], b[k]
        if x.dtype.kind in "fc":
            nanmis += int((np.isnan(x) != np.isnan(y)).sum())
            ok = ~(np.isnan(x) | np.isnan(y))
            if ok.any():
                worst = max(worst, float(np.abs(x[ok] - y[ok]).max()))
        else:
            worst = max(worst, 0.0 if np.array_equal(x, y) else math.inf)
        n += 1
    pc = {}
    for N in (20, 100):
        for h in ("A", "B", "full"):
            ga, aa = b[f"{h}{N}_gamma"], b[f"{h}{N}_auroc"]
            pc[f"{h}{N}"] = bool(np.array_equal(np.isnan(ga), np.isnan(aa)))
    return {"arrays_compared": n, "max_abs_diff": worst, "nan_pattern_mismatches": nanmis,
            "pass": bool(worst <= TOL and nanmis == 0),
            "positive_control_auroc_nan_equals_gamma_nan": pc,
            "positive_control_pass": bool(all(pc.values()))}


def main():
    out = {"tolerance": TOL}
    for s in ("t12", "t34"):
        blk = {}
        blk["curves_vs_baseline"] = cmp_curves(
            os.path.join(BASE, "results", f"reliability_curves_{s}.csv"),
            os.path.join(PRE, "results", f"reliability_curves_{s}.csv"))
        blk["minN_vs_baseline"] = cmp_json(
            os.path.join(BASE, "results", f"minN_{s}.json"),
            os.path.join(PRE, "results", f"minN_{s}.json"), is_auroc)
        blk["cache_vs_baseline"] = cmp_npz(
            os.path.join(BASE, "data", "cache", f"designA_{s}.npz"),
            os.path.join(PRE, "data", "cache", f"designA_{s}.npz"))
        cpath = os.path.join(PRE, "results", f"confound_selection_{s}.json")
        if os.path.exists(cpath):
            blk["confound_vs_submitted"] = cmp_json(
                os.path.join(PAPER, "results", f"confound_selection_{s}.json"),
                cpath, is_auroc)
        out[s] = blk
    cf = os.path.join(PRE, "results", "confirmatory_stats.json")
    if os.path.exists(cf):
        out["confirmatory_vs_submitted"] = cmp_json(
            os.path.join(PAPER, "results", "confirmatory_stats.json"), cf, is_auroc)
    passes = []
    for s in ("t12", "t34"):
        for v in out[s].values():
            passes.append(v["pass"])
            if "positive_control_pass" in v:
                passes.append(v["positive_control_pass"])
    if "confirmatory_vs_submitted" in out:
        passes.append(out["confirmatory_vs_submitted"]["pass"])
    out["all_pass"] = bool(all(passes))
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
