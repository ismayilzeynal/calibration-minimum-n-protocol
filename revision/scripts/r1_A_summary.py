"""Stream A (E5) assembly: every number the manuscript may quote about the
AUROC rule and about index thresholds by intended use, in one JSON.

Rules compared (same subsamples, same splits, same bootstrap draws):
  auroc_mw     = submitted code (scripts/indices.py line 46): Mann-Whitney
                 AUROC with 0.5 tie credit; undefined only if correctness is
                 constant; exactly .5 if confidence is constant and
                 correctness varies. Numbers = results/baseline (= submitted).
  auroc_prereg = rule stated in analysis_plan.md section 2 and in the
                 submitted manuscript text: undefined if confidence OR
                 correctness is constant. Numbers = results/r1_A/prereg.

Inputs: results/baseline/, results/r1_A/prereg/, results/r1_A/mw/ (downstream
rerun under the submitted rule), results/r1_A/light_*.json,
results/r1_A/indices_check.json, results/r1_A/gate.json.

Writes revision-2026-11/results/r1_A.json and prints a digest.
Usage: python -B r1_A_summary.py
"""
import json
import math
import os

import numpy as np
import pandas as pd

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
RES = os.path.join(REV, "results")
BASE = os.path.join(RES, "baseline")
PRE = os.path.join(RES, "r1_A", "prereg")
R1A = os.path.join(RES, "r1_A")
OUT = os.path.join(RES, "r1_A.json")

GRID = [5, 10, 12, 15, 20, 25, 30, 40, 50, 100, 200]
GRID_B = [5, 10, 20, 50, 100, 200]
B_BOOT = 1000
SAMPLES = {"t12": "S1", "t34": "S2"}
THR = (0.70, 0.80, 0.90)
ALL_IDX = ["bias", "brier", "gamma", "auroc_mw", "auroc_prereg",
           "conf", "acc", "bias3", "cal", "res"]
CONF_IDX = ["bias", "brier", "gamma", "auroc_mw", "auroc_prereg"]


def J(v):
    """JSON-safe number: inf -> '>200', nan -> None."""
    if v is None:
        return None
    if isinstance(v, (bool, str)):
        return v
    if isinstance(v, (int, np.integer)):
        return int(v)
    v = float(v)
    if math.isnan(v):
        return None
    if math.isinf(v):
        return ">200" if v > 0 else "-inf"
    return v


def Jtree(o):
    if isinstance(o, dict):
        return {k: Jtree(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [Jtree(v) for v in o]
    if isinstance(o, (int, float, np.floating, np.integer)) and not isinstance(o, bool):
        return J(o)
    return o


def cross(rvec, thr, grid):
    """Identical logic to reliability.min_n_cross / cross_on."""
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
    """Identical logic to reliability.ci_or_partial."""
    fin = draws[np.isfinite(draws)]
    if len(fin) == B_BOOT:
        return [float(v) for v in np.percentile(fin, [2.5, 97.5])]
    if len(fin):
        return [float(np.percentile(fin, 2.5)), None]
    return [None, None]


def load_curves(sample):
    mw = pd.read_csv(os.path.join(BASE, "results", f"reliability_curves_{sample}.csv"))
    pr = pd.read_csv(os.path.join(PRE, "results", f"reliability_curves_{sample}.csv"))
    pr = pr[pr["index"] == "auroc"].copy()
    pr["index"] = "auroc_prereg"
    mw.loc[mw["index"] == "auroc", "index"] = "auroc_mw"
    return pd.concat([mw, pr], ignore_index=True)


def curve(cv, idx, design, split):
    sub = cv[(cv["index"] == idx) & (cv.design == design)
             & (cv.split == split)].sort_values("N")
    return sub


def share_ge(light, v):
    if v is None or not np.isfinite(v):
        return None
    k = int(math.ceil(v - 1e-9))
    if k > 400:
        return None
    return {"k": k, "share_reporters_ge_k": light["share_reporters_ge_k"][str(k)]}


def main():
    out = {"_definitions": {
        "auroc_mw": "Submitted code rule (scripts/indices.py line 46): "
                    "AUROC = (C + 0.5*ties)/(n0*n1), Mann-Whitney with 0.5 tie "
                    "credit; undefined only when correctness is constant; = .5 "
                    "exactly when confidence is constant and correctness varies. "
                    "Numbers equal "
                    "the submitted results (baseline rerun).",
        "auroc_prereg": "Rule written in analysis_plan.md section 2 and in the "
                        "submitted manuscript text: AUROC undefined when confidence "
                        "OR correctness is constant (then identical to gamma's "
                        "undefined set). Sensitivity rerun, same seeds.",
        "gamma": "Goodman-Kruskal gamma, (C - D)/(C + D); undefined when C + D = 0 "
                 "(constant confidence or constant correctness). Same in both runs.",
        "forced_05": "Cells where confidence is constant and correctness varies: "
                     "AUROC = .5 under auroc_mw, undefined under auroc_prereg; "
                     "equals gamma_undefined - auroc_mw_undefined exactly.",
        "full": "share of students whose index is undefined on their full "
                "N-response set (practitioner-facing incomputability)",
        "halves": "share undefined in split half A (random: students x 50 "
                  "splits; temporal: earlier half), which governs exclusion",
        "minN": "N at which the reliability curve first reaches the threshold, "
                "log-linear interpolation between grid points; '>200' = not "
                "reached within the grid; CI = bootstrap percentile over students "
                "(B = 1000, joint across the Design A grid); a CI of [x, None] "
                "means some draws never cross (upper bound > 200)",
        "S1": "t12, primary sample", "S2": "t34, internal replication",
    }}

    ic = json.load(open(os.path.join(R1A, "indices_check.json")))
    out["code_treatment"] = {
        "conclusion": ic["conclusion"],
        "hand_cases": ic["hand_cases"],
        "random_vectors": ic["random_vectors"],
        "real_S2_full_sets": ic["real_S2_full_sets"],
        "manuscript_text_claim": "Gamma and AUROC are undefined whenever a student "
                                 "shows no variance in confidence or in correctness "
                                 "(Method, Indices and the Undefined-Value Protocol)",
        "text_vs_code": "text and plan state auroc_prereg; code computed auroc_mw; "
                        "all submitted AUROC numbers are auroc_mw numbers"}
    gate = json.load(open(os.path.join(R1A, "gate.json")))
    out["regression_gate"] = gate
    mwcmp = json.load(open(os.path.join(R1A, "mw", "compare_submitted.json")))
    out["regression_gate_mw_downstream_vs_submitted"] = mwcmp

    selfcheck = {"minN_point_max_abs_diff": 0.0, "minN_ci_max_abs_diff": 0.0,
                 "n_compared": 0}
    for s, S in SAMPLES.items():
        light = json.load(open(os.path.join(R1A, f"light_{s}.json")))
        cv = load_curves(s)
        mn = {"mw": json.load(open(os.path.join(BASE, "results", f"minN_{s}.json"))),
              "prereg": json.load(open(os.path.join(PRE, "results", f"minN_{s}.json")))}
        npz = {"mw": np.load(os.path.join(BASE, "data", "cache", f"designA_{s}.npz")),
               "prereg": np.load(os.path.join(PRE, "data", "cache", f"designA_{s}.npz"))}
        blk = {"n_poolA": mn["mw"]["_meta"]["n_poolA"], "reporters": light["reporters"]}

        # ---------- incomputability ----------
        inc = {"designA_full": {}, "designA_halves_random": {},
               "designA_halves_temporal": {}, "designB_full": {},
               "designB_halves_random": {}, "designB_halves_temporal": {}}
        for N, v in light["designA_full"].items():
            inc["designA_full"][N] = {
                "gamma": v["gamma_undefined"],
                "auroc_mw": v["auroc_undefined_submitted_rule"],
                "auroc_prereg": v["auroc_undefined_prereg_rule"],
                "forced_05": v["auroc_forced_05_share"]}
        for N, v in light["designB_full"].items():
            inc["designB_full"][N] = {
                "n_pool": v["n_students"],
                "gamma": v["gamma_undefined"],
                "auroc_mw": v["auroc_undefined_submitted_rule"],
                "auroc_prereg": v["auroc_undefined_prereg_rule"],
                "forced_05": v["auroc_forced_05_share"]}
        for design, split, key in (("A", "random", "designA_halves_random"),
                                   ("A", "temporal", "designA_halves_temporal"),
                                   ("B", "random", "designB_halves_random"),
                                   ("B", "temporal", "designB_halves_temporal")):
            g = curve(cv, "gamma", design, split).set_index("N")
            a = curve(cv, "auroc_mw", design, split).set_index("N")
            p = curve(cv, "auroc_prereg", design, split).set_index("N")
            for N in g.index:
                inc[key][str(N)] = {
                    "gamma": float(g.loc[N, "pct_undefined"]),
                    "auroc_mw": float(a.loc[N, "pct_undefined"]),
                    "auroc_prereg": float(p.loc[N, "pct_undefined"]),
                    "forced_05": float(g.loc[N, "pct_undefined"] - a.loc[N, "pct_undefined"])}
        blk["incomputability"] = inc
        blk["forced_05_persistence_designA"] = light["forced_persistence_designA"]

        # ---------- curves ----------
        cblk = {}
        for design, split in (("A", "random"), ("A", "temporal"),
                              ("B", "random"), ("B", "temporal")):
            d = {}
            for idx in ALL_IDX:
                sub = curve(cv, idx, design, split)
                d[idx] = {str(int(r.N)): {"r": r.r, "ci_lo": r.ci_lo, "ci_hi": r.ci_hi,
                                          "n_pool": int(r.n_pool)}
                          for r in sub.itertuples()}
            cblk[f"design{design}_{split}"] = d
        blk["curves"] = cblk

        # ---------- minN (all fields) for confirmatory indices ----------
        mblk = {}
        for idx in CONF_IDX:
            rule = "prereg" if idx == "auroc_prereg" else "mw"
            key = "auroc" if idx.startswith("auroc") else idx
            mblk[idx] = mn[rule][key]
        blk["minN"] = mblk

        # ---------- thresholds by index, incl. .90, with self-check ----------
        th = {}
        for idx in ALL_IDX:
            rule = "prereg" if idx == "auroc_prereg" else "mw"
            key = "auroc" if idx.startswith("auroc") else idx
            t = {}
            for design, split, grid in (("A", "random", GRID), ("A", "temporal", GRID),
                                        ("B", "random", GRID_B), ("B", "temporal", GRID_B)):
                sub = curve(cv, idx, design, split)
                assert list(sub.N) == grid, (idx, design, split)
                rv = sub.r.to_numpy(float)
                tt = {}
                for thr in THR:
                    pt = cross(rv, thr, grid)
                    e = {"point": pt, "at_or_below_first_grid_point": bool(pt == grid[0]),
                         "coverage": share_ge(light, pt)}
                    if design == "A" and key in ("bias", "brier", "gamma", "auroc"):
                        if split == "random":
                            draws = np.array([cross(npz[rule][f"boot_{key}"][b], thr, GRID)
                                              for b in range(B_BOOT)])
                            e["ci"] = ci_or_partial(draws)
                            e["share_draws_gt200"] = float(np.mean(~np.isfinite(draws)))
                            if thr in (0.70, 0.80):
                                ref = mn[rule][key][f"minN{int(thr * 100)}_ci"]
                                for u, w in zip(e["ci"], ref):
                                    if u is not None and w is not None:
                                        selfcheck["minN_ci_max_abs_diff"] = max(
                                            selfcheck["minN_ci_max_abs_diff"], abs(u - w))
                                    elif (u is None) != (w is None):
                                        selfcheck["minN_ci_max_abs_diff"] = math.inf
                        elif thr in (0.70, 0.80):
                            e["ci"] = mn[rule][key][f"minN{int(thr * 100)}_temporal_ci"]
                        else:
                            e["ci"] = "not computed (temporal bootstrap draws are not cached)"
                    if key in ("bias", "brier", "gamma", "auroc") and thr in (0.70, 0.80):
                        fld = {("A", "random"): f"minN{int(thr * 100)}_point",
                               ("A", "temporal"): f"minN{int(thr * 100)}_temporal_point",
                               ("B", "random"): f"minN{int(thr * 100)}_designB_point"}.get(
                            (design, split))
                        if fld:
                            ref = mn[rule][key][fld]
                            dd = 0.0 if (math.isinf(ref) and math.isinf(pt)) else abs(ref - pt)
                            selfcheck["minN_point_max_abs_diff"] = max(
                                selfcheck["minN_point_max_abs_diff"], dd)
                            selfcheck["n_compared"] += 1
                    tt[f"{thr:.2f}"] = e
                t[f"design{design}_{split}"] = tt
            th[idx] = t
        blk["thresholds_by_index"] = th

        # ---------- downstream analyses ----------
        cs_mw = json.load(open(os.path.join(R1A, "mw", "results", f"confound_selection_{s}.json")))
        cs_pr = json.load(open(os.path.join(PRE, "results", f"confound_selection_{s}.json")))
        dn = {}
        for N in ("N20", "N100"):
            dn[N] = {"auroc_mw": cs_mw["accuracy_confound"][N]["auroc"],
                     "auroc_prereg": cs_pr["accuracy_confound"][N]["auroc"],
                     "gamma": cs_mw["accuracy_confound"][N]["gamma"],
                     "brier": cs_mw["accuracy_confound"][N]["brier"],
                     "bias": cs_mw["accuracy_confound"][N]["bias"]}
        blk["accuracy_confound"] = dn
        blk["ipw_designB_N20"] = {"auroc_mw": cs_mw["selection"]["ipw_designB_N20"]["auroc"],
                                  "auroc_prereg": cs_pr["selection"]["ipw_designB_N20"]["auroc"],
                                  "gamma": cs_mw["selection"]["ipw_designB_N20"]["gamma"]}
        blk["test_retest"] = {"n_students": cs_mw["test_retest"]["n_students"],
                              "auroc_mw": cs_mw["test_retest"]["auroc"],
                              "auroc_prereg": cs_pr["test_retest"]["auroc"],
                              "gamma": cs_mw["test_retest"]["gamma"],
                              "brier": cs_mw["test_retest"]["brier"],
                              "bias": cs_mw["test_retest"]["bias"]}

        # ---------- qualitative checks ----------
        q = {}
        for rule_idx in ("auroc_mw", "auroc_prereg"):
            ra = curve(cv, rule_idx, "A", "random").r.to_numpy(float)
            rg = curve(cv, "gamma", "A", "random").r.to_numpy(float)
            q[f"{rule_idx}_r_above_gamma_at_all_designA_N"] = bool(np.all(ra > rg))
            q[f"{rule_idx}_minus_gamma_r_by_N"] = {str(N): float(a - g)
                                                  for N, a, g in zip(GRID, ra, rg)}
            m70 = {k: th[k]["designA_random"]["0.70"]["point"]
                   for k in ("bias", "brier", "gamma", rule_idx)}
            q[f"order_minN70_designA_with_{rule_idx}"] = sorted(
                m70, key=lambda k: (m70[k], k))
            q[f"{rule_idx}_reaches_070_within_200_designA"] = bool(
                np.isfinite(th[rule_idx]["designA_random"]["0.70"]["point"]))
            q[f"{rule_idx}_reaches_080_within_200_any_design"] = bool(any(
                np.isfinite(th[rule_idx][d]["0.80"]["point"]) for d in th[rule_idx]))
        q["gamma_reaches_070_within_200_designA"] = bool(
            np.isfinite(th["gamma"]["designA_random"]["0.70"]["point"]))
        blk["qualitative_checks"] = q
        out[S] = blk

    # ---------- confirmatory (H2 inner Holm includes AUROC) ----------
    cf_mw = json.load(open(os.path.join(R1A, "mw", "results", "confirmatory_stats.json")))
    cf_pr = json.load(open(os.path.join(PRE, "results", "confirmatory_stats.json")))
    out["confirmatory"] = {
        rule: {"H2": c["H2"], "holm_family": c["holm_family"], "supported": c["supported"],
               "replication_S2_b": c["replication_S2"]["b_minN70_gt_median31_for_3of4"]}
        for rule, c in (("auroc_mw", cf_mw), ("auroc_prereg", cf_pr))}
    out["threshold_selfcheck_vs_minN_json"] = selfcheck

    rt = {}
    for f in sorted(os.listdir(PRE)):
        if f.startswith("runtime_"):
            rt[f[8:-5]] = json.load(open(os.path.join(PRE, f)))["wall_seconds"]
    rt["indices_check"] = ic["wall_seconds"]
    rt["mw_downstream_confound_both_plus_confirmatory"] = mwcmp["wall_seconds"]
    for s in SAMPLES:
        rt[f"light_{s}"] = json.load(open(os.path.join(R1A, f"light_{s}.json")))["wall_seconds"]
    rt["baseline_reliability_t12_from_REPRO"] = 342.8
    rt["baseline_reliability_t34_from_REPRO"] = 56.6
    out["runtime_seconds"] = rt

    out = Jtree(out)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("selfcheck", selfcheck)
    print("gate all_pass", gate.get("all_pass"))


if __name__ == "__main__":
    main()
