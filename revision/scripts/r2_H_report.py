"""Assemble the Monte Carlo precision results (r2_E_parts) into
results/r2_E_mc.json and results/r2_E_mc.md, and write results/r2_F_grid.md
from results/r2_F_grid.json. No computation beyond selection, comparison
with the printed intervals and formatting.

Usage: python -B r2_H_report.py   (after r2_E_mc_designA/B, r2_E_mc_structured and r2_F_grid)
"""
import json
import os
import sys

import numpy as np

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_E_common as EC  # noqa: E402

RES = EC.RESULTS
SAMPLES = (("t12", "S1"), ("t34", "S2"))
NAMES = {"bias": "Bias", "brier": "Brier", "gamma": "Gamma", "auroc_pr": "AUROC, plan rule",
         "auroc": "AUROC, tie convention", "conf": "Mean confidence", "acc": "Accuracy",
         "bias_g": "Bias, guess-corrected coding", "brier_g": "Brier, guess-corrected coding",
         "pool_sensitivity_bias": "Bias, without near-perfect and constant raters",
         "pool_sensitivity_brier": "Brier, without near-perfect and constant raters"}
VNAMES = {"random": "Random halves (reference)", "strand": "Halves by content strand", "topic": "Halves by topic",
          "term": "Halves by approximate school term", "difficulty": "Harder versus easier items",
          "w_diff_random": "Weighted to the unrated difficulty mix",
          "w_null_diff_random": "Same weights shuffled within student (difficulty)",
          "w_diffstrand_random": "Weighted to the unrated difficulty-by-strand mix",
          "w_null_diffstrand_random": "Same weights shuffled within student (difficulty by strand)",
          "balanced_random": "Difficulty-balanced subset, equal weights"}


def load(p):
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def val(v):
    return EC.BIG if v == "inf" else (np.nan if v is None else float(v))


def recentred(p, cis, centre):
    """Printed CI shifted so that the draw centre takes the place of the printed
    point; True when the printed point lies inside the shifted interval."""
    if p in ("inf", None) or centre in ("inf", None):
        return None
    lo, hi = val(cis[0]), val(cis[1] if cis[1] is not None else "inf")
    c = float(centre)
    return bool(c - (float(p) - lo) <= float(p) <= c + (hi - float(p)))


def c1(v):
    if v is None:
        return "n/a"
    if v == "inf" or (isinstance(v, float) and v > 200):
        return "> 200"
    return f"{v:.1f}"


def c2(v):
    return "n/a" if v is None else f"{v:.2f}"


def r3(v):
    return "n/a" if v is None else f"{v:.3f}"


def r4(v):
    return "n/a" if v is None else f"{v:.4f}"


def ci(lst):
    return f"[{c1(lst[0])}, {c1(lst[1])}]"


def main():
    A = {lab: load(os.path.join(EC.PARTS, f"designA_{s}.json")) for s, lab in SAMPLES}
    Bd = {lab: load(os.path.join(EC.PARTS, f"designB_{s}.json")) for s, lab in SAMPLES}
    T6 = {lab: load(os.path.join(EC.PARTS, f"structured_{s}.json")) for s, lab in SAMPLES}
    R2A = load(os.path.join(RES, "r2_A.json"))
    VC = {lab: load(os.path.join(RES, "r1_B_parts", f"vc_{s}.json")) for s, lab in SAMPLES}
    for d in list(A.values()) + list(Bd.values()) + list(T6.values()):
        assert d["regression_pass"]
    K = A["S1"]["K_further_draws"]
    KB = Bd["S1"]["K_further_draws"]

    # ---------------- split error ----------------
    split = {}
    for lab in ("S1", "S2"):
        se = A[lab]["split_error"]
        blk = {"designA": {}, "designB": {}, "structured_halves_bias_brier": {}}
        for name, v in se.items():
            blk["designA"][name] = {"per_N": v["per_N_split_se"], "max_all_N": v["max_split_se_all_N"],
                                    "N_at_max": v["argmax_N"], "max_N_ge_20": v["max_split_se_N_ge_20"],
                                    "max_N_ge_50": v["max_split_se_N_ge_50"]}
        for name, v in Bd[lab]["split_error"].items():
            blk["designB"][name] = {"per_N": v["per_N_split_se"], "max_all_N": v["max_split_se_all_N"],
                                    "N_at_max": v["argmax_N"]}
        for vname, vv in T6[lab]["variants"].items():
            for k, b in vv.items():
                blk["structured_halves_bias_brier"][f"{vname}_{k}"] = {"max_all_N": b["max_split_se_all_N"],
                                                            "r20": b["r20"], "r20_split_se": b["r20_split_se"]}
        blk["no_split_component"] = ["temporal halves (one deterministic split per subsample)",
                                     "variance-components ICC(1,N) curve (no split)",
                                     "split-free all-answer estimator (no subsample, no split)"]
        split[lab] = blk

    # crossing rows (split SE and subsample draws together)
    rows = []
    flags_ci = []
    tails = []

    def drawinfo(dc):
        m = dc["mean"] if dc["mean"] is not None else None
        centre = m if m is not None else dc["median"]
        return m, centre

    for lab in ("S1", "S2"):
        a = A[lab]
        r2a = R2A[lab]
        tab = f"{lab} minimum-N table"
        for field in ("random", "temporal"):
            for k in ("bias", "brier", "gamma", "auroc_pr", "auroc", "bias_g", "brier_g"):
                for t in ("70", "80", "90"):
                    dc = a["draw_crossings"][f"{field}_{k}_full_grid_{t}"]
                    ref = r2a["indices"][k][f"crossings_{field}"][t]
                    if ref["point"] == "inf" and dc["share_reaching_within_200"] == 0:
                        continue
                    if field == "temporal" and t == "90":
                        continue
                    sp = a["split_error"][k]["crossings"][t] if field == "random" else None
                    m, centre = drawinfo(dc)
                    lo, hi = val(ref["ci_all_draws"][0]), val(ref["ci_all_draws"][1])
                    inside = bool(lo <= val(centre) <= hi)
                    where = tab if k not in ("bias_g", "brier_g") else "Table S1f"
                    if k in ("auroc_pr", "auroc", "gamma"):
                        where += "; Table S2b"
                    if field == "random" and k in ("bias", "brier"):
                        where += "; Table S1b; Table S1f"
                    if field == "temporal" and k == "bias":
                        where += "; bias sensitivity table"
                    printed_tab = ((field == "random" and k in ("bias", "brier", "auroc_pr", "auroc"))
                                   or (field == "temporal" and k in ("bias", "brier") and t in ("70", "80"))
                                   or (field == "random" and k in ("bias_g", "brier_g") and t in ("70", "80")))
                    if not printed_tab:
                        where = "not printed"
                    row = {"sample": lab, "where": where, "design": f"A, {field} halves", "index": k, "level": f".{t}",
                           "printed": ref["point"], "printed_ci": ref["ci_all_draws"],
                           "split_se_delta": sp["split_se_delta"] if sp else None,
                           "split_se_resample": sp["resample_sd"] if sp else None,
                           "draw_mean": m, "draw_sd": dc["sd"] if dc["sd"] is not None else dc.get("sd_of_draws_reaching"),
                           "draw_median": dc["median"], "draw_p2_5": dc["p2_5"], "draw_p97_5": dc["p97_5"],
                           "draw_share_reaching_within_200": dc["share_reaching_within_200"],
                           "crossing_of_mean_curve": dc["crossing_of_mean_curve"],
                           "share_draws_below_printed": dc["share_draws_below_draw0"],
                           "draw_centre_inside_printed_ci": inside,
                           "printed_inside_ci_recentred_on_draw_centre": recentred(ref["point"], ref["ci_all_draws"], centre),
                           "in_printed_table": bool(printed_tab)}
                    rows.append(row)
                    sb = dc["share_draws_below_draw0"]
                    if printed_tab and not inside:
                        flags_ci.append(row)
                    if printed_tab and sb is not None and (sb >= 0.975 or sb <= 0.025):
                        tails.append(row)
        # pool sensitivity (draw 0 split error only)
        for k in ("bias", "brier"):
            for t in ("70", "80"):
                sp = a["split_error"][f"pool_sensitivity_{k}"]["crossings"][t]
                ref = r2a["pool_sensitivity"][k]["crossings_random"][t]
                rows.append({"sample": lab, "where": "Table S1f", "design": "A, random halves, without near-perfect and constant raters",
                             "index": k, "level": f".{t}", "printed": ref["point"], "printed_ci": ref["ci_all_draws"],
                             "split_se_delta": sp["split_se_delta"], "split_se_resample": sp["resample_sd"],
                             "draw_mean": None, "draw_sd": None, "note": "subsample draws not rerun for this variant"})
        # ICC curve (no split error; subsample draws)
        for k in ("bias", "brier"):
            for t in ("70", "80"):
                dc = a["draw_crossings"][f"icc_{k}_full_grid_{t}"]
                cis = VC[lab]["summary"][k][f"icc_minN{t}_ci"]
                m, centre = drawinfo(dc)
                inside = bool(val(cis[0]) <= val(centre) <= val(cis[1] if cis[1] is not None else "inf"))
                row = {"sample": lab, "where": "Table S1b", "design": "A, variance-components ICC(1,N) curve",
                       "index": k, "level": f".{t}", "printed": dc["draw0"], "printed_ci": cis,
                       "split_se_delta": 0.0, "split_se_resample": 0.0, "draw_mean": m, "draw_sd": dc["sd"],
                       "draw_median": dc["median"], "draw_p2_5": dc["p2_5"], "draw_p97_5": dc["p97_5"],
                       "draw_share_reaching_within_200": dc["share_reaching_within_200"],
                       "crossing_of_mean_curve": dc["crossing_of_mean_curve"],
                       "share_draws_below_printed": dc["share_draws_below_draw0"],
                       "draw_centre_inside_printed_ci": inside,
                       "printed_inside_ci_recentred_on_draw_centre": recentred(dc["draw0"], cis, centre),
                       "in_printed_table": True}
                rows.append(row)
                if not inside:
                    flags_ci.append(row)
                sb = dc["share_draws_below_draw0"]
                if sb is not None and (sb >= 0.975 or sb <= 0.025):
                    tails.append(row)
        # split-free estimator: no Monte Carlo error
        for k in ("bias", "brier"):
            am = VC[lab]["all_record_moment_estimator"][k]
            for t in ("70", "80"):
                rows.append({"sample": lab, "where": "Table S1b", "design": "split-free estimator on all rated answers",
                             "index": k, "level": f".{t}", "printed": am[f"implied_n{t}"],
                             "printed_ci": am.get(f"implied_n{t}_ci"), "split_se_delta": 0.0, "split_se_resample": 0.0,
                             "draw_mean": None, "draw_sd": 0.0, "note": "no Monte Carlo error (no subsample, no split)"})
        # Design B
        b = Bd[lab]
        for k in ("bias", "brier", "gamma", "auroc_pr", "auroc"):
            for t in ("70", "80", "90"):
                sp = b["split_error"][k]["crossings"][t]
                dc = b["draw_crossings"][f"{k}_{t}"]
                if sp["point"] == "inf" and dc["share_reaching_within_200"] == 0:
                    continue
                m, centre = drawinfo(dc)
                rows.append({"sample": lab, "where": tab + ("; Table S2b" if k in ("auroc_pr", "auroc") else ""),
                             "design": "B, random halves (point estimate, no bootstrap)", "index": k, "level": f".{t}",
                             "printed": sp["point"], "printed_ci": None, "split_se_delta": sp["split_se_delta"],
                             "split_se_resample": sp["resample_sd"], "draw_mean": m,
                             "draw_sd": dc["sd"] if dc["sd"] is not None else dc.get("sd_of_draws_reaching"),
                             "draw_median": dc["median"], "draw_p2_5": dc["p2_5"], "draw_p97_5": dc["p97_5"],
                             "draw_share_reaching_within_200": dc["share_reaching_within_200"],
                             "crossing_of_mean_curve": dc["crossing_of_mean_curve"], "K": KB})
        # structured halves and weighting, with the target-subset note (split error only)
        for vname, vv in T6[lab]["variants"].items():
            for k in ("bias", "brier"):
                for t in ("70", "80"):
                    c = vv[k]["crossings"][t]
                    rows.append({"sample": lab, "where": "bias sensitivity table; Table S3b" if k == "bias" else "Results text (Brier, harder versus easier)",
                                 "design": f"A, {VNAMES[vname]}", "index": k, "level": f".{t}", "printed": c["point"],
                                 "split_se_delta": c["split_se_delta"], "split_se_resample": c["resample_sd"],
                                 "split_resample_share_reaching_within_200": c["resample_share_reaching_within_200"]})
        for vname, vv in T6[lab]["target_subset"].items():
            for t in ("70",):
                c = vv["bias"]["crossings"][t]
                rows.append({"sample": lab, "where": "bias sensitivity table note (students with a weighting target)",
                             "design": f"A, {VNAMES[vname]}, students with a weighting target", "index": "bias",
                             "level": f".{t}", "printed": c["point"], "split_se_delta": c["split_se_delta"],
                             "split_se_resample": c["resample_sd"]})

    # S2b quantities that are not crossings
    s2b = {}
    for lab in ("S1", "S2"):
        se = A[lab]["split_error"]
        s2b[lab] = {k: {"r15_interpolated": se[k]["r15_interpolated"],
                        "spearman_brown_projection_070_from_N20": se[k]["spearman_brown_projection_070_from_N20"]}
                    for k in ("bias", "brier", "gamma", "auroc_pr", "auroc")}

    # reliability values over draws
    rel = {}
    for lab in ("S1", "S2"):
        dcur = A[lab]["draw_curves"]
        rel[lab] = {name: {N: {"draw0": blk["draw0"][N], "mean_over_draws": blk["mean_over_draws"][N],
                               "sd_over_draws": blk["sd_over_draws"][N]} for N in blk["draw0"]}
                    for name, blk in dcur.items()}
        rel[lab]["designB"] = {k: {N: {"draw0": blk["draw0"][N], "mean_over_draws": blk["mean_over_draws"][N],
                                       "sd_over_draws": blk["sd_over_draws"][N]} for N in blk["draw0"]}
                               for k, blk in Bd[lab]["draw_curves"].items()}

    def maxsd(lab, name, Nmin=5):
        blk = A[lab]["draw_curves"][name]["sd_over_draws"]
        vals = [(float(v), int(N)) for N, v in blk.items() if v is not None and int(N) >= Nmin]
        return max(vals)

    summary = {}
    for lab in ("S1", "S2"):
        s = {}
        for k in ("bias", "brier", "gamma", "auroc_pr", "auroc", "conf", "acc"):
            v, n = maxsd(lab, f"random_{k}")
            v20, n20 = maxsd(lab, f"random_{k}", 20)
            sp = A[lab]["split_error"][k]
            s[k] = {"split_se_max_all_N": sp["max_split_se_all_N"], "split_se_N_at_max": sp["argmax_N"],
                    "split_se_max_N_ge_20": sp["max_split_se_N_ge_20"],
                    "draw_sd_random_max_all_N": v, "draw_sd_random_N_at_max": n, "draw_sd_random_max_N_ge_20": v20,
                    "draw_sd_random_N20": A[lab]["draw_curves"][f"random_{k}"]["sd_over_draws"]["20"],
                    "draw_sd_temporal_N20": A[lab]["draw_curves"][f"temporal_{k}"]["sd_over_draws"]["20"]}
        summary[lab] = s
    summary["max_split_se_crossing_bias_brier_designA"] = max(
        r["split_se_delta"] for r in rows if r.get("design") == "A, random halves" and r["index"] in ("bias", "brier")
        and r["split_se_delta"] is not None)
    summary["max_split_se_crossing_resolution_designA"] = max(
        r["split_se_delta"] for r in rows if r.get("design") == "A, random halves"
        and r["index"] in ("auroc_pr", "auroc") and r["split_se_delta"] is not None)
    summary["max_split_se_crossing_structured_halves"] = max(
        r["split_se_delta"] for r in rows if str(r.get("where", "")).startswith("bias sensitivity table")
        and r["split_se_delta"] is not None)
    summary["criteria"] = {"S2": A["S2"]["draw_quantities"], "S1": A["S1"]["draw_quantities"]}

    out = {"_definitions": {
        "split_se": "Monte Carlo SE of a printed reliability from the finite number of random splits: SD of the 50 per-split Spearman-Brown values divided by sqrt(50) (draw 0 = the submitted subsample). For a crossing: delta method from the split SEs of the two bracketing grid points (split_se_delta), checked by resampling the 50 per-split values within each N (split_se_resample, B = 2,000)",
        "draw": "one subsample of N rated answers per student at every grid point, with its own 50 random splits and its temporal split; draw 0 = the submitted seeds [42, 5, N, sid] (Design A) or [42, 6, N, sid] (Design B), further draws [42, 501, N, sid, k] (Design A, K = 100) and [42, 502, N, sid, k] (Design B, K = 50)",
        "draw_sd": "SD over the K further draws of a single-draw value: total Monte Carlo SD (subsample plus split component) of the printed value",
        "draw_mean, draw_median, draw_p2_5, draw_p97_5": "over the K further draws (draw 0 excluded); a draw whose curve never reaches the level within 200 is ranked above 200; the mean and SD are given only when every draw reaches the level",
        "crossing_of_mean_curve": "crossing of the reliability curve averaged over the K further draws",
        "share_draws_below_printed": "share of the K further draws whose crossing is smaller than the printed (draw 0) value; 1.0 = printed value above every further draw",
        "draw_centre_inside_printed_ci": "the draw mean (draw median when some draws do not reach the level) lies inside the printed bootstrap 95% CI, which holds the subsample and splits fixed",
        "printed_inside_ci_recentred_on_draw_centre": "the printed value lies inside the printed CI shifted so that the draw centre replaces the printed point (equivalently, the draw centre lies between 2 x printed - upper bound and 2 x printed - lower bound)",
        "no_split_component": "temporal halves, the ICC curve and the split-free estimator have no split error; the split-free estimator has no Monte Carlo error at all"},
        "K_further_draws_designA": K, "K_further_draws_designB": KB, "R_splits": 50,
        "regression_draw0": {lab: {"designA_max_abs_diff": max(A[lab]["regression_max_abs_diff"].values()),
                                   "designA_split_se_vs_r2_A_max_abs_diff": max(v for kk, v in A[lab]["regression_max_abs_diff"].items() if kk.startswith("split_se_")),
                                   "designB_max_abs_diff": max(Bd[lab]["regression_max_abs_diff"].values()),
                                   "structured_halves_max_abs_diff": max(T6[lab]["regression_max_abs_diff"].values())}
                             for lab in ("S1", "S2")},
        "flag_draw_mean_outside_printed_ci": flags_ci,
        "printed_value_in_outer_2_5_percent_of_draws": tails,
        "printed_outside_ci_recentred_on_draw_centre": [r for r in rows if r.get("in_printed_table") and r.get("printed_inside_ci_recentred_on_draw_centre") is False],
        "summary": summary, "split_error": split, "crossings": rows, "resolution_quantities_split_se": s2b,
        "reliability_over_draws": rel}
    with open(os.path.join(RES, "r2_E_mc.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    # ---------------- r2_E_mc.md ----------------
    L = []
    if flags_ci:
        L.append("FLAG: for " + ", ".join(f"{r['sample']} {NAMES.get(r['index'], r['index'])} {r['level']} ({r['design']})"
                                         for r in flags_ci)
                 + f", the mean over {K} further subsample draws lies outside the printed 95% CI; report the draw mean beside the printed value.")
    else:
        L.append(f"No flag: for every crossing with a printed 95% CI, the mean (or, where some draws never reach the level, the median) over {K} further subsample draws lies inside that CI.")
    L.append("")
    if tails:
        L.append("Printed values in the outer 2.5% of the further draws (the printed CI still covers the draw centre):")
        for r in tails:
            L.append(f"- {r['sample']} {NAMES.get(r['index'], r['index'])} {r['level']}, {r['design']}: printed {c1(r['printed'])} "
                     f"{ci(r['printed_ci'])}; {r['share_draws_below_printed']:.2f} of draws lie below it; draw mean {c1(r['draw_mean'])} "
                     f"(SD {c2(r['draw_sd'])}), 95% of draws {c1(r['draw_p2_5'])} to {c1(r['draw_p97_5'])}.")
        L.append("")
    rc_fail = [r for r in rows if r.get("in_printed_table") and r.get("printed_inside_ci_recentred_on_draw_centre") is False]
    L.append(("Second check (printed value inside its CI shifted to the draw centre): fails for " + "; ".join(
        f"{r['sample']} {NAMES.get(r['index'], r['index'])} {r['level']} ({r['design']})" for r in rc_fail) + ".")
             if rc_fail else "Second check (printed value inside its CI shifted to the draw centre): passes for every crossing with a printed CI.")
    L.append("")
    L.append("# Monte Carlo precision of the reliability curves and minimum N")
    L.append("")
    L.append("Two Monte Carlo components sit in every printed Design A or Design B value: the split component (each value averages 50 random splits of one subsample) and the subsample component (one seeded subsample of N rated answers per student). "
             f"The split component is measured from the 50 per-split values of the printed draw. The total of both is measured by rerunning the frozen estimator on {K} further subsample draws (Design A; {KB} for Design B) and taking the SD over draws. "
             "Printed values come from the pre-specified single draw and stay as printed; the bootstrap CIs hold the subsample and the splits fixed and cover student sampling only. "
             "Temporal halves, the variance-components (ICC) curve and the split-free estimator have no split component; the split-free estimator has no Monte Carlo error at all.")
    L.append("")
    L.append("## 1. Maximum Monte Carlo error of a reliability value, Design A random halves")
    L.append("")
    L.append("| Sample | Index | Split SE, max over all N (N) | Split SE, max over N >= 20 | Draw SD, max over all N (N) | Draw SD, max over N >= 20 | Draw SD at N = 20 (random / temporal) |")
    L.append("|---|---|---|---|---|---|---|")
    for lab in ("S1", "S2"):
        for k in ("bias", "brier", "gamma", "auroc_pr", "auroc", "conf", "acc"):
            s = summary[lab][k]
            L.append(f"| {lab} | {NAMES[k]} | {r4(s['split_se_max_all_N'])} ({s['split_se_N_at_max']}) | {r4(s['split_se_max_N_ge_20'])} | "
                     f"{r4(s['draw_sd_random_max_all_N'])} ({s['draw_sd_random_N_at_max']}) | {r4(s['draw_sd_random_max_N_ge_20'])} | "
                     f"{r4(s['draw_sd_random_N20'])} / {r4(s['draw_sd_temporal_N20'])} |")
    L.append("")
    L.append(f"Largest split SE of any printed Design A random-half crossing: bias and Brier {c2(summary['max_split_se_crossing_bias_brier_designA'])} answers; AUROC (either rule) {c2(summary['max_split_se_crossing_resolution_designA'])} answers; structured-halves and weighting table {c2(summary['max_split_se_crossing_structured_halves'])} answers.")
    L.append("")
    L.append("## 2. Per-N Monte Carlo error, Design A")
    L.append("")
    L.append("Cells: split SE / total draw SD. Temporal halves: no split SE; their draw SD is in r2_E_mc.json (reliability_over_draws).")
    L.append("")
    L.append("| Sample, N | Bias | Brier | Gamma | AUROC, plan rule | AUROC, tie convention |")
    L.append("|---|---|---|---|---|---|")
    for lab in ("S1", "S2"):
        for N in map(str, EC.FULL_GRID):
            cells = []
            for k in ("bias", "brier", "gamma", "auroc_pr", "auroc"):
                cells.append(f"{r4(split[lab]['designA'][k]['per_N'][N])} / {r4(rel[lab]['random_' + k][N]['sd_over_draws'])}")
            L.append(f"| {lab}, {N} | " + " | ".join(cells) + " |")
    L.append("")
    L.append("## 3. Minimum N: printed value, split SE and spread over further subsample draws")
    L.append("")
    L.append("| Sample | Where | Design | Index | Level | Printed [95% CI] | Split SE (delta / resample) | Draw mean (SD) | Draw median; 95% of draws | Reaching within 200 | Crossing of mean curve | Draws below printed | Draw centre in printed CI | Printed in CI recentred on draw centre |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        if "draw_mean" not in r and "printed_ci" not in r:
            continue
        if r["design"].startswith("A, Halves") or "weighting target" in r["design"]:
            continue
        pc = ci(r["printed_ci"]) if r.get("printed_ci") else ""
        dm = f"{c1(r.get('draw_mean'))} ({c2(r.get('draw_sd'))})" if r.get("draw_median") is not None else (r.get("note") or "")
        med = f"{c1(r.get('draw_median'))}; {c1(r.get('draw_p2_5'))} to {c1(r.get('draw_p97_5'))}" if r.get("draw_median") is not None else ""
        reach = f"{r['draw_share_reaching_within_200']:.2f}" if r.get("draw_share_reaching_within_200") is not None else ""
        below = f"{r['share_draws_below_printed']:.2f}" if r.get("share_draws_below_printed") is not None else ""
        inside = ("yes" if r["draw_centre_inside_printed_ci"] else "NO") if "draw_centre_inside_printed_ci" in r else ""
        rc = r.get("printed_inside_ci_recentred_on_draw_centre")
        rc = "" if rc is None else ("yes" if rc else "NO")
        L.append(f"| {r['sample']} | {r['where']} | {r['design']} | {NAMES.get(r['index'], r['index'])} | {r['level']} | {c1(r['printed'])} {pc} | "
                 f"{c2(r.get('split_se_delta'))} / {c2(r.get('split_se_resample'))} | {dm} | {med} | {reach} | {c1(r.get('crossing_of_mean_curve'))} | {below} | {inside} | {rc} |")
    L.append("")
    L.append("## 4. Structured halves and weighting (bias table and its note): split SE of the crossings")
    L.append("")
    L.append("These variants were not rerun on further subsample draws; they share the printed subsamples of the random-half curve.")
    L.append("")
    L.append("| Sample | Variant | Index | .70: printed, split SE | .80: printed, split SE | r(20) split SE |")
    L.append("|---|---|---|---|---|---|")
    for lab in ("S1", "S2"):
        for vname, vv in T6[lab]["variants"].items():
            for k in ("bias", "brier"):
                c70, c80 = vv[k]["crossings"]["70"], vv[k]["crossings"]["80"]
                L.append(f"| {lab} | {VNAMES[vname]} | {NAMES[k]} | {c1(c70['point'])}, {c2(c70['split_se_delta'])} | "
                         f"{c1(c80['point'])}, {c2(c80['split_se_delta'])} | {r4(vv[k]['r20_split_se'])} |")
        for vname, vv in T6[lab]["target_subset"].items():
            c70 = vv["bias"]["crossings"]["70"]
            L.append(f"| {lab} | {VNAMES[vname]}, students with a weighting target (n = {T6[lab]['n_target_subset']}) | Bias | {c1(c70['point'])}, {c2(c70['split_se_delta'])} | | |")
    L.append("")
    L.append("## 5. Confirmatory and replication quantities over the further draws (informational; the confirmatory record stays as computed)")
    L.append("")
    q1 = A["S1"]["draw_quantities"]
    q2 = A["S2"]["draw_quantities"]
    L.append(f"- S1 H2 quantity, bias r(15) interpolated between N = 10 and 20: printed {r3(q1['bias_r15_interpolated']['draw0'])}; over {K} further draws mean {r3(q1['bias_r15_interpolated']['mean'])} (SD {r3(q1['bias_r15_interpolated']['sd'])}, range {r3(q1['bias_r15_interpolated']['min'])} to {r3(q1['bias_r15_interpolated']['max'])}); at or above .70 in {q1['bias_r15_interpolated']['share_ge_070']:.2f} of draws.")
    L.append(f"- S1 H1 quantity, bias r(20): printed {r3(q1['bias_r20']['draw0'])}; draw mean {r3(q1['bias_r20']['mean'])} (SD {r3(q1['bias_r20']['sd'])}); below .70 in {q1['bias_r20']['share_below_070']:.2f} of draws.")
    cr = q2["replication_criteria_S2"]
    L.append(f"- S2 criterion (a), bias r(20) below .70: printed {r3(cr['criterion_a_bias_r20_below_070']['draw0'])} (not met); met in {cr['criterion_a_bias_r20_below_070']['share_further_draws_below_070']:.2f} of further draws.")
    for g in ("full_grid", "frozen_grid"):
        gl = "full 11-point grid (printed)" if g == "full_grid" else "pre-specified 6-point grid"
        b = cr[f"{g}_brier70"]
        L.append(f"- S2 criterion (b), {gl}: Brier .70 crossing {c1(b['draw0'])} in the printed draw; over {K} further draws mean {c1(b['mean'])}, median {c1(b['median'])}, above 31 in {b['share_further_draws_above_31']:.2f}; criterion met (three of four indices above 31) in {cr[g + '_plan_rule']['share_further_draws_met']:.2f} of further draws under the plan rule and {cr[g + '_tie_convention']['share_further_draws_met']:.2f} under the tie convention (printed draw: met under both).")
    L.append("")
    L.append("## 6. Files")
    L.append("")
    L.append("Per-sample detail, including per-split values, per-draw curves and crossings: r2_E_parts/designA_<sample>.json, designB_<sample>.json, structured_<sample>.json (per-split and per-draw aggregates; no per-student values). Scripts: scripts/r2_E_common.py, r2_E_mc_designA.py, r2_E_mc_designB.py, r2_E_mc_structured.py, r2_H_report.py.")
    with open(os.path.join(RES, "r2_E_mc.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")

    # ---------------- r2_F_grid.md ----------------
    F = load(os.path.join(RES, "r2_F_grid.json"))
    M = []
    g1, g2 = F["grid_S1"], F["grid_S2"]
    s1b = g1["comparisons"]
    s2c = g2["comparisons"]
    M.append(f"Pre-specified grid: S1 bias .70 crossing {c1(s1b['bias70_vs_S1_median_15_frozen_grid']['minN70'])} on the 6-point grid against {c1(s1b['bias70_vs_S1_median_15_full_grid']['minN70'])} on the full grid (S1 median volume 15: above 15 on the 6-point grid {s1b['bias70_vs_S1_median_15_frozen_grid']['above_15']}, on the full grid {s1b['bias70_vs_S1_median_15_full_grid']['above_15']}). "
             f"S2 Brier .70 crossing {c1(s2c['criterion_b_frozen_grid_plan_rule']['minN70']['brier'])} on the 6-point grid against {c1(s2c['criterion_b_full_grid_plan_rule']['minN70']['brier'])} on the full grid (S2 median 31); criterion (b) met on both grids in the printed draw.")
    M.append("")
    M.append("# Crossings on the pre-specified grid, alpha-type checks, AUROC diagnostics")
    M.append("")
    M.append("## 1. Design A crossings on the full grid and on the pre-specified grid")
    M.append("")
    M.append("The analysis plan fixed N = 5, 10, 20, 50, 100 and 200; addendum 1 (2026-07-09, before submission) added 12, 15, 25, 30 and 40 to Design A as an interpolation-sensitivity check. The printed crossings use all eleven points. Below, the same curves and the same crossing rule on the six pre-specified points; CIs from the same 1,000 joint bootstrap draws. Design B already uses the six pre-specified points. "
             f"Regression: regenerated bootstrap draws equal the cached draws and the r2_A intervals (pass: S1 {g1['regression_pass']}, S2 {g2['regression_pass']}); full-grid crossings and CIs reproduced (S1 {g1['full_grid_reproduces_r2_A']}, S2 {g2['full_grid_reproduces_r2_A']}).")
    M.append("")
    M.append("| Sample | Halves | Index | Level | Full grid (printed) [95% CI] | Pre-specified grid [95% CI] | Difference |")
    M.append("|---|---|---|---|---|---|---|")
    for lab, g in (("S1", g1), ("S2", g2)):
        for sp in ("random", "temporal"):
            for k in ("bias", "brier", "gamma", "auroc_pr", "auroc", "bias_g", "brier_g"):
                for t in ("70", "80", "90"):
                    c = g["crossings"][f"{sp}_{k}_{t}"]
                    if c["full_grid"]["point"] == "inf" and c["frozen_grid"]["point"] == "inf" and \
                            c["full_grid"]["ci_all_draws"][0] == "inf" and c["frozen_grid"]["ci_all_draws"][0] == "inf":
                        continue
                    dfz = c["frozen_minus_full"]
                    M.append(f"| {lab} | {sp} | {NAMES[k]} | .{t} | {c1(c['full_grid']['point'])} {ci(c['full_grid']['ci_all_draws'])} | "
                             f"{c1(c['frozen_grid']['point'])} {ci(c['frozen_grid']['ci_all_draws'])} | {'n/a' if dfz is None else f'{dfz:+.1f}'} |")
    M.append("")
    M.append("Rows where neither grid reaches the level and every draw stays above 200 are omitted (gamma at every level; the resolution indices at .80 and .90).")
    M.append("")
    M.append("S2 replication criterion (b) on both grids (printed draw):")
    M.append("")
    for key, v in s2c.items():
        M.append(f"- {key.replace('_', ' ')}: " + ", ".join(f"{NAMES.get(k, k)} {c1(x)}" for k, x in v["minN70"].items()) + f"; met {v['met']}")
    M.append("")
    al = F["alpha_subset_S2"]
    M.append("## 2. Crossed S2 subset: alpha and Phi against ICC(1,k)")
    M.append("")
    M.append(f"{al['n_quizzes']} quizzes, {al['k_range'][0]} to {al['k_range'][1]} items, {al['n_quiz_student_rows']:,} quiz-student rows. ICC(1,k), alpha and Phi recomputed from the data and matched to the stored values (max abs diff {max(al['checks_max_abs_diff_vs_r1_B_parts'][m] for m in ('icc1_k', 'alpha', 'phi')):.1e}; pass {al['checks_pass']}). Differences against ICC(1,k) are deterministic.")
    M.append("")
    M.append("| Index | Alpha minus ICC(1,k): median (range) | Phi minus ICC(1,k): median (range) | Printed alpha minus split-half median | Over 200 split reruns: mean (SD); 95% range | Reruns rounding to printed | Printed Phi minus split-half median | Over 200 reruns: mean (SD); 95% range | Reruns rounding to printed |")
    M.append("|---|---|---|---|---|---|---|---|---|")
    for key in ("bias", "brier"):
        b = al[key]
        am, pm = b["alpha_minus_sb_median_over_reruns"], b["phi_minus_sb_median_over_reruns"]
        M.append(f"| {NAMES[key]} | {r3(b['alpha_minus_icc1_k']['median'])} ({r3(b['alpha_minus_icc1_k']['min'])} to {r3(b['alpha_minus_icc1_k']['max'])}) | "
                 f"{r3(b['phi_minus_icc1_k']['median'])} ({r3(b['phi_minus_icc1_k']['min'])} to {r3(b['phi_minus_icc1_k']['max'])}) | "
                 f"{am['printed_value']:.3f} | {r3(am['mean'])} ({r3(am['sd'])}); {r3(am['p2_5'])} to {r3(am['p97_5'])} | {am['share_reruns_rounding_to_printed']:.2f} | "
                 f"{pm['printed_value']:.3f} | {r3(pm['mean'])} ({r3(pm['sd'])}); {r3(pm['p2_5'])} to {r3(pm['p97_5'])} | {pm['share_reruns_rounding_to_printed']:.2f} |")
    M.append("")
    for key in ("bias", "brier"):
        b = al[key]
        M.append(f"- {NAMES[key]}: alpha exceeded ICC(1,k) in {b['n_quizzes_alpha_gt_icc1_k']} of 14 quizzes, Phi in {b['n_quizzes_phi_gt_icc1_k']} of 14; medians alpha {r3(b['alpha_median'])}, Phi {r3(b['phi_median'])}, ICC(1,k) {r3(b['icc1_k_median'])}; split-half minus ICC(1,k), median over quizzes, mean over reruns {r3(b['sb_minus_icc1_k_median_over_reruns']['mean'])} (SD {r3(b['sb_minus_icc1_k_median_over_reruns']['sd'])}).")
    M.append("")
    M.append("## 3. AUROC minus gamma reliability (Design A, random halves, N <= 50)")
    M.append("")
    ag = F["auroc_gamma"]
    for lab in ("S1", "S2"):
        for rule in ("plan_rule", "tie_convention"):
            d = ag[lab][f"auroc_minus_gamma_{rule}"]
            od = d["over_further_draws"]
            M.append(f"- {lab}, {rule.replace('_', ' ')}: printed draw {d['min']:+.3f} (N = {d['N_at_min']}) to {d['max']:+.3f} (N = {d['N_at_max']}); "
                     f"mean over {od['K']} further draws {od['min_of_means']:+.3f} to {od['max_of_means']:+.3f} "
                     f"(SD of the difference per draw {min(od['sd_per_N'].values()):.3f} to {max(od['sd_per_N'].values()):.3f})")
    M.append("")
    M.append("## 4. Students entering the gamma and AUROC curves (Design A, random halves)")
    M.append("")
    M.append("A student enters the correlation at N only when both halves give a defined value. Shares are means over the 50 splits of the printed draw; one-half undefined = share with half A undefined (the definition behind the halves row of Table S2a).")
    M.append("")
    M.append("| Sample, N | Enter gamma and plan-rule AUROC curves | Enter tie-convention AUROC curve | One half undefined (gamma) | Enter temporal gamma curve |")
    M.append("|---|---|---|---|---|")
    for lab in ("S1", "S2"):
        ks = ag[lab]["students_entering_resolution_curves"]
        for N in map(str, EC.FULL_GRID):
            M.append(f"| {lab}, {N} | {100 * ks['gamma_and_auroc_plan_rule_draw0'][N]:.1f}% | {100 * ks['auroc_tie_convention_draw0'][N]:.1f}% | "
                     f"{100 * ks['gamma_one_half_undefined_draw0'][N]:.1f}% | {100 * ks['gamma_temporal_both_halves_defined_draw0'][N]:.1f}% |")
    M.append("")
    M.append("## 5. The AUROC curve between N = 15 and 20")
    M.append("")
    for lab in ("S1", "S2"):
        dg = ag[lab]["auroc_dip_diagnostic"]
        a15, a20 = dg["auroc_pr"]["15"], dg["auroc_pr"]["20"]
        ks = ag[lab]["students_entering_resolution_curves"]["gamma_and_auroc_plan_rule_draw0"]
        M.append(f"- {lab}: printed plan-rule AUROC reliability {r3(a15['draw0'])} at N = 15 and {r3(a20['draw0'])} at N = 20 (difference {dg['auroc_pr_draw0_r15_minus_r20']:+.3f}); "
                 f"over {K} further draws the means are {r3(a15['mean_over_draws'])} (SD {r3(a15['sd_over_draws'])}) and {r3(a20['mean_over_draws'])} (SD {r3(a20['sd_over_draws'])}), "
                 f"difference {dg['auroc_pr_mean_r15_minus_r20']:+.3f} (SD of the difference {r3(dg['auroc_pr_sd_of_r15_minus_r20_over_draws'])}); r(15) > r(20) in {dg['auroc_pr_share_draws_r15_gt_r20']:.2f} of draws; "
                 f"mean curve increasing at every step up to N = 50: {dg['auroc_pr_mean_curve_monotone_N_le_50']}. Students entering the curve: {100 * ks['15']:.1f}% at N = 15, {100 * ks['20']:.1f}% at N = 20. Halves at N = 15 hold 7 and 8 answers.")
    M.append("")
    M.append("Scripts: scripts/r2_F_grid.py (computation), scripts/r2_H_report.py (this file).")
    with open(os.path.join(RES, "r2_F_grid.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(M) + "\n")
    print("written r2_E_mc.json, r2_E_mc.md, r2_F_grid.md")


if __name__ == "__main__":
    main()
