"""Merge the per-sample stream D outputs into results/r1_D.json, add the
definitions, the written decision rule and the decision table, and write the
tables part of results/r1_D.md (between the TABLES markers).

Usage: python -B r1_D_report.py
"""
import json
import os

import pandas as pd

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
SAMPLES = {"t12": "S1", "t34": "S2"}
REQ_N = [10, 20, 30, 40, 60]
SHOW_N = [10, 15, 20, 30, 40, 50, 60, 100]

DEFS = {
    "pool": "Design A pool: students with >= 200 confidence-tagged responses (S1 n = 1,389; S2 n = 554). Same students at every N.",
    "N": "Number of confidence-tagged responses in ONE score. Split-half values: score on N responses, estimated from two N/2 halves with Spearman-Brown. Disjoint-set values: each of the two sets has N responses, no step-up.",
    "bias": "mean(c) - mean(x), c = Confidence/100 in {0, .25, .5, .75, 1}, x = correctness. Positive = overconfident.",
    "conf": "mean confidence of the N responses (0-1 scale).",
    "acc": "proportion correct of the N responses.",
    "splithalf_SB": "Paper machinery (reliability.subsample_matrix + half_indices + sb_mean): one N-subsample per student (seed [42, 5, N, sid]), 50 random splits (random) or first vs second half by date (temporal), Pearson r across students, Spearman-Brown 2r/(1+r), mean over splits. Values at N in the original grid equal the baseline curves (regression gate).",
    "disjoint_sets": "rev_round2 design: per repeat, 2N responses drawn without replacement per student (seed [42, 92, N, sid]); set A = first N, set B = last N. 50 repeats. Flags and the decomposition use these N-response scores.",
    "r_parallel_*": "Pearson r between set A and set B scores, mean over repeats: a direct parallel-forms reliability at N (no Spearman-Brown). Cross-check of splithalf_SB.",
    "var_true_*, cov_true_conf_acc": "Cross-set covariances (set A with set B), averaged over repeats. Because the sets share no response, these estimate between-student true-score (co)variances at volume N. Identity: var_true_bias = var_true_conf + var_true_acc - 2 cov_true_conf_acc (residual stored).",
    "true_share_conf / true_share_acc / true_share_cross": "var_true_conf / var_true_bias, var_true_acc / var_true_bias, -2 cov_true_conf_acc / var_true_bias. Sum to 1. Shares above 1 mean the component's stable variance exceeds that of bias and is cancelled by the positive conf-acc covariance.",
    "true_attrib_conf / true_attrib_acc": "Covariance attribution of the stable bias variance: cov_true(bias, conf) / var_true_bias = (var_true_conf - cov_true_conf_acc) / var_true_bias and cov_true(bias, -acc) / var_true_bias = (var_true_acc - cov_true_conf_acc) / var_true_bias. Sum to 1; each is the part of the stable between-student bias variance carried by that component.",
    "obs_*": "Same decompositions on observed (same-set) variances; include measurement error.",
    "rho_true_bias_conf, rho_true_bias_acc, rho_true_conf_acc": "Disattenuated (true-score) correlations from the cross-set covariances. rho_true_bias_conf^2 = share of stable bias variance linearly accounted for by stable confidence alone.",
    "r_obs_*_sameset": "Correlations within one set (part-whole: bias contains conf and acc, so these are not independent evidence).",
    "r_cross_bias_conf / r_cross_bias_acc": "Raw correlation of bias in one N-set with conf (acc) in the other N-set, symmetrized. Not the same quantity as the submitted .351 / -.279, which are 10-response half scores (rev_extras, N = 20 split).",
    "flag rules": {
        "decile": "top 10% of bias in the set (k = round(.10 n)), argsort order for ties; identical to rev_round2.",
        "thr15": "bias > +.15 (strict), computed on integers: 25 sum(code) - 100 sum(x) > 15 N.",
        "thr20": "bias > +.20 (strict), same integer rule with 20 N.",
        "thr15_ge / thr20_ge": "same with >= (tie sensitivity).",
        "lcb_z": "lower bound of the two-sided 95% CI of the student's bias above 0: mean(d) - 1.96 sd(d) / sqrt(N) > 0, d = c - x, sd with ddof 1. sd = 0 sets: bound = mean, flagged iff mean > 0 (counts in set_properties).",
        "lcb_t": "same with the t quantile, df = N - 1 (sensitivity).",
        "ci_three_way_z": "over = lcb_z > 0, under = upper bound < 0, else unclassified; agreement and Cohen's kappa across sets.",
    },
    "flag metrics": {
        "share_flagged": "share of students flagged in a set (mean of the two sets, mean over repeats).",
        "retention": "2 |A and B| / (|A| + |B|) (Dice), mean over repeats. Equals the submitted top-decile retention |A and B| / k when both sets have k flags.",
        "retention_directional": "P(flag in other set | flag in this set), both directions pooled.",
        "chance_retention": "expected retention if the two sets' flags were independent: 2 |A| |B| / (n (|A| + |B|)); about equal to share_flagged.",
        "kappa": "Cohen's kappa of the binary flag across the two sets (agreement beyond chance, both flagged and not flagged).",
        "share_flagged_with_other_set_bias_le_0": "among flags, share whose bias in the other set is <= 0 (not overconfident at all on fresh data).",
    },
    "bootstrap": "Percentile 95% CI, B = 1000 student resamples, the joint draws of reliability.py (seed [42, 77, sid]); the same resample is applied to every statistic and every N, with subsamples, splits and disjoint sets held fixed. Bias split-half CIs at original grid N equal the baseline curve CIs (regression gate).",
}

RULE = {
    "basis": "Reliability anchors already used in the submitted manuscript: .70 minimally acceptable for research use, .80 lower-stakes applied interpretation (Nunnally, 1978), about .90 for decisions with individual consequences (AERA et al., 2014). Applied to the Design A random-split bias reliability at N (point estimate); the labeling rule uses the temporal split because a label assumes stability over time.",
    "tiers": [
        {"tier": "group", "condition": "r_bias_random < .70",
         "guidance": "Individual scores not interpretable; only aggregate (class or cohort) use remains, and the reliability of aggregates was not estimated here. Do not report individual scores."},
        {"tier": "research", "condition": ".70 <= r_bias_random < .80",
         "guidance": "Research use of individual scores in correlational or group-comparison analyses; do not show individual scores or labels to students or teachers."},
        {"tier": "feedback", "condition": "r_bias_random >= .80",
         "guidance": "Low-stakes individual feedback: show the score with its uncertainty band as a provisional description; no 'overconfident' label."},
    ],
    "labeling": "An individual 'overconfident' or 'underconfident' label is treated as a classification with individual consequences and requires temporal-split bias reliability >= .90. The maximum temporal-split bias reliability over the whole Design A grid (N <= 200) is stored per sample; if it is below .90, no N in these data supports individual labeling.",
    "caveat": "Group-level reliability (of class or cohort means) was NOT computed in this stream; the 'group' tier rests on averaging over students, not on a computed coefficient.",
    "known_asymmetries": "Tiers use the random-split point estimate. Temporal-split bias r at N = 30 is .750 (S1) and .747 (S2), below the .80 feedback cut. Bootstrap intervals straddle a cut at S1 N = 15 (.70), S2 N = 15 (.70) and S2 N = 30 (.80). The choice of split and of point vs interval is left to the writer.",
}
DEFS["rho2_true_bias_conf / rho2_true_bias_acc"] = "Squares of rho_true_bias_conf / rho_true_bias_acc: share of stable bias variance linearly accounted for by stable confidence (accuracy) alone. Computed in r1_D_report.py from the stored point values (no CI)."


def load(s):
    with open(os.path.join(REV, "tmp", f"r1_D_{s}.json")) as fh:
        return json.load(fh)


def f3(v):
    if v is None:
        return "unknown"
    s = f"{v:.3f}"
    if s.startswith("0."):
        return s[1:]
    if s.startswith("-0."):
        return "-" + s[2:]
    return s


def f3ci(v, c):
    if c is None or c[0] is None:
        return f3(v)
    return f"{f3(v)} [{f3(c[0])}, {f3(c[1])}]"


def pct(v):
    return "unknown" if v is None else f"{100 * v:.1f}%"


def tier_of(r):
    if r < 0.70:
        return RULE["tiers"][0]
    if r < 0.80:
        return RULE["tiers"][1]
    return RULE["tiers"][2]


def main():
    data = {s: load(s) for s in SAMPLES}
    out = {
        "stream": "D (construct relevance and flag stability), revision round 1, post-review exploratory",
        "date": "2026-11-24",
        "scripts": ["revision-2026-11/scripts/r1_D_construct_flags.py",
                    "revision-2026-11/scripts/r1_D_report.py"],
        "definitions": DEFS,
        "decision_rule": RULE,
        "samples": {},
        "decision_table": {},
        "runtime_seconds": {},
    }
    md_tables = []
    for s, lab in SAMPLES.items():
        d = data[s]
        base = pd.read_csv(os.path.join(REV, "results", "baseline", "results",
                                        f"reliability_curves_{s}.csv"))
        bt = base[(base.design == "A") & (base.split == "temporal") & (base["index"] == "bias")]
        max_temp = float(bt.r.max())
        max_temp_N = int(bt.loc[bt.r.idxmax(), "N"])
        d["max_temporal_bias_reliability_designA"] = {"r": max_temp, "N": max_temp_N}
        for blk in d["per_N"].values():
            dc = blk["decomposition"]
            dc["rho2_true_bias_conf"] = dc["rho_true_bias_conf"] ** 2
            dc["rho2_true_bias_acc"] = dc["rho_true_bias_acc"] ** 2
        out["samples"][f"{lab}_{s}"] = d
        out["runtime_seconds"][f"{lab}_{s}"] = d["runtime_seconds"]
        label_ok = max_temp >= 0.90
        rows = []
        for N in SHOW_N:
            b = d["per_N"][str(N)]
            sh = b["splithalf_SB"]
            fl = b["flags"]
            t = tier_of(sh["bias"]["random"]["r"])
            row = {
                "N": N, "requested_N": N in REQ_N,
                "r_bias_random": sh["bias"]["random"]["r"], "r_bias_random_ci": sh["bias"]["random"]["ci"],
                "r_bias_temporal": sh["bias"]["temporal"]["r"], "r_bias_temporal_ci": sh["bias"]["temporal"]["ci"],
                "r_conf_random": sh["conf"]["random"]["r"], "r_acc_random": sh["acc"]["random"]["r"],
                "decile_retention": fl["decile"]["retention"], "decile_retention_ci": fl["decile"]["retention_ci"],
                "decile_kappa": fl["decile"]["kappa"], "decile_chance_retention": 0.10,
                "thr15_share_flagged": fl["thr15"]["share_flagged"],
                "thr15_retention": fl["thr15"]["retention"], "thr15_retention_ci": fl["thr15"]["retention_ci"],
                "thr15_chance_retention": fl["thr15"]["chance_retention"], "thr15_kappa": fl["thr15"]["kappa"],
                "thr20_share_flagged": fl["thr20"]["share_flagged"],
                "thr20_retention": fl["thr20"]["retention"], "thr20_retention_ci": fl["thr20"]["retention_ci"],
                "thr20_chance_retention": fl["thr20"]["chance_retention"], "thr20_kappa": fl["thr20"]["kappa"],
                "lcb_share_flagged": fl["lcb_z"]["share_flagged"],
                "lcb_retention": fl["lcb_z"]["retention"], "lcb_retention_ci": fl["lcb_z"]["retention_ci"],
                "lcb_chance_retention": fl["lcb_z"]["chance_retention"], "lcb_kappa": fl["lcb_z"]["kappa"],
                "tier": t["tier"], "guidance": t["guidance"],
                "labeling": ("supported" if label_ok else
                             f"not supported (temporal bias r max {max_temp:.3f} at N = {max_temp_N} < .90)"),
                "thr15_lost_on_retest": 1 - fl["thr15"]["retention_directional"],
                "source_r_bias": ("baseline curves (equal to recomputation)" if N != 60
                                  else "new grid point, seed [42, 5, 60, sid]"),
                "source_decile": ("rev_round2 (equal to recomputation)" if N in (10, 15, 20, 30, 50, 100)
                                  else "new grid point, seed [42, 92, N, sid]"),
            }
            rows.append(row)
        out["decision_table"][f"{lab}_{s}"] = rows

        # ---------- markdown tables for this sample
        g = d["regression_gates"]
        md_tables.append(f"### {lab} ({s}), Design A pool n = {d['n_pool']:,} (top decile k = {d['k_top_decile']})\n")
        md_tables.append(
            f"Regression gates: split-half r vs baseline curves max abs diff {g['splithalf_r_vs_baseline_curves_max_abs_diff']:.1e} "
            f"(N {g['splithalf_N_checked']}); bias bootstrap CI vs baseline {g['bias_bootstrap_ci_vs_baseline_curves_max_abs_diff']:.1e}; "
            f"top-decile flags vs rev_round2 {g['decile_flags_vs_rev_round2_max_abs_diff']:.1e} (N {g['decile_N_checked']}).\n")
        md_tables.append("**(1) What bias is made of** (split-half SB reliability; cross-set true-score attribution)\n")
        md_tables.append("| N | r bias | r conf | r acc | stable bias var carried by conf | by acc | rho_true(bias, conf) | rho_true(bias, acc) | rho_true(conf, acc) | r parallel bias (cross-check) |")
        md_tables.append("|---|---|---|---|---|---|---|---|---|---|")
        for N in SHOW_N:
            b = d["per_N"][str(N)]
            sh, dc = b["splithalf_SB"], b["decomposition"]
            md_tables.append(
                f"| {N} | {f3ci(sh['bias']['random']['r'], sh['bias']['random']['ci'])} | "
                f"{f3ci(sh['conf']['random']['r'], sh['conf']['random']['ci'])} | "
                f"{f3ci(sh['acc']['random']['r'], sh['acc']['random']['ci'])} | "
                f"{f3ci(dc['true_attrib_conf'], dc['true_attrib_conf_ci'])} | "
                f"{f3ci(dc['true_attrib_acc'], dc['true_attrib_acc_ci'])} | "
                f"{f3ci(dc['rho_true_bias_conf'], dc['rho_true_bias_conf_ci'])} | "
                f"{f3ci(dc['rho_true_bias_acc'], dc['rho_true_bias_acc_ci'])} | "
                f"{f3ci(dc['rho_true_conf_acc'], dc['rho_true_conf_acc_ci'])} | "
                f"{f3(dc['r_parallel_bias'])} |")
        md_tables.append("")
        md_tables.append("**(2) Flag stability across two disjoint N-sets** (share flagged / retention [95% CI] / chance retention / kappa)\n")
        md_tables.append("| N | top decile | bias > +.15 | bias > +.20 | 95% lower bound > 0 |")
        md_tables.append("|---|---|---|---|---|")
        for N in SHOW_N:
            fl = d["per_N"][str(N)]["flags"]
            cells = []
            for r in ("decile", "thr15", "thr20", "lcb_z"):
                x = fl[r]
                cells.append(f"{pct(x['share_flagged'])} / {f3ci(x['retention'], x['retention_ci'])} / "
                             f"{f3(x['chance_retention'])} / {f3(x['kappa'])}")
            md_tables.append(f"| {N} | " + " | ".join(cells) + " |")
        md_tables.append("")
        md_tables.append("**(3) Decision table** (rule in r1_D.json `decision_rule`)\n")
        md_tables.append("| N | r bias random | r bias temporal | top-decile retention (chance .10) | +.15 retention (chance) | lower-bound retention (chance) | guidance | individual labels |")
        md_tables.append("|---|---|---|---|---|---|---|---|")
        for row in rows:
            star = "" if row["requested_N"] else " (extra)"
            md_tables.append(
                f"| {row['N']}{star} | {f3(row['r_bias_random'])} | {f3(row['r_bias_temporal'])} | "
                f"{f3(row['decile_retention'])} | {f3(row['thr15_retention'])} ({f3(row['thr15_chance_retention'])}) | "
                f"{f3(row['lcb_retention'])} ({f3(row['lcb_chance_retention'])}) | {row['tier']} | "
                f"{'no' if not label_ok else 'yes'} |")
        md_tables.append(f"\nMax temporal-split bias reliability on the Design A grid: {max_temp:.3f} at N = {max_temp_N}.\n")
        rt = d["runtime_seconds"]
        md_tables.append(f"Runtime {lab}: total {rt['total']} s (load {rt['load']}, split-half {rt['splithalf']}, "
                         f"disjoint draws {rt['disjoint_draws']}, bootstrap {rt['bootstrap']}).\n")

    with open(os.path.join(REV, "results", "r1_D.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=1)

    md_path = os.path.join(REV, "results", "r1_D.md")
    tables = "\n".join(md_tables)
    start, end = "<!-- TABLES START -->", "<!-- TABLES END -->"
    if os.path.exists(md_path):
        with open(md_path, encoding="utf-8") as fh:
            txt = fh.read()
        pre, rest = txt.split(start, 1)
        _, post = rest.split(end, 1)
        txt = pre + start + "\n" + tables + "\n" + end + post
    else:
        txt = "# r1_D\n\n" + start + "\n" + tables + "\n" + end + "\n"
    with open(md_path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(txt)
    print("written", md_path)


if __name__ == "__main__":
    main()
