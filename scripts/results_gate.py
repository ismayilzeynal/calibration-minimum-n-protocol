"""Assemble RESULTS_GATE.md: every number destined for the manuscript,
pulled from results/ artifacts (no hand-typed values).
Run after confirmatory.py:  python scripts/results_gate.py
"""
import json

import numpy as np
import pandas as pd

IDX = ["bias", "brier", "gamma", "auroc"]
LAB = {"bias": "Bias", "brier": "Absolute accuracy (Brier)",
       "gamma": "Relative accuracy (gamma)", "auroc": "Relative accuracy (AUROC)"}


def f(x, d=3):
    if x is None or (isinstance(x, float) and (np.isnan(x))):
        return "NA"
    if isinstance(x, float) and np.isinf(x):
        return "> 200"
    return f"{x:.{d}f}" if isinstance(x, float) else str(x)


def fN(x):
    if x is None:
        return "> 200"
    if isinstance(x, float) and np.isinf(x):
        return "> 200"
    return f"{x:.1f}"


def main():
    out = ["# RESULTS GATE - single source of truth for manuscript numbers",
           "", "Generated from results/ artifacts. Recompute, never retype.", ""]
    conf = json.load(open("results/confirmatory_stats.json"))
    out += ["## Confirmatory", "```json",
            json.dumps({"supported": conf["supported"],
                        "holm_family": conf["holm_family"]}, indent=1), "```"]
    h1 = conf["H1"]
    out += [f"- H1 bias r(20) = {f(h1['r20_point'])} "
            f"[{f(h1['r20_ci'][0])}, {f(h1['r20_ci'][1])}], "
            f"one-sided upper {f(h1['r20_upper_1sided95'])}; "
            f"p(r20>=.70) = {f(h1['p_raw'])} -> H1 refuted-direction"]
    h2 = conf["H2"]
    out += ["- H2 r(15) points: " + ", ".join(
        f"{k}={f(v)}" for k, v in h2["r15_points"].items())]
    out += ["- H2 inner Holm p: " + ", ".join(
        f"{k}={f(v)}" for k, v in h2["p_inner_holm"].items())]
    h3 = conf["H3"]
    out += [f"- H3 rho raw {f(h3['rho_raw'])} -> EB {f(h3['rho_eb'])}, "
            f"diff {f(h3['diff_point'],4)} CI {[f(v,4) for v in h3['diff_ci95']]}, "
            f"lower1s {f(h3['diff_lower_1sided95'],4)}, p {f(h3['p_onesided'])}"]
    out += ["- Replication S2: " + json.dumps(conf["replication_S2"],
                                              default=str)]

    for s, tag in (("t12", "S1"), ("t34", "S2")):
        mn = json.load(open(f"results/minN_{s}.json"))

        def cis(m, key):
            ci = m.get(key)
            if ci and ci[0] is not None and ci[1] is not None:
                return f"[{fN(ci[0])}, {fN(ci[1])}]"
            if ci and ci[0] is not None:
                return f"[{fN(ci[0])}, >200]"
            return "-"
        out += ["", f"## Minimum-N decision table inputs ({tag})",
                "| Index | A minN70 | CI | A minN80 | CI | B minN70 | B minN80 | T minN70 | T CI | T minN80 | proph70 | r(20) | r(20) CI | r15interp | r15direct |",
                "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for k in IDX:
            m = mn[k]
            out += [f"| {LAB[k]} | {fN(m['minN70_point'])} | {cis(m,'minN70_ci')} | "
                    f"{fN(m['minN80_point'])} | {cis(m,'minN80_ci')} | "
                    f"{fN(m['minN70_designB_point'])} | {fN(m['minN80_designB_point'])} | "
                    f"{fN(m['minN70_temporal_point'])} | {cis(m,'minN70_temporal_ci')} | "
                    f"{fN(m['minN80_temporal_point'])} | "
                    f"{fN(m['prophecy_minN70_from20']) if m['prophecy_minN70_from20'] is not None else 'undef (r20<=0)'} | "
                    f"{f(m['r20_point'])} | [{f(m['r20_ci'][0])}, {f(m['r20_ci'][1])}] | "
                    f"{f(m['r15_point'])} | {f(m['r15_direct_point'])} |"]
        out += ["", f"### {tag} flag consistency HALF-BASED (addendum 5; released only)"]
        for fl in mn["flag_consistency_designA_bias"]:
            out += [f"- N={fl['N']}: top-decile Jaccard {f(fl['jaccard'])}, "
                    f"share of flags retained {f(fl['pct_flags_retained'])}"]
        r2 = json.load(open(f"results/rev_round2_{s}.json"))
        out += ["", f"### {tag} flag consistency DISJOINT FULL-N (addendum 8; manuscript)",
                f"(chance Jaccard {f(r2['flag_meta']['chance_jaccard'])}, "
                f"chance retention 0.100; pool n={r2['flag_meta']['n_pool']})"]
        for fl in r2["disjoint_flag_consistency"]:
            out += [f"- N={fl['N']}: Jaccard {f(fl['jaccard'])}, "
                    f"retained {f(fl['pct_flags_retained'])}"]
        nf = r2["null_flags"]
        out += [f"- NULL flags (addendum 9): raw-raw {f(nf['null_jaccard_raw_raw'])}, "
                f"EB-vs-raw-criterion {f(nf['null_jaccard_eb_vs_rawcriterion'])} "
                f"(chance .053); n={nf['n_students']}"]
        dsc = r2["descriptives"]
        out += [f"- descriptives (addendum 10): reporters {dsc['reporters']}, "
                f"median {dsc['median_tagged']:.0f} / mean {f(dsc['mean_tagged'],1)} tagged; "
                f"mean conf {f(dsc['mean_confidence'],1)}; acc tagged {f(dsc['accuracy_tagged'],4)}; "
                "shares >= " + ", ".join(f"{k}: {f(v)}" for k, v in
                                         dsc["share_ge"].items())]

        cv = pd.read_csv(f"results/reliability_curves_{s}.csv")
        out += ["", f"### {tag} Design A random curves (r; pct undefined half-based; full-set)"]
        a = cv[(cv.design == "A") & (cv.split == "random")
               & cv["index"].isin(IDX + ["conf", "acc"])]
        piv = a.pivot(index="N", columns="index", values="r").round(3)
        und = a.pivot(index="N", columns="index", values="pct_undefined").round(3)
        undf = a.pivot(index="N", columns="index",
                       values="pct_undefined_full").round(3)
        out += ["```", piv.to_string(), "", "undefined (half-based):",
                und[["gamma", "auroc"]].to_string(), "",
                "undefined (full-set):", undf[["gamma", "auroc"]].to_string(),
                "```"]
        out += [f"### {tag} Design B random curves"]
        b = cv[(cv.design == "B") & (cv.split == "random") & cv["index"].isin(IDX)]
        out += ["```", b.pivot(index="N", columns="index", values="r").round(3)
                .to_string(), "",
                "n_pool: " + str(b[b['index'] == 'bias'].sort_values('N')
                                 .n_pool.tolist()), "```"]
        t = cv[(cv.design == "A") & (cv.split == "temporal") & cv["index"].isin(IDX)]
        out += [f"### {tag} Design A temporal curves", "```",
                t.pivot(index="N", columns="index", values="r").round(3)
                .to_string(), "```"]

        cs = json.load(open(f"results/confound_selection_{s}.json"))
        out += ["", f"## {tag} confound / selection / retest"]
        for N in ("N20", "N100"):
            blk = cs["accuracy_confound"][N]
            out += [f"- {N}: " + "; ".join(
                f"{k}: r_acc={f(blk[k]['pearson_vs_acc'])}, "
                f"rel {f(blk[k]['reliability_raw'])}->"
                f"{f(blk[k]['reliability_acc_partialed'])}" for k in IDX)]
            out += [f"  top-decile Jaccard raw-vs-partialed bias: "
                    f"{f(blk['top_decile_jaccard_raw_vs_partialed_bias'])}"]
        sel = cs["selection"]
        w = sel["within_student"]
        out += [f"- within-student tagged-untagged acc: "
                f"{f(w['mean_acc_tagged_minus_untagged'],4)} "
                f"[{f(w['ci95'][0],4)}, {f(w['ci95'][1],4)}], n={w['n_students']}"]
        sl = sel["student_level"]
        out += [f"- student-level: acc {f(sl['acc_mean_rep'])} vs "
                f"{f(sl['acc_mean_non'])} (d={f(sl['d_accuracy'])}); "
                f"vol median {sl['vol_median_rep']:.0f} vs {sl['vol_median_non']:.0f} "
                f"(d_log={f(sl['logvol_d'])}); reporters {sl['n_reporters']}, "
                f"non {sl['n_nonreporters']}"]
        pm = sel["propensity_model"]
        out += ["- propensity (per-SD log odds): " + "; ".join(
            f"{k}: {f(pm['coef_std'][k])} (z={f(pm['z'][k],1)})"
            for k in ("acc", "log10_vol", "first_month"))]
        out += ["- IPW Design B N20: " + "; ".join(
            f"{k}: {f(v['unweighted'])} -> {f(v['ipw'])}"
            for k, v in sel["ipw_designB_N20"].items())]
        rt = cs["test_retest"]
        out += [f"- retest n={rt['n_students']}, median n "
                f"{rt['median_n_winA']:.0f}/{rt['median_n_winB']:.0f}: " +
                "; ".join(f"{k}: {f(rt[k]['r_obs'])} obs / "
                          f"{f(rt[k]['r_disattenuated'])} disatt"
                          for k in IDX if k in rt)]
        h3s = json.load(open(f"results/h3_shrinkage_{s}.json"))
        out += [f"- H3 {tag}: rho {f(h3s['rho_raw'])} -> {f(h3s['rho_eb'])} "
                f"(diff {f(h3s['diff_point'],4)}, p<={f(h3s['p_onesided_diff_le_0'])}); "
                f"Jaccard top-decile {f(h3s['jaccard_top_decile_raw'])} -> "
                f"{f(h3s['jaccard_top_decile_eb'])}; n={h3s['n_students']}"]
        out += [f"- H3 {tag} VALIDATION (addendum 4b): raw {f(h3s['rho_raw'])} -> "
                f"EB-val {f(h3s['rho_eb_validation'])} "
                f"(diff {f(h3s['diff_validation_point'],4)} "
                f"CI {[f(v,4) for v in h3s['diff_validation_ci95']]}, "
                f"p<={f(h3s['p_onesided_validation_le_0'])}); "
                f"flag validation Jaccard raw {f(h3s['jaccard_validation_raw'])} "
                f"-> EB {f(h3s['jaccard_validation_eb'])}"]
        h3n = json.load(open(f"results/h3_null_{s}.json"))
        out += [f"- H3 {tag} NULL: raw {f(h3n['null_rho_raw'])}, "
                f"EB-sym {f(h3n['null_rho_eb'])} (ARTIFACT), "
                f"EB-validation {f(h3n['null_rho_eb_validation'],4)} (~0 as expected)"]
        rx = json.load(open(f"results/rev_extras_{s}.json"))
        gaps = rx["temporal_gap_days"]
        out += ["- temporal gap days (median between half median-dates): " +
                "; ".join(f"N={k}: {v['median_gap_days']:.0f}"
                          for k, v in gaps.items())]
        for N in ("N20", "N100"):
            bd = rx["bias_decomposition"][N]
            out += [f"- bias decomposition {N}: sb_bias {f(bd['sb_bias'])}, "
                    f"sb_conf {f(bd['sb_conf'])}, sb_acc {f(bd['sb_acc'])}; "
                    f"cross-half r(biasA,biasB) {f(bd['r_biasA_biasB'])}, "
                    f"r(biasA,confB) {f(bd['r_biasA_confB'])}, "
                    f"r(biasA,accB) {f(bd['r_biasA_accB'])}; "
                    f"r(partialed bias, conf same-half) "
                    f"{f(bd['r_partialedBias_conf_samehalf'])}"]

    with open("results/RESULTS_GATE.md", "w", encoding="utf-8") as fh:
        fh.write("\n".join(out) + "\n")
    print("RESULTS_GATE.md written")


if __name__ == "__main__":
    main()
