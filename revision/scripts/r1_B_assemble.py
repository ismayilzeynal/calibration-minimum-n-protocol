"""Stream B: assemble r1_B.json from the part files and print the markdown
tables used in r1_B.md (numbers are copied from this output, not retyped).

Inputs:  revision-2026-11/results/r1_B_parts/{vc,alpha}_{t12,t34}.json
Outputs: revision-2026-11/results/r1_B.json
         revision-2026-11/tmp/r1_B_tables.md (scratch, pasted into r1_B.md)
"""
import json
import os
import platform
import sys

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(REV, "results", "r1_B_parts")
sys.dont_write_bytecode = True
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SAMPLES = {"t12": "S1", "t34": "S2"}
SEM_N = (15, 20, 30, 50)


def load(name):
    with open(os.path.join(PARTS, name)) as fh:
        return json.load(fh)


def f(v, d=3):
    if v is None:
        return "unknown"
    if v == "inf":
        return "> 200"
    return f"{v:.{d}f}"


def ci(v, d=1):
    if v is None or v[0] is None:
        return "unknown"
    hi = "> 200" if v[1] is None else f"{v[1]:.{d}f}"
    return f"[{v[0]:.{d}f}, {hi}]"


def main():
    vc = {s: load(f"vc_{s}.json") for s in SAMPLES}
    al = {s: load(f"alpha_{s}.json") for s in SAMPLES}
    lines = []

    # Table 1: split-half vs ICC at selected N
    lines.append("| Sample | Index | N | Split-half SB (baseline) | ICC(1,N) [95% CI] | ICC - split |")
    lines.append("|---|---|---|---|---|---|")
    for s, lab in SAMPLES.items():
        for k in ("bias", "brier"):
            for row in vc[s]["by_N"]:
                if row["N"] in (5, 10, 15, 20, 30, 50, 100, 200):
                    v = row[k]
                    lines.append(f"| {lab} | {k} | {row['N']} | {f(v['split_half_sb'])} | "
                                 f"{f(v['icc1_N'])} {ci(v['icc1_N_ci'], 3)} | "
                                 f"{v['gap_icc_minus_split']:+.4f} |")
    lines.append("")

    # Table 2: min-N comparison
    lines.append("| Sample | Index | Split-half minN .70 [CI] | ICC curve minN .70 [CI] | "
                 "Split-half minN .80 [CI] | ICC curve minN .80 [CI] | SB prophecy .70 from N=20 | "
                 "VC n*(.70) at N=20 [CI] | VC n*(.70) all responses [CI] | VC n*(.80) all responses [CI] | "
                 "max abs gap (at N) |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|")
    for s, lab in SAMPLES.items():
        for k in ("bias", "brier"):
            sm = vc[s]["summary"][k]
            rc = vc[s]["all_record_moment_estimator"][k]
            lines.append(
                f"| {lab} | {k} | {f(sm['split_minN70_point'], 1)} {ci(sm['split_minN70_ci'])} | "
                f"{f(sm['icc_minN70_point'], 1)} {ci(sm['icc_minN70_ci'])} | "
                f"{f(sm['split_minN80_point'], 1)} {ci(sm['split_minN80_ci'])} | "
                f"{f(sm['icc_minN80_point'], 1)} {ci(sm['icc_minN80_ci'])} | "
                f"{f(sm['split_prophecy_n70_from20'], 1)} | "
                f"{f(sm['vc_implied_n70_at20'], 1)} {ci(sm['vc_implied_n70_at20_ci'])} | "
                f"{f(rc['implied_n70'], 1)} {ci(rc['implied_n70_ci'])} | "
                f"{f(rc['implied_n80'], 1)} {ci(rc['implied_n80_ci'])} | "
                f"{sm['max_abs_gap_icc_vs_split']:.4f} (N={sm['N_at_max_abs_gap']}) |")
    lines.append("")

    # Table 3: variance components
    lines.append("| Sample | Index | sigma_p^2 (all responses) | sigma_e^2 (all responses) | "
                 "single-response ICC(1,1) | pool n | median responses per student |")
    lines.append("|---|---|---|---|---|---|---|")
    for s, lab in SAMPLES.items():
        for k in ("bias", "brier"):
            rc = vc[s]["all_record_moment_estimator"][k]
            lines.append(f"| {lab} | {k} | {rc['sigma_p2']:.5f} | {rc['sigma_e2']:.5f} | "
                         f"{rc['sigma_p2'] / (rc['sigma_p2'] + rc['sigma_e2']):.4f} | "
                         f"{rc['n_students']} | {rc['median_responses']:.0f} |")
    lines.append("")

    # Table 4: SEM for bias (confidence-scale points, 0-100)
    sem = {}
    lines.append("| Sample | N | SD of observed bias | SEM from split-half | SEM from VC | "
                 "median per-student 95% CI half-width | share with zero within-student SD |")
    lines.append("|---|---|---|---|---|---|---|")
    for s, lab in SAMPLES.items():
        sem[s] = {}
        for row in vc[s]["by_N"]:
            if row["N"] in SEM_N:
                v = row["bias"]
                sem[s][str(row["N"])] = dict(
                    sd_obs_points=100 * v["sd_obs"], sem_split_points=100 * v["sem_split"],
                    sem_vc_points=100 * v["sem_vc"],
                    median_ci_halfwidth_points=100 * v["median_ci_halfwidth"],
                    share_sd_zero=v["share_sd_zero"])
                lines.append(f"| {lab} | {row['N']} | {100 * v['sd_obs']:.1f} | "
                             f"{100 * v['sem_split']:.1f} | {100 * v['sem_vc']:.1f} | "
                             f"{100 * v['median_ci_halfwidth']:.1f} | {v['share_sd_zero']:.3f} |")
    lines.append("")

    # Table 5: alpha subset (S2)
    lines.append("| Set | Index | quizzes | (quiz, student) rows | items k (median, range) | alpha | "
                 "two-way Phi | one-way ICC(1,k) | split-half, per-student partitions (paper estimator) | "
                 "split-half, common partition | alpha - per-student split (median, range) | "
                 "Phi - per-student split (median) | "
                 "quizzes alpha > per-student split | CI of difference excludes 0 | "
                 "item share of error variance | n70 from alpha | n70 from Phi | n70 from per-student split |")
    lines.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    a2 = al["t34"]
    for setkey, setlab in (("summary_primary_ge300", ">= 300 students"),
                           ("summary_all_ge100", ">= 100 students")):
        for k in ("bias", "brier"):
            sm = a2[setkey][k]
            lines.append(
                f"| S2 {setlab} | {k} | {sm['n_quizzes']} | {sm['n_students_total']} | "
                f"{sm['k_median']:.0f} ({sm['k_range'][0]}-{sm['k_range'][1]}) | "
                f"{sm['alpha_median']:.3f} | {sm['phi_median']:.3f} | {sm['icc1_k_median']:.3f} | "
                f"{sm['sb_per_student_median']:.3f} | {sm['sb_common_median']:.3f} | "
                f"{sm['alpha_minus_sb_per_student_median']:+.3f} ({sm['alpha_minus_sb_per_student_min']:+.3f} to "
                f"{sm['alpha_minus_sb_per_student_max']:+.3f}) | "
                f"{sm['phi_minus_sb_per_student_median']:+.3f} | "
                f"{sm['n_quizzes_alpha_gt_sb_per_student']} of {sm['n_quizzes']} | "
                f"{sm['n_quizzes_diff_ci_excludes_0']} of {sm['n_quizzes']} | "
                f"{sm['share_item_var_of_error_median']:.3f} | {sm['n70_from_alpha_median']:.1f} | "
                f"{sm['n70_from_phi_median']:.1f} | "
                f"{sm['n70_from_sb_per_student_median']:.1f} |")
    lines.append("")

    headline = {}
    for s, lab in SAMPLES.items():
        for k in ("bias", "brier"):
            sm = vc[s]["summary"][k]
            rc = vc[s]["all_record_moment_estimator"][k]
            headline[f"{lab}_{k}"] = dict(
                split_minN70=sm["split_minN70_point"], icc_minN70=sm["icc_minN70_point"],
                icc_minN70_ci=sm["icc_minN70_ci"],
                split_minN80=sm["split_minN80_point"], icc_minN80=sm["icc_minN80_point"],
                icc_minN80_ci=sm["icc_minN80_ci"],
                vc_all_responses_n70=rc["implied_n70"], vc_all_responses_n70_ci=rc["implied_n70_ci"],
                vc_all_responses_n80=rc["implied_n80"], vc_all_responses_n80_ci=rc["implied_n80_ci"],
                max_abs_gap_icc_vs_split=sm["max_abs_gap_icc_vs_split"],
                N_at_max_abs_gap=sm["N_at_max_abs_gap"],
                mean_abs_gap_odd_N=sm["mean_abs_gap_odd_N"],
                mean_abs_gap_even_N=sm["mean_abs_gap_even_N"],
                n_grid_points=len(vc[s]["by_N"]),
                n_grid_points_split_inside_icc_ci=int(sum(
                    1 for row in vc[s]["by_N"]
                    if row[k]["icc1_N_ci"][0] <= row[k]["split_half_sb"] <= row[k]["icc1_N_ci"][1])),
                single_response_icc=rc["sigma_p2"] / (rc["sigma_p2"] + rc["sigma_e2"]))
    for k in ("bias", "brier"):
        sm = a2["summary_primary_ge300"][k]
        headline[f"S2_alpha_subset_{k}"] = {m: sm[m] for m in (
            "n_quizzes", "n_students_total", "k_median", "k_range", "alpha_median", "phi_median",
            "icc1_k_median", "sb_per_student_median", "sb_common_median",
            "alpha_minus_sb_per_student_median", "alpha_minus_sb_per_student_min",
            "alpha_minus_sb_per_student_max", "n_quizzes_alpha_gt_sb_per_student",
            "n_quizzes_diff_ci_excludes_0", "sb_per_student_minus_icc1_k_median",
            "sb_common_minus_alpha_median", "share_item_var_of_error_median",
            "n70_from_alpha_median", "n70_from_phi_median", "n70_from_sb_per_student_median",
            "phi_minus_sb_per_student_median")}

    gates = {s: vc[s]["gates"] for s in SAMPLES}
    alpha_checks = {k: dict(
        max_abs_alpha_minus_icc_C_k=a2["summary_all_ge100"][k]["max_abs_alpha_minus_icc_C_k"],
        n_quizzes_rulon_checked=a2["summary_all_ge100"][k]["n_quizzes_rulon_checked"],
        max_abs_rulon_minus_alpha=a2["summary_all_ge100"][k]["max_abs_rulon_minus_alpha"])
        for k in ("bias", "brier")}

    runtime = {f"vc_{s}": vc[s]["runtime_seconds"]["total"] for s in SAMPLES}
    runtime.update({f"alpha_{s}": al[s]["runtime_seconds"]["total"] for s in SAMPLES})
    runtime["total_compute_seconds"] = float(sum(runtime.values()))

    out = dict(
        stream="B (E4, R1-2): reliability framework justification",
        paper="PARE article 4262, major revision 2026-11-23",
        environment=dict(python=platform.python_version(), numpy=np.__version__,
                         pandas=pd.__version__),
        scripts=["revision-2026-11/scripts/r1_B_vc.py", "revision-2026-11/scripts/r1_B_alpha.py",
                 "revision-2026-11/scripts/r1_B_assemble.py"],
        inputs=dict(parquet_t12_sha256="c27099e0ccd4072918085bfdea984a769753cdbc27c26ab38fada58759559ce5",
                    parquet_t34_sha256="49971f533a40f2a8b7fee769c22c931ca3085693f9eeed5d9d01c6d8dfd1aaab",
                    baseline_curves="revision-2026-11/results/baseline/results/reliability_curves_{t12,t34}.csv",
                    baseline_minN="revision-2026-11/results/baseline/results/minN_{t12,t34}.json",
                    sha256_source="revision-2026-11/results/baseline/REPRO.md (not recomputed here)"),
        definitions=dict(
            item_scores="d = c - x and q = (c - x)^2 per tagged response, c = LEVELS[round(Confidence/25)] in {0,.25,.5,.75,1}, x = IsCorrect; bias = mean(d), Brier = mean(q)",
            design="Design A pool (>= 200 tagged responses); same N-subsamples as the submitted curves (default_rng([42, 5, N, sid]), sid 1 = S1/t12, 2 = S2/t34)",
            split_half_sb="mean over 50 random per-student splits of the Spearman-Brown corrected Pearson r (reliability.py), recomputed here on the same generator and equal to the baseline curves (gate G1)",
            icc1_N="one-way random-effects ICC(1,N) = (MSB - MSW) / MSB on the n x N subsample matrix; = G coefficient of an r:p (responses nested in persons) design; equals sigma_p^2 / (sigma_p^2 + sigma_e^2 / N)",
            sigma_p2="person variance (MSB - MSW) / N (per-N) or var(m_i) - mean(s_i^2 / n_i) (all responses of each pool student, unweighted moment estimator, as in the paper's EB model)",
            sigma_e2="within-person residual variance of the item-level score (MSW per-N; mean s_i^2 for all responses)",
            implied_n="n*(t) = t / (1 - t) * sigma_e^2 / sigma_p^2, the volume at which sigma_p^2 / (sigma_p^2 + sigma_e^2 / n) = t",
            icc_curve_minN="crossing of the ICC(1,N) grid curve with the same rule as the paper (first crossing, linear in log N); CI from a joint student bootstrap over the grid (B = 1000, default_rng([42, 101, sid]))",
            vc_all_responses_ci="student bootstrap B = 1000, default_rng([42, 102, sid])",
            sem_split="SD of observed full-N bias across pool students x sqrt(1 - split-half SB)",
            sem_vc="sqrt(MSW / N) = SD_obs x sqrt(1 - ICC(1,N)) exactly",
            median_ci_halfwidth="median over pool students of 1.96 x s_i / sqrt(N), s_i = within-student SD of d on the N-subsample (ddof 1)",
            points="SEM tables are in confidence-scale points (x 100; bias of 0.10 = 10 points)",
            alpha_subset="S2 quizzes: students who rated EVERY question of the quiz (selected subset); complete crossed students x items matrix; QuizId -1 excluded; >= 100 complete students and >= 4 items (primary set >= 300)",
            alpha="Cronbach alpha on item-level d (bias) or q (Brier); self-checked equal to two-way ICC(C,k) and to the mean Flanagan-Rulon coefficient over all equal splits (even k <= 14)",
            phi="G-theory absolute coefficient for the crossed p x i design, item main effect counted as error",
            icc1_k="one-way ICC(1,k) on the quiz matrix, items treated as nested (item variance in error)",
            sb_per_student="the paper's estimator on the quiz matrix: per-student random partitions (reliability.half_indices), 50 splits",
            sb_common="one common random item partition for all students, 50 partitions, mean SB",
            share_item_var_of_error="sigma_i^2 / (sigma_i^2 + sigma_pi,e^2) in the crossed G-study",
            n70_from_m="k x (.7/.3) x (1 - m) / m",
            quiz_student_rows="n_students_total in the alpha summaries = number of (quiz, student) rows summed over quizzes; one student can appear in several quizzes, so this is NOT a count of distinct students",
            relative_vs_absolute="alpha (= E rho^2) is the coefficient for relative (rank-order) use on a common form; Phi is the coefficient for absolute use (bias against zero, e.g. calling a student overconfident), where item effects do not cancel even on a crossed form",
        ),
        not_computable=dict(
            gamma_auroc_alpha="Not defined. Gamma and AUROC are pair-count ratios over a student's correct-incorrect pairs, not sums or means of item-level scores, so there is no item-level score whose covariance matrix alpha or a variance-components model could use. Not approximated.",
            gamma_auroc_vc="Not defined for the same reason.",
            s1_alpha="Not computed. No S1 quiz has >= 100 students who rated every item; the largest has 70 students with 5 items (feasibility block).",
            temporal_split="Not part of stream B: the variance-components model assumes exchangeable responses, which temporal halves violate by design.",
        ),
        auroc_rule_independence="bias and Brier do not use the AUROC/gamma tie rule; every number here is unaffected by stream A.",
        gates=gates,
        alpha_self_checks=alpha_checks,
        headline=headline,
        sem_bias_points=sem,
        variance_components=vc,
        alpha_subset=al,
        runtime_seconds=runtime,
    )
    with open(os.path.join(REV, "results", "r1_B.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    os.makedirs(os.path.join(REV, "tmp"), exist_ok=True)
    with open(os.path.join(REV, "tmp", "r1_B_tables.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines))
    print("\n".join(lines))
    print(json.dumps(runtime))


if __name__ == "__main__":
    main()
