"""Stream C: assemble the part files into results/r1_C.json (all numbers the
manuscript might quote, with definitions) and print compact tables.

Usage: python r1_C_assemble.py
Reads:  revision-2026-11/results/C/r1_C_selection_{t12,t34}.json,
        revision-2026-11/results/C/r1_C_splits_{t12,t34}.json
Writes: revision-2026-11/results/r1_C.json
"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_C_common as RC  # noqa: E402

SAMPLES = {"t12": "S1", "t34": "S2"}
OUT = os.path.join(RC.REV, "results", "r1_C.json")

DEFINITIONS = {
    "S1, S2": "Sample 1 = Eedi task 1-2 files (t12), Sample 2 = task 3-4 files (t34), as in the submitted paper.",
    "tagged answer": "an answer with a confidence rating, Confidence non-missing and in [0, 100] (common.tagged).",
    "reporter": "a student with at least one tagged answer.",
    "item p-value (p_all)": "proportion correct for the QuestionId over ALL answers in the sample.",
    "item p-value (p_loo)": "leave-one-out version: proportion correct for the QuestionId over all OTHER answers "
                            "in the sample (removes the answer's own correctness). Primary difficulty measure.",
    "relative accuracy (rel)": "IsCorrect - p_loo per answer; mean = accuracy beyond what item difficulty predicts.",
    "strand": "level-1 subject of the question in subject_metadata.csv (Number, Algebra, Geometry and Measure, "
              "Data and Statistics, ...). A question listing several level-1 subjects gets the combined label "
              "(e.g. 'Number+Algebra') as its own category.",
    "topic": "level-2 subject of the question (69 level-2 subjects in the tree), combination rule as for strand.",
    "term": "approximate English school term from the calendar month of DateAnswered: Autumn Sep-Dec, Spring Jan-Mar, "
            "Summer Apr-Aug, within academic year. Approximation; actual term dates are not in the data.",
    "difficulty strata": "global answer-level quintiles of p_loo over all answers of the sample (edges in "
                         "weighting.difficulty_strata_edges_p_loo); answers with undefined p_loo form a sixth stratum.",
    "paired difference": "per student: mean over tagged answers minus mean over untagged answers; summary over "
                         "students (mean, median, SD, IQR, share > 0) with a student-bootstrap 95% CI of the mean "
                         "(B = 1000).",
    "basis all_untagged": "students with >= 5 tagged and >= 5 untagged answers; all their untagged answers.",
    "basis same_quiz": "only answers in student x QuizId pairs that contain both tagged and untagged answers; "
                       "students with >= 5 tagged and >= 5 untagged such answers. Isolates which answers a student "
                       "rates when rating was evidently available in that quiz.",
    "TVD": "total variation distance between the student's strand (or topic) distribution on tagged vs untagged "
           "answers, 0 = identical, 1 = disjoint. The permutation baseline shuffles the tag labels within student "
           "once (same counts), so it is the TVD expected from sampling noise alone.",
    "eta2": "share of the variance of the 0/1 tag indicator that lies between groups (student, quiz, question, "
            "student x quiz, calendar month).",
    "Design A": "pool of students with >= 200 tagged answers (S1 1,389; S2 554); the same per-N subsamples as the "
                "submitted analysis (seed [42, 5, N, sid]).",
    "reliability": "Spearman-Brown corrected split-half correlation across students, averaged over 50 splits "
                   "(temporal: 1 split), exactly as reliability.py.",
    "minN70 / minN80": "first crossing of r = .70 / .80 on the Design A grid, log-linear interpolation, as "
                       "reliability.py; 'inf' = not reached by N = 200. CI = percentile over the 1000 joint "
                       "bootstrap draws of the submitted analysis; an upper bound of null = some draws never cross.",
    "minus_random": "paired difference to the random-split curve of the same subsamples, per bootstrap draw "
                    "(CI over draws where both curves cross).",
    "w_diff": "each tagged answer weighted by P_untagged(stratum) / P_tagged(stratum) within student, difficulty "
              "strata; weights renormalised to mean 1 per student. Students with < 20 untagged answers keep "
              "weight 1 (no reliable target). Indices are weighted means of the same responses.",
    "w_diffstrand": "as w_diff with strata = difficulty quintile x strand.",
    "balanced": "difficulty-balanced subset: the N responses are drawn without replacement with probability "
                "proportional to the w_diff weight (balanced_ds: w_diffstrand weight), then 50 random splits, "
                "unweighted indices. Equal-weight alternative to w_diff.",
    "strand / topic / term split": "half A = the first N//2 responses after ordering the student's categories at "
                                   "random (whole categories go to one half; only one category can straddle the "
                                   "cut); 50 random orders.",
    "difficulty split": "half A = the N//2 responses on the hardest items (lowest p_loo), half B = the rest; ties "
                        "broken at random, 50 repetitions.",
    "shared_category_response_share": "share of a student's N responses whose category occurs in both halves; "
                                      "measures how content-disjoint the halves are (keyed vs random split).",
    "restricted_largest_category_le_half": "point estimates on the students whose largest category holds <= 50% "
                                           "of their N responses, so that the halves are clearly content-separated; "
                                           "random, keyed and temporal on the same students.",
    "target_subset": "point estimates on the Design A students that have a weighting target for both weight "
                     "schemes (>= 20 untagged answers), raw vs weighted vs balanced on the same students.",
    "fullset_bias_shift": "bias on ALL of a student's tagged answers, raw vs difficulty-reweighted, for students "
                          "with a weighting target; mean shift, Pearson r and top-decile Jaccard.",
    "kish_ess_ratio": "(sum w)^2 / sum w^2 / N per student: effective share of the N responses left after weighting.",
    "w_null_diff / w_null_diffstrand": "null weights: each student's own w_diff (w_diffstrand) weight vector permuted "
                                       "at random across that student's rated answers. Same weight distribution and "
                                       "Kish ESS, no link to difficulty or strand. The gap between a real-weight result "
                                       "and its null is the composition effect; the gap between the null and the raw "
                                       "result is the cost of unequal weights alone.",
}

NOT_COMPUTED = [
    "Gamma and AUROC under structured splits and reweighting: not computed in this stream, because the AUROC "
    "constant-confidence rule is being changed in stream A, and both rank indices are below .30 at N = 20 under "
    "random splits in the submitted analysis (headline: baseline_r20_gamma, baseline_r20_auroc).",
    "Selection on unobservables (for example, a student skipping the rating when unsure of the answer) cannot "
    "be separated from difficulty or content in these data: confidence on unrated answers is not observed.",
    "Monitoring on untagged answers is unobservable; the reweighting assumes that, within a difficulty (x strand) "
    "stratum, a student's tagged answers stand for the untagged ones.",
    "Whether a quiz offered the rating prompt to a student is not recorded; the same-quiz basis uses quiz "
    "attempts in which the student rated at least one answer as the evidence of availability.",
    "Gender codes are reported raw; their labels are not verified (unknown).",
    "Term split uses calendar-month approximations of English school terms; actual term dates are unknown.",
]


def load(name):
    with open(os.path.join(RC.OUT_DIR, name)) as fh:
        return json.load(fh)


def headline(X, s):
    import pandas as pd
    h = {}
    T = X["part1_tagging_structure"]["all_answers"]
    h["eta2_tag_student"] = T["eta2_student"]
    h["eta2_tag_student_x_quiz"] = T["eta2_student_x_quiz"]
    h["share_student_quiz_pairs_mixed"] = T["share_student_quiz_pairs_mixed"]
    TR = X["part1_tagging_structure"]["reporters_answers"]
    for k in ["student", "student_x_quiz", "quiz", "question", "month"]:
        h[f"reporters_eta2_tag_{k}"] = TR[f"eta2_{k}"]
    h["reporters_within_student_x_quiz_share"] = 1.0 - TR["eta2_student_x_quiz"]
    h["reporters_share_student_quiz_pairs_mixed"] = TR["share_student_quiz_pairs_mixed"]
    for b in ["all_untagged", "same_quiz"]:
        P = X[f"part1_paired_{b}"]
        h[f"{b}_n_students"] = P["n_students"]
        for k, nm in [("acc", "accuracy_gap"), ("p", "item_p_gap"), ("rel", "beyond_difficulty_gap"),
                      ("days", "days_gap"), ("wkend", "weekend_share_gap")]:
            h[f"{b}_{nm}"] = P[f"{k}_tagged_minus_untagged"]["mean"]
            h[f"{b}_{nm}_ci95"] = P[f"{k}_tagged_minus_untagged"]["mean_ci95"]
        for c in ["strand", "topic"]:
            h[f"{b}_tvd_{c}"] = P[f"tvd_{c}"]["mean"]
            h[f"{b}_tvd_{c}_permutation_baseline"] = P[f"tvd_{c}_permutation_baseline"]["mean"]
    h["same_quiz_relpos_gap"] = X["part1_paired_same_quiz"]["relpos_tagged_minus_untagged"]["mean"]
    h["same_quiz_relpos_gap_ci95"] = X["part1_paired_same_quiz"]["relpos_tagged_minus_untagged"]["mean_ci95"]
    RW = X["part2_reweighting"]
    for v in ["random", "w_diff_random", "w_diffstrand_random", "w_null_diff_random", "w_null_diffstrand_random",
              "balanced_random", "balanced_ds_random"]:
        b = RW["variants"][v]["bias"]
        for t in ["70", "80"]:
            h[f"bias_minN{t}_{v}"] = b[f"minN{t}_point"]
            h[f"bias_minN{t}_{v}_ci"] = b.get(f"minN{t}_ci")
            if v != "random":
                h[f"bias_minN{t}_change_{v}"] = b[f"minN{t}_minus_random"]["point"]
                h[f"bias_minN{t}_change_{v}_ci"] = b[f"minN{t}_minus_random"]["ci_finite_draws"]
        h[f"bias_r20_{v}"] = b["r20_point"]
    w20 = [x for x in RW["weights_diagnostics_by_N"] if x["N"] == 20][0]
    h["kish_ratio_N20_w_diff_target_students"] = w20["kish_ess_ratio_mean_target_students"]
    h["kish_ratio_N20_w_diffstrand_target_students"] = w20["diffstrand_kish_ess_ratio_mean_target_students"]
    ts = RW["target_subset_minN"]
    h["target_subset_n_students"] = [x for x in RW["target_subset_curves"] if x["N"] == 20][0]["n_students"]
    for v in ["random", "w_diff_random", "w_diffstrand_random", "w_null_diff_random", "w_null_diffstrand_random",
              "balanced_random", "balanced_ds_random"]:
        h[f"target_subset_bias_minN70_{v}"] = ts[f"{v}_bias"]["minN70_point"]
    b20 = [x for x in RW["balanced_diagnostics_by_N"] if x["N"] == 20][0]
    h["N20_tvd_to_unrated_difficulty_mix_plain"] = b20["tvd_to_untagged_mix_plain_mean"]
    h["N20_tvd_to_unrated_difficulty_mix_balanced"] = b20["tvd_to_untagged_mix_balanced_mean"]
    h["N20_tvd_to_unrated_diffstrand_mix_plain"] = b20["ds_tvd_to_untagged_mix_plain_mean"]
    h["N20_tvd_to_unrated_diffstrand_mix_balanced_ds"] = b20["ds_tvd_to_untagged_mix_balanced_ds_mean"]
    for pool in ["ge5_tagged", "ge20_tagged", "designA_pool"]:
        for wl in ["diff", "diffstrand", "null_diff", "null_diffstrand"]:
            f = RW["fullset_bias_shift"][pool][wl]
            h[f"fullset_{pool}_{wl}_n"] = f["n_students"]
            h[f"fullset_{pool}_{wl}_bias_shift"] = f["weighted_minus_raw"]["mean"]
            h[f"fullset_{pool}_{wl}_bias_shift_ci95"] = f["weighted_minus_raw"]["mean_ci95"]
            h[f"fullset_{pool}_{wl}_pearson"] = f["pearson_raw_weighted"]
            h[f"fullset_{pool}_{wl}_top_decile_jaccard"] = f["top_decile_jaccard_raw_vs_weighted"]
    SP = X["part3_structured_splits"]
    for v in ["random", "temporal", "strand", "topic", "term", "difficulty"]:
        for k in ["bias", "brier"]:
            b = SP["variants"][v][k]
            h[f"{k}_r20_{v}"] = b["r20_point"]
            h[f"{k}_r20_{v}_ci"] = b["r20_ci"]
            h[f"{k}_minN70_{v}"] = b["minN70_point"]
            h[f"{k}_minN80_{v}"] = b["minN80_point"]
            if v != "random":
                h[f"{k}_r20_change_{v}"] = b["r20_minus_random"]["point"]
                h[f"{k}_r20_change_{v}_ci"] = b["r20_minus_random"]["ci"]
                h[f"{k}_minN70_change_{v}"] = b["minN70_minus_random"]["point"]
                h[f"{k}_minN70_change_{v}_ci"] = b["minN70_minus_random"]["ci_finite_draws"]
                h[f"{k}_minN80_change_{v}"] = b["minN80_minus_random"]["point"]
                h[f"{k}_minN80_change_{v}_ci"] = b["minN80_minus_random"]["ci_finite_draws"]
        for k in ["conf", "acc"]:
            h[f"{k}_r20_{v}"] = SP["variants"][v][k]["r20_point"]
    base = pd.read_csv(os.path.join(RC.REV, "results", "baseline", "results", f"reliability_curves_{s}.csv"))
    for k in ["gamma", "auroc"]:
        h[f"baseline_r20_{k}"] = float(base[(base.design == "A") & (base.split == "random")
                                            & (base["index"] == k) & (base.N == 20)].r.iloc[0])
    return h


def main():
    out = {"stream": "C - selection and content heterogeneity (E1b, E3c, E3d, R1-3, R1-5)",
           "definitions": DEFINITIONS, "not_computed": NOT_COMPUTED,
           "scripts": ["revision-2026-11/scripts/r1_C_common.py", "revision-2026-11/scripts/r1_C_selection.py",
                       "revision-2026-11/scripts/r1_C_splits.py", "revision-2026-11/scripts/r1_C_assemble.py"],
           "inputs": {"parquet": "data/resp_{t12,t34}.parquet (paper folder, read-only)",
                      "metadata": "revision-2026-11/data/question_metadata_task_{1_2,3_4}.csv, subject_metadata.csv, "
                                  "student_metadata_task_{1_2,3_4}.csv (hashes in METADATA_SOURCE.json)",
                      "baseline": "revision-2026-11/results/baseline (regression check)"},
           "runtime_seconds": {}}
    for s, lab in SAMPLES.items():
        sel = load(f"r1_C_selection_{s}.json")
        spl = load(f"r1_C_splits_{s}.json")
        out[lab] = {"sample_code": s,
                    "part4_sample_descriptives": sel["descriptives"],
                    "part1_tagging_structure": sel["tagging_structure"],
                    "part1_paired_all_untagged": sel["paired_all_untagged"],
                    "part1_paired_same_quiz": sel["paired_same_quiz"],
                    "part1_designA_pool_difficulty": sel["designA_pool_difficulty"],
                    "part2_reweighting": {
                        "variants": {v: spl["summary"][v] for v in
                                     ["random", "w_diff_random", "w_diff_temporal", "w_diffstrand_random",
                                      "w_null_diff_random", "w_null_diffstrand_random",
                                      "balanced_random", "balanced_ds_random"]},
                        "target_subset_minN": spl["target_subset_minN"],
                        "target_subset_curves": spl["target_subset_curves"],
                        "weighting": spl["weighting"],
                        "weights_diagnostics_by_N": spl["diagnostics"]["weights"],
                        "balanced_diagnostics_by_N": spl["diagnostics"]["balanced"],
                        "fullset_bias_shift": spl["fullset_bias_shift"]},
                    "part3_structured_splits": {
                        "variants": {v: spl["summary"][v] for v in
                                     ["random", "temporal", "strand", "topic", "term", "difficulty"]},
                        "diagnostics_by_N": {k: spl["diagnostics"][k] for k in ["strand", "topic", "term", "difficulty"]},
                        "restricted_largest_category_le_half": spl["restricted_largest_category_le_half"],
                        "category_labels": spl["category_labels"]},
                    "curves_all_variants": spl["curves"], "grid": spl["grid"], "n_poolA": spl["n_poolA"],
                    "seeds": spl["seeds"],
                    "regression_max_abs_diff_vs_baseline": spl["regression_max_abs_diff_vs_baseline"]}
        out["runtime_seconds"][f"{lab}_selection"] = sel["runtime_seconds"]["total"]
        out["runtime_seconds"][f"{lab}_splits"] = spl["runtime_seconds"]["total"]
    out["runtime_seconds"]["total_compute"] = sum(v for v in out["runtime_seconds"].values())
    out["headline"] = {lab: headline(out[lab], s) for s, lab in SAMPLES.items()}
    out["definitions"]["headline"] = (
        "flat copies of the most quotable numbers. Differences in accuracy / item p / shares are on the 0-1 scale "
        "(x 100 = percentage points); bias is on the 0-1 confidence scale. reporters_* = computed on the answers "
        "of reporters only (students with >= 1 rated answer). baseline_r20_gamma / _auroc = submitted Design A random-split values at N = 20 (AUROC under the submitted "
        "constant-confidence rule), read from results/baseline.")
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print("written", OUT)


if __name__ == "__main__":
    main()
