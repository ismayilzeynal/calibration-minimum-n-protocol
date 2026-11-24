"""Stream C: write results/r1_C.md from results/r1_C.json. Every number in
the summary is formatted from the JSON, nothing is typed by hand.

Usage: python r1_C_report.py
"""
import json
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_C_common as RC  # noqa: E402

J = os.path.join(RC.REV, "results", "r1_C.json")
MD = os.path.join(RC.REV, "results", "r1_C.md")


def f3(v):
    if v is None:
        return "n/a"
    s = f"{v:.3f}"
    if s.startswith("0."):
        return s[1:]
    if s.startswith("-0."):
        return "-" + s[2:]
    return s


def f1(v):
    if v is None:
        return "n/a"
    if v == "inf":
        return "> 200"
    return f"{v:.1f}"


def pct(v, d=1):
    return "n/a" if v is None else f"{100 * v:.{d}f}%"


def pts(v, d=1):
    """proportion difference -> percentage points"""
    return "n/a" if v is None else f"{100 * v:+.{d}f}"


def ci1(c):
    if c is None:
        return "[n/a]"
    lo = f1(c[0]) if c[0] is not None else "n/a"
    hi = f1(c[1]) if c[1] is not None else "> 200"
    return f"[{lo}, {hi}]"


def cipts(c, d=1):
    return f"[{100 * c[0]:+.{d}f}, {100 * c[1]:+.{d}f}]"


def ci3(c):
    return f"[{f3(c[0])}, {f3(c[1])}]"


def main():
    d = json.load(open(J))
    S = {"S1": d["S1"], "S2": d["S2"]}
    S["S1"]["H"], S["S2"]["H"] = d["headline"]["S1"], d["headline"]["S2"]
    L = []
    w = L.append
    w("# Stream C - selection and content heterogeneity (revision analyses)")
    w("")
    w("Machine-readable results with definitions: `results/r1_C.json`. Scripts: "
      "`scripts/r1_C_common.py`, `r1_C_selection.py`, `r1_C_splits.py`, `r1_C_assemble.py`, "
      "`r1_C_report.py` (this file is generated; no number is typed by hand).")
    rt = d["runtime_seconds"]
    w("")
    w(f"Runtime (single process each): S1 selection {rt['S1_selection']:.0f} s, S1 splits {rt['S1_splits']:.0f} s, "
      f"S2 selection {rt['S2_selection']:.0f} s, S2 splits {rt['S2_splits']:.0f} s; total {rt['total_compute']:.0f} s.")
    w("")
    w("## What was computed, and how")
    w("")
    w("1. **Who rates which answers.** Within each student, tagged answers were compared with untagged answers on item "
      "difficulty (leave-one-out item p-value), accuracy, accuracy beyond item difficulty (IsCorrect - p), strand and "
      "topic mix (total variation distance, TVD, against a within-student permutation baseline) and timing (days since "
      "the student's first answer; weekend share; position within the quiz). Two bases: all untagged answers "
      "(students with >= 5 of each, the basis of the submitted +4.0-point result), and the same quiz attempts only "
      "(student x quiz pairs that mix rated and unrated answers). Also: where the variation in tagging sits (eta-squared).")
    w("2. **Difficulty reweighting sensitivity of the bias minimum N.** On the Design A pool and the identical "
      "subsamples and splits of the submitted analysis, each tagged answer was reweighted so that the student's "
      "difficulty-stratum mix (and, as a second variant, difficulty x strand mix) matches that student's untagged "
      "answers; an equal-weight alternative draws difficulty-balanced subsets instead, and a permuted-weight null "
      "(each student's own weights shuffled across their rated answers: same effective sample size, no link to "
      "difficulty) separates composition effects from the cost of unequal weights. Same joint bootstrap draws.")
    w("3. **Structured splits.** Instead of random halves, halves were formed by whole strands, whole topics, whole "
      "(approximate) school terms, or harder vs easier items, on the same subsamples; compared with the random and "
      "temporal splits already reported.")
    w("4. **Sample 1 vs Sample 2 table.**")
    w("")
    w("Regression check: the regenerated random and temporal curves match the baseline reproduction to at most "
      f"{max(max(S[s]['regression_max_abs_diff_vs_baseline'][k] for k in S[s]['regression_max_abs_diff_vs_baseline']) for s in S):.1e} "
      "(summation order); the random-split bootstrap draws match the baseline cache exactly (0).")
    w("")

    # ---------------- Table 4 ----------------
    w("## (4) Sample 1 vs Sample 2")
    w("")
    w("| Quantity | Sample 1 | Sample 2 |")
    w("|---|---|---|")
    D = {k: S[k]["part4_sample_descriptives"] for k in S}

    def row(lab, fn):
        w(f"| {lab} | {fn(D['S1'], 'S1')} | {fn(D['S2'], 'S2')} |")
    row("Students; answers", lambda x, s: f"{x['n_students']:,}; {x['n_answers']:,}")
    row("Questions; quizzes; groups", lambda x, s: f"{x['n_questions']:,}; {x['n_quizzes']:,}; {x['n_groups']:,}")
    row("Answers with a confidence rating", lambda x, s: f"{x['n_tagged']:,} ({pct(x['share_answers_tagged'])})")
    row("Students with any rating (reporters)", lambda x, s: f"{x['n_reporters']:,} ({pct(x['share_students_reporting'])})")
    row("Median (mean) rated answers per reporter",
        lambda x, s: f"{x['median_tagged_per_reporter']:.0f} ({x['mean_tagged_per_reporter']:.1f})")
    row("Reporters with >= 20 / >= 200 rated", lambda x, s: f"{x['reporters_with_ge']['20']:,} / {x['reporters_with_ge']['200']:,}")
    row("Median answers per student", lambda x, s: f"{x['median_answers_per_student']:.0f}")
    row("Accuracy: all / rated / unrated answers",
        lambda x, s: f"{pct(x['accuracy_all'])} / {pct(x['accuracy_tagged'])} / {pct(x['accuracy_untagged'])}")
    row("Mean confidence (rated); share at 100",
        lambda x, s: f"{x['mean_confidence_tagged']:.1f}; {pct(x['confidence_level_shares']['100'])}")
    row("Answer window", lambda x, s: f"{x['first_answer'][:7]} to {x['last_answer'][:7]}")
    row("Rated-answer window", lambda x, s: f"{x['first_tagged'][:10]} to {x['last_tagged'][:10]}")
    row("Item p-value, question level: median [IQR]",
        lambda x, s: f"{f3(x['item_p_question_level']['median'])} [{f3(x['item_p_question_level']['q25'])}, {f3(x['item_p_question_level']['q75'])}]")
    row("Strand labels; topic labels (incl. combined labels for multi-subject questions)",
        lambda x, s: f"{x['n_strand_categories']}; {x['n_topic_categories']}")
    for st in ["Number", "Algebra", "Geometry and Measure", "Data and Statistics"]:
        row(f"Share of answers: {st} (all / rated)",
            lambda x, s, st=st: (f"{pct(x['strand_share_all_answers'].get(st))} / {pct(x['strand_share_tagged_answers'].get(st))}"
                                 if st in x['strand_share_all_answers'] else "0 (not in pool)"))
    row("Median age at first answer, years [IQR] (birth date known)",
        lambda x, s: (f"{x['demographics_all_students']['median_age_years_at_first_answer']:.1f} "
                      f"[{x['demographics_all_students']['age_q25_q75'][0]:.1f}, {x['demographics_all_students']['age_q25_q75'][1]:.1f}] "
                      f"({pct(x['demographics_all_students']['share_dob_known'], 0)} known)"))
    row("PremiumPupil = 1 among known (share known)",
        lambda x, s: (f"{pct(x['demographics_all_students']['share_premiumpupil_1_among_known'])} "
                      f"({pct(x['demographics_all_students']['share_premiumpupil_known'], 0)})"))
    T = {k: S[k]["part1_tagging_structure"]["reporters_answers"] for k in S}
    w(f"| Reporters' answers: variance of the rating indicator between students / between student x quiz pairs / "
      f"between months (eta-squared) | {f3(T['S1']['eta2_student'])} / {f3(T['S1']['eta2_student_x_quiz'])} / "
      f"{f3(T['S1']['eta2_month'])} | {f3(T['S2']['eta2_student'])} / {f3(T['S2']['eta2_student_x_quiz'])} / "
      f"{f3(T['S2']['eta2_month'])} |")
    w(f"| Reporters' student x quiz pairs with some but not all answers rated | {pct(T['S1']['share_student_quiz_pairs_mixed'])} | "
      f"{pct(T['S2']['share_student_quiz_pairs_mixed'])} |")
    w("")
    w("Month-level rating share (share of answers rated): S1 "
      + ", ".join(f"{m} {pct(v, 0)}" for m, v in D["S1"]["tag_share_by_month"].items() if m in ("2018-12", "2019-01", "2019-02", "2019-09", "2020-04"))
      + "; S2 " + ", ".join(f"{m} {pct(v, 0)}" for m, v in D["S2"]["tag_share_by_month"].items() if m in ("2019-09", "2019-10", "2020-01", "2020-04"))
      + ". Gender is stored as raw codes with unverified labels and is not tabulated here.")
    w("")

    # ---------------- Table 1 ----------------
    w("## (1) Rated vs unrated answers within student")
    w("")
    w("Paired differences, rated minus unrated, mean over students [student-bootstrap 95% CI]; percentage points "
      "unless stated.")
    w("")
    w("| Quantity | S1 all unrated | S1 same quiz | S2 all unrated | S2 same quiz |")
    w("|---|---|---|---|---|")
    P = {(s, b): S[s][f"part1_paired_{b}"] for s in S for b in ["all_untagged", "same_quiz"]}
    cols = [("S1", "all_untagged"), ("S1", "same_quiz"), ("S2", "all_untagged"), ("S2", "same_quiz")]
    w("| Students | " + " | ".join(f"{P[c]['n_students']:,}" for c in cols) + " |")
    for key, lab in [("acc", "Accuracy"), ("p", "Item p-value (leave-one-out; + = rated items easier)"),
                     ("rel", "Accuracy beyond item difficulty")]:
        w(f"| {lab} | " + " | ".join(f"{pts(P[c][key + '_tagged_minus_untagged']['mean'])} {cipts(P[c][key + '_tagged_minus_untagged']['mean_ci95'])}" for c in cols) + " |")
    w("| Days since the student's first answer (days) | " + " | ".join(
        f"{P[c]['days_tagged_minus_untagged']['mean']:+.1f} [{P[c]['days_tagged_minus_untagged']['mean_ci95'][0]:+.1f}, {P[c]['days_tagged_minus_untagged']['mean_ci95'][1]:+.1f}]" for c in cols) + " |")
    w("| Weekend share | " + " | ".join(f"{pts(P[c]['wkend_tagged_minus_untagged']['mean'])} {cipts(P[c]['wkend_tagged_minus_untagged']['mean_ci95'])}" for c in cols) + " |")
    w("| Relative position within the quiz (0 = first, 1 = last) | n/a | "
      f"{s3(P[('S1', 'same_quiz')]['relpos_tagged_minus_untagged']['mean'])} {ci3(P[('S1', 'same_quiz')]['relpos_tagged_minus_untagged']['mean_ci95'])} | n/a | "
      f"{s3(P[('S2', 'same_quiz')]['relpos_tagged_minus_untagged']['mean'])} {ci3(P[('S2', 'same_quiz')]['relpos_tagged_minus_untagged']['mean_ci95'])} |")
    for cat in ["strand", "topic"]:
        w(f"| {cat.capitalize()} mix: TVD observed / permutation baseline | " + " | ".join(
            f"{f3(P[c]['tvd_' + cat]['mean'])} / {f3(P[c]['tvd_' + cat + '_permutation_baseline']['mean'])}" for c in cols) + " |")
    w("| Share of students with rated items easier (p difference > 0) | " + " | ".join(
        pct(P[c]["p_tagged_minus_untagged"]["share_positive"]) for c in cols) + " |")
    w("")
    for s in S:
        pa = S[s]["part1_paired_all_untagged"]
        w(f"{s}: mean strand share, rated vs unrated (all-unrated basis): " + "; ".join(
            f"{k} {pct(pa['strand_share_tagged_mean'][k])} vs {pct(pa['strand_share_untagged_mean'][k])}"
            for k in list(pa["strand_share_tagged_mean"])[:4]) + ".")
    w("")

    # ---------------- Table 2 ----------------
    w("## (2) Bias minimum N under difficulty reweighting (Design A)")
    w("")
    w("| Variant | S1 minN .70 [CI] | S1 minN .80 [CI] | S1 r(20) | S2 minN .70 [CI] | S2 minN .80 [CI] | S2 r(20) |")
    w("|---|---|---|---|---|---|---|")
    labels = [("random", "As submitted (random splits)"), ("w_diff_random", "Reweighted to unrated difficulty mix"),
              ("w_null_diff_random", "Null: same difficulty weights permuted within student"),
              ("w_diffstrand_random", "Reweighted to unrated difficulty x strand mix"),
              ("w_null_diffstrand_random", "Null: same difficulty x strand weights permuted within student"),
              ("balanced_random", "Difficulty-balanced subset (equal weights)"),
              ("balanced_ds_random", "Difficulty x strand balanced subset (equal weights)")]
    for v, lab in labels:
        cells = []
        for s in S:
            b = S[s]["part2_reweighting"]["variants"][v]["bias"]
            cells += [f"{f1(b['minN70_point'])} {ci1(b.get('minN70_ci'))}", f"{f1(b['minN80_point'])} {ci1(b.get('minN80_ci'))}",
                      f3(b["r20_point"])]
        w(f"| {lab} | " + " | ".join(cells) + " |")
    w("")
    w("Paired change in minN .70 against the submitted random-split curve (same students, same bootstrap draws): "
      + "; ".join(f"{s} {lab.lower()} {S[s]['part2_reweighting']['variants'][v]['bias']['minN70_minus_random']['point']:+.1f} "
                  f"{ci1(S[s]['part2_reweighting']['variants'][v]['bias']['minN70_minus_random']['ci_finite_draws'])}"
                  for s in S for v, lab in labels[1:]) + ".")
    w("")
    for s in S:
        rw = S[s]["part2_reweighting"]
        wd = [x for x in rw["weights_diagnostics_by_N"] if x["N"] == 20][0]
        tsm = rw["target_subset_minN"]
        ts20 = [x for x in rw["target_subset_curves"] if x["N"] == 20][0]
        w(f"{s}: {rw['weighting']['pool_students_with_target']} of {S[s]['n_poolA']:,} Design A students have >= 20 unrated "
          f"answers and so a weighting target (the rest keep weight 1). Effective share of responses kept by the weights at N = 20 "
          f"(Kish): difficulty {pct(wd['kish_ess_ratio_mean_target_students'])}, difficulty x strand "
          f"{pct(wd['diffstrand_kish_ess_ratio_mean_target_students'])} (students with a target). On those "
          f"{ts20['n_students']} students only: minN .70 raw {f1(tsm['random_bias']['minN70_point'])}, difficulty-weighted "
          f"{f1(tsm['w_diff_random_bias']['minN70_point'])} (null {f1(tsm['w_null_diff_random_bias']['minN70_point'])}), "
          f"difficulty x strand weighted {f1(tsm['w_diffstrand_random_bias']['minN70_point'])} "
          f"(null {f1(tsm['w_null_diffstrand_random_bias']['minN70_point'])}), "
          f"difficulty-balanced {f1(tsm['balanced_random_bias']['minN70_point'])}, difficulty x strand balanced "
          f"{f1(tsm['balanced_ds_random_bias']['minN70_point'])}.")
        w("")
    w("Bias level on all of a student's rated answers, raw vs reweighted to the unrated mix, with the "
      "permuted-weight null (same weights, no link to difficulty / strand):")
    w("")
    w("| Sample, pool, strata | Students | Mean raw bias | Shift [CI] | Null shift [CI] | r raw vs reweighted (null) | "
      "Top-decile Jaccard (null) |")
    w("|---|---|---|---|---|---|---|")
    for s in S:
        for pool, plab in [("ge5_tagged", ">= 5 rated"), ("ge20_tagged", ">= 20 rated"), ("designA_pool", ">= 200 rated")]:
            for wl, wlab in [("diff", "difficulty"), ("diffstrand", "difficulty x strand")]:
                b = S[s]["part2_reweighting"]["fullset_bias_shift"][pool][wl]
                bn = S[s]["part2_reweighting"]["fullset_bias_shift"][pool]["null_" + wl]
                w(f"| {s}, {plab}, {wlab} | {b['n_students']:,} | {f3(b['mean_raw_bias'])} | "
                  f"{s3(b['weighted_minus_raw']['mean'])} {ci3(b['weighted_minus_raw']['mean_ci95'])} | "
                  f"{s3(bn['weighted_minus_raw']['mean'])} {ci3(bn['weighted_minus_raw']['mean_ci95'])} | "
                  f"{f3(b['pearson_raw_weighted'])} ({f3(bn['pearson_raw_weighted'])}) | "
                  f"{f3(b['top_decile_jaccard_raw_vs_weighted'])} ({f3(bn['top_decile_jaccard_raw_vs_weighted'])}) |")
    w("")

    # ---------------- Table 3 ----------------
    w("## (3) Structured splits vs random and temporal splits (Design A)")
    w("")
    w("| Split | S1 bias r(20) | S1 bias minN .70 / .80 | S1 Brier minN .70 / .80 | S2 bias r(20) | S2 bias minN .70 / .80 | S2 Brier minN .70 / .80 |")
    w("|---|---|---|---|---|---|---|")
    for v, lab in [("random", "Random (submitted)"), ("temporal", "Temporal (submitted)"), ("strand", "By strand"),
                   ("topic", "By topic"), ("term", "By school term (approx.)"), ("difficulty", "Harder vs easier items")]:
        cells = []
        for s in S:
            vv = S[s]["part3_structured_splits"]["variants"][v]
            b, br = vv["bias"], vv["brier"]
            cells += [f"{f3(b['r20_point'])} {ci3(b['r20_ci'])}", f"{f1(b['minN70_point'])} / {f1(b['minN80_point'])}",
                      f"{f1(br['minN70_point'])} / {f1(br['minN80_point'])}"]
        w(f"| {lab} | " + " | ".join(cells) + " |")
    w("")
    w("Paired change against the random split (same students and bootstrap draws), bias: r(20) change [CI]; "
      "minN .70 change [CI]; minN .80 change [CI]:")
    w("")
    for s in S:
        H = S[s]["H"]
        w(f"- {s}: " + "; ".join(
            f"{v} {s3(H[f'bias_r20_change_{v}'])} {ci3(H[f'bias_r20_change_{v}_ci'])}, "
            f"{H[f'bias_minN70_change_{v}']:+.1f} {ci1(H[f'bias_minN70_change_{v}_ci'])}, "
            f"{H[f'bias_minN80_change_{v}']:+.1f} {ci1(H[f'bias_minN80_change_{v}_ci'])}"
            for v in ["temporal", "strand", "topic", "term", "difficulty"]) + ".")
    w("")
    w("Reliability at N = 20 of the components (random / strand / topic / term / difficulty split): "
      + "; ".join(f"{s} confidence " + " / ".join(f3(S[s]['part3_structured_splits']['variants'][v]['conf']['r20_point'])
                                                   for v in ['random', 'strand', 'topic', 'term', 'difficulty'])
                  + ", accuracy " + " / ".join(f3(S[s]['part3_structured_splits']['variants'][v]['acc']['r20_point'])
                                               for v in ['random', 'strand', 'topic', 'term', 'difficulty'])
                  for s in S) + ".")
    w("")
    w("How content-separated the halves are at N = 20 (share of responses whose category falls in both halves, "
      "keyed vs random split; share of students with a single category):")
    w("")
    for s in S:
        dg = S[s]["part3_structured_splits"]["diagnostics_by_N"]
        parts = []
        for k in ["strand", "topic", "term"]:
            x = [y for y in dg[k] if y["N"] == 20][0]
            parts.append(f"{k} {pct(x['shared_category_response_share_keyed'])} vs {pct(x['shared_category_response_share_random'])} "
                         f"(single category {pct(x['share_students_single_category'])})")
        dd = [y for y in dg["difficulty"] if y["N"] == 20][0]
        parts.append(f"difficulty halves mean p {f3(dd['mean_p_harder_half'])} vs {f3(dd['mean_p_easier_half'])}")
        w(f"- {s}: " + "; ".join(parts) + ".")
    w("")
    w("Students whose largest category holds <= 50% of their N responses (clearly separated halves), bias at N = 20, "
      "random / keyed / temporal on the same students:")
    w("")
    for s in S:
        rr = [x for x in S[s]["part3_structured_splits"]["restricted_largest_category_le_half"] if x["N"] == 20]
        w(f"- {s}: " + "; ".join(
            f"{x['key']} (n = {x['n_students']}{', descriptive, n < 300' if x['n_students'] < 300 else ''}) {f3(x.get('bias_random'))} / {f3(x.get('bias_keyed'))} / {f3(x.get('bias_temporal'))}"
            for x in rr) + ".")
    w("")
    w("## What could not be computed")
    w("")
    for t in d["not_computed"]:
        w(f"- {t}")
    w("")
    bullets, para = interpretation(S)
    w("## Key findings")
    w("")
    for b in bullets:
        w(f"- {b}")
    w("")
    w("## Interpretation for the manuscript (one paragraph)")
    w("")
    w(para)
    w("")
    txt = "\n".join(L)
    assert chr(0x2014) not in txt and chr(0x2013) not in txt
    with open(MD, "w", encoding="utf-8") as fh:
        fh.write(txt)
    print("written", MD)


def s3(v):
    """signed 3-decimal without leading zero"""
    t = f3(abs(v))
    return ("+" if v >= 0 else "-") + t


def rng_(vals, fmt):
    lo, hi = min(vals), max(vals)
    return f"{fmt(lo)} to {fmt(hi)}"


def interpretation(S):
    H1, H2 = S["S1"]["H"], S["S2"]["H"]
    D1, D2 = S["S1"]["part4_sample_descriptives"], S["S2"]["part4_sample_descriptives"]
    keyed = ["strand", "topic", "term", "difficulty"]
    r1 = S["S1"]["part3_structured_splits"]["restricted_largest_category_le_half"]
    rt1 = [x for x in r1 if x["N"] == 20 and x["key"] == "topic"][0]
    chg = [H[f"bias_minN70_change_{v}"] for H in (H1, H2) for v in keyed]
    chg += [H[f"bias_minN70_change_{v}"] for H in (H1, H2)
            for v in ["w_diff_random", "balanced_random", "balanced_ds_random"]]
    assert all(c < min(H1["bias_minN70_change_temporal"], H2["bias_minN70_change_temporal"]) for c in chg)
    fs = "fullset_ge20_tagged_"
    p = []
    p.append(
        "Among reporters' answers, whether an answer carries a rating varies mostly between students (eta-squared "
        f"{f3(H1['reporters_eta2_tag_student'])} in S1, {f3(H2['reporters_eta2_tag_student'])} in S2) and between a "
        f"student's quiz attempts (student x quiz {f3(H1['reporters_eta2_tag_student_x_quiz'])} and "
        f"{f3(H2['reporters_eta2_tag_student_x_quiz'])}); about {pct(H1['reporters_within_student_x_quiz_share'])} and "
        f"{pct(H2['reporters_within_student_x_quiz_share'])} of the variance lies within an attempt, and calendar month "
        f"accounts for {f3(H1['reporters_eta2_tag_month'])} and {f3(H2['reporters_eta2_tag_month'])}. Of reporters' student x "
        f"quiz pairs, {pct(H1['reporters_share_student_quiz_pairs_mixed'])} (S1) and "
        f"{pct(H2['reporters_share_student_quiz_pairs_mixed'])} (S2) mix rated and unrated answers. Rating also depends on "
        f"the period: in S1 no answer is rated before {D1['first_tagged'][:7]}, and in S2 the rated share falls from "
        f"{pct(D2['tag_share_by_month']['2019-09'], 0)} (2019-09) to {pct(D2['tag_share_by_month']['2020-04'], 0)} (2020-04). "
        f"Within students, rated items are only slightly easier than unrated ones ({pts(H1['all_untagged_item_p_gap'])} "
        f"points of item p-value in S1, {pts(H2['all_untagged_item_p_gap'])} in S2), so most of the known accuracy gap "
        f"({pts(H1['all_untagged_accuracy_gap'])} and {pts(H2['all_untagged_accuracy_gap'])} points) is accuracy beyond item "
        f"difficulty ({pts(H1['all_untagged_beyond_difficulty_gap'])} and {pts(H2['all_untagged_beyond_difficulty_gap'])}); the "
        "submitted wording that students rate their easier answers therefore needs narrowing. Within the same quiz "
        f"attempt, rated answers sit earlier in the quiz (relative position {s3(H1['same_quiz_relpos_gap'])} and "
        f"{s3(H2['same_quiz_relpos_gap'])}) and on slightly easier items ({pts(H1['same_quiz_item_p_gap'])} and "
        f"{pts(H2['same_quiz_item_p_gap'])} points); whether this reflects student choice or the platform's prompting is "
        "not recorded. The topic mix of rated and unrated answers differs beyond chance (topic TVD "
        f"{f3(H1['all_untagged_tvd_topic'])} vs {f3(H1['all_untagged_tvd_topic_permutation_baseline'])} expected from sampling "
        f"in S1; {f3(H2['all_untagged_tvd_topic'])} vs {f3(H2['all_untagged_tvd_topic_permutation_baseline'])} in S2).")
    p.append(
        "These selections move the precision boundary of the bias score little, and what moves it is unequal "
        "weighting, not composition. Reweighting each student's rated answers to the difficulty mix of their unrated "
        f"answers moves the .70 crossing from {f1(H1['bias_minN70_random'])} to {f1(H1['bias_minN70_w_diff_random'])} "
        f"responses in S1 (paired change {H1['bias_minN70_change_w_diff_random']:+.1f} "
        f"{ci1(H1['bias_minN70_change_w_diff_random_ci'])}) and from {f1(H2['bias_minN70_random'])} to "
        f"{f1(H2['bias_minN70_w_diff_random'])} in S2; the same weights permuted at random within student (same "
        f"effective sample size, no link to difficulty) give {f1(H1['bias_minN70_w_null_diff_random'])} and "
        f"{f1(H2['bias_minN70_w_null_diff_random'])}. The difficulty x strand weights raise the crossing to "
        f"{f1(H1['bias_minN70_w_diffstrand_random'])} and {f1(H2['bias_minN70_w_diffstrand_random'])}, but their permuted "
        f"null does the same ({f1(H1['bias_minN70_w_null_diffstrand_random'])} and "
        f"{f1(H2['bias_minN70_w_null_diffstrand_random'])}): the rise is the cost of very unequal weights (effective sample "
        f"size {pct(H1['kish_ratio_N20_w_diffstrand_target_students'])} and {pct(H2['kish_ratio_N20_w_diffstrand_target_students'])} "
        "of N at N = 20, students with a target). Equal-weight subsets drawn toward the same mixes leave the crossing "
        f"within bootstrap error of the submitted value ({f1(H1['bias_minN70_balanced_random'])} in S1, "
        f"{f1(H2['bias_minN70_balanced_random'])} in S2; difficulty x strand: {f1(H1['bias_minN70_balanced_ds_random'])} and "
        f"{f1(H2['bias_minN70_balanced_ds_random'])}); these subsets move the difficulty mix only partway toward the unrated "
        f"mix (TVD to the unrated mix at N = 20 from {f3(H1['N20_tvd_to_unrated_difficulty_mix_plain'])} to "
        f"{f3(H1['N20_tvd_to_unrated_difficulty_mix_balanced'])} in S1, {f3(H2['N20_tvd_to_unrated_difficulty_mix_plain'])} to "
        f"{f3(H2['N20_tvd_to_unrated_difficulty_mix_balanced'])} in S2, part of it sampling noise).")
    p.append(
        "Content, difficulty and period heterogeneity cost some reliability, but less than the temporal split already "
        "reported. Splitting by whole strands, topics, approximate school terms or harder vs easier items instead of at "
        f"random changes bias reliability at N = 20 by {rng_([H1[f'bias_r20_change_{v}'] for v in keyed], s3)} in S1 and "
        f"{rng_([H2[f'bias_r20_change_{v}'] for v in keyed], s3)} in S2 (temporal: {s3(H1['bias_r20_change_temporal'])} and "
        f"{s3(H2['bias_r20_change_temporal'])}); the .70 crossing moves to {rng_([H1[f'bias_minN70_{v}'] for v in keyed], f1)} "
        f"responses in S1 and {rng_([H2[f'bias_minN70_{v}'] for v in keyed], f1)} in S2, against {f1(H1['bias_minN70_temporal'])} "
        f"and {f1(H2['bias_minN70_temporal'])} under the temporal split, and the .80 crossing to "
        f"{rng_([H1[f'bias_minN80_{v}'] for v in keyed], f1)} (S1) and {rng_([H2[f'bias_minN80_{v}'] for v in keyed], f1)} (S2), "
        f"against {f1(H1['bias_minN80_temporal'])} and {f1(H2['bias_minN80_temporal'])}. Among S1 students whose halves are "
        f"clearly topic-separated (n = {rt1['n_students']:,}), bias reliability at N = 20 is {f3(rt1['bias_random'])} random, "
        f"{f3(rt1['bias_keyed'])} by topic and {f3(rt1['bias_temporal'])} temporal. Mean confidence stays the most stable "
        f"component across all splits (r(20) {rng_([H1[f'conf_r20_{v}'] for v in ['random'] + keyed], f3)} in S1). The Brier "
        "score is more sensitive to item difficulty: under harder vs easier halves its .80 crossing is "
        f"{f1(H1['brier_minN80_difficulty'])} responses in S1 (random: {f1(H1['brier_minN80_random'])}) and it does not reach .70 "
        "by 200 responses in S2.")
    p.append(
        "What the observable correction does change is the level of the score and, to a lesser degree, who is "
        "flagged. On all of a student's rated answers (students with >= 20 rated answers and a weighting target), "
        f"reweighting to the unrated difficulty mix raises mean bias by {s3(H1[fs + 'diff_bias_shift'])} in S1 and "
        f"{s3(H2[fs + 'diff_bias_shift'])} in S2 on the 0-1 scale, against {s3(H1[fs + 'null_diff_bias_shift'])} and "
        f"{s3(H2[fs + 'null_diff_bias_shift'])} under the permuted-weight null; the reweighted top decile overlaps the raw "
        f"top decile by Jaccard {f3(H1[fs + 'diff_top_decile_jaccard'])} and {f3(H2[fs + 'diff_top_decile_jaccard'])}, against "
        f"{f3(H1[fs + 'null_diff_top_decile_jaccard'])} and {f3(H2[fs + 'null_diff_top_decile_jaccard'])} under the null.")
    p.append(
        "Sample 2 differs from Sample 1 in ways that bear on generalization: a far smaller item pool "
        f"({D2['n_questions']:,} vs {D1['n_questions']:,} questions, {D2['n_strand_categories']} strands), harder items "
        f"(median item p {f3(D2['item_p_question_level']['median'])} vs {f3(D1['item_p_question_level']['median'])}), younger "
        f"students (median age at first answer {D2['demographics_all_students']['median_age_years_at_first_answer']:.1f} vs "
        f"{D1['demographics_all_students']['median_age_years_at_first_answer']:.1f} years, where birth date is recorded), far "
        f"more common rating ({pct(D2['share_answers_tagged'])} vs {pct(D1['share_answers_tagged'])} of answers; "
        f"{pct(D2['share_students_reporting'])} vs {pct(D1['share_students_reporting'])} of students) and stronger item-level "
        f"selection (accuracy gap {pts(H2['all_untagged_accuracy_gap'])} vs {pts(H1['all_untagged_accuracy_gap'])} points).")
    para = (
        "Optional rating selects answers mainly by student and quiz attempt, and only modestly by item difficulty: "
        f"among reporters' answers, {pct(H1['reporters_eta2_tag_student_x_quiz'])} (S1) and "
        f"{pct(H2['reporters_eta2_tag_student_x_quiz'])} (S2) of the variance in whether an answer is rated lies between "
        f"student x quiz pairs (about {pct(H1['reporters_within_student_x_quiz_share'])} and "
        f"{pct(H2['reporters_within_student_x_quiz_share'])} within an attempt), rated items are only slightly easier "
        f"({pts(H1['all_untagged_item_p_gap'])} and {pts(H2['all_untagged_item_p_gap'])} points of item p-value), and most "
        "of the accuracy advantage of rated answers lies beyond item difficulty. The minimum-response estimates for "
        "bias are nonetheless stable under every observable correction tried here: difficulty reweighting, "
        "difficulty-balanced subsets, and halves formed by whole strands, topics, school terms or by item difficulty "
        f"shift the .70 crossing by {rng_(chg, lambda v: f'{v:+.1f}')} responses, less than the temporal split already "
        f"reported ({H1['bias_minN70_change_temporal']:+.1f} in S1, {H2['bias_minN70_change_temporal']:+.1f} in S2); the one "
        "larger shift, under difficulty x strand weights, is reproduced by the same weights permuted at random within "
        "student and so reflects unequal weighting, not composition. What more responses cannot buy is "
        "representativeness: reweighting to the unrated difficulty mix raises mean bias "
        f"({s3(H1[fs + 'diff_bias_shift'])} and {s3(H2[fs + 'diff_bias_shift'])} on the 0-1 scale; permuted-weight null "
        f"{s3(H1[fs + 'null_diff_bias_shift'])} and {s3(H2[fs + 'null_diff_bias_shift'])}) and changes part of the top-decile "
        f"flag set beyond the null (Jaccard {f3(H1[fs + 'diff_top_decile_jaccard'])} vs {f3(H1[fs + 'null_diff_top_decile_jaccard'])} "
        f"in S1, {f3(H2[fs + 'diff_top_decile_jaccard'])} vs {f3(H2[fs + 'null_diff_top_decile_jaccard'])} in S2; students with "
        ">= 20 rated answers), and calibration on unrated answers, including any tendency to rate only when sure, is "
        "not observable. The thresholds should therefore be stated as conditional on rated answers, in quiz attempts "
        "and months where rating was in use, on this platform, and the agreement of the two samples (bias .70 crossing "
        f"{f1(H1['bias_minN70_random'])} and {f1(H2['bias_minN70_random'])}) read as stability across two reporting regimes "
        "and item pools of one platform, not across populations.")
    return p, para


if __name__ == "__main__":
    main()
