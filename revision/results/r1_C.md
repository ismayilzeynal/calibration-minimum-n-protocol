# Stream C - selection and content heterogeneity (revision analyses)

Machine-readable results with definitions: `results/r1_C.json`. Scripts: `scripts/r1_C_common.py`, `r1_C_selection.py`, `r1_C_splits.py`, `r1_C_assemble.py`, `r1_C_report.py` (this file is generated; no number is typed by hand).

Runtime (single process each): S1 selection 31 s, S1 splits 205 s, S2 selection 4 s, S2 splits 52 s; total 293 s.

## What was computed, and how

1. **Who rates which answers.** Within each student, tagged answers were compared with untagged answers on item difficulty (leave-one-out item p-value), accuracy, accuracy beyond item difficulty (IsCorrect - p), strand and topic mix (total variation distance, TVD, against a within-student permutation baseline) and timing (days since the student's first answer; weekend share; position within the quiz). Two bases: all untagged answers (students with >= 5 of each, the basis of the submitted +4.0-point result), and the same quiz attempts only (student x quiz pairs that mix rated and unrated answers). Also: where the variation in tagging sits (eta-squared).
2. **Difficulty reweighting sensitivity of the bias minimum N.** On the Design A pool and the identical subsamples and splits of the submitted analysis, each tagged answer was reweighted so that the student's difficulty-stratum mix (and, as a second variant, difficulty x strand mix) matches that student's untagged answers; an equal-weight alternative draws difficulty-balanced subsets instead, and a permuted-weight null (each student's own weights shuffled across their rated answers: same effective sample size, no link to difficulty) separates composition effects from the cost of unequal weights. Same joint bootstrap draws.
3. **Structured splits.** Instead of random halves, halves were formed by whole strands, whole topics, whole (approximate) school terms, or harder vs easier items, on the same subsamples; compared with the random and temporal splits already reported.
4. **Sample 1 vs Sample 2 table.**

Regression check: the regenerated random and temporal curves match the baseline reproduction to at most 1.1e-16 (summation order); the random-split bootstrap draws match the baseline cache exactly (0).

## (4) Sample 1 vs Sample 2

| Quantity | Sample 1 | Sample 2 |
|---|---|---|
| Students; answers | 118,971; 15,867,850 | 4,918; 1,382,727 |
| Questions; quizzes; groups | 27,613; 17,227; 11,838 | 948; 120; 370 |
| Answers with a confidence rating | 1,506,052 (9.5%) | 346,428 (25.1%) |
| Students with any rating (reporters) | 39,842 (33.5%) | 4,284 (87.1%) |
| Median (mean) rated answers per reporter | 15 (37.8) | 31 (80.9) |
| Reporters with >= 20 / >= 200 rated | 17,156 / 1,389 | 2,617 / 554 |
| Median answers per student | 88 | 239 |
| Accuracy: all / rated / unrated answers | 64.3% / 67.6% / 63.9% | 53.7% / 62.2% / 50.9% |
| Mean confidence (rated); share at 100 | 75.0; 48.9% | 73.9; 46.0% |
| Answer window | 2018-09 to 2020-05 | 2019-09 to 2020-05 |
| Rated-answer window | 2019-01-24 to 2020-04-29 | 2019-09-06 to 2020-04-29 |
| Item p-value, question level: median [IQR] | .685 [.558, .797] | .495 [.387, .618] |
| Strand labels; topic labels (incl. combined labels for multi-subject questions) | 29; 354 | 3; 26 |
| Share of answers: Number (all / rated) | 46.5% / 53.8% | 49.6% / 60.3% |
| Share of answers: Algebra (all / rated) | 21.8% / 19.7% | 20.3% / 16.5% |
| Share of answers: Geometry and Measure (all / rated) | 21.4% / 18.5% | 30.1% / 23.2% |
| Share of answers: Data and Statistics (all / rated) | 8.6% / 6.8% | 0 (not in pool) |
| Median age at first answer, years [IQR] (birth date known) | 13.4 [11.9, 14.8] (68% known) | 11.7 [11.4, 12.0] (47% known) |
| PremiumPupil = 1 among known (share known) | 25.8% (29%) | 28.4% (43%) |
| Reporters' answers: variance of the rating indicator between students / between student x quiz pairs / between months (eta-squared) | .501 / .840 / .069 | .540 / .829 / .052 |
| Reporters' student x quiz pairs with some but not all answers rated | 15.7% | 25.2% |

Month-level rating share (share of answers rated): S1 2018-12 0%, 2019-01 2%, 2019-02 14%, 2019-09 16%, 2020-04 13%; S2 2019-09 59%, 2019-10 40%, 2020-01 20%, 2020-04 16%. Gender is stored as raw codes with unverified labels and is not tabulated here.

## (1) Rated vs unrated answers within student

Paired differences, rated minus unrated, mean over students [student-bootstrap 95% CI]; percentage points unless stated.

| Quantity | S1 all unrated | S1 same quiz | S2 all unrated | S2 same quiz |
|---|---|---|---|---|
| Students | 27,625 | 20,748 | 3,539 | 3,315 |
| Accuracy | +4.0 [+3.8, +4.2] | +4.1 [+3.8, +4.4] | +7.8 [+7.2, +8.3] | +4.8 [+4.2, +5.4] |
| Item p-value (leave-one-out; + = rated items easier) | +1.2 [+1.1, +1.2] | +2.0 [+1.9, +2.1] | +3.6 [+3.4, +3.8] | +2.5 [+2.3, +2.7] |
| Accuracy beyond item difficulty | +2.8 [+2.6, +3.0] | +2.1 [+1.8, +2.3] | +4.2 [+3.7, +4.7] | +2.3 [+1.8, +2.8] |
| Days since the student's first answer (days) | -2.6 [-3.4, -1.7] | -9.3 [-9.8, -8.9] | -38.4 [-39.7, -37.1] | -15.2 [-16.0, -14.3] |
| Weekend share | -0.3 [-0.5, +0.0] | +0.2 [-0.0, +0.4] | +1.3 [+0.4, +2.1] | +0.5 [-0.1, +1.1] |
| Relative position within the quiz (0 = first, 1 = last) | n/a | -.122 [-.124, -.120] | n/a | -.110 [-.115, -.106] |
| Strand mix: TVD observed / permutation baseline | .337 / .151 | .212 / .150 | .343 / .110 | .199 / .093 |
| Topic mix: TVD observed / permutation baseline | .659 / .391 | .447 / .353 | .589 / .272 | .382 / .223 |
| Share of students with rated items easier (p difference > 0) | 59.1% | 62.8% | 76.3% | 67.5% |

S1: mean strand share, rated vs unrated (all-unrated basis): Number 50.0% vs 49.2%; Algebra 21.4% vs 21.5%; Geometry and Measure 19.5% vs 20.2%; Data and Statistics 7.8% vs 7.7%.
S2: mean strand share, rated vs unrated (all-unrated basis): Number 66.9% vs 48.0%; Geometry and Measure 21.3% vs 30.6%; Algebra 11.8% vs 21.4%.

## (2) Bias minimum N under difficulty reweighting (Design A)

| Variant | S1 minN .70 [CI] | S1 minN .80 [CI] | S1 r(20) | S2 minN .70 [CI] | S2 minN .80 [CI] | S2 r(20) |
|---|---|---|---|---|---|---|
| As submitted (random splits) | 14.7 [13.1, 16.9] | 26.0 [23.4, 28.4] | .754 | 16.2 [14.4, 19.2] | 28.7 [24.2, 33.8] | .742 |
| Reweighted to unrated difficulty mix | 15.9 [13.9, 18.1] | 27.2 [24.6, 30.3] | .744 | 17.8 [14.9, 21.2] | 32.4 [27.4, 36.8] | .727 |
| Null: same difficulty weights permuted within student | 15.5 [13.9, 17.9] | 27.5 [25.3, 30.0] | .745 | 17.6 [14.8, 21.0] | 31.2 [26.9, 36.6] | .727 |
| Reweighted to unrated difficulty x strand mix | 24.6 [20.1, 27.1] | 43.8 [37.0, 50.6] | .667 | 34.6 [24.8, 39.1] | 65.0 [52.2, 77.4] | .623 |
| Null: same difficulty x strand weights permuted within student | 25.3 [19.7, 27.7] | 42.9 [36.8, 47.1] | .669 | 33.8 [24.6, 39.1] | 61.4 [48.5, 73.8] | .616 |
| Difficulty-balanced subset (equal weights) | 13.9 [12.5, 15.6] | 24.3 [22.0, 27.6] | .764 | 15.0 [13.2, 17.9] | 27.6 [23.4, 32.4] | .747 |
| Difficulty x strand balanced subset (equal weights) | 14.7 [12.9, 16.6] | 25.2 [22.3, 29.0] | .764 | 16.7 [13.4, 19.1] | 27.3 [24.0, 31.2] | .741 |

Paired change in minN .70 against the submitted random-split curve (same students, same bootstrap draws): S1 reweighted to unrated difficulty mix +1.3 [0.5, 2.0]; S1 null: same difficulty weights permuted within student +0.8 [0.4, 1.5]; S1 reweighted to unrated difficulty x strand mix +9.9 [6.4, 11.9]; S1 null: same difficulty x strand weights permuted within student +10.7 [5.8, 12.3]; S1 difficulty-balanced subset (equal weights) -0.8 [-2.3, 0.4]; S1 difficulty x strand balanced subset (equal weights) +0.0 [-1.6, 1.5]; S2 reweighted to unrated difficulty mix +1.5 [0.5, 3.0]; S2 null: same difficulty weights permuted within student +1.3 [0.3, 3.0]; S2 reweighted to unrated difficulty x strand mix +18.4 [9.6, 22.7]; S2 null: same difficulty x strand weights permuted within student +17.6 [9.4, 22.1]; S2 difficulty-balanced subset (equal weights) -1.2 [-3.5, 1.2]; S2 difficulty x strand balanced subset (equal weights) +0.4 [-2.8, 2.7].

S1: 849 of 1,389 Design A students have >= 20 unrated answers and so a weighting target (the rest keep weight 1). Effective share of responses kept by the weights at N = 20 (Kish): difficulty 88.7%, difficulty x strand 61.2% (students with a target). On those 849 students only: minN .70 raw 14.6, difficulty-weighted 16.6 (null 15.8), difficulty x strand weighted 27.5 (null 28.0), difficulty-balanced 13.2, difficulty x strand balanced 14.4.

S2: 352 of 554 Design A students have >= 20 unrated answers and so a weighting target (the rest keep weight 1). Effective share of responses kept by the weights at N = 20 (Kish): difficulty 86.8%, difficulty x strand 49.9% (students with a target). On those 352 students only: minN .70 raw 15.7, difficulty-weighted 18.2 (null 17.7), difficulty x strand weighted 38.6 (null 38.8), difficulty-balanced 14.3, difficulty x strand balanced 15.6.

Bias level on all of a student's rated answers, raw vs reweighted to the unrated mix, with the permuted-weight null (same weights, no link to difficulty / strand):

| Sample, pool, strata | Students | Mean raw bias | Shift [CI] | Null shift [CI] | r raw vs reweighted (null) | Top-decile Jaccard (null) |
|---|---|---|---|---|---|---|
| S1, >= 5 rated, difficulty | 24,642 | .063 | +.005 [.004, .006] | +.000 [-.001, .001] | .935 (.945) | .600 (.646) |
| S1, >= 5 rated, difficulty x strand | 24,549 | .063 | +.002 [.001, .004] | +.001 [-.001, .002] | .887 (.900) | .526 (.553) |
| S1, >= 20 rated, difficulty | 11,929 | .069 | +.006 [.005, .007] | +.001 [.000, .002] | .961 (.974) | .660 (.729) |
| S1, >= 20 rated, difficulty x strand | 11,909 | .069 | +.004 [.002, .005] | +.001 [-.000, .003] | .905 (.925) | .529 (.581) |
| S1, >= 200 rated, difficulty | 849 | .089 | +.005 [.003, .007] | +.000 [-.000, .001] | .986 (.998) | .809 (.954) |
| S1, >= 200 rated, difficulty x strand | 849 | .089 | +.005 [.001, .008] | +.002 [-.001, .004] | .959 (.977) | .700 (.735) |
| S2, >= 5 rated, difficulty | 3,230 | .100 | +.019 [.017, .022] | -.001 [-.003, .001] | .946 (.959) | .603 (.691) |
| S2, >= 5 rated, difficulty x strand | 3,228 | .100 | +.017 [.013, .020] | +.001 [-.002, .004] | .899 (.924) | .499 (.603) |
| S2, >= 20 rated, difficulty | 2,154 | .106 | +.021 [.019, .023] | +.001 [-.001, .003] | .962 (.975) | .648 (.784) |
| S2, >= 20 rated, difficulty x strand | 2,152 | .106 | +.023 [.019, .026] | +.001 [-.002, .004] | .902 (.937) | .488 (.629) |
| S2, >= 200 rated, difficulty | 352 | .130 | +.019 [.014, .023] | +.001 [-.000, .002] | .979 (.998) | .628 (.842) |
| S2, >= 200 rated, difficulty x strand | 352 | .130 | +.023 [.015, .031] | -.001 [-.005, .003] | .929 (.975) | .458 (.795) |

## (3) Structured splits vs random and temporal splits (Design A)

| Split | S1 bias r(20) | S1 bias minN .70 / .80 | S1 Brier minN .70 / .80 | S2 bias r(20) | S2 bias minN .70 / .80 | S2 Brier minN .70 / .80 |
|---|---|---|---|---|---|---|
| Random (submitted) | .754 [.728, .777] | 14.7 / 26.0 | 21.3 / 34.3 | .742 [.707, .771] | 16.2 / 28.7 | 31.1 / 58.8 |
| Temporal (submitted) | .698 [.654, .736] | 20.4 / 39.4 | 30.4 / 56.0 | .697 [.640, .745] | 20.4 / 48.1 | 48.1 / 95.1 |
| By strand | .748 [.719, .775] | 15.5 / 29.0 | 23.1 / 39.7 | .725 [.687, .760] | 18.2 / 35.3 | 36.6 / 81.3 |
| By topic | .739 [.710, .764] | 16.1 / 28.8 | 23.8 / 40.6 | .724 [.685, .757] | 17.7 / 32.8 | 33.2 / 67.8 |
| By school term (approx.) | .727 [.693, .756] | 17.4 / 31.5 | 26.2 / 42.7 | .719 [.676, .754] | 18.3 / 35.9 | 34.9 / 76.5 |
| Harder vs easier items | .736 [.702, .765] | 14.6 / 35.2 | 33.1 / 103.7 | .721 [.670, .766] | 18.4 / 35.3 | > 200 / > 200 |

Paired change against the random split (same students and bootstrap draws), bias: r(20) change [CI]; minN .70 change [CI]; minN .80 change [CI]:

- S1: temporal -.056 [-.084, -.030], +5.8 [0.7, 11.1], +13.4 [10.2, 25.7]; strand -.006 [-.017, .005], +0.9 [-0.2, 1.8], +3.0 [1.8, 5.4]; topic -.014 [-.022, -.007], +1.5 [0.4, 2.5], +2.8 [1.8, 5.1]; term -.026 [-.042, -.013], +2.8 [1.2, 4.8], +5.5 [3.1, 9.9]; difficulty -.018 [-.038, .002], -0.1 [-1.4, 2.2], +9.2 [3.6, 15.5].
- S2: temporal -.045 [-.084, -.007], +4.2 [0.2, 8.8], +19.4 [7.8, 36.6]; strand -.017 [-.033, -.001], +1.9 [0.6, 3.6], +6.6 [2.8, 11.1]; topic -.018 [-.031, -.004], +1.5 [0.2, 3.4], +4.1 [1.4, 7.9]; term -.023 [-.042, -.008], +2.0 [0.4, 4.6], +7.2 [2.7, 12.6]; difficulty -.021 [-.052, .014], +2.2 [-0.0, 6.5], +6.6 [-0.4, 22.6].

Reliability at N = 20 of the components (random / strand / topic / term / difficulty split): S1 confidence .929 / .913 / .916 / .894 / .909, accuracy .777 / .758 / .761 / .764 / .747; S2 confidence .926 / .906 / .913 / .901 / .920, accuracy .710 / .697 / .694 / .695 / .711.

How content-separated the halves are at N = 20 (share of responses whose category falls in both halves, keyed vs random split; share of students with a single category):

- S1: strand 59.0% vs 93.4% (single category 3.0%); topic 17.3% vs 59.1% (single category 0.0%); term 58.2% vs 97.0% (single category 6.3%); difficulty halves mean p .537 vs .753.
- S2: strand 54.2% vs 96.0% (single category 0.5%); topic 16.8% vs 68.6% (single category 0.0%); term 70.7% vs 98.4% (single category 19.3%); difficulty halves mean p .428 vs .663.

Students whose largest category holds <= 50% of their N responses (clearly separated halves), bias at N = 20, random / keyed / temporal on the same students:

- S1: strand (n = 454) .736 / .716 / .696; topic (n = 1255) .757 / .742 / .704; term (n = 341) .780 / .769 / .745.
- S2: strand (n = 213, descriptive, n < 300) .782 / .757 / .751; topic (n = 543) .745 / .727 / .698; term (n = 61, descriptive, n < 300) .792 / .766 / .779.

## What could not be computed

- Gamma and AUROC under structured splits and reweighting: not computed in this stream, because the AUROC constant-confidence rule is being changed in stream A, and both rank indices are below .30 at N = 20 under random splits in the submitted analysis (headline: baseline_r20_gamma, baseline_r20_auroc).
- Selection on unobservables (for example, a student skipping the rating when unsure of the answer) cannot be separated from difficulty or content in these data: confidence on unrated answers is not observed.
- Monitoring on untagged answers is unobservable; the reweighting assumes that, within a difficulty (x strand) stratum, a student's tagged answers stand for the untagged ones.
- Whether a quiz offered the rating prompt to a student is not recorded; the same-quiz basis uses quiz attempts in which the student rated at least one answer as the evidence of availability.
- Gender codes are reported raw; their labels are not verified (unknown).
- Term split uses calendar-month approximations of English school terms; actual term dates are unknown.

## Key findings

- Among reporters' answers, whether an answer carries a rating varies mostly between students (eta-squared .501 in S1, .540 in S2) and between a student's quiz attempts (student x quiz .840 and .829); about 16.0% and 17.1% of the variance lies within an attempt, and calendar month accounts for .069 and .052. Of reporters' student x quiz pairs, 15.7% (S1) and 25.2% (S2) mix rated and unrated answers. Rating also depends on the period: in S1 no answer is rated before 2019-01, and in S2 the rated share falls from 59% (2019-09) to 16% (2020-04). Within students, rated items are only slightly easier than unrated ones (+1.2 points of item p-value in S1, +3.6 in S2), so most of the known accuracy gap (+4.0 and +7.8 points) is accuracy beyond item difficulty (+2.8 and +4.2); the submitted wording that students rate their easier answers therefore needs narrowing. Within the same quiz attempt, rated answers sit earlier in the quiz (relative position -.122 and -.110) and on slightly easier items (+2.0 and +2.5 points); whether this reflects student choice or the platform's prompting is not recorded. The topic mix of rated and unrated answers differs beyond chance (topic TVD .659 vs .391 expected from sampling in S1; .589 vs .272 in S2).
- These selections move the precision boundary of the bias score little, and what moves it is unequal weighting, not composition. Reweighting each student's rated answers to the difficulty mix of their unrated answers moves the .70 crossing from 14.7 to 15.9 responses in S1 (paired change +1.3 [0.5, 2.0]) and from 16.2 to 17.8 in S2; the same weights permuted at random within student (same effective sample size, no link to difficulty) give 15.5 and 17.6. The difficulty x strand weights raise the crossing to 24.6 and 34.6, but their permuted null does the same (25.3 and 33.8): the rise is the cost of very unequal weights (effective sample size 61.2% and 49.9% of N at N = 20, students with a target). Equal-weight subsets drawn toward the same mixes leave the crossing within bootstrap error of the submitted value (13.9 in S1, 15.0 in S2; difficulty x strand: 14.7 and 16.7); these subsets move the difficulty mix only partway toward the unrated mix (TVD to the unrated mix at N = 20 from .226 to .166 in S1, .238 to .169 in S2, part of it sampling noise).
- Content, difficulty and period heterogeneity cost some reliability, but less than the temporal split already reported. Splitting by whole strands, topics, approximate school terms or harder vs easier items instead of at random changes bias reliability at N = 20 by -.026 to -.006 in S1 and -.023 to -.017 in S2 (temporal: -.056 and -.045); the .70 crossing moves to 14.6 to 17.4 responses in S1 and 17.7 to 18.4 in S2, against 20.4 and 20.4 under the temporal split, and the .80 crossing to 28.8 to 35.2 (S1) and 32.8 to 35.9 (S2), against 39.4 and 48.1. Among S1 students whose halves are clearly topic-separated (n = 1,255), bias reliability at N = 20 is .757 random, .742 by topic and .704 temporal. Mean confidence stays the most stable component across all splits (r(20) .894 to .929 in S1). The Brier score is more sensitive to item difficulty: under harder vs easier halves its .80 crossing is 103.7 responses in S1 (random: 34.3) and it does not reach .70 by 200 responses in S2.
- What the observable correction does change is the level of the score and, to a lesser degree, who is flagged. On all of a student's rated answers (students with >= 20 rated answers and a weighting target), reweighting to the unrated difficulty mix raises mean bias by +.006 in S1 and +.021 in S2 on the 0-1 scale, against +.001 and +.001 under the permuted-weight null; the reweighted top decile overlaps the raw top decile by Jaccard .660 and .648, against .729 and .784 under the null.
- Sample 2 differs from Sample 1 in ways that bear on generalization: a far smaller item pool (948 vs 27,613 questions, 3 strands), harder items (median item p .495 vs .685), younger students (median age at first answer 11.7 vs 13.4 years, where birth date is recorded), far more common rating (25.1% vs 9.5% of answers; 87.1% vs 33.5% of students) and stronger item-level selection (accuracy gap +7.8 vs +4.0 points).

## Interpretation for the manuscript (one paragraph)

Optional rating selects answers mainly by student and quiz attempt, and only modestly by item difficulty: among reporters' answers, 84.0% (S1) and 82.9% (S2) of the variance in whether an answer is rated lies between student x quiz pairs (about 16.0% and 17.1% within an attempt), rated items are only slightly easier (+1.2 and +3.6 points of item p-value), and most of the accuracy advantage of rated answers lies beyond item difficulty. The minimum-response estimates for bias are nonetheless stable under every observable correction tried here: difficulty reweighting, difficulty-balanced subsets, and halves formed by whole strands, topics, school terms or by item difficulty shift the .70 crossing by -1.2 to +2.8 responses, less than the temporal split already reported (+5.8 in S1, +4.2 in S2); the one larger shift, under difficulty x strand weights, is reproduced by the same weights permuted at random within student and so reflects unequal weighting, not composition. What more responses cannot buy is representativeness: reweighting to the unrated difficulty mix raises mean bias (+.006 and +.021 on the 0-1 scale; permuted-weight null +.001 and +.001) and changes part of the top-decile flag set beyond the null (Jaccard .660 vs .729 in S1, .648 vs .784 in S2; students with >= 20 rated answers), and calibration on unrated answers, including any tendency to rate only when sure, is not observable. The thresholds should therefore be stated as conditional on rated answers, in quiz attempts and months where rating was in use, on this platform, and the agreement of the two samples (bias .70 crossing 14.7 and 16.2) read as stability across two reporting regimes and item pools of one platform, not across populations.
