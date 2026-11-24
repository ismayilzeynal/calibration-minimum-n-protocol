Pre-specified grid: S1 bias .70 crossing 15.4 on the 6-point grid against 14.7 on the full grid (S1 median volume 15: above 15 on the 6-point grid True, on the full grid False). S2 Brier .70 crossing 35.9 on the 6-point grid against 31.1 on the full grid (S2 median 31); criterion (b) met on both grids in the printed draw.

# Crossings on the pre-specified grid, alpha-type checks, AUROC diagnostics

## 1. Design A crossings on the full grid and on the pre-specified grid

The analysis plan fixed N = 5, 10, 20, 50, 100 and 200; addendum 1 (2026-07-09, before submission) added 12, 15, 25, 30 and 40 to Design A as an interpolation-sensitivity check. The printed crossings use all eleven points. Below, the same curves and the same crossing rule on the six pre-specified points; CIs from the same 1,000 joint bootstrap draws. Design B already uses the six pre-specified points. Regression: regenerated bootstrap draws equal the cached draws and the r2_A intervals (pass: S1 True, S2 True); full-grid crossings and CIs reproduced (S1 True, S2 True).

| Sample | Halves | Index | Level | Full grid (printed) [95% CI] | Pre-specified grid [95% CI] | Difference |
|---|---|---|---|---|---|---|
| S1 | random | Bias | .70 | 14.7 [13.1, 16.9] | 15.4 [13.5, 17.4] | +0.7 |
| S1 | random | Bias | .80 | 26.0 [23.4, 28.4] | 27.4 [23.8, 30.9] | +1.4 |
| S1 | random | Bias | .90 | 58.0 [50.2, 66.3] | 58.0 [50.2, 66.3] | +0.0 |
| S1 | random | Brier | .70 | 21.3 [18.2, 23.9] | 21.6 [18.4, 25.1] | +0.3 |
| S1 | random | Brier | .80 | 34.3 [30.3, 38.4] | 37.4 [33.1, 41.4] | +3.1 |
| S1 | random | Brier | .90 | 80.2 [71.1, 88.5] | 80.2 [71.1, 88.5] | +0.0 |
| S1 | random | AUROC, plan rule | .70 | 172.9 [154.6, 199.0] | 172.9 [154.6, 199.0] | +0.0 |
| S1 | random | AUROC, tie convention | .70 | 151.1 [133.1, 172.0] | 151.1 [133.1, 172.0] | +0.0 |
| S1 | random | Bias, guess-corrected coding | .70 | 17.5 [15.1, 20.1] | 17.9 [15.9, 20.1] | +0.4 |
| S1 | random | Bias, guess-corrected coding | .80 | 29.0 [26.8, 32.6] | 31.9 [28.5, 35.1] | +2.8 |
| S1 | random | Bias, guess-corrected coding | .90 | 68.6 [61.1, 76.1] | 68.6 [61.1, 76.1] | +0.0 |
| S1 | random | Brier, guess-corrected coding | .70 | 19.7 [16.6, 22.0] | 19.8 [17.8, 22.7] | +0.1 |
| S1 | random | Brier, guess-corrected coding | .80 | 32.2 [29.1, 35.5] | 35.0 [31.5, 38.2] | +2.8 |
| S1 | random | Brier, guess-corrected coding | .90 | 75.0 [67.6, 82.0] | 75.0 [67.6, 82.0] | +0.0 |
| S1 | temporal | Bias | .70 | 20.4 [14.5, 26.7] | 20.3 [16.9, 26.5] | -0.1 |
| S1 | temporal | Bias | .80 | 39.4 [34.6, 53.1] | 43.6 [36.3, 53.5] | +4.2 |
| S1 | temporal | Bias | .90 | > 200 [172.7, > 200] | > 200 [172.7, > 200] | n/a |
| S1 | temporal | Brier | .70 | 30.4 [22.8, 33.9] | 27.7 [20.9, 33.8] | -2.7 |
| S1 | temporal | Brier | .80 | 56.0 [39.2, 69.2] | 56.0 [45.0, 69.2] | +0.0 |
| S1 | temporal | Brier | .90 | > 200 [174.7, > 200] | > 200 [174.7, > 200] | n/a |
| S1 | temporal | Bias, guess-corrected coding | .70 | 24.4 [18.4, 30.7] | 23.5 [18.8, 29.2] | -0.9 |
| S1 | temporal | Bias, guess-corrected coding | .80 | 46.2 [37.5, 61.1] | 48.2 [40.8, 61.2] | +2.0 |
| S1 | temporal | Bias, guess-corrected coding | .90 | > 200 [181.4, > 200] | > 200 [181.4, > 200] | n/a |
| S1 | temporal | Brier, guess-corrected coding | .70 | 25.7 [19.4, 32.4] | 25.1 [19.6, 30.6] | -0.6 |
| S1 | temporal | Brier, guess-corrected coding | .80 | 50.9 [37.7, 62.3] | 50.9 [43.1, 62.4] | +0.0 |
| S1 | temporal | Brier, guess-corrected coding | .90 | > 200 [152.3, > 200] | > 200 [152.3, > 200] | n/a |
| S2 | random | Bias | .70 | 16.2 [14.4, 19.2] | 16.5 [14.2, 19.3] | +0.3 |
| S2 | random | Bias | .80 | 28.7 [24.2, 33.8] | 29.9 [25.2, 34.7] | +1.3 |
| S2 | random | Bias | .90 | 68.0 [58.8, 77.5] | 68.0 [58.8, 77.5] | +0.0 |
| S2 | random | Brier | .70 | 31.1 [24.6, 35.5] | 35.9 [30.0, 42.4] | +4.8 |
| S2 | random | Brier | .80 | 58.8 [46.0, 69.2] | 58.8 [49.1, 69.2] | +0.0 |
| S2 | random | Brier | .90 | 121.6 [101.9, 142.5] | 121.6 [101.9, 142.5] | +0.0 |
| S2 | random | AUROC, plan rule | .70 | 183.8 [159.5, > 200] | 183.8 [159.5, > 200] | +0.0 |
| S2 | random | AUROC, tie convention | .70 | 168.7 [144.6, 194.2] | 168.7 [144.6, 194.2] | +0.0 |
| S2 | random | Bias, guess-corrected coding | .70 | 20.7 [17.8, 24.6] | 20.8 [17.8, 25.6] | +0.2 |
| S2 | random | Bias, guess-corrected coding | .80 | 36.7 [31.8, 42.4] | 38.7 [33.7, 44.3] | +2.0 |
| S2 | random | Bias, guess-corrected coding | .90 | 86.3 [77.1, 96.0] | 86.3 [77.1, 96.0] | +0.0 |
| S2 | random | Brier, guess-corrected coding | .70 | 26.2 [22.8, 31.6] | 30.0 [24.3, 35.7] | +3.8 |
| S2 | random | Brier, guess-corrected coding | .80 | 46.4 [37.2, 58.4] | 49.2 [43.4, 58.4] | +2.8 |
| S2 | random | Brier, guess-corrected coding | .90 | 104.9 [92.8, 123.7] | 104.9 [92.8, 123.7] | +0.0 |
| S2 | temporal | Bias | .70 | 20.4 [14.8, 26.1] | 20.6 [16.6, 29.0] | +0.2 |
| S2 | temporal | Bias | .80 | 48.1 [35.0, 68.3] | 48.6 [37.5, 68.3] | +0.6 |
| S2 | temporal | Bias | .90 | > 200 [169.0, > 200] | > 200 [169.0, > 200] | n/a |
| S2 | temporal | Brier | .70 | 48.1 [28.6, 66.3] | 49.4 [40.1, 66.3] | +1.2 |
| S2 | temporal | Brier | .80 | 95.1 [75.1, 162.0] | 95.1 [75.1, 162.0] | +0.0 |
| S2 | temporal | Bias, guess-corrected coding | .70 | 26.3 [19.3, 36.1] | 28.2 [19.2, 38.0] | +1.9 |
| S2 | temporal | Bias, guess-corrected coding | .80 | 66.8 [47.1, 88.8] | 66.8 [47.5, 88.8] | +0.0 |
| S2 | temporal | Brier, guess-corrected coding | .70 | 29.9 [26.0, 48.9] | 41.8 [34.0, 52.7] | +11.9 |
| S2 | temporal | Brier, guess-corrected coding | .80 | 84.8 [59.8, 125.3] | 84.8 [59.8, 125.3] | +0.0 |

Rows where neither grid reaches the level and every draw stays above 200 are omitted (gamma at every level; the resolution indices at .80 and .90).

S2 replication criterion (b) on both grids (printed draw):

- criterion b full grid plan rule: Bias 16.2, Brier 31.1, Gamma > 200, AUROC, plan rule 183.8; met True
- criterion b full grid tie convention: Bias 16.2, Brier 31.1, Gamma > 200, AUROC, tie convention 168.7; met True
- criterion b frozen grid plan rule: Bias 16.5, Brier 35.9, Gamma > 200, AUROC, plan rule 183.8; met True
- criterion b frozen grid tie convention: Bias 16.5, Brier 35.9, Gamma > 200, AUROC, tie convention 168.7; met True

## 2. Crossed S2 subset: alpha and Phi against ICC(1,k)

14 quizzes, 12 to 18 items, 4,973 quiz-student rows. ICC(1,k), alpha and Phi recomputed from the data and matched to the stored values (max abs diff 0.0e+00; pass True). Differences against ICC(1,k) are deterministic.

| Index | Alpha minus ICC(1,k): median (range) | Phi minus ICC(1,k): median (range) | Printed alpha minus split-half median | Over 200 split reruns: mean (SD); 95% range | Reruns rounding to printed | Printed Phi minus split-half median | Over 200 reruns: mean (SD); 95% range | Reruns rounding to printed |
|---|---|---|---|---|---|---|---|---|
| Bias | 0.019 (0.010 to 0.064) | 0.006 (0.003 to 0.021) | 0.022 | 0.020 (0.001); 0.017 to 0.023 | 0.07 | 0.007 | 0.006 (0.001); 0.004 to 0.009 | 0.28 |
| Brier | 0.036 (0.015 to 0.084) | 0.017 (0.005 to 0.041) | 0.039 | 0.035 (0.003); 0.029 to 0.039 | 0.04 | 0.018 | 0.015 (0.002); 0.011 to 0.019 | 0.07 |

- Bias: alpha exceeded ICC(1,k) in 14 of 14 quizzes, Phi in 14 of 14; medians alpha 0.718, Phi 0.699, ICC(1,k) 0.690; split-half minus ICC(1,k), median over quizzes, mean over reruns 0.000 (SD 0.001).
- Brier: alpha exceeded ICC(1,k) in 14 of 14 quizzes, Phi in 14 of 14; medians alpha 0.572, Phi 0.562, ICC(1,k) 0.554; split-half minus ICC(1,k), median over quizzes, mean over reruns 0.000 (SD 0.002).

## 3. AUROC minus gamma reliability (Design A, random halves, N <= 50)

- S1, plan rule: printed draw -0.017 (N = 12) to +0.063 (N = 50); mean over 100 further draws +0.013 to +0.061 (SD of the difference per draw 0.028 to 0.036)
- S1, tie convention: printed draw +0.058 (N = 25) to +0.149 (N = 5); mean over 100 further draws +0.069 to +0.129 (SD of the difference per draw 0.031 to 0.066)
- S2, plan rule: printed draw +0.009 (N = 20) to +0.167 (N = 40); mean over 100 further draws +0.010 to +0.146 (SD of the difference per draw 0.034 to 0.059)
- S2, tie convention: printed draw +0.049 (N = 12) to +0.210 (N = 40); mean over 100 further draws +0.043 to +0.202 (SD of the difference per draw 0.046 to 0.071)

## 4. Students entering the gamma and AUROC curves (Design A, random halves)

A student enters the correlation at N only when both halves give a defined value. Shares are means over the 50 splits of the printed draw; one-half undefined = share with half A undefined (the definition behind the halves row of Table S2a).

| Sample, N | Enter gamma and plan-rule AUROC curves | Enter tie-convention AUROC curve | One half undefined (gamma) | Enter temporal gamma curve |
|---|---|---|---|---|
| S1, 5 | 8.7% | 20.7% | 82.7% | 7.7% |
| S1, 10 | 34.9% | 53.1% | 49.4% | 30.6% |
| S1, 12 | 41.0% | 59.8% | 45.1% | 37.2% |
| S1, 15 | 49.0% | 67.5% | 40.0% | 44.4% |
| S1, 20 | 59.2% | 76.4% | 31.2% | 53.3% |
| S1, 25 | 63.6% | 79.8% | 28.3% | 57.5% |
| S1, 30 | 67.2% | 83.0% | 24.6% | 59.5% |
| S1, 40 | 73.2% | 86.5% | 20.3% | 66.2% |
| S1, 50 | 76.7% | 88.1% | 17.8% | 70.4% |
| S1, 100 | 84.5% | 92.8% | 12.5% | 76.4% |
| S1, 200 | 89.4% | 94.4% | 8.8% | 81.8% |
| S2, 5 | 12.3% | 28.0% | 78.4% | 9.7% |
| S2, 10 | 47.4% | 68.0% | 37.3% | 43.3% |
| S2, 12 | 59.0% | 77.0% | 27.5% | 49.8% |
| S2, 15 | 68.5% | 84.5% | 23.5% | 64.6% |
| S2, 20 | 74.7% | 90.1% | 17.9% | 69.9% |
| S2, 25 | 81.4% | 93.5% | 14.0% | 74.9% |
| S2, 30 | 84.2% | 94.7% | 11.4% | 75.1% |
| S2, 40 | 87.4% | 96.4% | 9.7% | 82.3% |
| S2, 50 | 89.9% | 97.1% | 7.7% | 84.3% |
| S2, 100 | 93.6% | 97.8% | 5.0% | 87.5% |
| S2, 200 | 95.0% | 97.8% | 4.1% | 89.4% |

## 5. The AUROC curve between N = 15 and 20

- S1: printed plan-rule AUROC reliability 0.236 at N = 15 and 0.172 at N = 20 (difference +0.064); over 100 further draws the means are 0.134 (SD 0.044) and 0.179 (SD 0.038), difference -0.045 (SD of the difference 0.060); r(15) > r(20) in 0.22 of draws; mean curve increasing at every step up to N = 50: True. Students entering the curve: 49.0% at N = 15, 59.2% at N = 20. Halves at N = 15 hold 7 and 8 answers.
- S2: printed plan-rule AUROC reliability 0.003 at N = 15 and 0.056 at N = 20 (difference -0.053); over 100 further draws the means are 0.098 (SD 0.062) and 0.129 (SD 0.053), difference -0.031 (SD of the difference 0.078); r(15) > r(20) in 0.35 of draws; mean curve increasing at every step up to N = 50: True. Students entering the curve: 68.5% at N = 15, 74.7% at N = 20. Halves at N = 15 hold 7 and 8 answers.

Scripts: scripts/r2_F_grid.py (computation), scripts/r2_H_report.py (this file).
