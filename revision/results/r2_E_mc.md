No flag: for every crossing with a printed 95% CI, the mean (or, where some draws never reach the level, the median) over 100 further subsample draws lies inside that CI.

Printed values in the outer 2.5% of the further draws (the printed CI still covers the draw centre):
- S1 Brier .70, A, random halves: printed 21.3 [18.2, 23.9]; 1.00 of draws lie below it; draw mean 19.8 (SD 0.74), 95% of draws 18.6 to 21.2.
- S1 Brier .70, A, temporal halves: printed 30.4 [22.8, 33.9]; 0.99 of draws lie below it; draw mean 25.4 (SD 1.89), 95% of draws 22.3 to 28.9.
- S1 Brier .70, A, variance-components ICC(1,N) curve: printed 21.3 [18.2, 24.0]; 0.99 of draws lie below it; draw mean 19.8 (SD 0.73), 95% of draws 18.5 to 21.1.
- S2 Brier .80, A, random halves: printed 58.8 [46.0, 69.2]; 1.00 of draws lie below it; draw mean 51.7 (SD 2.71), 95% of draws 47.7 to 57.3.
- S2 Brier .80, A, variance-components ICC(1,N) curve: printed 58.3 [48.1, 69.2]; 1.00 of draws lie below it; draw mean 51.7 (SD 2.61), 95% of draws 47.8 to 57.1.

Second check (printed value inside its CI shifted to the draw centre): fails for S1 Brier .70 (A, temporal halves).

# Monte Carlo precision of the reliability curves and minimum N

Two Monte Carlo components sit in every printed Design A or Design B value: the split component (each value averages 50 random splits of one subsample) and the subsample component (one seeded subsample of N rated answers per student). The split component is measured from the 50 per-split values of the printed draw. The total of both is measured by rerunning the frozen estimator on 100 further subsample draws (Design A; 50 for Design B) and taking the SD over draws. Printed values come from the pre-specified single draw and stay as printed; the bootstrap CIs hold the subsample and the splits fixed and cover student sampling only. Temporal halves, the variance-components (ICC) curve and the split-free estimator have no split component; the split-free estimator has no Monte Carlo error at all.

## 1. Maximum Monte Carlo error of a reliability value, Design A random halves

| Sample | Index | Split SE, max over all N (N) | Split SE, max over N >= 20 | Draw SD, max over all N (N) | Draw SD, max over N >= 20 | Draw SD at N = 20 (random / temporal) |
|---|---|---|---|---|---|---|
| S1 | Bias | 0.0029 (5) | 0.0016 | 0.0223 (5) | 0.0056 | 0.0056 / 0.0130 |
| S1 | Brier | 0.0041 (5) | 0.0024 | 0.0243 (5) | 0.0095 | 0.0095 / 0.0204 |
| S1 | Gamma | 0.0235 (5) | 0.0068 | 0.1222 (5) | 0.0474 | 0.0474 / 0.0721 |
| S1 | AUROC, plan rule | 0.0205 (5) | 0.0065 | 0.1221 (5) | 0.0382 | 0.0382 / 0.0667 |
| S1 | AUROC, tie convention | 0.0104 (5) | 0.0051 | 0.0734 (5) | 0.0314 | 0.0314 / 0.0510 |
| S1 | Mean confidence | 0.0013 (5) | 0.0005 | 0.0093 (5) | 0.0013 | 0.0013 / 0.0061 |
| S1 | Accuracy | 0.0032 (5) | 0.0013 | 0.0213 (5) | 0.0056 | 0.0056 / 0.0126 |
| S2 | Bias | 0.0051 (5) | 0.0021 | 0.0423 (5) | 0.0123 | 0.0123 / 0.0231 |
| S2 | Brier | 0.0062 (5) | 0.0038 | 0.0578 (5) | 0.0244 | 0.0244 / 0.0373 |
| S2 | Gamma | 0.0287 (5) | 0.0125 | 0.1555 (5) | 0.0755 | 0.0601 / 0.0999 |
| S2 | AUROC, plan rule | 0.0314 (5) | 0.0103 | 0.1590 (5) | 0.0599 | 0.0531 / 0.0924 |
| S2 | AUROC, tie convention | 0.0216 (5) | 0.0091 | 0.1057 (5) | 0.0532 | 0.0477 / 0.0787 |
| S2 | Mean confidence | 0.0025 (5) | 0.0009 | 0.0152 (5) | 0.0030 | 0.0030 / 0.0101 |
| S2 | Accuracy | 0.0053 (5) | 0.0026 | 0.0405 (5) | 0.0128 | 0.0128 / 0.0235 |

Largest split SE of any printed Design A random-half crossing: bias and Brier 1.04 answers; AUROC (either rule) 1.70 answers; structured-halves and weighting table 0.64 answers.

## 2. Per-N Monte Carlo error, Design A

Cells: split SE / total draw SD. Temporal halves: no split SE; their draw SD is in r2_E_mc.json (reliability_over_draws).

| Sample, N | Bias | Brier | Gamma | AUROC, plan rule | AUROC, tie convention |
|---|---|---|---|---|---|
| S1, 5 | 0.0029 / 0.0223 | 0.0041 / 0.0243 | 0.0235 / 0.1222 | 0.0205 / 0.1221 | 0.0104 / 0.0734 |
| S1, 10 | 0.0020 / 0.0124 | 0.0028 / 0.0172 | 0.0097 / 0.0567 | 0.0087 / 0.0523 | 0.0061 / 0.0396 |
| S1, 12 | 0.0020 / 0.0105 | 0.0027 / 0.0151 | 0.0094 / 0.0452 | 0.0088 / 0.0487 | 0.0066 / 0.0384 |
| S1, 15 | 0.0020 / 0.0091 | 0.0024 / 0.0137 | 0.0071 / 0.0435 | 0.0057 / 0.0443 | 0.0046 / 0.0363 |
| S1, 20 | 0.0016 / 0.0056 | 0.0024 / 0.0095 | 0.0064 / 0.0474 | 0.0062 / 0.0382 | 0.0051 / 0.0314 |
| S1, 25 | 0.0011 / 0.0051 | 0.0016 / 0.0081 | 0.0068 / 0.0437 | 0.0065 / 0.0374 | 0.0050 / 0.0293 |
| S1, 30 | 0.0009 / 0.0039 | 0.0012 / 0.0059 | 0.0059 / 0.0440 | 0.0054 / 0.0379 | 0.0044 / 0.0307 |
| S1, 40 | 0.0007 / 0.0027 | 0.0009 / 0.0038 | 0.0042 / 0.0366 | 0.0041 / 0.0281 | 0.0033 / 0.0233 |
| S1, 50 | 0.0006 / 0.0020 | 0.0009 / 0.0029 | 0.0053 / 0.0411 | 0.0040 / 0.0264 | 0.0034 / 0.0217 |
| S1, 100 | 0.0003 / 0.0007 | 0.0005 / 0.0011 | 0.0036 / 0.0374 | 0.0034 / 0.0156 | 0.0030 / 0.0137 |
| S1, 200 | 0.0002 / 0.0002 | 0.0002 / 0.0004 | 0.0041 / 0.0217 | 0.0022 / 0.0086 | 0.0021 / 0.0075 |
| S2, 5 | 0.0051 / 0.0423 | 0.0062 / 0.0578 | 0.0287 / 0.1555 | 0.0314 / 0.1590 | 0.0216 / 0.1057 |
| S2, 10 | 0.0033 / 0.0228 | 0.0053 / 0.0393 | 0.0124 / 0.0846 | 0.0124 / 0.0843 | 0.0096 / 0.0707 |
| S2, 12 | 0.0033 / 0.0185 | 0.0051 / 0.0366 | 0.0106 / 0.0700 | 0.0101 / 0.0697 | 0.0083 / 0.0620 |
| S2, 15 | 0.0025 / 0.0157 | 0.0051 / 0.0262 | 0.0127 / 0.0695 | 0.0116 / 0.0617 | 0.0100 / 0.0537 |
| S2, 20 | 0.0021 / 0.0123 | 0.0038 / 0.0244 | 0.0104 / 0.0601 | 0.0092 / 0.0531 | 0.0073 / 0.0477 |
| S2, 25 | 0.0017 / 0.0088 | 0.0031 / 0.0154 | 0.0102 / 0.0737 | 0.0103 / 0.0599 | 0.0089 / 0.0532 |
| S2, 30 | 0.0014 / 0.0069 | 0.0028 / 0.0126 | 0.0125 / 0.0726 | 0.0103 / 0.0548 | 0.0091 / 0.0482 |
| S2, 40 | 0.0016 / 0.0055 | 0.0025 / 0.0110 | 0.0089 / 0.0704 | 0.0068 / 0.0455 | 0.0064 / 0.0396 |
| S2, 50 | 0.0009 / 0.0042 | 0.0022 / 0.0086 | 0.0093 / 0.0755 | 0.0075 / 0.0458 | 0.0066 / 0.0403 |
| S2, 100 | 0.0006 / 0.0012 | 0.0010 / 0.0031 | 0.0087 / 0.0643 | 0.0040 / 0.0277 | 0.0038 / 0.0240 |
| S2, 200 | 0.0003 / 0.0005 | 0.0006 / 0.0012 | 0.0077 / 0.0459 | 0.0028 / 0.0096 | 0.0026 / 0.0083 |

## 3. Minimum N: printed value, split SE and spread over further subsample draws

| Sample | Where | Design | Index | Level | Printed [95% CI] | Split SE (delta / resample) | Draw mean (SD) | Draw median; 95% of draws | Reaching within 200 | Crossing of mean curve | Draws below printed | Draw centre in printed CI | Printed in CI recentred on draw centre |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .70 | 14.7 [13.1, 16.9] | 0.11 / 0.11 | 14.8 (0.54) | 14.8; 14.0 to 15.9 | 1.00 | 14.8 | 0.45 | yes | yes |
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .80 | 26.0 [23.4, 28.4] | 0.11 / 0.11 | 25.1 (0.67) | 25.1; 24.0 to 26.4 | 1.00 | 25.1 | 0.86 | yes | yes |
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .90 | 58.0 [50.2, 66.3] | 0.35 / 0.34 | 58.5 (1.21) | 58.5; 56.1 to 60.8 | 1.00 | 58.5 | 0.34 | yes | yes |
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .70 | 21.3 [18.2, 23.9] | 0.17 / 0.17 | 19.8 (0.74) | 19.7; 18.6 to 21.2 | 1.00 | 19.7 | 1.00 | yes | yes |
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .80 | 34.3 [30.3, 38.4] | 0.16 / 0.16 | 34.1 (0.84) | 34.1; 32.5 to 35.7 | 1.00 | 34.1 | 0.60 | yes | yes |
| S1 | S1 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .90 | 80.2 [71.1, 88.5] | 0.36 / 0.35 | 80.0 (1.04) | 79.9; 77.9 to 81.6 | 1.00 | 80.0 | 0.58 | yes | yes |
| S1 | S1 minimum-N table; Table S2b | A, random halves | AUROC, plan rule | .70 | 172.9 [154.6, 199.0] | 1.43 / 1.40 | 171.7 (5.72) | 171.5; 161.2 to 181.4 | 1.00 | 171.6 | 0.56 | yes | yes |
| S1 | S1 minimum-N table; Table S2b | A, random halves | AUROC, tie convention | .70 | 151.1 [133.1, 172.0] | 1.32 / 1.29 | 149.1 (5.33) | 149.4; 137.7 to 157.3 | 1.00 | 149.2 | 0.59 | yes | yes |
| S1 | Table S1f | A, random halves | Bias, guess-corrected coding | .70 | 17.5 [15.1, 20.1] | 0.13 / 0.14 | 17.1 (0.51) | 17.1; 16.0 to 17.9 | 1.00 | 17.1 | 0.82 | yes | yes |
| S1 | Table S1f | A, random halves | Bias, guess-corrected coding | .80 | 29.0 [26.8, 32.6] | 0.12 / 0.12 | 29.4 (0.78) | 29.3; 27.9 to 31.0 | 1.00 | 29.3 | 0.37 | yes | yes |
| S1 | not printed | A, random halves | Bias, guess-corrected coding | .90 | 68.6 [61.1, 76.1] | 0.32 / 0.31 | 68.9 (1.12) | 69.0; 66.6 to 70.9 | 1.00 | 69.0 | 0.42 | yes | yes |
| S1 | Table S1f | A, random halves | Brier, guess-corrected coding | .70 | 19.7 [16.6, 22.0] | 0.31 / 0.28 | 18.5 (0.67) | 18.4; 17.4 to 19.9 | 1.00 | 18.5 | 0.93 | yes | yes |
| S1 | Table S1f | A, random halves | Brier, guess-corrected coding | .80 | 32.2 [29.1, 35.5] | 0.18 / 0.19 | 31.5 (0.96) | 31.6; 29.9 to 33.5 | 1.00 | 31.6 | 0.74 | yes | yes |
| S1 | not printed | A, random halves | Brier, guess-corrected coding | .90 | 75.0 [67.6, 82.0] | 0.36 / 0.36 | 74.5 (1.10) | 74.6; 72.5 to 76.6 | 1.00 | 74.6 | 0.65 | yes | yes |
| S1 | S1 minimum-N table; bias sensitivity table | A, temporal halves | Bias | .70 | 20.4 [14.5, 26.7] | n/a / n/a | 19.4 (1.09) | 19.3; 17.6 to 21.7 | 1.00 | 19.4 | 0.80 | yes | yes |
| S1 | S1 minimum-N table; bias sensitivity table | A, temporal halves | Bias | .80 | 39.4 [34.6, 53.1] | n/a / n/a | 41.6 (2.68) | 41.7; 36.8 to 45.4 | 1.00 | 41.7 | 0.22 | yes | yes |
| S1 | S1 minimum-N table | A, temporal halves | Brier | .70 | 30.4 [22.8, 33.9] | n/a / n/a | 25.4 (1.89) | 25.2; 22.3 to 28.9 | 1.00 | 25.6 | 0.99 | yes | NO |
| S1 | S1 minimum-N table | A, temporal halves | Brier | .80 | 56.0 [39.2, 69.2] | n/a / n/a | 52.1 (4.02) | 51.0; 46.7 to 59.6 | 1.00 | 51.6 | 0.79 | yes | yes |
| S1 | not printed | A, temporal halves | Bias, guess-corrected coding | .70 | 24.4 [18.4, 30.7] | n/a / n/a | 22.3 (1.66) | 22.4; 19.6 to 25.7 | 1.00 | 22.4 | 0.89 | yes | yes |
| S1 | not printed | A, temporal halves | Bias, guess-corrected coding | .80 | 46.2 [37.5, 61.1] | n/a / n/a | 47.5 (3.32) | 47.1; 41.2 to 56.3 | 1.00 | 47.2 | 0.32 | yes | yes |
| S1 | not printed | A, temporal halves | Brier, guess-corrected coding | .70 | 25.7 [19.4, 32.4] | n/a / n/a | 23.1 (1.94) | 22.9; 19.5 to 26.1 | 1.00 | 23.1 | 0.93 | yes | yes |
| S1 | not printed | A, temporal halves | Brier, guess-corrected coding | .80 | 50.9 [37.7, 62.3] | n/a / n/a | 47.1 (3.37) | 46.6; 41.6 to 55.7 | 1.00 | 46.8 | 0.89 | yes | yes |
| S1 | Table S1f | A, random halves, without near-perfect and constant raters | Bias | .70 | 16.0 [14.0, 18.3] | 0.15 / 0.15 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S1 | Table S1f | A, random halves, without near-perfect and constant raters | Bias | .80 | 27.5 [24.8, 30.2] | 0.10 / 0.10 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S1 | Table S1f | A, random halves, without near-perfect and constant raters | Brier | .70 | 25.6 [22.7, 28.9] | 0.18 / 0.17 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S1 | Table S1f | A, random halves, without near-perfect and constant raters | Brier | .80 | 41.4 [36.6, 47.4] | 0.30 / 0.30 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S1 | Table S1b | A, variance-components ICC(1,N) curve | Bias | .70 | 14.6 [13.2, 16.5] | 0.00 / 0.00 | 14.8 (0.52) | 14.7; 13.9 to 15.9 | 1.00 | 14.7 | 0.38 | yes | yes |
| S1 | Table S1b | A, variance-components ICC(1,N) curve | Bias | .80 | 26.1 [23.3, 28.4] | 0.00 / 0.00 | 25.2 (0.66) | 25.1; 24.0 to 26.3 | 1.00 | 25.2 | 0.91 | yes | yes |
| S1 | Table S1b | A, variance-components ICC(1,N) curve | Brier | .70 | 21.3 [18.2, 24.0] | 0.00 / 0.00 | 19.8 (0.73) | 19.7; 18.5 to 21.1 | 1.00 | 19.7 | 0.99 | yes | yes |
| S1 | Table S1b | A, variance-components ICC(1,N) curve | Brier | .80 | 34.7 [31.1, 38.7] | 0.00 / 0.00 | 34.1 (0.80) | 34.1; 32.7 to 35.6 | 1.00 | 34.2 | 0.74 | yes | yes |
| S1 | Table S1b | split-free estimator on all rated answers | Bias | .70 | 14.7 [13.2, 16.4] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S1 | Table S1b | split-free estimator on all rated answers | Bias | .80 | 25.2 [22.7, 28.2] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S1 | Table S1b | split-free estimator on all rated answers | Brier | .70 | 19.8 [17.9, 22.1] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S1 | Table S1b | split-free estimator on all rated answers | Brier | .80 | 33.9 [30.7, 37.8] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .70 | 17.7  | 0.03 / 0.03 | 17.8 (0.13) | 17.8; 17.5 to 18.0 | 1.00 | 17.8 |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .80 | 31.7  | 0.05 / 0.05 | 31.7 (0.20) | 31.7; 31.3 to 32.1 | 1.00 | 31.7 |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .90 | 66.4  | 0.13 / 0.13 | 66.4 (0.33) | 66.4; 65.8 to 67.0 | 1.00 | 66.4 |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .70 | 27.9  | 0.05 / 0.05 | 28.2 (0.25) | 28.2; 27.7 to 28.6 | 1.00 | 28.2 |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .80 | 45.0  | 0.08 / 0.08 | 45.3 (0.23) | 45.3; 44.9 to 45.7 | 1.00 | 45.3 |  |  |  |
| S1 | S1 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .90 | 91.2  | 0.22 / 0.23 | 90.7 (0.47) | 90.7; 89.9 to 91.6 | 1.00 | 90.7 |  |  |  |
| S1 | S1 minimum-N table; Table S2b | B, random halves (point estimate, no bootstrap) | AUROC, plan rule | .70 | 170.8  | 1.14 / 1.12 | 172.8 (3.98) | 172.6; 166.4 to 180.2 | 1.00 | 172.7 |  |  |  |
| S1 | S1 minimum-N table; Table S2b | B, random halves (point estimate, no bootstrap) | AUROC, tie convention | .70 | 151.1  | 0.89 / 0.89 | 153.3 (3.30) | 153.4; 147.6 to 159.3 | 1.00 | 153.3 |  |  |  |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .70 | 16.2 [14.4, 19.2] | 0.15 / 0.15 | 16.8 (0.93) | 16.9; 15.0 to 18.5 | 1.00 | 16.8 | 0.27 | yes | yes |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .80 | 28.7 [24.2, 33.8] | 0.24 / 0.24 | 28.4 (1.12) | 28.4; 26.3 to 30.8 | 1.00 | 28.3 | 0.69 | yes | yes |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Bias | .90 | 68.0 [58.8, 77.5] | 0.47 / 0.46 | 66.9 (1.94) | 66.9; 63.0 to 70.4 | 1.00 | 67.0 | 0.70 | yes | yes |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .70 | 31.1 [24.6, 35.5] | 0.28 / 0.28 | 29.9 (1.66) | 29.7; 27.5 to 33.0 | 1.00 | 29.9 | 0.74 | yes | yes |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .80 | 58.8 [46.0, 69.2] | 0.64 / 0.64 | 51.7 (2.71) | 51.7; 47.7 to 57.3 | 1.00 | 51.5 | 1.00 | yes | yes |
| S2 | S2 minimum-N table; Table S1b; Table S1f | A, random halves | Brier | .90 | 121.6 [101.9, 142.5] | 1.04 / 1.01 | 119.0 (3.58) | 118.9; 112.4 to 125.5 | 1.00 | 119.1 | 0.73 | yes | yes |
| S2 | S2 minimum-N table; Table S2b | A, random halves | AUROC, plan rule | .70 | 183.8 [159.5, > 200] | 1.70 / 1.70 | n/a (5.31) | 192.5; 179.5 to > 200 | 0.85 | 192.7 | 0.09 | yes | yes |
| S2 | S2 minimum-N table; Table S2b | A, random halves | AUROC, tie convention | .70 | 168.7 [144.6, 194.2] | 1.48 / 1.48 | 175.6 (6.64) | 176.3; 160.2 to 187.8 | 1.00 | 175.8 | 0.14 | yes | yes |
| S2 | Table S1f | A, random halves | Bias, guess-corrected coding | .70 | 20.7 [17.8, 24.6] | 0.23 / 0.23 | 21.0 (1.08) | 21.1; 18.8 to 22.8 | 1.00 | 21.1 | 0.34 | yes | yes |
| S2 | Table S1f | A, random halves | Bias, guess-corrected coding | .80 | 36.7 [31.8, 42.4] | 0.30 / 0.29 | 36.3 (1.33) | 36.3; 34.1 to 38.8 | 1.00 | 36.3 | 0.62 | yes | yes |
| S2 | not printed | A, random halves | Bias, guess-corrected coding | .90 | 86.3 [77.1, 96.0] | 0.54 / 0.52 | 84.3 (1.45) | 84.2; 81.5 to 86.7 | 1.00 | 84.3 | 0.94 | yes | yes |
| S2 | Table S1f | A, random halves | Brier, guess-corrected coding | .70 | 26.2 [22.8, 31.6] | 0.30 / 0.30 | 25.6 (1.43) | 25.7; 23.1 to 28.0 | 1.00 | 25.7 | 0.63 | yes | yes |
| S2 | Table S1f | A, random halves | Brier, guess-corrected coding | .80 | 46.4 [37.2, 58.4] | 1.63 / 1.83 | 44.5 (1.92) | 44.6; 39.7 to 47.7 | 1.00 | 44.6 | 0.88 | yes | yes |
| S2 | not printed | A, random halves | Brier, guess-corrected coding | .90 | 104.9 [92.8, 123.7] | 1.10 / 1.10 | 100.0 (2.79) | 99.6; 95.9 to 106.0 | 1.00 | 99.5 | 0.93 | yes | yes |
| S2 | S2 minimum-N table; bias sensitivity table | A, temporal halves | Bias | .70 | 20.4 [14.8, 26.1] | n/a / n/a | 22.1 (2.52) | 22.0; 18.4 to 26.2 | 1.00 | 22.4 | 0.32 | yes | yes |
| S2 | S2 minimum-N table; bias sensitivity table | A, temporal halves | Bias | .80 | 48.1 [35.0, 68.3] | n/a / n/a | 47.7 (5.16) | 47.1; 38.7 to 57.1 | 1.00 | 47.7 | 0.55 | yes | yes |
| S2 | S2 minimum-N table | A, temporal halves | Brier | .70 | 48.1 [28.6, 66.3] | n/a / n/a | 41.2 (4.61) | 40.7; 30.9 to 50.5 | 1.00 | 42.1 | 0.94 | yes | yes |
| S2 | S2 minimum-N table | A, temporal halves | Brier | .80 | 95.1 [75.1, 162.0] | n/a / n/a | 97.8 (9.41) | 95.9; 82.3 to 116.6 | 1.00 | 96.2 | 0.47 | yes | yes |
| S2 | not printed | A, temporal halves | Bias, guess-corrected coding | .70 | 26.3 [19.3, 36.1] | n/a / n/a | 28.3 (2.76) | 28.4; 23.6 to 33.2 | 1.00 | 28.6 | 0.20 | yes | yes |
| S2 | not printed | A, temporal halves | Bias, guess-corrected coding | .80 | 66.8 [47.1, 88.8] | n/a / n/a | 62.4 (6.76) | 63.2; 48.1 to 73.3 | 1.00 | 63.4 | 0.75 | yes | yes |
| S2 | not printed | A, temporal halves | Brier, guess-corrected coding | .70 | 29.9 [26.0, 48.9] | n/a / n/a | 35.4 (4.26) | 35.8; 28.1 to 44.1 | 1.00 | 35.7 | 0.12 | yes | NO |
| S2 | not printed | A, temporal halves | Brier, guess-corrected coding | .80 | 84.8 [59.8, 125.3] | n/a / n/a | 80.5 (6.86) | 80.2; 69.4 to 95.4 | 1.00 | 80.4 | 0.76 | yes | yes |
| S2 | Table S1f | A, random halves, without near-perfect and constant raters | Bias | .70 | 17.0 [14.7, 20.5] | 0.15 / 0.14 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S2 | Table S1f | A, random halves, without near-perfect and constant raters | Bias | .80 | 30.5 [25.8, 35.6] | 0.26 / 0.26 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S2 | Table S1f | A, random halves, without near-perfect and constant raters | Brier | .70 | 33.9 [29.0, 39.0] | 0.26 / 0.26 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S2 | Table S1f | A, random halves, without near-perfect and constant raters | Brier | .80 | 66.2 [55.6, 77.7] | 0.59 / 0.59 | subsample draws not rerun for this variant |  |  | n/a |  |  |  |
| S2 | Table S1b | A, variance-components ICC(1,N) curve | Bias | .70 | 16.3 [14.4, 19.6] | 0.00 / 0.00 | 16.7 (0.93) | 16.9; 14.9 to 18.4 | 1.00 | 16.8 | 0.30 | yes | yes |
| S2 | Table S1b | A, variance-components ICC(1,N) curve | Bias | .80 | 28.8 [24.4, 33.5] | 0.00 / 0.00 | 28.4 (1.12) | 28.4; 26.5 to 30.7 | 1.00 | 28.4 | 0.68 | yes | yes |
| S2 | Table S1b | A, variance-components ICC(1,N) curve | Brier | .70 | 31.1 [24.9, 35.6] | 0.00 / 0.00 | 30.0 (1.61) | 29.8; 27.6 to 33.0 | 1.00 | 30.0 | 0.75 | yes | yes |
| S2 | Table S1b | A, variance-components ICC(1,N) curve | Brier | .80 | 58.3 [48.1, 69.2] | 0.00 / 0.00 | 51.7 (2.61) | 51.6; 47.8 to 57.1 | 1.00 | 51.6 | 1.00 | yes | yes |
| S2 | Table S1b | split-free estimator on all rated answers | Bias | .70 | 16.6 [14.5, 19.1] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S2 | Table S1b | split-free estimator on all rated answers | Bias | .80 | 28.4 [24.9, 32.8] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S2 | Table S1b | split-free estimator on all rated answers | Brier | .70 | 29.9 [26.0, 35.6] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S2 | Table S1b | split-free estimator on all rated answers | Brier | .80 | 51.3 [44.5, 61.0] | 0.00 / 0.00 | no Monte Carlo error (no subsample, no split) |  |  | n/a |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .70 | 18.4  | 0.08 / 0.08 | 18.6 (0.29) | 18.6; 18.2 to 19.2 | 1.00 | 18.6 |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .80 | 32.7  | 0.12 / 0.12 | 32.9 (0.42) | 32.9; 32.2 to 33.6 | 1.00 | 32.9 |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Bias | .90 | 70.6  | 0.30 / 0.30 | 70.1 (0.83) | 70.0; 68.7 to 71.5 | 1.00 | 70.1 |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .70 | 32.2  | 0.14 / 0.14 | 31.8 (0.59) | 32.0; 30.7 to 32.8 | 1.00 | 31.8 |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .80 | 49.7  | 0.22 / 0.25 | 50.5 (1.24) | 50.1; 48.8 to 53.0 | 1.00 | 50.4 |  |  |  |
| S2 | S2 minimum-N table | B, random halves (point estimate, no bootstrap) | Brier | .90 | 122.0  | 1.03 / 0.99 | 120.5 (2.72) | 120.6; 114.9 to 124.6 | 1.00 | 120.5 |  |  |  |
| S2 | S2 minimum-N table; Table S2b | B, random halves (point estimate, no bootstrap) | AUROC, plan rule | .70 | 198.4  | 1.68 / 1.31 | n/a (4.75) | 190.1; 181.8 to > 200 | 0.90 | 190.5 |  |  |  |
| S2 | S2 minimum-N table; Table S2b | B, random halves (point estimate, no bootstrap) | AUROC, tie convention | .70 | 182.2  | 1.39 / 1.36 | 175.4 (4.93) | 175.4; 167.2 to 184.8 | 1.00 | 175.2 |  |  |  |

## 4. Structured halves and weighting (bias table and its note): split SE of the crossings

These variants were not rerun on further subsample draws; they share the printed subsamples of the random-half curve.

| Sample | Variant | Index | .70: printed, split SE | .80: printed, split SE | r(20) split SE |
|---|---|---|---|---|---|
| S1 | Random halves (reference) | Bias | 14.7, 0.11 | 26.0, 0.11 | 0.0016 |
| S1 | Random halves (reference) | Brier | 21.3, 0.17 | 34.3, 0.16 | 0.0024 |
| S1 | Halves by content strand | Bias | 15.5, 0.11 | 29.0, 0.12 | 0.0012 |
| S1 | Halves by content strand | Brier | 23.1, 0.17 | 39.7, 0.20 | 0.0015 |
| S1 | Halves by topic | Bias | 16.1, 0.12 | 28.8, 0.13 | 0.0014 |
| S1 | Halves by topic | Brier | 23.8, 0.16 | 40.6, 0.35 | 0.0021 |
| S1 | Halves by approximate school term | Bias | 17.4, 0.08 | 31.5, 0.21 | 0.0013 |
| S1 | Halves by approximate school term | Brier | 26.2, 0.09 | 42.7, 0.33 | 0.0018 |
| S1 | Harder versus easier items | Bias | 14.6, 0.00 | 35.2, 0.00 | 0.0000 |
| S1 | Harder versus easier items | Brier | 33.1, 0.00 | 103.7, 0.00 | 0.0000 |
| S1 | Weighted to the unrated difficulty mix | Bias | 15.9, 0.14 | 27.2, 0.11 | 0.0018 |
| S1 | Weighted to the unrated difficulty mix | Brier | 22.7, 0.14 | 37.0, 0.18 | 0.0024 |
| S1 | Same weights shuffled within student (difficulty) | Bias | 15.5, 0.16 | 27.5, 0.09 | 0.0018 |
| S1 | Same weights shuffled within student (difficulty) | Brier | 22.7, 0.21 | 37.0, 0.18 | 0.0023 |
| S1 | Weighted to the unrated difficulty-by-strand mix | Bias | 24.6, 0.23 | 43.8, 0.22 | 0.0020 |
| S1 | Weighted to the unrated difficulty-by-strand mix | Brier | 29.2, 0.10 | 65.9, 0.47 | 0.0028 |
| S1 | Same weights shuffled within student (difficulty by strand) | Bias | 25.3, 0.14 | 42.9, 0.19 | 0.0024 |
| S1 | Same weights shuffled within student (difficulty by strand) | Brier | 30.0, 0.15 | 59.8, 0.61 | 0.0027 |
| S1 | Difficulty-balanced subset, equal weights | Bias | 13.9, 0.06 | 24.3, 0.15 | 0.0013 |
| S1 | Difficulty-balanced subset, equal weights | Brier | 18.9, 0.13 | 32.1, 0.22 | 0.0018 |
| S1 | Random halves (reference), students with a weighting target (n = 849) | Bias | 14.6, 0.09 | | |
| S1 | Weighted to the unrated difficulty mix, students with a weighting target (n = 849) | Bias | 16.6, 0.18 | | |
| S1 | Weighted to the unrated difficulty-by-strand mix, students with a weighting target (n = 849) | Bias | 27.5, 0.16 | | |
| S1 | Same weights shuffled within student (difficulty), students with a weighting target (n = 849) | Bias | 15.8, 0.23 | | |
| S1 | Same weights shuffled within student (difficulty by strand), students with a weighting target (n = 849) | Bias | 28.0, 0.14 | | |
| S1 | Difficulty-balanced subset, equal weights, students with a weighting target (n = 849) | Bias | 13.2, 0.08 | | |
| S2 | Random halves (reference) | Bias | 16.2, 0.15 | 28.7, 0.24 | 0.0021 |
| S2 | Random halves (reference) | Brier | 31.1, 0.28 | 58.8, 0.64 | 0.0038 |
| S2 | Halves by content strand | Bias | 18.2, 0.12 | 35.3, 0.24 | 0.0022 |
| S2 | Halves by content strand | Brier | 36.6, 0.29 | 81.3, 0.55 | 0.0035 |
| S2 | Halves by topic | Bias | 17.7, 0.15 | 32.8, 0.30 | 0.0021 |
| S2 | Halves by topic | Brier | 33.2, 0.30 | 67.8, 0.52 | 0.0037 |
| S2 | Halves by approximate school term | Bias | 18.3, 0.15 | 35.9, 0.27 | 0.0022 |
| S2 | Halves by approximate school term | Brier | 34.9, 0.31 | 76.5, 0.63 | 0.0031 |
| S2 | Harder versus easier items | Bias | 18.4, 0.00 | 35.3, 0.00 | 0.0000 |
| S2 | Harder versus easier items | Brier | > 200, n/a | > 200, n/a | 0.0000 |
| S2 | Weighted to the unrated difficulty mix | Bias | 17.8, 0.13 | 32.4, 0.23 | 0.0021 |
| S2 | Weighted to the unrated difficulty mix | Brier | 33.9, 0.21 | 60.3, 0.72 | 0.0040 |
| S2 | Same weights shuffled within student (difficulty) | Bias | 17.6, 0.15 | 31.2, 0.31 | 0.0022 |
| S2 | Same weights shuffled within student (difficulty) | Brier | 33.6, 0.24 | 65.6, 0.61 | 0.0037 |
| S2 | Weighted to the unrated difficulty-by-strand mix | Bias | 34.6, 0.22 | 65.0, 0.63 | 0.0028 |
| S2 | Weighted to the unrated difficulty-by-strand mix | Brier | 72.1, 0.73 | 115.2, 1.74 | 0.0044 |
| S2 | Same weights shuffled within student (difficulty by strand) | Bias | 33.8, 0.27 | 61.4, 0.64 | 0.0029 |
| S2 | Same weights shuffled within student (difficulty by strand) | Brier | 71.1, 0.78 | 125.5, 1.88 | 0.0060 |
| S2 | Difficulty-balanced subset, equal weights | Bias | 15.0, 0.13 | 27.6, 0.34 | 0.0020 |
| S2 | Difficulty-balanced subset, equal weights | Brier | 32.3, 0.28 | 48.3, 0.32 | 0.0032 |
| S2 | Random halves (reference), students with a weighting target (n = 352) | Bias | 15.7, 0.22 | | |
| S2 | Weighted to the unrated difficulty mix, students with a weighting target (n = 352) | Bias | 18.2, 0.17 | | |
| S2 | Weighted to the unrated difficulty-by-strand mix, students with a weighting target (n = 352) | Bias | 38.6, 0.30 | | |
| S2 | Same weights shuffled within student (difficulty), students with a weighting target (n = 352) | Bias | 17.7, 0.20 | | |
| S2 | Same weights shuffled within student (difficulty by strand), students with a weighting target (n = 352) | Bias | 38.8, 0.39 | | |
| S2 | Difficulty-balanced subset, equal weights, students with a weighting target (n = 352) | Bias | 14.3, 0.13 | | |

## 5. Confirmatory and replication quantities over the further draws (informational; the confirmatory record stays as computed)

- S1 H2 quantity, bias r(15) interpolated between N = 10 and 20: printed 0.694; over 100 further draws mean 0.701 (SD 0.006, range 0.687 to 0.718); at or above .70 in 0.55 of draws.
- S1 H1 quantity, bias r(20): printed 0.754; draw mean 0.762 (SD 0.006); below .70 in 0.00 of draws.
- S2 criterion (a), bias r(20) below .70: printed 0.742 (not met); met in 0.00 of further draws.
- S2 criterion (b), full 11-point grid (printed): Brier .70 crossing 31.1 in the printed draw; over 100 further draws mean 29.9, median 29.7, above 31 in 0.26; criterion met (three of four indices above 31) in 0.26 of further draws under the plan rule and 0.26 under the tie convention (printed draw: met under both).
- S2 criterion (b), pre-specified 6-point grid: Brier .70 crossing 35.9 in the printed draw; over 100 further draws mean 31.2, median 31.2, above 31 in 0.54; criterion met (three of four indices above 31) in 0.54 of further draws under the plan rule and 0.54 under the tie convention (printed draw: met under both).

## 6. Files

Per-sample detail, including per-split values, per-draw curves and crossings: r2_E_parts/designA_<sample>.json, designB_<sample>.json, structured_<sample>.json (per-split and per-draw aggregates; no per-student values). Scripts: scripts/r2_E_common.py, r2_E_mc_designA.py, r2_E_mc_designB.py, r2_E_mc_structured.py, r2_H_report.py.
