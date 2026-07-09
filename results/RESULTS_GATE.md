# RESULTS GATE - single source of truth for manuscript numbers

Generated from results/ artifacts. Recompute, never retype.

## Confirmatory
```json
{
 "supported": {
  "H1": false,
  "H2": false,
  "H3": true
 },
 "holm_family": {
  "H3": 0.003,
  "H2": 0.66,
  "H1": 1.0
 }
}
```
- H1 bias r(20) = 0.754 [0.728, 0.777], one-sided upper 0.774; p(r20>=.70) = 1.000 -> H1 refuted-direction
- H2 r(15) points: bias=0.694, brier=0.620, gamma=0.146, auroc=0.225
- H2 inner Holm p: brier=0.000, gamma=0.000, auroc=0.000, bias=0.330
- H3 rho raw 0.620 -> EB 0.659, diff 0.0389 CI ['0.0365', '0.0412'], lower1s 0.0369, p 0.000
- Replication S2: {"a_r20_bias_lt_070": {"value": 0.7420257807209322, "met": false}, "b_minN70_gt_median31_for_3of4": {"minN70": {"bias": 16.235298688721766, "brier": 31.101049019896752, "gamma": null, "auroc": 168.65235482586777}, "hits": {"bias": false, "brier": true, "gamma": true, "auroc": true}, "met": true}, "c_eb_improves": {"diff": 0.05179055199495808, "met": true}}

## Minimum-N decision table inputs (S1)
| Index | A minN70 | CI | A minN80 | CI | B minN70 | B minN80 | T minN70 | T CI | T minN80 | proph70 | r(20) | r(20) CI | r15interp | r15direct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Bias | 14.7 | [13.1, 16.9] | 26.0 | [23.4, 28.4] | 17.7 | 31.7 | 20.4 | [14.5, 26.7] | 39.4 | 15.3 | 0.754 | [0.728, 0.777] | 0.694 | 0.705 |
| Absolute accuracy (Brier) | 21.3 | [18.2, 23.9] | 34.3 | [30.3, 38.4] | 27.9 | 45.0 | 30.4 | [22.8, 33.9] | 56.0 | 21.4 | 0.686 | [0.654, 0.717] | 0.620 | 0.644 |
| Relative accuracy (gamma) | > 200 | - | > 200 | - | > 200 | > 200 | > 200 | - | > 200 | 204.8 | 0.186 | [0.089, 0.272] | 0.146 | 0.180 |
| Relative accuracy (AUROC) | 151.1 | [133.1, 172.0] | > 200 | - | 151.1 | > 200 | > 200 | [196.5, >200] | > 200 | 135.2 | 0.257 | [0.189, 0.321] | 0.225 | 0.271 |

### S1 flag consistency HALF-BASED (addendum 5; released only)
- N=10: top-decile Jaccard 0.151, share of flags retained 0.262
- N=15: top-decile Jaccard 0.205, share of flags retained 0.340
- N=20: top-decile Jaccard 0.250, share of flags retained 0.400
- N=30: top-decile Jaccard 0.302, share of flags retained 0.464

### S1 flag consistency DISJOINT FULL-N (addendum 8; manuscript)
(chance Jaccard 0.053, chance retention 0.100; pool n=1389)
- N=10: Jaccard 0.232, retained 0.376
- N=15: Jaccard 0.290, retained 0.449
- N=20: Jaccard 0.335, retained 0.501
- N=30: Jaccard 0.386, retained 0.557
- N=50: Jaccard 0.462, retained 0.632
- N=100: Jaccard 0.558, retained 0.716
- NULL flags (addendum 9): raw-raw 0.084, EB-vs-raw-criterion 0.036 (chance .053); n=23825
- descriptives (addendum 10): reporters 39842, median 15 / mean 37.8 tagged; mean conf 75.0; acc tagged 0.6762; shares >= 15: 0.502, 18: 0.456, 20: 0.431, 26: 0.369, 33: 0.312, 50: 0.211, 100: 0.094, 200: 0.035

### S1 Design A random curves (r; pct undefined half-based; full-set)
```
index    acc  auroc   bias  brier   conf  gamma
N                                              
5      0.433  0.125  0.420  0.379  0.751 -0.024
10     0.635  0.180  0.611  0.527  0.869  0.091
12     0.693  0.201  0.653  0.584  0.897  0.127
15     0.723  0.271  0.705  0.644  0.912  0.180
20     0.777  0.257  0.754  0.686  0.929  0.186
25     0.815  0.285  0.792  0.735  0.944  0.227
30     0.844  0.342  0.830  0.777  0.954  0.239
40     0.878  0.421  0.863  0.826  0.964  0.340
50     0.901  0.477  0.889  0.853  0.971  0.347
100    0.948  0.617  0.940  0.922  0.985  0.440
200    0.973  0.756  0.970  0.960  0.993  0.585

undefined (half-based):
index  gamma  auroc
N                  
5      0.827  0.660
10     0.494  0.321
12     0.451  0.277
15     0.400  0.240
20     0.312  0.167
25     0.283  0.146
30     0.246  0.123
40     0.203  0.099
50     0.178  0.090
100    0.125  0.061
200    0.088  0.049

undefined (full-set):
index  gamma  auroc
N                  
5      0.482  0.314
10     0.308  0.171
12     0.294  0.153
15     0.239  0.127
20     0.204  0.099
25     0.174  0.081
30     0.158  0.076
40     0.131  0.065
50     0.120  0.063
100    0.091  0.047
200    0.068  0.041
```
### S1 Design B random curves
```
index  auroc   bias  brier  gamma
N                                
5      0.088  0.400  0.300  0.068
10     0.127  0.580  0.460  0.075
20     0.223  0.726  0.631  0.135
50     0.429  0.874  0.822  0.308
100    0.614  0.938  0.912  0.464
200    0.759  0.970  0.959  0.588

n_pool: [29889, 23825, 17156, 8387, 3738, 1389]
```
### S1 Design A temporal curves
```
index  auroc   bias  brier  gamma
N                                
5      0.125  0.387  0.371 -0.055
10     0.204  0.560  0.491  0.143
12     0.150  0.624  0.565  0.045
15     0.255  0.671  0.618  0.143
20     0.223  0.698  0.651  0.133
25     0.216  0.720  0.670  0.174
30     0.322  0.750  0.696  0.207
40     0.260  0.803  0.778  0.129
50     0.373  0.818  0.789  0.238
100    0.456  0.867  0.859  0.225
200    0.643  0.890  0.890  0.377
```

## S1 confound / selection / retest
- N20: bias: r_acc=-0.547, rel 0.754->0.883; brier: r_acc=-0.744, rel 0.686->0.564; gamma: r_acc=-0.124, rel 0.186->0.186; auroc: r_acc=-0.039, rel 0.257->0.257
  top-decile Jaccard raw-vs-partialed bias: 0.553
- N100: bias: r_acc=-0.481, rel 0.940->0.969; brier: r_acc=-0.791, rel 0.922->0.871; gamma: r_acc=-0.119, rel 0.440->0.435; auroc: r_acc=-0.082, rel 0.617->0.616
  top-decile Jaccard raw-vs-partialed bias: 0.626
- within-student tagged-untagged acc: 0.0397 [0.0375, 0.0418], n=27625
- student-level: acc 0.602 vs 0.671 (d=-0.377); vol median 121 vs 79 (d_log=0.511); reporters 39842, non 79129
- propensity (per-SD log odds): acc: -0.424 (z=-64.0); log10_vol: 0.623 (z=90.1); first_month: 0.416 (z=60.0)
- IPW Design B N20: bias: 0.727 -> 0.727; brier: 0.629 -> 0.634; gamma: 0.141 -> 0.150; auroc: 0.222 -> 0.223
- retest n=3796, median n 39/67: bias: 0.623 obs / 0.719 disatt; brier: 0.599 obs / 0.736 disatt; gamma: 0.124 obs / 0.397 disatt; auroc: 0.247 obs / 0.568 disatt
- H3 S1: rho 0.620 -> 0.659 (diff 0.0389, p<=0.000); Jaccard top-decile 0.283 -> 0.316; n=23825
- H3 S1 VALIDATION (addendum 4b): raw 0.620 -> EB-val 0.633 (diff 0.0134 CI ['0.0121', '0.0146'], p<=0.000); flag validation Jaccard raw 0.283 -> EB 0.258
- H3 S1 NULL: raw 0.005, EB-sym 0.134 (ARTIFACT), EB-validation -0.0047 (~0 as expected)
- temporal gap days (median between half median-dates): N=5: 68; N=10: 75; N=12: 78; N=15: 82; N=20: 82; N=25: 82; N=30: 84; N=40: 86; N=50: 86; N=100: 88; N=200: 89
- bias decomposition N20: sb_bias 0.754, sb_conf 0.929, sb_acc 0.777; cross-half r(biasA,biasB) 0.605, r(biasA,confB) 0.351, r(biasA,accB) -0.279; r(partialed bias, conf same-half) 0.880
- bias decomposition N100: sb_bias 0.940, sb_conf 0.985, sb_acc 0.948; cross-half r(biasA,biasB) 0.888, r(biasA,confB) 0.444, r(biasA,accB) -0.412; r(partialed bias, conf same-half) 0.843

## Minimum-N decision table inputs (S2)
| Index | A minN70 | CI | A minN80 | CI | B minN70 | B minN80 | T minN70 | T CI | T minN80 | proph70 | r(20) | r(20) CI | r15interp | r15direct |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Bias | 16.2 | [14.4, 19.2] | 28.7 | [24.2, 33.8] | 18.4 | 32.7 | 20.4 | [14.8, 26.1] | 48.1 | 16.2 | 0.742 | [0.707, 0.771] | 0.679 | 0.684 |
| Absolute accuracy (Brier) | 31.1 | [24.6, 35.5] | 58.8 | [46.0, 69.2] | 32.2 | 49.7 | 48.1 | [28.6, 66.3] | 95.1 | 35.5 | 0.568 | [0.495, 0.623] | 0.531 | 0.543 |
| Relative accuracy (gamma) | > 200 | [175.8, >200] | > 200 | - | > 200 | > 200 | > 200 | [5.0, >200] | > 200 | 945.0 | 0.047 | [-0.092, 0.175] | 0.057 | -0.034 |
| Relative accuracy (AUROC) | 168.7 | [144.6, >200] | > 200 | - | 182.2 | > 200 | > 200 | - | > 200 | 272.5 | 0.146 | [0.036, 0.244] | 0.149 | 0.060 |

### S2 flag consistency HALF-BASED (addendum 5; released only)
- N=10: top-decile Jaccard 0.160, share of flags retained 0.275
- N=15: top-decile Jaccard 0.198, share of flags retained 0.329
- N=20: top-decile Jaccard 0.216, share of flags retained 0.353
- N=30: top-decile Jaccard 0.264, share of flags retained 0.417

### S2 flag consistency DISJOINT FULL-N (addendum 8; manuscript)
(chance Jaccard 0.053, chance retention 0.100; pool n=554)
- N=10: Jaccard 0.227, retained 0.368
- N=15: Jaccard 0.273, retained 0.427
- N=20: Jaccard 0.306, retained 0.467
- N=30: Jaccard 0.381, retained 0.551
- N=50: Jaccard 0.472, retained 0.640
- N=100: Jaccard 0.596, retained 0.745
- NULL flags (addendum 9): raw-raw 0.095, EB-vs-raw-criterion 0.042 (chance .053); n=3230
- descriptives (addendum 10): reporters 4284, median 31 / mean 80.9 tagged; mean conf 73.9; acc tagged 0.6220; shares >= 15: 0.673, 18: 0.631, 20: 0.611, 26: 0.549, 33: 0.490, 50: 0.393, 100: 0.231, 200: 0.129

### S2 Design A random curves (r; pct undefined half-based; full-set)
```
index    acc  auroc   bias  brier   conf  gamma
N                                              
5      0.432 -0.044  0.454  0.330  0.741 -0.104
10     0.578  0.153  0.591  0.478  0.871  0.070
12     0.571  0.176  0.587  0.429  0.871  0.126
15     0.646  0.060  0.684  0.543  0.903 -0.034
20     0.710  0.146  0.742  0.568  0.926  0.047
25     0.748  0.231  0.781  0.664  0.938  0.105
30     0.789  0.154  0.806  0.690  0.952 -0.032
40     0.829  0.323  0.848  0.769  0.961  0.113
50     0.846  0.347  0.874  0.774  0.967  0.252
100    0.920  0.572  0.933  0.884  0.984  0.414
200    0.961  0.742  0.967  0.941  0.992  0.587

undefined (half-based):
index  gamma  auroc
N                  
5      0.784  0.586
10     0.373  0.202
12     0.275  0.138
15     0.235  0.105
20     0.179  0.065
25     0.140  0.049
30     0.114  0.039
40     0.097  0.030
50     0.077  0.025
100    0.050  0.022
200    0.041  0.021

undefined (full-set):
index  gamma  auroc
N                  
5      0.368  0.204
10     0.190  0.081
12     0.123  0.043
15     0.114  0.038
20     0.097  0.029
25     0.083  0.029
30     0.069  0.025
40     0.069  0.025
50     0.054  0.022
100    0.036  0.022
200    0.031  0.020
```
### S2 Design B random curves
```
index  auroc   bias  brier  gamma
N                                
5      0.086  0.379  0.272  0.070
10     0.102  0.556  0.419  0.068
20     0.208  0.720  0.591  0.106
50     0.369  0.869  0.801  0.196
100    0.553  0.931  0.884  0.329
200    0.723  0.965  0.939  0.477

n_pool: [3693, 3230, 2617, 1685, 988, 554]
```
### S2 Design A temporal curves
```
index  auroc   bias  brier  gamma
N                                
5      0.117  0.369  0.312  0.292
10     0.160  0.550  0.363  0.074
12     0.042  0.495  0.338 -0.008
15    -0.048  0.650  0.488 -0.147
20     0.025  0.697  0.477 -0.100
25     0.193  0.734  0.598  0.025
30     0.130  0.747  0.667  0.010
40     0.274  0.785  0.685  0.096
50     0.293  0.803  0.703  0.149
100    0.459  0.860  0.807  0.229
200    0.576  0.888  0.847  0.151
```

## S2 confound / selection / retest
- N20: bias: r_acc=-0.508, rel 0.742->0.891; brier: r_acc=-0.692, rel 0.568->0.444; gamma: r_acc=-0.032, rel 0.047->0.049; auroc: r_acc=0.045, rel 0.146->0.147
  top-decile Jaccard raw-vs-partialed bias: 0.642
- N100: bias: r_acc=-0.400, rel 0.933->0.967; brier: r_acc=-0.739, rel 0.884->0.828; gamma: r_acc=0.182, rel 0.414->0.402; auroc: r_acc=0.087, rel 0.572->0.569
  top-decile Jaccard raw-vs-partialed bias: 0.692
- within-student tagged-untagged acc: 0.0776 [0.0720, 0.0830], n=3539
- student-level: acc 0.541 vs 0.578 (d=-0.205); vol median 259 vs 103 (d_log=0.942); reporters 4284, non 634
- propensity (per-SD log odds): acc: -0.101 (z=-2.1); log10_vol: 0.672 (z=12.6); first_month: -0.599 (z=-14.6)
- IPW Design B N20: bias: 0.715 -> 0.714; brier: 0.595 -> 0.594; gamma: 0.128 -> 0.124; auroc: 0.193 -> 0.193
- retest n=823, median n 63/105: bias: 0.631 obs / 0.692 disatt; brier: 0.589 obs / 0.687 disatt; gamma: 0.091 obs / 0.320 disatt; auroc: 0.276 obs / 0.559 disatt
- H3 S2: rho 0.674 -> 0.725 (diff 0.0518, p<=0.000); Jaccard top-decile 0.322 -> 0.370; n=3230
- H3 S2 VALIDATION (addendum 4b): raw 0.674 -> EB-val 0.693 (diff 0.0194 CI ['0.0160', '0.0227'], p<=0.000); flag validation Jaccard raw 0.322 -> EB 0.302
- H3 S2 NULL: raw -0.023, EB-sym 0.472 (ARTIFACT), EB-validation -0.0170 (~0 as expected)
- temporal gap days (median between half median-dates): N=5: 49; N=10: 53; N=12: 53; N=15: 52; N=20: 55; N=25: 55; N=30: 53; N=40: 55; N=50: 56; N=100: 58; N=200: 56
- bias decomposition N20: sb_bias 0.742, sb_conf 0.926, sb_acc 0.710; cross-half r(biasA,biasB) 0.590, r(biasA,confB) 0.432, r(biasA,accB) -0.219; r(partialed bias, conf same-half) 0.912
- bias decomposition N100: sb_bias 0.933, sb_conf 0.984, sb_acc 0.920; cross-half r(biasA,biasB) 0.875, r(biasA,confB) 0.581, r(biasA,accB) -0.314; r(partialed bias, conf same-half) 0.888
