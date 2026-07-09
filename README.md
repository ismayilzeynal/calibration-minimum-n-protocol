# Minimum-Response Reliability Protocol for Confidence-Calibration Indices

Code, pre-registered analysis plan, and derived aggregate results for the
manuscript "How Many Answers Before You Can Call a Student Overconfident?
Reliability Boundaries and Minimum-Response Guidelines for
Confidence-Calibration Indices in Large-Scale Mathematics Practice Data"
(under review). Author identity withheld for double-anonymous review.

## What this is

A transferable estimation protocol answering one question for any learning
platform that collects item-level confidence: how many confidence-tagged
responses does a student need before person-level calibration scores
(bias, Brier, gamma, AUROC) are reliable enough to use?

Headline defaults from the two samples analyzed (Eedi mathematics data,
1.85M confidence-tagged answers; random-split, Designs A and B):

| Index | N for r = .70 | N for r = .80 |
|---|---|---|
| Bias (over/underconfidence) | ~15-18 | ~26-33 |
| Absolute accuracy (Brier) | ~21-32 | ~34-59 |
| Relative accuracy (AUROC) | ~151-182 | > 200 |
| Relative accuracy (gamma) | > 200 | > 200 |

Also in the released aggregates: temporal-split (conservative) thresholds,
disjoint full-N flag-stability curves, permutation-null certifications of
the shrinkage evaluation metrics, selection/IPW checks, and the
accuracy-confound decomposition.

## Data

Row-level data are NOT redistributed here (dataset license). Obtain the
NeurIPS 2020 Education Challenge data ("Instructions and Guide for
Diagnostic Questions", Wang et al. 2020, arXiv:2007.12061) from the
official distribution, place `data.zip` where `scripts/common.py` expects
it (see ZIP_PATH), and run the pipeline.

## Reproduce

```
pip install -r scripts/requirements.txt
python scripts/reliability.py t12      # sample 1 curves + minimum-N
python scripts/reliability.py t34      # sample 2 (replication)
python scripts/rev_h3null.py t12       # permutation-null certification
python scripts/rev_h3null.py t34
python scripts/h3_shrinkage.py t12     # shrinkage + validation metric
python scripts/h3_shrinkage.py t34
python scripts/rev_extras.py t12       # temporal gaps + decomposition
python scripts/rev_extras.py t34
python scripts/confound_selection.py t12
python scripts/confound_selection.py t34
python scripts/rev_round2.py t12       # disjoint flags + null flags + descriptives
python scripts/rev_round2.py t34
python scripts/confirmatory.py
python scripts/results_gate.py
python scripts/figures.py
```

Master seed 42; 50 Monte Carlo splits; 1,000 bootstrap draws. The frozen
pre-registration is `analysis_plan.md`; `results/` ships the derived
aggregates behind every number in the manuscript.

## License

Code: MIT (see LICENSE). Derived aggregate statistics: CC BY 4.0.
Source data remain governed by the original Eedi/NeurIPS challenge terms.
