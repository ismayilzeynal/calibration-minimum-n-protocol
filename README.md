# Minimum-Response Reliability Protocol for Confidence-Calibration Indices

Code, the pre-specified analysis plan with its dated addenda, and derived
aggregate results for the manuscript "How Many Rated Answers Does a
Confidence-Calibration Score Need? Conditional Minimum-Response Estimates and
a Reusable Protocol from Mathematics Practice Data" (under review, revised
version). Author identity withheld for double-anonymous review.

Version 1.0.0 of this repository accompanied the original submission. The
`revision/` folder accompanies the revised manuscript; the original
`scripts/` and `results/` are unchanged.

## What this is

A reusable estimation protocol for one question that any learning platform
collecting item-level confidence ratings can ask: how many rated answers does
a student need before a person-level calibration score (signed bias, Brier
score, gamma, AUROC) is reliable enough for a given use?

## Headline estimates (conditional on these data)

Number of rated answers per student at which split-half reliability
(Spearman-Brown corrected, mean over 50 random splits) first reaches .70 or
.80. Random-half ranges span both samples (Eedi mathematics practice data,
1.85 million rated answers) and both sampling designs (A: the same students
at every N; B: every student with at least N rated answers). Temporal halves
(a student's earlier versus later answers) are Design A, both samples. AUROC
under the plan rule (undefined whenever gamma is undefined) is primary; the
Mann-Whitney tie convention, which gives .50 to a student whose answers are
mixed in correctness but all carry the same rating, is a sensitivity
analysis.

| Index | .70, random halves | .80, random halves | .70, temporal halves | .80, temporal halves |
|---|---|---|---|---|
| Signed bias (over- or underconfidence) | about 15-18 | about 26-33 | about 20 | about 39-48 |
| Brier score | about 21-32 | about 34-59 | about 30-48 | about 56-95 |
| AUROC, plan rule | about 171-198 | > 200 | > 200 | > 200 |
| AUROC, tie convention | about 151-182 | > 200 | > 200 | > 200 |
| Gamma | > 200 | > 200 | > 200 | > 200 |

Sources: `revision/results/baseline/results/minN_{t12,t34}.json` (bias,
Brier, gamma, AUROC tie convention) and
`revision/results/r1_A/prereg/results/minN_{t12,t34}.json` (AUROC plan
rule). The estimates hold for this platform, subject, rating instrument and
student pool; the protocol is meant to be rerun locally.

## Layout

- `analysis_plan.md`: the pre-specified analysis plan, frozen on 2026-07-09
  before any inferential analysis, with dated, append-only addenda.
  Addenda 1 and 2 were part of version 1.0.0; Addendum 3 (2026-11-24)
  records the analyses added after the editor's decision; Addendum 4
  (2026-11-24) records the Monte Carlo and pre-specified-grid checks, the
  crossed-subset comparison with ICC(1,k), the AUROC minus gamma
  differences, and two further exploratory analyses (where rating occurs;
  the flag table at 60 answers). The first 15,027 bytes of the file are
  byte-identical to version 1.0.0, and each addendum was appended without
  changing the text before it. The title
  line keeps the term "pre-registered" used at the original submission;
  because the plan was first made public with that submission, the revised
  manuscript calls it pre-specified.
- `scripts/`, `results/`: original pipeline and aggregates (version 1.0.0).
- `revision/scripts/`: revision analyses (exploratory, recorded in
  Addenda 3 and 4) and reruns of the original code under the plan's AUROC
  rule.
- `revision/results/`: aggregate outputs (the summary files `r1_*.json`
  and `r2_*.json` carry a definitions block; CSV reliability curves;
  figures). `figures_r2/` holds Figures 1 to 3 of the revised manuscript;
  `figures/` holds Figure S1 of the Supplementary Notes and earlier
  versions of Figures 1 to 3. `r2_E_parts/` holds the Monte Carlo
  components per sample and design (values per random split and per
  further subsample draw, each already aggregated over students).
- `revision/Supplementary_Notes.md`: Supplementary Notes S1 to S4 of the
  revised manuscript.

Script docstrings call the revision folder `revision-2026-11/`; here it is
`revision/`. Every revision script derives its paths from its own location,
with the repository root in the role of the paper folder, so no path needs
editing. Analysis groups: A (AUROC rule and the reruns), B (variance
components, coefficient alpha, standard error of measurement), C (rated
versus unrated answers, reweighting, halves by content, difficulty and term),
D (what the bias score is made of; flag stability), and `r2_*` (bootstrap
intervals over all draws, guess-corrected coding, pool sensitivity, level of
bias, flag rules, figures). Within `r2_*`, the Addendum 4 checks are
`r2_E_*` (Monte Carlo error from the random splits and from the single
subsample draw per student), `r2_F_grid.py` (crossings on the
pre-specified grid, the crossed-subset comparison with ICC(1,k), AUROC
minus gamma) and `r2_H_report.py` (assembles them); `r2_G_figures.py` and
`r3_figures.py` produce the final Figures 1 and 3. The two exploratory
analyses of Addendum 4 come from `r1_C_selection.py` (where rating occurs)
and `r1_D_construct_flags.py` (the flag table at 60 answers).

Not included: row-level data and anything derived per student. The
per-student cache arrays (`data/cache/designA_*.npz`: split-half values per
student and split) are regenerated by the run order below.

## Data

Row-level data are not redistributed here (dataset terms). Obtain the
NeurIPS 2020 Education Challenge data ("Instructions and Guide for
Diagnostic Questions", Wang et al. 2020, arXiv:2007.12061) from the
official distribution and place `data.zip` at `data/raw/eedi/data.zip`, as
`scripts/common.py` expects.

## Reproduce: original pipeline

Run from the repository root:

```
pip install -r scripts/requirements.txt
python scripts/reliability.py t12      # sample 1 curves + minimum N
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

Master seed 42; 50 random splits; 1,000 bootstrap draws.

## Reproduce: revision

Run from the repository root, after `scripts/reliability.py` has built
`data/resp_t12.parquet`, `data/resp_t34.parquet` and `data/cache/`. The
revision scripts read the shipped `results/` files of version 1.0.0.

```
pip install -r revision/scripts/requirements.txt

# baseline: rerun the original curves into revision/results/baseline
python revision/scripts/extract_metadata.py data/raw/eedi/data.zip
python -B revision/scripts/repro_baseline.py t12
python -B revision/scripts/repro_baseline.py t34
python -B revision/scripts/compare_baseline.py
python -B revision/scripts/auroc_tie_diagnostic.py

# A: AUROC under the plan rule and under the tie convention
python -B revision/scripts/r1_A_indices_check.py
python -B revision/scripts/r1_A_prereg_run.py reliability t12
python -B revision/scripts/r1_A_prereg_run.py reliability t34
python -B revision/scripts/r1_A_prereg_run.py confound t12
python -B revision/scripts/r1_A_prereg_run.py confound t34
python -B revision/scripts/r1_A_prereg_run.py confirmatory
python -B revision/scripts/r1_A_mw_run.py
python -B revision/scripts/r1_A_gate.py
python -B revision/scripts/r1_A_light.py t12
python -B revision/scripts/r1_A_light.py t34
python -B revision/scripts/r1_A_summary.py

# B: variance components, coefficient alpha, SEM
python -B revision/scripts/r1_B_vc.py t12
python -B revision/scripts/r1_B_vc.py t34
python -B revision/scripts/r1_B_alpha.py t34
python -B revision/scripts/r1_B_alpha.py t12
python -B revision/scripts/r1_B_assemble.py

# C: selection, reweighting, structured halves
python -B revision/scripts/r1_C_selection.py t12
python -B revision/scripts/r1_C_selection.py t34
python -B revision/scripts/r1_C_splits.py t12
python -B revision/scripts/r1_C_splits.py t34
python -B revision/scripts/r1_C_assemble.py
python -B revision/scripts/r1_C_report.py

# D: decomposition and flag stability (writes revision/tmp/r1_D_*.json)
python -B revision/scripts/r1_D_construct_flags.py t12
python -B revision/scripts/r1_D_construct_flags.py t34
python -B revision/scripts/r1_D_report.py

# accuracy at each confidence level
python -B revision/scripts/r1_W_calibration_table.py

# r2: intervals over all draws, coding and pool sensitivity, levels, flags, figures
python -B revision/scripts/r2_A_rerun.py
python -B revision/scripts/r2_C_levels.py
python -B revision/scripts/r2_D_flags.py
python -B revision/scripts/r2_summary.py
python -B revision/scripts/r2_figures.py

# Addendum 4: Monte Carlo error, pre-specified grid, crossed subset, AUROC minus gamma
python -B revision/scripts/r2_E_mc_designA.py t12
python -B revision/scripts/r2_E_mc_designA.py t34
python -B revision/scripts/r2_E_mc_designB.py t12
python -B revision/scripts/r2_E_mc_designB.py t34
python -B revision/scripts/r2_E_mc_structured.py t12
python -B revision/scripts/r2_E_mc_structured.py t34
python -B revision/scripts/r2_F_grid.py        # reads r2_E_parts/designA_*.json
python -B revision/scripts/r2_H_report.py

# final Figures 1 and 3 (revision/results/figures_r2)
python -B revision/scripts/r2_G_figures.py
python -B revision/scripts/r3_figures.py
```

`r2_G_figures.py` must run before `r3_figures.py`, and neither
`r2_G_figures.py` nor `r2_figures.py` may be rerun after `r3_figures.py`:
they restore the earlier y-axis limits (and `r2_G_figures.py` copies the
earlier Figure 3 back).

The revision reruns used Python 3.12.10 with
`revision/scripts/requirements.txt` (numpy 2.4.6, scipy 1.18.1; the original
pins were numpy 2.2.6 and scipy 1.16.1). Under these versions the original
minimum-N files reproduce exactly and the reliability curves to at most
3.0e-15 (`revision/results/baseline/compare_baseline.json`).

## License

Code: MIT (see LICENSE). Derived aggregate statistics: CC BY 4.0.
Source data remain governed by the original Eedi/NeurIPS challenge terms.
