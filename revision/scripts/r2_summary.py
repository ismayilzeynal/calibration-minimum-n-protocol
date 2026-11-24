"""Revision, gate fix round 1: small derived quantities quoted in the revised
manuscript that are not single keys of another results file.

  - Design A pool students whose rated answers are all correct: their unrated
    answers (count, pooled accuracy, students perfect on unrated answers);
  - S1 strand shares of all answers outside the four main strands (combined
    labels, minor single strands, no strand);
  - the inverse-probability-weighting change in Design B bias reliability at
    N = 20 (results/confound_selection_t12.json, its own subsample);
  - maxima of the Monte Carlo SE over the 50 splits (results/r2_A.json);
  - ranges of the level of bias across strands and terms (results/r2_C.json).

Usage: python -B r2_summary.py
Writes: revision-2026-11/results/r2_summary.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
PAPER = os.path.dirname(REV)
sys.path.insert(0, HERE)
import r1_C_common as RC  # noqa: E402

MAIN = ["Number", "Algebra", "Geometry and Measure", "Data and Statistics"]
out = {"_definitions": {
    "all_correct_unrated": "Design A pool (>= 200 rated answers) students whose rated answers are all correct; their answers without a rating",
    "strand_other_S1": "share of all S1 answers whose question carries more than one level-1 strand (combined), one level-1 strand outside the four main strands (minor), or none",
    "ipw_change": "IPW minus unweighted Spearman-Brown reliability, Design B pool N = 20, on the subsample of confound_selection.py (seed [42, 81, sid])",
    "mc_se_max": "maximum Monte Carlo SE of the mean over 50 splits",
    "level_ranges": "min and max of answer-level mean bias across the four main strands (or terms), all rated answers"}}

for s, lab in (("t12", "S1"), ("t34", "S2")):
    df = pd.read_parquet(os.path.join(PAPER, "data", f"resp_{s}.parquet"))
    c = df["Confidence"]
    tag = (c.notna() & (c >= 0) & (c <= 100)).to_numpy()
    g = df[tag].groupby("UserId")["IsCorrect"].agg(["size", "mean"])
    ids = g[(g["size"] >= 200) & (g["mean"] == 1)].index
    un = df[(~tag) & df.UserId.isin(ids)].groupby("UserId")["IsCorrect"]
    out.setdefault("all_correct_unrated", {})[lab] = {
        "n_students_all_rated_correct": int(len(ids)),
        "n_with_unrated": int(un.size().shape[0]),
        "unrated_answers": int(un.size().sum()),
        "pooled_unrated_accuracy": float(un.sum().sum() / un.size().sum()),
        "students_perfect_on_unrated": int((un.mean() == 1).sum())}
    if s == "t12":
        qk = RC.question_keys(s)
        m = df[["QuestionId"]].merge(qk, on="QuestionId", how="left")
        n1 = m["n_strand"].fillna(0).to_numpy()
        lab1 = m["strand"].fillna("none").to_numpy()
        out["strand_other_S1"] = {
            "combined": float(np.mean(n1 > 1)),
            "minor_single": float(np.mean((n1 == 1) & ~np.isin(lab1, MAIN))),
            "none": float(np.mean(n1 == 0)),
            "main_four_total": float(np.mean((n1 == 1) & np.isin(lab1, MAIN)))}
    del df

cs = json.load(open(os.path.join(PAPER, "results", "confound_selection_t12.json")))["selection"]["ipw_designB_N20"]
out["ipw_change_S1"] = {k: cs[k]["ipw"] - cs[k]["unweighted"] for k in cs}

A = json.load(open(os.path.join(REV, "results", "r2_A.json")))
mc = {}
for lab in ("S1", "S2"):
    I = A[lab]["indices"]
    mc[lab] = {
        "resolution_N100_200": max(I[k]["mc_se_over_splits"][n] for k in ("gamma", "auroc", "auroc_pr") for n in ("100", "200")),
        "bias_brier_N_ge_20": max(I[k]["mc_se_over_splits"][n] for k in ("bias", "brier")
                                  for n in ("20", "25", "30", "40", "50", "100", "200")),
        "bias_brier_all_N": max(I[k]["mc_se_over_splits"][n] for k in ("bias", "brier") for n in map(str, A[lab]["grid"]))}
out["mc_se_max"] = mc

C = json.load(open(os.path.join(REV, "results", "r2_C.json")))
lr = {}
for lab in ("S1", "S2"):
    L = C[lab]["levels"]["all_rated"]
    st = [v["mean_bias"] for k, v in L["strand"].items() if k in MAIN]
    tm = [v["mean_bias"] for v in L["term"].values()]
    lr[lab] = {"strand_min": min(st), "strand_max": max(st), "term_min": min(tm), "term_max": max(tm)}
out["level_ranges"] = lr

with open(os.path.join(REV, "results", "r2_summary.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
