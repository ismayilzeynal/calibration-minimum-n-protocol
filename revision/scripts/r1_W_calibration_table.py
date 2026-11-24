"""Aggregate calibration table (plan row E3b / D4a): accuracy of confidence-rated
answers at each of the five confidence levels, Sample 1 (t12) and Sample 2 (t34).

Reads the paper's answer-level parquet caches read-only; writes
revision-2026-11/results/r1_W_calibration.json. Tagged = Confidence non-missing
and in [0, 100] (common.tagged). Level = 25 * round(Confidence / 25).
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
OUT = os.path.join(REV, "results", "r1_W_calibration.json")


def table(code):
    df = pd.read_parquet(os.path.join(PAPER, "data", f"resp_{code}.parquet"),
                         columns=["UserId", "IsCorrect", "Confidence"])
    c = df["Confidence"]
    t = df[c.notna() & (c >= 0) & (c <= 100)].copy()
    t["level"] = (np.round(t["Confidence"].astype(float) / 25.0) * 25).astype(int)
    out = {"n_tagged": int(len(t)),
           "n_values_off_grid": int((t["Confidence"].astype(float) % 25 != 0).sum()),
           "levels": {}}
    for lv, g in t.groupby("level"):
        acc = float(g["IsCorrect"].mean())
        n = int(len(g))
        se = float(np.sqrt(acc * (1 - acc) / n))
        out["levels"][str(lv)] = {
            "n_answers": n,
            "share_of_tagged": n / len(t),
            "n_students": int(g["UserId"].nunique()),
            "accuracy": acc,
            "accuracy_minus_confidence": acc - lv / 100.0,
            "answer_level_se_naive": se,
        }
    accs = [out["levels"][k]["accuracy"] for k in sorted(out["levels"], key=int)]
    out["accuracy_monotone_in_level"] = bool(all(a < b for a, b in zip(accs, accs[1:])))
    out["accuracy_range_0_to_100"] = [accs[0], accs[-1]]
    return out


def main():
    res = {
        "what": "accuracy of confidence-rated answers at each confidence level (answer level, pooled over students)",
        "definitions": {
            "tagged": "Confidence non-missing and in [0, 100] (common.tagged)",
            "level": "25 * round(Confidence / 25); the release stores 0, 25, 50, 75, 100",
            "accuracy": "share of answers at that level that are correct (IsCorrect = 1)",
            "accuracy_minus_confidence": "accuracy - level/100; negative = overconfident at that level, on the coded scale",
            "answer_level_se_naive": "binomial SE ignoring clustering within students; descriptive only",
        },
        "S1": table("t12"),
        "S2": table("t34"),
    }
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
