"""Diagnostic for the AUROC constant-confidence rule (editor point / stream A).

indices.py sets AUROC = (C + 0.5*ties)/(n0*n1) whenever n0*n1 > 0, so a
student whose confidence is constant but whose correctness varies gets
AUROC = 0.5 exactly (all correct-incorrect pairs tied), while gamma is NaN
for the same student. This script counts those cases in the ORIGINAL
Design A cache (read-only) and writes a JSON under revision results.

Usage: python auroc_tie_diagnostic.py
"""
import json
import os

import numpy as np

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUT = os.path.join(REV, "results", "baseline", "auroc_tie_diagnostic.json")


def block(g, a):
    g = np.asarray(g, float)
    a = np.asarray(a, float)
    tot = g.size
    g_nan = np.isnan(g)
    a_nan = np.isnan(a)
    imputed_half = (~a_nan) & g_nan & (a == 0.5)
    return {
        "cells": int(tot),
        "share_gamma_undefined": float(g_nan.mean()),
        "share_auroc_undefined": float(a_nan.mean()),
        "share_auroc_exactly_0.5_where_gamma_undefined": float(imputed_half.mean()),
        "share_gamma_nan_auroc_defined_not_0.5": float(((~a_nan) & g_nan & (a != 0.5)).mean()),
        "share_auroc_nan_gamma_defined": float((a_nan & ~g_nan).mean()),
    }


def main():
    out = {"rule_source": "scripts/indices.py lines 45-46",
           "note": "cells = students (full set) or students x 50 splits (halves); "
                   "Design A pools from the submitted run's cache"}
    for s in ("t12", "t34"):
        c = np.load(os.path.join(PAPER, "data", "cache", f"designA_{s}.npz"))
        blk = {}
        for N in (20, 100):
            blk[f"full{N}"] = block(c[f"full{N}_gamma"], c[f"full{N}_auroc"])
            blk[f"halfA{N}"] = block(c[f"A{N}_gamma"], c[f"A{N}_auroc"])
            blk[f"halfB{N}"] = block(c[f"B{N}_gamma"], c[f"B{N}_auroc"])
        out[s] = blk
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
