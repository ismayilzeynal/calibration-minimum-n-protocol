"""Revised Figure 1 (revision round 2): the label of the median-reporter line
is moved to the empty band above the curves; data, styles and layout are
those of r2_figures.py (data from results/r2_A.json). Figures 2 and 3 are
copied unchanged.

Checks written to results/figures_r2/figures_r2_check.json:
  - every plotted series equals the r2_A.json values (max abs diff);
  - the PNG differs from results/figures/fig1_curves_s1.png only inside the
    old and new label areas (bounding boxes of the differing pixels);
  - Figures 2 and 3 are byte-identical copies (sha256).

Superseded for Figures 1 and 3 by r3_figures.py (fix round 5, lower y-limits); do not
rerun this script after r3_figures.py, because it restores the old y-limits and copies
the old Figure 3 back.

Usage: python -B r2_G_figures.py
Writes: revision-2026-11/results/figures_r2/
"""
import hashlib
import json
import os
import shutil
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_figures as F  # noqa: E402  (palette, labels, data loaders; its main is not run)
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REV = F.REV
OLD = os.path.join(REV, "results", "figures")
NEW = os.path.join(REV, "results", "figures_r2")
GRID = F.GRID


def fig1():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    plotted = {}
    for k in ("bias", "brier", "gamma", "auroc_pr"):
        y = F.series("S1", k, "curve_random_point")
        ax.plot(GRID, y, "-o", color=F.COL[k], label=F.LAB[k], ms=4.5, lw=1.8)
        lo, hi = F.band("S1", k, "curve_random_ci")
        ax.fill_between(GRID, lo, hi, color=F.COL[k], alpha=0.15, lw=0)
        plotted[k] = (y, lo, hi)
    y = F.series("S1", "auroc", "curve_random_point")
    ax.plot(GRID, y, "--", color=F.COL["auroc"], lw=1.2, label=F.LAB["auroc"] + " (sensitivity)")
    plotted["auroc"] = (y, None, None)
    for thr, txt in ((0.70, ".70"), (0.80, ".80")):
        ax.axhline(thr, color="#555555", ls=":", lw=1)
        ax.text(214, thr, txt, fontsize=8, color="#555555", va="center")
    ax.axvline(15, color="#888888", ls="--", lw=1)
    ax.text(14.3, 0.985, "median reporter\n(15 rated answers)", fontsize=8, color="#666666",
            ha="right", va="top")
    F.logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlim(right=230)
    ax.legend(fontsize=8, loc="lower right", framealpha=1.0, edgecolor="none")
    fig.tight_layout()
    path = os.path.join(NEW, "fig1_curves_s1.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    J = F.J["S1"]["indices"]
    diff = 0.0
    for k, (yy, lo, hi) in plotted.items():
        ref = np.array([J[k]["curve_random_point"][str(n)] for n in GRID])
        diff = max(diff, float(np.max(np.abs(yy - ref))))
        if lo is not None:
            rl = np.array([J[k]["curve_random_ci"][str(n)][0] for n in GRID])
            rh = np.array([J[k]["curve_random_ci"][str(n)][1] for n in GRID])
            diff = max(diff, float(np.max(np.abs(lo - rl))), float(np.max(np.abs(hi - rh))))
    return path, diff


def pixel_diff(a_path, b_path):
    a = mpimg.imread(a_path)
    b = mpimg.imread(b_path)
    if a.shape != b.shape:
        return {"same_shape": False}
    d = np.any(np.abs(a - b) > 1e-6, axis=2)
    out = {"same_shape": True, "shape": list(a.shape), "n_pixels_differing": int(d.sum())}
    cols = np.where(d.any(axis=0))[0]
    # split the differing columns into separate horizontal groups
    groups = []
    if len(cols):
        start = prev = cols[0]
        for c in cols[1:]:
            if c > prev + 5:
                groups.append((start, prev))
                start = c
            prev = c
        groups.append((start, prev))
    boxes = []
    for c0, c1 in groups:
        rows = np.where(d[:, c0:c1 + 1].any(axis=1))[0]
        boxes.append({"x0": int(c0), "x1": int(c1), "y0": int(rows.min()), "y1": int(rows.max())})
    out["boxes_of_differing_pixels"] = boxes
    return out


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


if __name__ == "__main__":
    os.makedirs(NEW, exist_ok=True)
    p1, d1 = fig1()
    chk = {"_definitions": {
               "fig1_max_abs_diff_plotted_vs_r2_A": "largest absolute difference between any plotted value (curves and CI bands) and results/r2_A.json",
               "fig1_pixel_diff_vs_round1": "pixels differing from results/figures/fig1_curves_s1.png, grouped into boxes (pixel columns x0-x1, rows y0-y1 from the top); expected only at the old and the new label",
               "byte_identical": "sha256 of the copy equals sha256 of results/figures/<file>"},
           "fig1_max_abs_diff_plotted_vs_r2_A": d1,
           "fig1_pixel_diff_vs_round1": pixel_diff(os.path.join(OLD, "fig1_curves_s1.png"), p1)}
    for name in ("fig2_undefined_temporal.png", "fig3_replication.png"):
        shutil.copyfile(os.path.join(OLD, name), os.path.join(NEW, name))
        chk[f"{name}_byte_identical"] = sha(os.path.join(OLD, name)) == sha(os.path.join(NEW, name))
    with open(os.path.join(NEW, "figures_r2_check.json"), "w", encoding="utf-8") as fh:
        json.dump(chk, fh, indent=1)
    print(json.dumps(chk, indent=1))
