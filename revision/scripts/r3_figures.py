"""Figures 1 and 3, fix round 5: the y-axis floor is lowered so that no
plotted point or confidence band is clipped (Figure 1: the gamma and AUROC
bands at N = 5 reach -.25; Figure 3: the S2 gamma point at N = 5 is -.104).
Data, palette, markers, labels and layout are those of r2_G_figures.py
(Figure 1) and r2_figures.py (Figure 3); data from results/r2_A.json.
Figure 2 is not touched. This script supersedes r2_G_figures.py for
Figures 1 and 3; do not rerun r2_G_figures.py or r2_figures.py afterwards
(they restore the old axis limits).

Checks written to results/figures_r2/figures_r2_check.json:
  - every plotted series equals the r2_A.json values (max abs diff);
  - every plotted point and band lies inside the y-limits;
  - Figure 2 is unchanged (sha256 against results/figures).

Usage: python -B r3_figures.py
Writes: revision-2026-11/results/figures_r2/fig1_curves_s1.png,
        fig3_replication.png, figures_r2_check.json
"""
import hashlib
import json
import os
import sys

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import r2_figures as F  # noqa: E402  (palette, labels, data loaders; its main is not run)
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REV = F.REV
OLD = os.path.join(REV, "results", "figures")
NEW = os.path.join(REV, "results", "figures_r2")
GRID = F.GRID
YLIM1 = (-0.3, 1.0)
YLIM3 = (-0.2, 1.22)
YTICKS = [-0.2, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


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
    ax.set_ylim(*YLIM1)
    ax.set_yticks(YTICKS)
    ax.set_xlim(right=230)
    ax.legend(fontsize=8, loc="lower right", framealpha=1.0, edgecolor="none")
    fig.tight_layout()
    path = os.path.join(NEW, "fig1_curves_s1.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    return path, check("S1", plotted, YLIM1)


def fig3():
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    plotted = {}
    for s, ls, mfc, nm in (("S1", "-", None, "Sample 1"), ("S2", "--", "white", "Sample 2")):
        for k in ("bias", "gamma"):
            lbl = ("Signed bias" if k == "bias" else "Gamma (resolution)") + f" - {nm}"
            y = F.series(s, k, "curve_random_point")
            ax.plot(GRID, y, ls + "o", color=F.COL[k], ms=4.5, lw=1.8, mfc=mfc, label=lbl)
            plotted[(s, k)] = (y, None, None)
    ax.axhline(0.70, color="#555555", ls=":", lw=1)
    F.logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(*YLIM3)
    ax.set_yticks(YTICKS)
    ax.legend(fontsize=8, ncol=2, loc="upper left", framealpha=1.0, edgecolor="none")
    fig.tight_layout()
    path = os.path.join(NEW, "fig3_replication.png")
    fig.savefig(path, dpi=300)
    plt.close(fig)
    out = {}
    for s in ("S1", "S2"):
        out[s] = check(s, {k: v for (ss, k), v in plotted.items() if ss == s}, YLIM3)
    return path, out


def check(s, plotted, ylim):
    J = F.J[s]["indices"]
    diff = 0.0
    lo_min, hi_max = np.inf, -np.inf
    for k, (yy, lo, hi) in plotted.items():
        ref = np.array([J[k]["curve_random_point"][str(n)] for n in GRID])
        diff = max(diff, float(np.max(np.abs(yy - ref))))
        lo_min, hi_max = min(lo_min, float(np.min(yy))), max(hi_max, float(np.max(yy)))
        if lo is not None:
            rl = np.array([J[k]["curve_random_ci"][str(n)][0] for n in GRID])
            rh = np.array([J[k]["curve_random_ci"][str(n)][1] for n in GRID])
            diff = max(diff, float(np.max(np.abs(lo - rl))), float(np.max(np.abs(hi - rh))))
            lo_min, hi_max = min(lo_min, float(np.min(lo))), max(hi_max, float(np.max(hi)))
    return {"max_abs_diff_vs_r2_A": diff, "lowest_plotted": lo_min, "highest_plotted": hi_max,
            "ylim": list(ylim), "inside_ylim": bool(ylim[0] < lo_min and hi_max < ylim[1])}


def sha(p):
    with open(p, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


if __name__ == "__main__":
    p1, c1 = fig1()
    p3, c3 = fig3()
    f2 = "fig2_undefined_temporal.png"
    chk = {"_definitions": {
               "written_by": "scripts/r3_figures.py (fix round 5); supersedes the check written by r2_G_figures.py",
               "max_abs_diff_vs_r2_A": "largest absolute difference between any plotted value (points and CI bands) and results/r2_A.json",
               "inside_ylim": "every plotted point and band lies strictly inside the y-limits (nothing clipped)",
               "fig2_byte_identical": "sha256 of figures_r2/fig2 equals that of results/figures/fig2 (Figure 2 unchanged)"},
           "fig1": c1,
           "fig3": c3,
           "fig2_byte_identical": sha(os.path.join(OLD, f2)) == sha(os.path.join(NEW, f2))}
    ok = c1["inside_ylim"] and all(v["inside_ylim"] for v in c3.values()) and chk["fig2_byte_identical"] \
        and c1["max_abs_diff_vs_r2_A"] == 0 and all(v["max_abs_diff_vs_r2_A"] == 0 for v in c3.values())
    chk["all_pass"] = bool(ok)
    with open(os.path.join(NEW, "figures_r2_check.json"), "w", encoding="utf-8", newline="\n") as fh:
        json.dump(chk, fh, indent=1)
    print(json.dumps(chk, indent=1))
    sys.exit(0 if ok else 1)
