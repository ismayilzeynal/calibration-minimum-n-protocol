"""Revised manuscript figures (gate fix round 1). Same palette and layout as
the submitted scripts/figures.py; data from results/r2_A.json (Design A
curves regenerated with the submitted seeds, both AUROC rules).

  fig1_curves_s1.png        Figure 1: S1 Design A curves; AUROC under the plan
                            rule (primary) with CI band, tie convention dashed
  fig2_undefined_temporal.png  Figure 2: (a) students without a usable
                            resolution value; (b) bias, random vs temporal halves
  fig3_replication.png      Figure 3: bias and gamma, S1 vs S2
  figS1_shrinkage.png       Supplementary Figure S1: shrinkage (H3), revised labels

Usage: python -B r2_figures.py
Writes: revision-2026-11/results/figures/
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
REV = os.path.dirname(HERE)
FIGDIR = os.path.join(REV, "results", "figures")
COL = {"bias": "#2a78d6", "brier": "#1baf7a", "gamma": "#eda100", "auroc_pr": "#a34fc9", "auroc": "#a34fc9"}
LAB = {"bias": "Signed bias", "brier": "Brier score", "gamma": "Gamma (resolution)",
       "auroc_pr": "AUROC (resolution), plan rule", "auroc": "AUROC, tie convention"}
TICKS = [5, 10, 20, 50, 100, 200]
J = json.load(open(os.path.join(REV, "results", "r2_A.json"), encoding="utf-8"))
GRID = J["S1"]["grid"]


def series(s, k, field):
    blk = J[s]["indices"][k][field]
    return np.array([blk[str(n)] for n in GRID], dtype=float)


def band(s, k, field):
    blk = J[s]["indices"][k][field]
    return (np.array([blk[str(n)][0] for n in GRID]), np.array([blk[str(n)][1] for n in GRID]))


def logx(ax, label):
    ax.set_xscale("log")
    ax.set_xticks(TICKS)
    ax.set_xticklabels([str(n) for n in TICKS])
    ax.set_xlabel(label)
    ax.spines[["top", "right"]].set_visible(False)


def fig1():
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for k in ("bias", "brier", "gamma", "auroc_pr"):
        ax.plot(GRID, series("S1", k, "curve_random_point"), "-o", color=COL[k], label=LAB[k], ms=4.5, lw=1.8)
        lo, hi = band("S1", k, "curve_random_ci")
        ax.fill_between(GRID, lo, hi, color=COL[k], alpha=0.15, lw=0)
    ax.plot(GRID, series("S1", "auroc", "curve_random_point"), "--", color=COL["auroc"], lw=1.2,
            label=LAB["auroc"] + " (sensitivity)")
    for thr, txt in ((0.70, ".70"), (0.80, ".80")):
        ax.axhline(thr, color="#555555", ls=":", lw=1)
        ax.text(214, thr, txt, fontsize=8, color="#555555", va="center")
    ax.axvline(15, color="#888888", ls="--", lw=1)
    ax.text(15.6, 0.03, "median reporter (15 rated answers)", rotation=90, fontsize=8, color="#666666", va="bottom")
    logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlim(right=230)
    ax.legend(fontsize=8, loc="lower right", framealpha=1.0, edgecolor="none")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_curves_s1.png"), dpi=300)
    plt.close(fig)


def fig2():
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    g = 100 * series("S1", "gamma", "undefined_full_set")
    a = 100 * series("S1", "auroc", "undefined_full_set")
    ax.plot(GRID, g, "-o", color=COL["gamma"], ms=4.5, lw=1.8,
            label="Gamma, and AUROC under the plan rule: undefined")
    ax.plot(GRID, a, "--o", color=COL["auroc"], ms=4, lw=1.4, mfc="white",
            label="AUROC, tie convention: undefined")
    ax.fill_between(GRID, a, g, color=COL["auroc"], alpha=0.15, lw=0,
                    label="AUROC = .50 placeholder (identical ratings)")
    logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Students without a usable value (%)")
    ax.set_title("(a) No usable resolution value", fontsize=10)
    ax.legend(frameon=False, fontsize=7.5)
    ax = axes[1]
    for sp, ls, mfc, lab in (("curve_random_point", "-", None, "Random halves"),
                             ("curve_temporal_point", "--", "white", "Earlier versus later halves")):
        ax.plot(GRID, series("S1", "bias", sp), ls + "o", color=COL["bias"], ms=4.5, lw=1.8, mfc=mfc, label=lab)
    ax.axhline(0.70, color="#555555", ls=":", lw=1)
    logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Split-half reliability")
    ax.set_title("(b) Signed bias: random versus temporal halves", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_undefined_temporal.png"), dpi=300)
    plt.close(fig)


def fig3():
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for s, ls, mfc, nm in (("S1", "-", None, "Sample 1"), ("S2", "--", "white", "Sample 2")):
        for k in ("bias", "gamma"):
            lbl = ("Signed bias" if k == "bias" else "Gamma (resolution)") + f" - {nm}"
            ax.plot(GRID, series(s, k, "curve_random_point"), ls + "o", color=COL[k], ms=4.5, lw=1.8,
                    mfc=mfc, label=lbl)
    ax.axhline(0.70, color="#555555", ls=":", lw=1)
    logx(ax, "Number of rated answers per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(-0.1, 1.22)
    ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.legend(fontsize=8, ncol=2, loc="upper left", framealpha=1.0, edgecolor="none")
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig3_replication.png"), dpi=300)
    plt.close(fig)


def figS1():
    """Supplementary Figure S1: the submitted shrinkage figure with revised labels
    (data: results/h3_shrinkage_*.json and results/h3_null_*.json of the paper)."""
    pres = os.path.join(os.path.dirname(REV), "results")
    vals = {s: json.load(open(os.path.join(pres, f"h3_shrinkage_{s}.json"))) for s in ("t12", "t34")}
    nulls = {s: json.load(open(os.path.join(pres, f"h3_null_{s}.json"))) for s in ("t12", "t34")}
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xs = np.arange(2)
    w = 0.26
    raw = [vals[s]["rho_raw"] for s in ("t12", "t34")]
    eb = [vals[s]["rho_eb"] for s in ("t12", "t34")]
    ebv = [vals[s]["rho_eb_validation"] for s in ("t12", "t34")]
    nul = [nulls[s]["null_rho_eb"] for s in ("t12", "t34")]
    ax.bar(xs - w, raw, w, color="#9dbfe3", label="Raw bias")
    ax.bar(xs, eb, w, color=COL["bias"], label="Shrunken, pre-registered metric (inflated; see null)")
    ax.bar(xs + w, ebv, w, color="#1baf7a", label="Shrunken, validation metric (other half raw)")
    for i in range(2):
        ax.text(i - w, raw[i] + 0.012, f"{raw[i]:.3f}", ha="center", fontsize=8)
        ax.text(i, eb[i] + 0.012, f"{eb[i]:.3f}", ha="center", fontsize=8)
        ax.text(i + w, ebv[i] + 0.012, f"{ebv[i]:.3f}", ha="center", fontsize=8)
        ax.hlines(nul[i], i - w / 2, i + w / 2, color="#c0392b", ls="--", lw=1.6)
        ax.text(i, nul[i] + 0.015, f"value under the null {nul[i]:.2f}", ha="center", fontsize=7.5,
                color="#c0392b", bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1))
    ax.set_xticks(xs)
    ax.set_xticklabels(["Sample 1 (>= 10 rated answers)", "Sample 2 (>= 10 rated answers)"])
    ax.set_ylabel("Cross-half Spearman rank stability")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "figS1_shrinkage.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    fig1()
    fig2()
    fig3()
    figS1()
    print("figures written to results/figures")
