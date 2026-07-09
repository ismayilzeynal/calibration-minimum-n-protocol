"""Manuscript figures. Run after reliability + h3 for both samples:
    python scripts/figures.py
Palette kept colorblind-safe and consistent across figures.
"""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COL = {"bias": "#2a78d6", "brier": "#1baf7a", "gamma": "#eda100",
       "auroc": "#a34fc9"}
LAB = {"bias": "Bias (over/underconfidence)", "brier": "Absolute accuracy (Brier)",
       "gamma": "Relative accuracy (gamma)", "auroc": "Relative accuracy (AUROC)"}
TICKS = [5, 10, 20, 50, 100, 200]
FIGDIR = os.path.join("results", "figures")


def curves(sample):
    return pd.read_csv(os.path.join("results",
                                    f"reliability_curves_{sample}.csv"))


def fig1():
    """Design A reliability-vs-N curves, S1, 4 indices, CI bands."""
    cv = curves("t12")
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for k in COL:
        sub = cv[(cv.design == "A") & (cv.split == "random")
                 & (cv["index"] == k)].sort_values("N")
        ax.plot(sub.N, sub.r, "-o", color=COL[k], label=LAB[k], ms=4.5, lw=1.8)
        ax.fill_between(sub.N, sub.ci_lo, sub.ci_hi, color=COL[k], alpha=0.15,
                        lw=0)
    for thr, txt in ((0.70, ".70"), (0.80, ".80")):
        ax.axhline(thr, color="#555555", ls=":", lw=1)
        ax.text(214, thr, txt, fontsize=8, color="#555555", va="center")
    ax.axvline(15, color="#888888", ls="--", lw=1)
    ax.text(15.6, 0.03, "median student volume (15)", rotation=90, fontsize=8,
            color="#666666", va="bottom")
    ax.set_xscale("log")
    ax.set_xticks(TICKS)
    ax.set_xticklabels([str(n) for n in TICKS])
    ax.set_xlabel("Number of confidence-tagged responses per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(-0.05, 1.0)
    ax.set_xlim(right=230)
    ax.legend(fontsize=9, loc="upper left", framealpha=1.0, edgecolor="none")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig1_curves_s1.png"), dpi=300)
    plt.close(fig)


def fig2():
    """(a) undefined share for gamma/AUROC; (b) temporal vs random (bias)."""
    cv = curves("t12")
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.9))
    ax = axes[0]
    for k in ("gamma", "auroc"):
        sub = cv[(cv.design == "A") & (cv.split == "random")
                 & (cv["index"] == k)].sort_values("N")
        ax.plot(sub.N, 100 * sub.pct_undefined_full, "-o", color=COL[k],
                label=LAB[k], ms=4.5, lw=1.8)
    ax.set_xscale("log"); ax.set_xticks(TICKS)
    ax.set_xticklabels([str(n) for n in TICKS])
    ax.set_xlabel("N per student"); ax.set_ylabel("Students with undefined index (%)")
    ax.set_title("(a) Incomputable at small N", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    for split, ls, mfc in (("random", "-", None), ("temporal", "--", "white")):
        sub = cv[(cv.design == "A") & (cv.split == split)
                 & (cv["index"] == "bias")].sort_values("N")
        ax.plot(sub.N, sub.r, ls + "o", color=COL["bias"], ms=4.5, lw=1.8,
                mfc=mfc, label=f"{split} split")
    ax.axhline(0.70, color="#555555", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_xticks(TICKS)
    ax.set_xticklabels([str(n) for n in TICKS])
    ax.set_xlabel("N per student"); ax.set_ylabel("Split-half reliability")
    ax.set_title("(b) Bias: temporal vs random split", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig2_undefined_temporal.png"), dpi=300)
    plt.close(fig)


def fig3():
    """S1 vs S2 replication overlay (bias + gamma, Design A)."""
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for sample, ls, mfc in (("t12", "-", None), ("t34", "--", "white")):
        cv = curves(sample)
        for k in ("bias", "gamma"):
            sub = cv[(cv.design == "A") & (cv.split == "random")
                     & (cv["index"] == k)].sort_values("N")
            nm = "Bias" if k == "bias" else "Gamma (relative accuracy)"
            lbl = f"{nm} - " + ("sample 1" if sample == "t12" else "sample 2")
            ax.plot(sub.N, sub.r, ls + "o", color=COL[k], ms=4.5, lw=1.8,
                    mfc=mfc, label=lbl)
    ax.axhline(0.70, color="#555555", ls=":", lw=1)
    ax.set_xscale("log"); ax.set_xticks(TICKS)
    ax.set_xticklabels([str(n) for n in TICKS])
    ax.set_xlabel("Number of confidence-tagged responses per student (N)")
    ax.set_ylabel("Split-half reliability (Spearman-Brown corrected)")
    ax.set_ylim(-0.1, 1.0)
    ax.legend(fontsize=8, ncol=2, loc="lower right", framealpha=1.0,
              edgecolor="none")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig4_replication.png"), dpi=300)
    plt.close(fig)


def fig4():
    """EB shrinkage: pre-registered symmetric metric vs artifact-immune
    validation metric, with the permutation-null artifact shown."""
    vals, nulls = {}, {}
    for s in ("t12", "t34"):
        vals[s] = json.load(open(os.path.join("results",
                                              f"h3_shrinkage_{s}.json")))
        nulls[s] = json.load(open(os.path.join("results",
                                               f"h3_null_{s}.json")))
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    xs = np.arange(2)
    w = 0.26
    raw = [vals["t12"]["rho_raw"], vals["t34"]["rho_raw"]]
    eb = [vals["t12"]["rho_eb"], vals["t34"]["rho_eb"]]
    ebv = [vals["t12"]["rho_eb_validation"], vals["t34"]["rho_eb_validation"]]
    nul = [nulls["t12"]["null_rho_eb"], nulls["t34"]["null_rho_eb"]]
    ax.bar(xs - w, raw, w, color="#9dbfe3", label="Raw bias")
    ax.bar(xs, eb, w, color=COL["bias"],
           label="EB, same-metric (inflated; see null)")
    ax.bar(xs + w, ebv, w, color="#1baf7a",
           label="EB, validation vs raw criterion")
    for i in range(2):
        ax.text(i - w, raw[i] + 0.012, f"{raw[i]:.3f}", ha="center", fontsize=8)
        ax.text(i, eb[i] + 0.012, f"{eb[i]:.3f}", ha="center", fontsize=8)
        ax.text(i + w, ebv[i] + 0.012, f"{ebv[i]:.3f}", ha="center", fontsize=8)
        ax.hlines(nul[i], i - w / 2, i + w / 2, color="#c0392b", ls="--", lw=1.6)
        ax.text(i, nul[i] + 0.015, f"null artifact {nul[i]:.2f}", ha="center",
                fontsize=7.5, color="#c0392b",
                bbox=dict(facecolor="white", alpha=0.85, edgecolor="none",
                          pad=1))
    ax.set_xticks(xs)
    ax.set_xticklabels(["Sample 1 (>=10 tagged)", "Sample 2 (>=10 tagged)"])
    ax.set_ylabel("Cross-half Spearman rank stability")
    ax.set_ylim(0, 1.0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIGDIR, "fig3_shrinkage.png"), dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIGDIR, exist_ok=True)
    fig1(); fig2(); fig3(); fig4()
    print("figures written to", FIGDIR)
