"""Confirmatory H1-H3 with Holm correction, per the frozen plan Sec 4.
Also evaluates the S2 replication criteria (Sec 7).

Run after reliability.py + h3_shrinkage.py for both samples:
    python scripts/confirmatory.py
"""
import json
import os

import numpy as np

B = 1000


def holm(pvals):
    """Holm step-down adjusted p-values (dict name -> p)."""
    items = sorted(pvals.items(), key=lambda kv: kv[1])
    m = len(items)
    adj, running = {}, 0.0
    for i, (k, p) in enumerate(items):
        running = max(running, (m - i) * p)
        adj[k] = min(running, 1.0)
    return adj


def main():
    mn1 = json.load(open("results/minN_t12.json"))
    mn2 = json.load(open("results/minN_t34.json"))
    h3_1 = json.load(open("results/h3_shrinkage_t12.json"))
    h3_2 = json.load(open("results/h3_shrinkage_t34.json"))
    cache1 = np.load(os.path.join("data", "cache", "designA_t12.npz"))
    GRID = mn1["_meta"]["grid"]

    # H1: bias r(20) < .70 in S1 Design A
    r20_draws = cache1["boot_bias"][:, GRID.index(20)]
    p_h1 = float(np.mean(r20_draws >= 0.70))
    h1 = {"r20_point": mn1["bias"]["r20_point"],
          "r20_ci": mn1["bias"]["r20_ci"],
          "r20_upper_1sided95": mn1["bias"]["r20_upper_1sided95"],
          "p_onesided": max(p_h1, 1.0 / B if p_h1 == 0 else p_h1),
          "p_raw": p_h1}

    # H2: r(15) < .70 for all four indices, Holm inner family
    inner = {k: mn1[k]["p_r15_ge_070"] for k in
             ["bias", "brier", "gamma", "auroc"]}
    inner_adj = holm(inner)
    p_h2 = max(inner_adj.values())
    h2 = {"r15_points": {k: mn1[k]["r15_point"] for k in inner},
          "p_inner_raw": inner, "p_inner_holm": inner_adj,
          "p_headline": p_h2,
          "minN70_points": {k: mn1[k]["minN70_point"] for k in inner},
          "median_volume_anchor": 15}

    # H3: EB shrinkage improves rank stability
    h3 = {"rho_raw": h3_1["rho_raw"], "rho_eb": h3_1["rho_eb"],
          "diff_point": h3_1["diff_point"], "diff_ci95": h3_1["diff_ci95"],
          "diff_lower_1sided95": h3_1["diff_lower_1sided95"],
          "p_onesided": h3_1["p_onesided_diff_le_0"]}

    fam = {"H1": max(h1["p_raw"], 1.0 / B),
           "H2": max(p_h2, 1.0 / B),
           "H3": max(h3["p_onesided"], 1.0 / B)}
    fam_adj = holm(fam)
    support = {
        "H1": bool(h1["r20_upper_1sided95"] < 0.70 and fam_adj["H1"] < 0.05),
        "H2": bool(all(v < 0.05 for v in inner_adj.values())
                   and fam_adj["H2"] < 0.05),
        "H3": bool(h3["diff_lower_1sided95"] > 0 and fam_adj["H3"] < 0.05),
    }

    # S2 replication criteria (qualitative, Sec 7)
    minn70_s2 = {k: mn2[k]["minN70_point"] for k in inner}
    crit_a = bool(mn2["bias"]["r20_point"] < 0.70)
    crit_b_hits = {k: bool((np.isinf(v) or v > 31)) for k, v in minn70_s2.items()}
    crit_b = sum(crit_b_hits.values()) >= 3
    crit_c = bool(h3_2["diff_point"] > 0)
    rep = {"a_r20_bias_lt_070": {"value": mn2["bias"]["r20_point"],
                                 "met": crit_a},
           "b_minN70_gt_median31_for_3of4": {"minN70": {
               k: (None if np.isinf(v) else v) for k, v in minn70_s2.items()},
               "hits": crit_b_hits, "met": bool(crit_b)},
           "c_eb_improves": {"diff": h3_2["diff_point"], "met": crit_c}}

    out = {"H1": h1, "H2": h2, "H3": h3,
           "holm_family": fam_adj, "supported": support,
           "replication_S2": rep,
           "note_p_granularity": "B=1000 bootstrap draws; smallest reportable "
                                 "one-sided p is 0.001 (reported as < .001)"}
    with open("results/confirmatory_stats.json", "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(json.dumps({"supported": support, "holm": fam_adj,
                      "replication": {k: v["met"] if isinstance(v, dict) and "met" in v
                                      else v for k, v in rep.items()}}, indent=1,
                     default=float))


if __name__ == "__main__":
    main()
