"""Print the key tables of results/r1_A.json (read-only) for writing r1_A.md.
Usage: python -B r1_A_digest.py
"""
import json
import os

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = json.load(open(os.path.join(REV, "results", "r1_A.json")))


def f(v, d=3):
    if v is None:
        return "None"
    if isinstance(v, str):
        return v
    return f"{v:.{d}f}"


def ci(c):
    if isinstance(c, str):
        return c
    if c is None:
        return ""
    return "[" + ", ".join("None" if v is None else f(v, 1) for v in c) + "]"


for S in ("S1", "S2"):
    b = R[S]
    print(f"===== {S} poolA={b['n_poolA']} reporters={b['reporters']}")
    print("-- incomputability (gamma | auroc_mw | auroc_prereg | forced)")
    for key in ("designA_full", "designA_halves_random", "designA_halves_temporal",
                "designB_full", "designB_halves_random"):
        print(" ", key)
        for N, v in b["incomputability"][key].items():
            extra = f" n={v['n_pool']}" if "n_pool" in v else ""
            print(f"    N={N:>3}: {f(v['gamma'])} | {f(v['auroc_mw'])} | "
                  f"{f(v['auroc_prereg'])} | {f(v['forced_05'])}{extra}")
    print("-- persistence", json.dumps(b["forced_05_persistence_designA"]))
    print("-- curves r (A random): N: gamma | auroc_mw [ci] | auroc_prereg [ci] | brier | bias")
    cA = b["curves"]["designA_random"]
    for N in cA["gamma"]:
        am, ap = cA["auroc_mw"][N], cA["auroc_prereg"][N]
        print(f"    N={N:>3}: {f(cA['gamma'][N]['r'])} | {f(am['r'])} [{f(am['ci_lo'])}, {f(am['ci_hi'])}]"
              f" | {f(ap['r'])} [{f(ap['ci_lo'])}, {f(ap['ci_hi'])}] | {f(cA['brier'][N]['r'])}"
              f" | {f(cA['bias'][N]['r'])} | conf {f(cA['conf'][N]['r'])} | acc {f(cA['acc'][N]['r'])}")
    for dk in ("designA_temporal", "designB_random", "designB_temporal"):
        cc = b["curves"][dk]
        print(f"-- curves r ({dk}): N: gamma | auroc_mw | auroc_prereg | brier | bias")
        for N in cc["gamma"]:
            print(f"    N={N:>3}: {f(cc['gamma'][N]['r'])} | {f(cc['auroc_mw'][N]['r'])} | "
                  f"{f(cc['auroc_prereg'][N]['r'])} | {f(cc['brier'][N]['r'])} | {f(cc['bias'][N]['r'])}")
    print("-- thresholds (point, ci, coverage)")
    for idx, t in b["thresholds_by_index"].items():
        for dk, tt in t.items():
            row = []
            for thr, e in tt.items():
                cov = e["coverage"]
                covs = f"{cov['k']}:{f(cov['share_reporters_ge_k'])}" if cov else "-"
                le = "<=" if e["at_or_below_first_grid_point"] else ""
                row.append(f"{thr}={le}{f(e['point'], 1)} {ci(e.get('ci'))} cov {covs}")
            print(f"  {idx:13s} {dk:17s} " + " ; ".join(row))
    print("-- minN extras")
    for idx in ("gamma", "auroc_mw", "auroc_prereg", "bias", "brier"):
        m = b["minN"][idx]
        print(f"  {idx}: prophecy70={f(m['prophecy_minN70_from20'], 1)} r15={f(m['r15_point'])} "
              f"r15ci={ci(m['r15_ci'])} p_r15={m['p_r15_ge_070']} r20={f(m['r20_point'])} "
              f"r20ci=[{f(m['r20_ci'][0])}, {f(m['r20_ci'][1])}] share70gt200={m['share_minN70_gt200']} "
              f"share80gt200={m['share_minN80_gt200']} temporal70ci={ci(m['minN70_temporal_ci'])}")
    print("-- confound", json.dumps(b["accuracy_confound"]))
    print("-- ipw", json.dumps(b["ipw_designB_N20"]))
    print("-- retest", json.dumps(b["test_retest"]))
    print("-- qualitative", json.dumps(b["qualitative_checks"]))
print("===== confirmatory")
print(json.dumps(R["confirmatory"], indent=0))
print("runtime", json.dumps(R["runtime_seconds"]))
