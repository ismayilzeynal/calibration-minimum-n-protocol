"""Stream A (E5): rerun the downstream AUROC-dependent analyses under the
SUBMITTED Mann-Whitney rule (unchanged code), so that every auroc_mw number
in r1_A.json comes from a revision rerun and not only from the submitted files.

Inputs: the baseline rerun (results/baseline: curves CSV and Design A cache,
which reproduce the submitted files) copied into results/r1_A/mw/, plus the
submitted h3_shrinkage JSONs (read-only). Runs confound_selection.py for both
samples and confirmatory.py, unchanged, with the working directory set to
results/r1_A/mw. Then compares every field (AUROC included) with the
submitted paper results/ and writes results/r1_A/mw/compare_submitted.json.

Usage: python -B r1_A_mw_run.py
"""
import json
import os
import shutil
import sys
import time

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
BASE = os.path.join(REV, "results", "baseline")
OUTDIR = os.path.join(REV, "results", "r1_A", "mw")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")


def main():
    os.makedirs(os.path.join(OUTDIR, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, "data", "cache"), exist_ok=True)
    for s in ("t12", "t34"):
        shutil.copyfile(os.path.join(BASE, "results", f"reliability_curves_{s}.csv"),
                        os.path.join(OUTDIR, "results", f"reliability_curves_{s}.csv"))
        shutil.copyfile(os.path.join(BASE, "results", f"minN_{s}.json"),
                        os.path.join(OUTDIR, "results", f"minN_{s}.json"))
        shutil.copyfile(os.path.join(BASE, "data", "cache", f"designA_{s}.npz"),
                        os.path.join(OUTDIR, "data", "cache", f"designA_{s}.npz"))
        shutil.copyfile(os.path.join(PAPER, "results", f"h3_shrinkage_{s}.json"),
                        os.path.join(OUTDIR, "results", f"h3_shrinkage_{s}.json"))
    os.chdir(OUTDIR)
    t0 = time.time()
    import confound_selection
    import confirmatory
    for s in ("t12", "t34"):
        confound_selection.run(s)
    confirmatory.main()
    wall = time.time() - t0

    import r1_A_gate as G
    cmp = {"wall_seconds": wall}
    for s in ("t12", "t34"):
        cmp[f"confound_selection_{s}"] = G.cmp_json(
            os.path.join(PAPER, "results", f"confound_selection_{s}.json"),
            os.path.join(OUTDIR, "results", f"confound_selection_{s}.json"),
            lambda k: False)
    cmp["confirmatory_stats"] = G.cmp_json(
        os.path.join(PAPER, "results", "confirmatory_stats.json"),
        os.path.join(OUTDIR, "results", "confirmatory_stats.json"), lambda k: False)
    cmp["all_pass"] = bool(all(v["pass"] for k, v in cmp.items()
                               if isinstance(v, dict)))
    with open(os.path.join(OUTDIR, "compare_submitted.json"), "w") as fh:
        json.dump(cmp, fh, indent=1)
    print(json.dumps(cmp, indent=1))


if __name__ == "__main__":
    main()
