"""Stream A (E5): rerun the submitted AUROC-dependent chain under the
PRE-REGISTERED AUROC rule, as a side-by-side sensitivity to the submitted rule.

Submitted rule (scripts/indices.py line 46, unchanged here): Mann-Whitney
AUROC = (C + 0.5 * ties) / (n0 * n1), undefined only when n0 * n1 = 0. A
student with constant confidence and mixed correctness gets exactly .5.

Pre-registered rule (analysis_plan.md section 2, manuscript "Indices and the
Undefined-Value Protocol"): AUROC undefined when confidence OR correctness has
zero variance. Implemented by wrapping indices.indices_from_counts and setting
AUROC to NaN where all responses sit at one confidence level. Nothing else
changes and RNG consumption is identical, so every non-AUROC output must
reproduce results/baseline exactly (checked by r1_A_gate.py).

The original modules are imported read-only; the working directory is set to
revision-2026-11/results/r1_A/prereg so every relative write lands there.

Usage (python -B):
  r1_A_prereg_run.py reliability t12|t34
  r1_A_prereg_run.py confound t12|t34
  r1_A_prereg_run.py confirmatory
"""
import json
import os
import shutil
import sys
import time

import numpy as np

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUTDIR = os.path.join(REV, "results", "r1_A", "prereg")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import indices as IX  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")
_orig_indices = IX.indices_from_counts


def indices_prereg(cnt, min_n=1):
    out = _orig_indices(cnt, min_n)
    n_levels_used = (cnt.sum(axis=2) > 0).sum(axis=1)
    const_conf = n_levels_used <= 1
    out["auroc"] = np.where(const_conf, np.nan, out["auroc"])
    return out


IX.indices_from_counts = indices_prereg


def main(stage, sample=None):
    os.makedirs(os.path.join(OUTDIR, "results"), exist_ok=True)
    os.makedirs(os.path.join(OUTDIR, "data", "cache"), exist_ok=True)
    os.chdir(OUTDIR)
    t0 = time.time()
    if stage == "reliability":
        import reliability
        reliability.run(sample)
    elif stage == "confound":
        import confound_selection
        confound_selection.run(sample)
    elif stage == "confirmatory":
        for s in ("t12", "t34"):
            src = os.path.join(PAPER, "results", f"h3_shrinkage_{s}.json")
            shutil.copyfile(src, os.path.join(OUTDIR, "results",
                                              f"h3_shrinkage_{s}.json"))
        import confirmatory
        confirmatory.main()
    else:
        raise SystemExit(f"unknown stage {stage}")
    wall = time.time() - t0
    tag = f"{stage}_{sample}" if sample else stage
    with open(os.path.join(OUTDIR, f"runtime_{tag}.json"), "w") as fh:
        json.dump({"stage": stage, "sample": sample, "wall_seconds": wall,
                   "python": sys.version.split()[0],
                   "numpy": np.__version__}, fh, indent=1)
    print(f"[{tag}] wall seconds: {wall:.1f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
