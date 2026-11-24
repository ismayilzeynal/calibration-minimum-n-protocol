"""Baseline reproduction of the submitted reliability curves + minN JSONs.

Runs the UNCHANGED original scripts/reliability.py, but with the working
directory set to revision-2026-11/results/baseline so that every relative
write ("results/...", "data/cache/...") lands there and the canonical
paper results/ folder is never touched. Input parquet is read from the
paper data/ folder via an absolute common.DATA_DIR.

Usage:  python repro_baseline.py t12   (or t34)
"""
import os
import sys
import time

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUTDIR = os.path.join(REV, "results", "baseline")

sys.dont_write_bytecode = True  # never write __pycache__ into the paper scripts/
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402
import reliability  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")


def main(sample):
    os.makedirs(OUTDIR, exist_ok=True)
    os.chdir(OUTDIR)
    t0 = time.time()
    reliability.run(sample)
    print(f"[{sample}] wall seconds: {time.time() - t0:.1f}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "t12")
