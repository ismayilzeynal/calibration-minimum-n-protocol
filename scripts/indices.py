"""Calibration indices computed from per-student 5x2 count tables.

Every index in the plan is a function of the (confidence level x correctness)
count table, which makes half-split computation fully vectorizable.
Levels: 0=0, 1=25, 2=50, 3=75, 4=100 -> c in {0, .25, .5, .75, 1}.
"""
import numpy as np

LEVELS = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
LEVELS3 = np.array([0.0, 0.0, 0.5, 1.0, 1.0])  # coarse 3-category recode


def count_tables(codes, x, mask):
    """codes,x,mask: (n, N). Returns counts (n, 5, 2) over masked entries."""
    n = codes.shape[0]
    row = np.broadcast_to(np.arange(n)[:, None], codes.shape)
    flat = (row * 10 + codes.astype(np.int64) * 2 + x).ravel()
    w = mask.ravel().astype(np.float64)
    cnt = np.bincount(flat, weights=w, minlength=n * 10)
    return cnt.reshape(n, 5, 2)


def indices_from_counts(cnt, min_n=1):
    """cnt: (n, 5, 2). Returns dict of per-student index arrays (NaN where undefined)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ntot = cnt.sum(axis=(1, 2))
        n1 = cnt[:, :, 1].sum(axis=1)
        n0 = cnt[:, :, 0].sum(axis=1)
        nk = cnt.sum(axis=2)
        acc = n1 / ntot
        conf = (nk * LEVELS).sum(axis=1) / ntot
        bias = conf - acc
        conf3 = (nk * LEVELS3).sum(axis=1) / ntot
        bias3 = conf3 - acc
        brier = ((cnt[:, :, 0] * LEVELS ** 2).sum(axis=1)
                 + (cnt[:, :, 1] * (LEVELS - 1.0) ** 2).sum(axis=1)) / ntot
        # pair counts for gamma / AUROC
        c0 = np.cumsum(cnt[:, :, 0], axis=1)
        c0s = np.concatenate([np.zeros((len(cnt), 1)), c0[:, :-1]], axis=1)
        C = (cnt[:, :, 1] * c0s).sum(axis=1)          # correct above incorrect
        c1 = np.cumsum(cnt[:, :, 1], axis=1)
        c1s = np.concatenate([np.zeros((len(cnt), 1)), c1[:, :-1]], axis=1)
        D = (cnt[:, :, 0] * c1s).sum(axis=1)          # incorrect above correct
        ties = (cnt[:, :, 0] * cnt[:, :, 1]).sum(axis=1)
        gamma = np.where(C + D > 0, (C - D) / (C + D), np.nan)
        auroc = np.where(n0 * n1 > 0, (C + 0.5 * ties) / (n0 * n1), np.nan)
        # Murphy decomposition over the 5 native bins
        obar_k = np.where(nk > 0, cnt[:, :, 1] / np.maximum(nk, 1), 0.0)
        cal = (nk * (LEVELS - obar_k) ** 2).sum(axis=1) / ntot
        res = (nk * (obar_k - acc[:, None]) ** 2).sum(axis=1) / ntot
        bad = ntot < min_n
        out = {"bias": bias, "brier": brier, "gamma": gamma, "auroc": auroc,
               "bias3": bias3, "cal": cal, "res": res, "acc": acc, "conf": conf,
               "n": ntot}
        for k in out:
            if k != "n":
                out[k] = np.where(bad, np.nan, out[k])
        return out


CONFIRMATORY = ["bias", "brier", "gamma", "auroc"]
EXPLORATORY = ["bias3", "cal", "res"]
