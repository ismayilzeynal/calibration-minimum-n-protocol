"""Stream C (selection and content heterogeneity): shared loading helpers.

Reads the paper parquet caches (read-only) and the Eedi question/subject/
student metadata copied into revision-2026-11/data/. Adds per-answer item
difficulty (leave-one-out p-value), strand / topic / term keys and the
tagged flag. Row order of the parquet is kept unchanged, so filtering on
`tag` yields exactly the order used by the original reliability.student_arrays.
"""
import ast
import os
import sys

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import common  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")
META = os.path.join(REV, "data")
OUT_DIR = os.path.join(REV, "results", "C")
SID = {"t12": 1, "t34": 2}
QFILE = {"t12": "question_metadata_task_1_2.csv", "t34": "question_metadata_task_3_4.csv"}
SFILE = {"t12": "student_metadata_task_1_2.csv", "t34": "student_metadata_task_3_4.csv"}
MASTER_SEED = 42
N_DIFF_STRATA = 5


def subject_table():
    sm = pd.read_csv(os.path.join(META, "subject_metadata.csv"), encoding="utf-8-sig")
    return sm


def question_keys(sample):
    """QuestionId -> strand (level-1) and topic (level-2) category labels.

    A question tagged with several level-k subjects gets the combined label
    (sorted ids joined by '+'), so no arbitrary choice is made; a question
    with no level-k subject gets 'none'.
    """
    sm = subject_table()
    lvl = dict(zip(sm.SubjectId, sm.Level))
    name = dict(zip(sm.SubjectId, sm.Name))
    qm = pd.read_csv(os.path.join(META, QFILE[sample]))
    lists = qm.SubjectId.apply(ast.literal_eval)

    def lab(lst, k):
        ids = sorted(s for s in lst if lvl.get(s) == k)
        if not ids:
            return "none"
        return "+".join(name[s] for s in ids)
    out = pd.DataFrame({"QuestionId": qm.QuestionId.astype(np.int32),
                        "strand": lists.apply(lambda l: lab(l, 1)),
                        "topic": lists.apply(lambda l: lab(l, 2)),
                        "n_strand": lists.apply(lambda l: sum(lvl.get(s) == 1 for s in l)),
                        "n_topic": lists.apply(lambda l: sum(lvl.get(s) == 2 for s in l))})
    return out


def term_label(dates):
    """Approximate English school term from the calendar month:
    Autumn = Sep-Dec, Spring = Jan-Mar, Summer = Apr-Aug, with academic year."""
    m = dates.dt.month.to_numpy()
    y = dates.dt.year.to_numpy()
    ay = np.where(m >= 9, y, y - 1)
    t = np.where(m >= 9, 0, np.where(m <= 3, 1, 2))
    code = ay.astype(np.int64) * 10 + t
    names = {0: "Autumn", 1: "Spring", 2: "Summer"}
    uniq = np.unique(code)
    lab = {c: f"{c // 10}/{(c // 10 + 1) % 100:02d} {names[c % 10]}" for c in uniq}
    return pd.Series(pd.Categorical.from_codes(np.searchsorted(uniq, code),
                                               [lab[c] for c in uniq]),
                     index=dates.index)


def load(sample):
    df = common.build_responses(sample)  # parquet read, original row order
    c = df["Confidence"]
    df["tag"] = (c.notna() & (c >= 0) & (c <= 100)).to_numpy()
    g = df.groupby("QuestionId")["IsCorrect"]
    n_q = g.transform("size").to_numpy().astype(np.float64)
    s_q = g.transform("sum").to_numpy().astype(np.float64)
    x = df["IsCorrect"].to_numpy().astype(np.float64)
    df["p_all"] = s_q / n_q
    with np.errstate(divide="ignore", invalid="ignore"):
        df["p_loo"] = np.where(n_q > 1, (s_q - x) / (n_q - 1), np.nan)
    qk = question_keys(sample)
    df = df.merge(qk[["QuestionId", "strand", "topic"]], on="QuestionId",
                  how="left", sort=False, validate="many_to_one")
    df["strand"] = df["strand"].fillna("unknown").astype("category")
    df["topic"] = df["topic"].fillna("unknown").astype("category")
    df["term"] = term_label(df["DateAnswered"])
    # global answer-level difficulty quintiles of the leave-one-out p-value
    edges = np.nanquantile(df["p_loo"].to_numpy(), np.linspace(0, 1, N_DIFF_STRATA + 1)[1:-1])
    p = df["p_loo"].to_numpy()
    st = np.searchsorted(edges, p, side="right").astype(np.int16)
    st[np.isnan(p)] = N_DIFF_STRATA  # undefined p (question answered once)
    df["dstrat"] = st
    df.attrs["dstrat_edges"] = [float(e) for e in edges]
    return df


def boot_mean_ci(v, seed, B=1000):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    if len(v) == 0:
        return [None, None]
    rng = np.random.default_rng(seed)
    bs = np.array([v[rng.integers(0, len(v), len(v))].mean() for _ in range(B)])
    return [float(q) for q in np.percentile(bs, [2.5, 97.5])]


def summ(v, seed=None):
    v = np.asarray(v, float)
    v = v[~np.isnan(v)]
    out = {"n": int(len(v)), "mean": float(v.mean()) if len(v) else None,
           "median": float(np.median(v)) if len(v) else None,
           "sd": float(v.std(ddof=1)) if len(v) > 1 else None,
           "q25": float(np.percentile(v, 25)) if len(v) else None,
           "q75": float(np.percentile(v, 75)) if len(v) else None,
           "share_positive": float((v > 0).mean()) if len(v) else None}
    if seed is not None:
        out["mean_ci95"] = boot_mean_ci(v, seed)
    return out
