"""Canonical data loading for the confidence-calibration reliability study.

All paths are relative to the paper folder root (run scripts from the paper root).
Place the official challenge data.zip at data/raw/eedi/data.zip.
"""
import io
import os
import zipfile

import numpy as np
import pandas as pd

ZIP_PATH = os.path.join("data", "raw", "eedi", "data.zip")
DATA_DIR = "data"

MEMBERS = {
    "train_t12": "data/train_data/train_task_1_2.csv",
    "train_t34": "data/train_data/train_task_3_4.csv",
    "ans_t12": "data/metadata/answer_metadata_task_1_2.csv",
    "ans_t34": "data/metadata/answer_metadata_task_3_4.csv",
    "stu_t12": "data/metadata/student_metadata_task_1_2.csv",
    "stu_t34": "data/metadata/student_metadata_task_3_4.csv",
    "q_t12": "data/metadata/question_metadata_task_1_2.csv",
    "q_t34": "data/metadata/question_metadata_task_3_4.csv",
}


def _read_member(z, member, **kw):
    with z.open(member) as fh:
        return pd.read_csv(io.BytesIO(fh.read()), **kw)


def build_responses(task="t12", force=False):
    """Answer-level table: UserId, QuestionId, IsCorrect, Confidence, DateAnswered, GroupId, QuizId.

    ALL rows kept (Confidence NaN retained) so the reporter-selection analysis
    can compare confidence-tagged vs untagged responses.
    """
    out = os.path.join(DATA_DIR, f"resp_{task}.parquet")
    if os.path.exists(out) and not force:
        return pd.read_parquet(out)
    z = zipfile.ZipFile(ZIP_PATH)
    tr = _read_member(z, MEMBERS[f"train_{task}"],
                      usecols=["QuestionId", "UserId", "AnswerId", "IsCorrect"],
                      dtype={"QuestionId": np.int32, "UserId": np.int32,
                             "AnswerId": np.int64, "IsCorrect": np.int8})
    am = _read_member(z, MEMBERS[f"ans_{task}"],
                      usecols=["AnswerId", "DateAnswered", "Confidence",
                               "GroupId", "QuizId"],
                      dtype={"Confidence": np.float32})
    am = am[am["AnswerId"].notna()].copy()
    am["AnswerId"] = am["AnswerId"].astype(np.float64).astype(np.int64)
    am = am.drop_duplicates(subset="AnswerId", keep="first")
    am["DateAnswered"] = pd.to_datetime(am["DateAnswered"])
    am["GroupId"] = am["GroupId"].fillna(-1).astype(np.int32)
    am["QuizId"] = am["QuizId"].fillna(-1).astype(np.int32)
    df = tr.merge(am, on="AnswerId", how="left", validate="one_to_one")
    df = df.drop(columns=["AnswerId"]).sort_values(
        ["UserId", "DateAnswered"], kind="stable").reset_index(drop=True)
    os.makedirs(DATA_DIR, exist_ok=True)
    df.to_parquet(out)
    return df


def tagged(df):
    """Rows with a usable confidence rating (non-missing, in [0, 100])."""
    c = df["Confidence"]
    return df[c.notna() & (c >= 0) & (c <= 100)]
