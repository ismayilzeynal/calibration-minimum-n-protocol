"""Stream C parts (1) and (4): who rates which answers, and the S1 vs S2 table.

(1) Within-student comparison of confidence-tagged vs untagged answers:
    item difficulty (leave-one-out item p-value), accuracy and accuracy
    relative to item difficulty, strand / topic mix (total variation
    distance with a within-student permutation baseline), timing (elapsed
    days since the student's first answer, weekend share), on two bases:
    (a) all untagged answers of the student, (b) only answers in quiz
    attempts (UserId x QuizId) where the student rated some but not all
    answers ("same-quiz" basis, isolates item-level choice).
    Also where the variation in tagging sits (eta-squared of the tag
    indicator by student, quiz, question, student x quiz) and tag share by
    calendar month.
(4) Sample descriptives table (S1 vs S2).

Usage: python r1_C_selection.py t12   (or t34)
Writes: revision-2026-11/results/C/r1_C_selection_{sample}.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r1_C_common as RC  # noqa: E402

MIN_EACH = 5


def eta2(y, keys):
    y = y.astype(np.float64)
    p = y.mean()
    sst = len(y) * p * (1 - p)
    g = pd.DataFrame({"y": y, **keys}).groupby(list(keys), sort=False)["y"].agg(["mean", "size"])
    ssb = float((g["size"] * (g["mean"] - p) ** 2).sum())
    return ssb / sst, int(len(g))


def tvd_by_student(uid, tag, cat):
    """Per-student total variation distance between the category
    distribution of tagged and of untagged answers."""
    d = pd.DataFrame({"u": uid, "t": tag, "c": cat})
    ct = d.groupby(["u", "t", "c"], sort=False).size().unstack("c", fill_value=0)
    ct = ct.div(ct.sum(axis=1), axis=0)
    tt = ct.xs(True, level="t")
    uu = ct.xs(False, level="t").reindex(tt.index)
    return 0.5 * (tt - uu).abs().sum(axis=1)


def within_perm(uid_sorted, vals, rng):
    """Permute vals within contiguous uid blocks."""
    r = rng.random(len(vals))
    idx = np.lexsort((r, uid_sorted))
    return vals[idx]


def paired_block(df, sid, base_tag, label, seed_base):
    """df: rows of eligible students only, contiguous by UserId."""
    t = df["tag"].to_numpy()
    out = {"basis": label, "n_students": int(df.UserId.nunique()),
           "n_answers": int(len(df)), "n_tagged": int(t.sum())}
    g = df.groupby(["UserId", "tag"], sort=False)
    agg = g.agg(acc=("IsCorrect", "mean"), p=("p_loo", "mean"), rel=("rel", "mean"),
                days=("days", "mean"), wkend=("wkend", "mean"),
                pall=("p_all", "mean"), n=("IsCorrect", "size")).unstack("tag")
    k = 0
    for name in ["acc", "p", "pall", "rel", "days", "wkend"]:
        diff = (agg[(name, True)] - agg[(name, False)]).to_numpy()
        out[f"{name}_tagged_minus_untagged"] = RC.summ(diff, seed=[RC.MASTER_SEED, seed_base + k, sid])
        out[f"{name}_mean_tagged"] = float(np.nanmean(agg[(name, True)]))
        out[f"{name}_mean_untagged"] = float(np.nanmean(agg[(name, False)]))
        k += 1
    if "pos" in df:
        agg2 = g["pos"].mean().unstack("tag")
        diff = (agg2[True] - agg2[False]).to_numpy()
        out["relpos_tagged_minus_untagged"] = RC.summ(diff, seed=[RC.MASTER_SEED, seed_base + k, sid])
        out["relpos_mean_tagged"] = float(np.nanmean(agg2[True]))
        out["relpos_mean_untagged"] = float(np.nanmean(agg2[False]))
        k += 1
    # strand / topic mix
    rng = np.random.default_rng([RC.MASTER_SEED, seed_base + 50, sid])
    uid = df.UserId.to_numpy()
    for cat in ["strand", "topic"]:
        tv = tvd_by_student(uid, t, df[cat].to_numpy())
        tp = within_perm(uid, t, rng)
        tvp = tvd_by_student(uid, tp, df[cat].to_numpy())
        out[f"tvd_{cat}"] = RC.summ(tv.to_numpy(), seed=[RC.MASTER_SEED, seed_base + k, sid])
        out[f"tvd_{cat}_permutation_baseline"] = RC.summ(tvp.to_numpy())
        out[f"tvd_{cat}_minus_baseline_mean"] = float(tv.mean() - tvp.mean())
        k += 1
    # per-strand share differences (mean over students of tagged share - untagged share)
    d = pd.DataFrame({"u": uid, "t": t, "c": df["strand"].to_numpy()})
    ct = d.groupby(["u", "t", "c"], sort=False).size().unstack("c", fill_value=0)
    ct = ct.div(ct.sum(axis=1), axis=0)
    tt = ct.xs(True, level="t")
    uu = ct.xs(False, level="t").reindex(tt.index)
    sd = (tt - uu)
    top = d["c"].value_counts().index[:6]
    out["strand_share_diff_mean"] = {c: float(sd[c].mean()) for c in top}
    out["strand_share_tagged_mean"] = {c: float(tt[c].mean()) for c in top}
    out["strand_share_untagged_mean"] = {c: float(uu[c].mean()) for c in top}
    return out


def run(sample):
    t0 = time.time()
    sid = RC.SID[sample]
    df = RC.load(sample)
    t_load = time.time() - t0
    res = {"sample": sample, "dstrat_edges_p_loo": df.attrs["dstrat_edges"]}
    tag = df["tag"].to_numpy()
    x = df["IsCorrect"].to_numpy()

    # ---------------- (4) sample descriptives ----------------
    per_u = df.groupby("UserId").agg(n=("tag", "size"), ntag=("tag", "sum"),
                                     first=("DateAnswered", "min"),
                                     last=("DateAnswered", "max"))
    rep = per_u[per_u.ntag > 0]
    tg = df[tag]
    lv = (tg["Confidence"] / 25.0).round().astype(int)
    q = df.groupby("QuestionId")["IsCorrect"].agg(["mean", "size"])
    month = df["DateAnswered"].dt.to_period("M").astype(str)
    mshare = df.assign(m=month).groupby("m")["tag"].agg(["mean", "size", "sum"])
    strand_all = df["strand"].value_counts(normalize=True)
    strand_tag = tg["strand"].value_counts(normalize=True)
    cover = df.groupby("strand", observed=True)["tag"].mean()
    assert df.QuizId.max() < 100000 and df.QuizId.min() >= 0
    st = pd.read_csv(os.path.join(RC.META, RC.SFILE[sample]))
    st = st[st.UserId.isin(per_u.index)].set_index("UserId").reindex(per_u.index)
    dob = pd.to_datetime(st["DateOfBirth"], errors="coerce")
    age = (per_u["first"] - dob).dt.days / 365.25

    def demo(mask):
        a = age[mask]
        pp = st["PremiumPupil"][mask]
        gd = st["Gender"][mask]
        return {"n_students": int(mask.sum()),
                "share_dob_known": float(dob[mask].notna().mean()),
                "median_birth_year": float(dob[mask].dt.year.median()),
                "median_age_years_at_first_answer": float(a.median()),
                "age_q25_q75": [float(a.quantile(.25)), float(a.quantile(.75))],
                "share_premiumpupil_known": float(pp.notna().mean()),
                "share_premiumpupil_1_among_known": float((pp[pp.notna()] == 1).mean()),
                "gender_code_shares_raw": {str(k): float(v) for k, v in
                                           gd.value_counts(normalize=True).sort_index().items()}}
    span_days = (per_u["last"] - per_u["first"]).dt.days
    desc = {
        "n_answers": int(len(df)), "n_students": int(len(per_u)),
        "n_questions": int(df.QuestionId.nunique()), "n_quizzes": int(df.QuizId.nunique()),
        "n_groups": int(df.GroupId.nunique()),
        "n_tagged": int(tag.sum()), "share_answers_tagged": float(tag.mean()),
        "n_reporters": int(len(rep)), "share_students_reporting": float(len(rep) / len(per_u)),
        "median_tagged_per_reporter": float(rep.ntag.median()),
        "mean_tagged_per_reporter": float(rep.ntag.mean()),
        "reporters_with_ge": {str(k): int((rep.ntag >= k).sum()) for k in (5, 10, 20, 50, 100, 200)},
        "median_answers_per_student": float(per_u.n.median()),
        "median_answers_per_reporter": float(rep.n.median()),
        "median_active_span_days": float(span_days.median()),
        "accuracy_all": float(x.mean()), "accuracy_tagged": float(x[tag].mean()),
        "accuracy_untagged": float(x[~tag].mean()),
        "mean_confidence_tagged": float(tg["Confidence"].mean()),
        "confidence_level_shares": {str(int(k * 25)): float(v) for k, v in
                                    lv.value_counts(normalize=True).sort_index().items()},
        "first_answer": str(df.DateAnswered.min()), "last_answer": str(df.DateAnswered.max()),
        "first_tagged": str(tg.DateAnswered.min()), "last_tagged": str(tg.DateAnswered.max()),
        "tag_share_by_month": {m: float(v) for m, v in mshare["mean"].items()},
        "answers_by_month": {m: int(v) for m, v in mshare["size"].items()},
        "item_p_question_level": {"mean": float(q["mean"].mean()),
                                  "q25": float(q["mean"].quantile(.25)),
                                  "median": float(q["mean"].median()),
                                  "q75": float(q["mean"].quantile(.75))},
        "answers_per_question_median": float(q["size"].median()),
        "mean_item_p_all_answers": float(df["p_all"].mean()),
        "mean_item_p_tagged_answers": float(tg["p_all"].mean()),
        "mean_item_p_loo_tagged_answers": float(np.nanmean(tg["p_loo"])),
        "mean_item_p_loo_untagged_answers": float(np.nanmean(df.loc[~tag, "p_loo"])),
        "n_strand_categories": int(df.strand.nunique()), "n_topic_categories": int(df.topic.nunique()),
        "strand_share_all_answers": {k: float(v) for k, v in strand_all.head(8).items()},
        "strand_share_tagged_answers": {k: float(v) for k, v in strand_tag.head(8).items()},
        "tag_coverage_by_strand": {k: float(cover[k]) for k in strand_all.head(8).index},
        "demographics_all_students": demo(np.ones(len(per_u), bool)),
        "demographics_reporters": demo((per_u.ntag > 0).to_numpy()),
        "student_metadata_note": "Gender is reported as raw codes; the code-to-label mapping is "
                                 "not verified here (unknown). DateOfBirth is month resolution.",
    }
    res["descriptives"] = desc

    # ---------------- (1a) where tagging varies ----------------
    uq = df["UserId"].astype(np.int64) * 100000 + df["QuizId"].astype(np.int64)
    struct = {}
    for lab, m in [("all_answers", np.ones(len(df), bool)),
                   ("reporters_answers", df.UserId.isin(rep.index).to_numpy())]:
        y = tag[m]
        blk = {"n_answers": int(m.sum()), "tag_share": float(y.mean())}
        for kname, kv in [("student", {"u": df.UserId.to_numpy()[m]}),
                          ("quiz", {"q": df.QuizId.to_numpy()[m]}),
                          ("question", {"i": df.QuestionId.to_numpy()[m]}),
                          ("student_x_quiz", {"uq": uq.to_numpy()[m]}),
                          ("month", {"mo": month.to_numpy()[m]})]:
            e, ng = eta2(y, kv)
            blk[f"eta2_{kname}"] = e
            blk[f"n_groups_{kname}"] = ng
        gq = pd.DataFrame({"uq": uq.to_numpy()[m], "t": y}).groupby("uq")["t"].agg(["mean", "size"])
        allu = gq["mean"] == 0
        allt = gq["mean"] == 1
        mixed = ~(allu | allt)
        blk["share_student_quiz_pairs_all_untagged"] = float(allu.mean())
        blk["share_student_quiz_pairs_all_tagged"] = float(allt.mean())
        blk["share_student_quiz_pairs_mixed"] = float(mixed.mean())
        blk["share_tagged_answers_in_mixed_pairs"] = float(
            (gq["mean"] * gq["size"])[mixed].sum() / (gq["mean"] * gq["size"]).sum())
        struct[lab] = blk
    res["tagging_structure"] = struct

    # ---------------- (1b) within-student paired comparisons ----------------
    first = df.groupby("UserId")["DateAnswered"].transform("min")
    df["days"] = (df["DateAnswered"] - first).dt.total_seconds() / 86400.0
    df["wkend"] = (df["DateAnswered"].dt.dayofweek >= 5).astype(np.float64)
    df["rel"] = df["IsCorrect"].astype(np.float64) - df["p_loo"]
    nt = per_u["ntag"]
    nu = per_u["n"] - per_u["ntag"]
    elig = per_u.index[(nt >= MIN_EACH) & (nu >= MIN_EACH)]
    sub = df[df.UserId.isin(elig)]
    res["paired_all_untagged"] = paired_block(sub, sid, True, "all untagged answers", 310)
    res["paired_all_untagged"]["eligibility"] = f">= {MIN_EACH} tagged and >= {MIN_EACH} untagged answers"

    # same-quiz basis: student x quiz pairs with mixed tagging
    gq = pd.DataFrame({"uq": uq, "t": tag}).groupby("uq")["t"].transform("mean").to_numpy()
    mixed_rows = (gq > 0) & (gq < 1)
    sq = df[mixed_rows].copy()
    sq = sq.sort_values(["UserId", "QuizId", "DateAnswered"], kind="stable")
    pos = sq.groupby(["UserId", "QuizId"]).cumcount().to_numpy().astype(np.float64)
    size = sq.groupby(["UserId", "QuizId"])["tag"].transform("size").to_numpy().astype(np.float64)
    sq["pos"] = pos / (size - 1)
    cnt = sq.groupby("UserId")["tag"].agg(["sum", "size"])
    el2 = cnt.index[(cnt["sum"] >= MIN_EACH) & ((cnt["size"] - cnt["sum"]) >= MIN_EACH)]
    sq = sq[sq.UserId.isin(el2)]
    res["paired_same_quiz"] = paired_block(sq, sid, True, "untagged answers in the same quiz attempts", 360)
    res["paired_same_quiz"]["eligibility"] = (f">= {MIN_EACH} tagged and >= {MIN_EACH} untagged answers "
                                              "within student x quiz pairs that mix tagged and untagged")
    res["paired_same_quiz"]["relpos_note"] = ("relative position within the student x quiz pair, 0 = first "
                                              "answer, 1 = last; order by DateAnswered (minute resolution), "
                                              "ties in original file order")

    # Design A pool (>= 200 tagged): difficulty composition
    poolA = per_u.index[per_u.ntag >= 200]
    pa = df[df.UserId.isin(poolA)]
    ga = pa.groupby(["UserId", "tag"])["p_loo"].mean().unstack("tag")
    res["designA_pool_difficulty"] = {
        "n_students": int(len(poolA)),
        "n_with_ge5_untagged": int(((per_u.loc[poolA, "n"] - per_u.loc[poolA, "ntag"]) >= MIN_EACH).sum()),
        "p_loo_tagged_minus_untagged": RC.summ((ga[True] - ga[False]).to_numpy(),
                                               seed=[RC.MASTER_SEED, 399, sid])}
    res["runtime_seconds"] = {"load": t_load, "total": time.time() - t0}
    os.makedirs(RC.OUT_DIR, exist_ok=True)
    with open(os.path.join(RC.OUT_DIR, f"r1_C_selection_{sample}.json"), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(sample, "selection done", round(time.time() - t0, 1), "s")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t34")
