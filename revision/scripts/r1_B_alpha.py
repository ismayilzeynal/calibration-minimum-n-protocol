"""Stream B, part 3: coefficient alpha on a balanced (crossed) subset.

Balanced subset = for one quiz (QuizId != -1, the fill value in common.py),
the students who gave a usable confidence rating on EVERY question of that
quiz (question set = distinct QuestionId answered in that quiz across all
answers of the sample). These students form a complete students x items
matrix, so a crossed p x i design exists and alpha is defined. The subset is
selected: only students who rated every item of the quiz are in it.

Item scores as in r1_B_vc.py: d = c - x (bias), q = (c - x)^2 (Brier).
Per quiz and index:
  alpha             Cronbach alpha (ddof 1)
  icc_C_k           two-way consistency ICC(C,k) = (MSp - MSres) / MSp (== alpha; self-check)
  phi               G-theory absolute coefficient sigma_p^2 / (sigma_p^2 + (sigma_i^2 + sigma_pi,e^2) / k)
  icc1_k            one-way ICC(1,k) (items nested in students; item variance counted as error)
  sb_per_student    the paper's estimator: per-student random partitions (reliability.half_indices), 50 splits, mean SB
  sb_common         one common random item partition for all students, 50 partitions, mean SB
  rulon_all_mean    mean Flanagan-Rulon over ALL equal splits (even k <= 14 only; equals alpha exactly)
  n70_from_<m>      Spearman-Brown projection of coefficient m from k items to .70: k * (.7/.3) * (1 - m) / m
Student bootstrap B = 1000 per quiz (splits and partitions held fixed).
Eligibility: >= 100 complete students and k >= 4 items. Primary set: >= 300.

Usage:  python -B r1_B_alpha.py t34   (or t12)
Writes: revision-2026-11/results/r1_B_parts/alpha_<sample>.json
"""
import itertools
import json
import os
import sys
import time

REV = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER = os.path.dirname(REV)
OUTDIR = os.path.join(REV, "results", "r1_B_parts")

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(PAPER, "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import common  # noqa: E402
import indices as IX  # noqa: E402
import reliability as REL  # noqa: E402

common.DATA_DIR = os.path.join(PAPER, "data")

B_BOOT = 1000
R = REL.R_SPLITS
MIN_STUD = 100
PRIMARY_STUD = 300
MIN_K = 4
KEYS = ("bias", "brier")


def alpha_of(Y):
    k = Y.shape[1]
    return k / (k - 1) * (1 - Y.var(axis=0, ddof=1).sum() / Y.sum(axis=1).var(ddof=1))


def anova2(Y):
    n, k = Y.shape
    M = Y.mean()
    rm = Y.mean(axis=1)
    cm = Y.mean(axis=0)
    ssp = k * ((rm - M) ** 2).sum()
    ssi = n * ((cm - M) ** 2).sum()
    sst = ((Y - M) ** 2).sum()
    ssr = sst - ssp - ssi
    msp = ssp / (n - 1)
    msi = ssi / (k - 1)
    msr = ssr / ((n - 1) * (k - 1))
    msw = (ssi + ssr) / (n * (k - 1))
    sp2 = (msp - msr) / k
    si2 = (msi - msr) / n
    return dict(icc_C_k=(msp - msr) / msp, icc1_k=(msp - msw) / msp,
                phi=sp2 / (sp2 + (si2 + msr) / k) if sp2 > 0 else np.nan,
                sigma_p2=sp2, sigma_i2=si2, sigma_pie2=msr)


def sb_cols(A, Bm):
    """Column-wise Pearson + SB (A, B: n x S, no NaN) -> mean SB."""
    Az = A - A.mean(0)
    Bz = Bm - Bm.mean(0)
    r = (Az * Bz).sum(0) / np.sqrt((Az ** 2).sum(0) * (Bz ** 2).sum(0))
    return float(np.mean(2 * r / (1 + r)))


def rulon_all(Y):
    n, k = Y.shape
    if k % 2 or k > 14:
        return None
    vals = []
    for comb in itertools.combinations(range(1, k), k // 2 - 1):
        a = [0] + list(comb)
        b = [j for j in range(k) if j not in a]
        sa, sb_ = Y[:, a].sum(1), Y[:, b].sum(1)
        vals.append(2 * (1 - (sa.var(ddof=1) + sb_.var(ddof=1)) / (sa + sb_).var(ddof=1)))
    return float(np.mean(vals))


def run(sample):
    t0 = time.time()
    sid = {"t12": 1, "t34": 2}[sample]
    df = common.build_responses(sample)
    df = df[df.QuizId != -1]
    qpq = df.groupby("QuizId").QuestionId.nunique().rename("k")
    tg = common.tagged(df)
    g = tg.groupby(["QuizId", "UserId"]).agg(nq=("QuestionId", "nunique"),
                                             n=("QuestionId", "size")).reset_index()
    g = g.merge(qpq.reset_index(), on="QuizId")
    dup_share = float((g.n > g.nq).mean())
    comp = g[(g.nq == g.k) & (g.n == g.nq)]
    per_quiz = comp.groupby("QuizId").agg(nstud=("UserId", "size"),
                                          k=("k", "first")).sort_values("nstud", ascending=False)
    feas = dict(
        n_quizzes=int(len(qpq)),
        quiz_user_pairs_with_duplicate_answers_share=dup_share,
        max_complete_students=int(per_quiz.nstud.max()) if len(per_quiz) else 0,
        k_at_max=int(per_quiz.k.iloc[0]) if len(per_quiz) else None,
        n_quizzes_ge100=int(((per_quiz.nstud >= 100) & (per_quiz.k >= MIN_K)).sum()),
        n_quizzes_ge300=int(((per_quiz.nstud >= 300) & (per_quiz.k >= MIN_K)).sum()),
        top10=[dict(QuizId=int(q), nstud=int(r.nstud), k=int(r.k))
               for q, r in per_quiz.head(10).iterrows()],
    )
    elig = per_quiz[(per_quiz.nstud >= MIN_STUD) & (per_quiz.k >= MIN_K)]
    tgc = tg.merge(comp[["QuizId", "UserId"]], on=["QuizId", "UserId"])
    quizzes = []
    for qid, info in elig.iterrows():
        sub = tgc[tgc.QuizId == qid]
        items = np.sort(sub.QuestionId.unique())
        uids = np.sort(sub.UserId.unique())
        piv_c = sub.pivot(index="UserId", columns="QuestionId", values="Confidence").loc[uids, items]
        piv_x = sub.pivot(index="UserId", columns="QuestionId", values="IsCorrect").loc[uids, items]
        C = (piv_c.to_numpy() / 25.0).round().astype(np.int8)
        X = piv_x.to_numpy().astype(np.int8)
        n, k = C.shape
        c = IX.LEVELS[C.astype(np.int64)]
        Ys = {"bias": c - X, "brier": (c - X) ** 2}
        rng = np.random.default_rng([REL.MASTER_SEED, 103, sid, int(qid)])
        A, Bh, _, _, full = REL.half_indices(C, X, rng)
        crng = np.random.default_rng([REL.MASTER_SEED, 104, sid, int(qid)])
        parts = [crng.permutation(k) for _ in range(R)]
        row = dict(QuizId=int(qid), n_students=int(n), k_items=int(k))
        brng = np.random.default_rng([REL.MASTER_SEED, 105, sid, int(qid)])
        bidx = [brng.integers(0, n, n) for _ in range(B_BOOT)]
        for key in KEYS:
            Y = Ys[key]
            if np.max(np.abs(Y.mean(1) - full[key])) > 1e-12:
                raise RuntimeError("item scores do not reproduce indices.py full-set value")
            cA = np.stack([Y[:, p[:k // 2]].mean(1) for p in parts], axis=1)
            cB = np.stack([Y[:, p[k // 2:]].mean(1) for p in parts], axis=1)
            psA, psB = A[key], Bh[key]
            if np.isnan(psA).any() or np.isnan(psB).any():
                raise RuntimeError("NaN in additive half index")
            av = anova2(Y)
            est = dict(alpha=float(alpha_of(Y)), icc_C_k=float(av["icc_C_k"]),
                       phi=float(av["phi"]), icc1_k=float(av["icc1_k"]),
                       sb_per_student=float(REL.sb_mean(psA, psB)),
                       sb_common=sb_cols(cA, cB),
                       rulon_all_mean=rulon_all(Y),
                       sigma_p2=float(av["sigma_p2"]), sigma_i2=float(av["sigma_i2"]),
                       sigma_pie2=float(av["sigma_pie2"]),
                       share_item_var_of_error=float(av["sigma_i2"] / (av["sigma_i2"] + av["sigma_pie2"])))
            for m in ("alpha", "phi", "icc1_k", "sb_per_student"):
                rv = est[m]
                est[f"n70_from_{m}"] = float(k * (0.70 / 0.30) * (1 - rv) / rv) if rv > 0 else None
            est["alpha_minus_sb_per_student"] = est["alpha"] - est["sb_per_student"]
            est["phi_minus_sb_per_student"] = est["phi"] - est["sb_per_student"]
            est["alpha_minus_icc1_k"] = est["alpha"] - est["icc1_k"]
            est["sb_per_student_minus_icc1_k"] = est["sb_per_student"] - est["icc1_k"]
            est["sb_common_minus_alpha"] = est["sb_common"] - est["alpha"]
            bt = {m: np.empty(B_BOOT) for m in ("alpha", "icc1_k", "sb_per_student",
                                                  "sb_common", "alpha_minus_sb_per_student")}
            for b, ix in enumerate(bidx):
                Yb = Y[ix]
                a = alpha_of(Yb)
                av_b = anova2(Yb)
                sp = REL.sb_mean(psA[ix], psB[ix])
                bt["alpha"][b] = a
                bt["icc1_k"][b] = av_b["icc1_k"]
                bt["sb_per_student"][b] = sp
                bt["sb_common"][b] = sb_cols(cA[ix], cB[ix])
                bt["alpha_minus_sb_per_student"][b] = a - sp
            for m, v in bt.items():
                est[f"{m}_ci"] = [float(z) for z in np.percentile(v, [2.5, 97.5])]
            row[key] = est
        quizzes.append(row)
    t_q = time.time() - t0

    def summ(rows, key):
        if not rows:
            return None
        f = lambda m: np.array([r[key][m] for r in rows])  # noqa: E731
        w = np.array([r["n_students"] for r in rows], float)
        out = dict(n_quizzes=len(rows),
                   n_students_total=int(w.sum()),
                   k_median=float(np.median([r["k_items"] for r in rows])),
                   k_range=[int(min(r["k_items"] for r in rows)), int(max(r["k_items"] for r in rows))])
        for m in ("alpha", "phi", "icc1_k", "sb_per_student", "sb_common",
                  "alpha_minus_sb_per_student", "alpha_minus_icc1_k",
                  "sb_per_student_minus_icc1_k", "sb_common_minus_alpha",
                  "share_item_var_of_error", "n70_from_alpha", "n70_from_phi",
                  "n70_from_icc1_k", "n70_from_sb_per_student",
                  "phi_minus_sb_per_student"):
            v = f(m)
            out[f"{m}_median"] = float(np.median(v))
            out[f"{m}_min"] = float(v.min())
            out[f"{m}_max"] = float(v.max())
            out[f"{m}_nweighted_mean"] = float(np.sum(w * v) / w.sum())
        out["n_quizzes_alpha_gt_sb_per_student"] = int(np.sum(f("alpha") > f("sb_per_student")))
        out["n_quizzes_diff_ci_excludes_0"] = int(sum(
            1 for r in rows if r[key]["alpha_minus_sb_per_student_ci"][0] > 0
            or r[key]["alpha_minus_sb_per_student_ci"][1] < 0))
        out["max_abs_alpha_minus_icc_C_k"] = float(np.max(np.abs(f("alpha") - f("icc_C_k"))))
        ru = [(r[key]["rulon_all_mean"], r[key]["alpha"]) for r in rows
              if r[key]["rulon_all_mean"] is not None]
        out["n_quizzes_rulon_checked"] = len(ru)
        out["max_abs_rulon_minus_alpha"] = float(max(abs(a - b) for a, b in ru)) if ru else None
        return out

    prim = [r for r in quizzes if r["n_students"] >= PRIMARY_STUD]
    res = dict(
        sample=sample, feasibility=feas,
        eligibility=dict(min_complete_students=MIN_STUD, primary_min_students=PRIMARY_STUD,
                         min_items=MIN_K),
        seeds=dict(per_student_splits=f"default_rng([42, 103, {sid}, QuizId])",
                   common_partitions=f"default_rng([42, 104, {sid}, QuizId])",
                   bootstrap=f"default_rng([42, 105, {sid}, QuizId])"),
        R_splits=R, B=B_BOOT,
        summary_primary_ge300={k: summ(prim, k) for k in KEYS},
        summary_all_ge100={k: summ(quizzes, k) for k in KEYS},
        quizzes=quizzes,
        runtime_seconds=dict(total=time.time() - t0, quizzes=t_q),
    )
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, f"alpha_{sample}.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(f"[{sample}] eligible quizzes {len(quizzes)} (primary {len(prim)}); "
          f"{res['runtime_seconds']['total']:.1f}s")


if __name__ == "__main__":
    run(sys.argv[1] if len(sys.argv) > 1 else "t34")
