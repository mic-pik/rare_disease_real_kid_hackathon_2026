"""Fit one Laplacian-modified naive Bayes model per target (one-vs-rest).

Each target gets its own model over only its own actives/inactives, so the
feature space stays small and the 700-odd fits parallelise cleanly.

Usage:
    python -I src/train_targets.py --data data/processed --fp data/processed/fp_r2.npz \
        --outdir models/lnb_r2 --n-jobs 64
"""

import os
import sys

# python -I drops the script's own directory from sys.path; this project's
# modules live beside this file, so put that directory back explicitly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import argparse
import json

import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from laplaciannb import LaplacianNB
from sklearn.metrics import average_precision_score, roc_auc_score

import model_store


def enrichment_factor(y_true: np.ndarray, score: np.ndarray, frac: float = 0.01) -> float:
    """Fold-enrichment of actives in the top `frac` of the ranking."""
    n_top = max(1, int(round(len(y_true) * frac)))
    hit_rate = y_true[np.argsort(-score)[:n_top]].mean()
    base = y_true.mean()
    return float(hit_rate / base) if base > 0 else float("nan")


def fit_one(tid: str, rows_tr: np.ndarray, y_tr: np.ndarray,
            rows_te: np.ndarray, y_te: np.ndarray, rows_bg: np.ndarray,
            fp_path: str, outdir: str, alpha: float, verify: bool) -> dict:
    X = sp.load_npz(fp_path)                      # memory-mapped by the OS page cache
    Xtr, Xte, Xbg = X[rows_tr], X[rows_te], X[rows_bg]

    model = LaplacianNB(alpha=alpha)
    model.fit(Xtr, y_tr)

    path = f"{outdir}/{tid}.npz"
    meta = {"target": tid, "n_train": int(len(y_tr)),
            "n_train_active": int(y_tr.sum()), "alpha": alpha}
    model_store.save(path, model, meta=meta)
    bg = model_store.LaplacianModel.load(path).decision_function(Xbg)
    model_store.save(path, model, meta=meta,
                     background=(float(bg.mean()), float(bg.std())),
                     background_scores=bg)
    saved = model_store.LaplacianModel.load(path)

    res = {"target_chembl_id": tid, "n_train": int(len(y_tr)),
           "n_train_active": int(y_tr.sum()), "n_test": int(len(y_te)),
           "n_test_active": int(y_te.sum()), "n_features": int(len(saved.feature_ids))}

    if verify and len(y_te):
        model_store.check_matches_library(model, saved, Xte[:min(50, len(y_te))])

    if len(y_te) and 0 < y_te.sum() < len(y_te):
        s = saved.decision_function(Xte)
        res |= {"roc_auc": float(roc_auc_score(y_te, s)),
                "pr_auc": float(average_precision_score(y_te, s)),
                "ef1": enrichment_factor(y_te, s, 0.01),
                "ef5": enrichment_factor(y_te, s, 0.05)}

    # Support tier, so a caller can tell a well-evidenced model from one backed
    # by 30 actives. Held-out AUC matters as much as raw count: a target with
    # 200 actives that still cannot separate its own test set is not reliable.
    n_act = int(y_tr.sum())
    auc = res.get("roc_auc")
    if n_act >= 100 and (auc is None or auc >= 0.85):
        tier = "high"
    elif n_act >= 50 and (auc is None or auc >= 0.75):
        tier = "medium"
    else:
        tier = "low"
    res["support_tier"] = tier

    # Record the tier inside the model file too, so the serving tool does not
    # depend on the metrics CSV being present and in sync.
    meta |= {"support_tier": tier, "roc_auc": auc}
    model_store.save(path, model, meta=meta,
                     background=(float(bg.mean()), float(bg.std())),
                     background_scores=bg)
    return res


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data", required=True, help="dir with pairs/molecules parquet")
    p.add_argument("--fp", required=True, help="fingerprint .npz")
    p.add_argument("--outdir", required=True)
    p.add_argument("--alpha", type=float, default=1.0)
    p.add_argument("--n-jobs", type=int, default=64)
    p.add_argument("--limit-targets", type=int, default=0, help="0 = all")
    p.add_argument("--max-negatives", type=int, default=50000,
                   help="cap negatives per target; PubChem screens push some "
                        "targets past 450k pairs and the fit loops in Python")
    p.add_argument("--background", type=int, default=3000,
                   help="random molecules used to calibrate cross-target scores")
    p.add_argument("--no-verify", action="store_true",
                   help="skip checking the compact scorer against laplaciannb")
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    pairs = pd.read_parquet(f"{args.data}/pairs.parquet")
    # The merged (ChEMBL+ExCAPE) build keys targets by gene symbol and
    # molecules by an internal id; accept either naming.
    pairs = pairs.rename(columns={"gene": "target_chembl_id",
                                  "mol_id": "mol_chembl_id"})
    fp_ids = np.load(args.fp.replace(".npz", "_ids.npy"), allow_pickle=True)
    row_of = pd.Series(np.arange(len(fp_ids)), index=fp_ids)
    pairs = pairs[pairs.mol_chembl_id.isin(row_of.index)].copy()
    pairs["row"] = row_of.loc[pairs.mol_chembl_id].to_numpy()

    targets = sorted(pairs.target_chembl_id.unique())
    if args.limit_targets:
        counts = pairs[pairs.label == 1].groupby("target_chembl_id").size()
        targets = counts.sort_values(ascending=False).index[:args.limit_targets].tolist()
    print(f"training {len(targets):,} targets with {args.n_jobs} workers", flush=True)

    rng = np.random.default_rng(0)
    rows_bg = rng.choice(len(fp_ids), size=min(args.background, len(fp_ids)),
                         replace=False)

    by_target = {t: g for t, g in pairs.groupby("target_chembl_id")}
    jobs = []
    sub_rng = np.random.default_rng(1)
    for t in targets:
        g = by_target[t]
        neg = g[g.label == 0]
        if len(neg) > args.max_negatives:
            drop = sub_rng.choice(neg.index, size=len(neg) - args.max_negatives,
                                  replace=False)
            g = g.drop(index=drop)
        tr, te = g[g.split == "train"], g[g.split == "test"]
        if tr.label.sum() < 10 or tr.label.nunique() < 2:
            continue
        jobs.append(delayed(fit_one)(
            t, tr.row.to_numpy(), tr.label.to_numpy(),
            te.row.to_numpy(), te.label.to_numpy(), rows_bg,
            args.fp, args.outdir, args.alpha, not args.no_verify))

    results = Parallel(n_jobs=args.n_jobs, verbose=5)(jobs)
    res = pd.DataFrame(results)
    res.to_csv(f"{args.outdir}/per_target_metrics.csv", index=False)

    with open(f"{args.outdir}/config.json", "w") as fh:
        json.dump(vars(args), fh, indent=2)

    print(f"\ntrained {len(res):,} models -> {args.outdir}", flush=True)
    if "roc_auc" in res:
        ev = res.dropna(subset=["roc_auc"])
        print(f"scaffold-split test metrics over {len(ev):,} targets (median):", flush=True)
        for m in ["roc_auc", "pr_auc", "ef1", "ef5"]:
            print(f"  {m:8s} {ev[m].median():.3f}   "
                  f"(IQR {ev[m].quantile(.25):.3f}-{ev[m].quantile(.75):.3f})", flush=True)


if __name__ == "__main__":
    main()
