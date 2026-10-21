"""Calibrate the z bars for the extended list of the tiered tool.

Extended models are smaller and noisier, so z = 4 does not mean what it means
for a core model. On held-out compounds (merged build test split, with at least
one measured active) this measures how often the top core hit and the top
extended hit are a measured active, binned by z.

    show       lowest z where extended precision >= half of core precision at z >= 4
    confident  lowest z where it matches core@4, or reaches 50% if it never does
               (it does not: extended models top out around 60-65%)

Precision is a lower bound: an untested compound-target pair counts as wrong.

Usage:
    python -I src/calibrate_extended.py --serving models/serving_v2 \
        --data data/processed_merged3
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import argparse
import json

import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed

import model_store

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _score(files, X):
    return [model_store.LaplacianModel.load(f).zscore(X) for f in files]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--serving", default="models/serving_v2")
    p.add_argument("--data", default="data/processed_merged3")
    p.add_argument("--n-queries", type=int, default=20000)
    p.add_argument("--n-jobs", type=int, default=32)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    os.chdir(ROOT)

    man = pd.read_csv(f"{args.serving}/manifest.csv")
    man = man[man.list.isin(["core", "extended"])].reset_index(drop=True)
    pairs = pd.read_parquet(f"{args.data}/pairs.parquet",
                            columns=["gene", "mol_id", "label", "split"])
    act = pairs[(pairs.label == 1)]
    truth = act[act.split == "test"].groupby("mol_id").gene.apply(set)
    rng = np.random.default_rng(args.seed)
    q = rng.choice(np.array(sorted(truth.index)),
                   size=min(args.n_queries, len(truth)), replace=False)
    ids = np.load(f"{args.data}/fp_r2_ids.npy", allow_pickle=True)
    row = pd.Series(np.arange(len(ids)), index=ids)
    X = sp.load_npz(f"{args.data}/fp_r2.npz")[row.loc[q].to_numpy()]

    files = [os.path.join(ROOT, f) for f in man.model_file]
    shards = np.array_split(np.arange(len(files)), args.n_jobs)
    out = Parallel(n_jobs=args.n_jobs)(
        delayed(_score)([files[i] for i in s], X) for s in shards)
    Z = np.column_stack([z for o in out for z in o])
    genes = man.gene.to_numpy()

    res = {}
    for lst in ("core", "extended"):
        cols = np.flatnonzero(man.list.to_numpy() == lst)
        best = cols[np.argmax(Z[:, cols], axis=1)]
        z = Z[np.arange(len(q)), best]
        ok = np.array([genes[b] in truth[m] for b, m in zip(best, q)])
        res[lst] = (z, ok)
        print(f"\n{lst}: top-1 hit is a measured active, by z of that hit")
        for lo, hi in [(-99, 2.5), (2.5, 4), (4, 5), (5, 6), (6, 7), (7, 8),
                       (8, 10), (10, 99)]:
            s = (z >= lo) & (z < hi)
            if s.sum():
                print(f"  z {lo:>5g}-{hi:<4g} n={s.sum():6d}  precision {ok[s].mean():6.1%}")

    zc, okc = res["core"]
    target = okc[zc >= 4].mean()
    ze, oke = res["extended"]
    grid = np.arange(2.5, 15.01, 0.25)
    prec = [(t, oke[ze >= t].mean(), int((ze >= t).sum())) for t in grid
            if (ze >= t).sum() >= 50]
    conf = next((t for t, pr, n in prec if pr >= target), None)
    if conf is None:
        conf = next((t for t, pr, n in prec if pr >= 0.5), None)
    show = next((t for t, pr, n in prec if pr >= 0.5 * target), None)
    print(f"\ncore precision at z >= 4: {target:.1%}")
    for t, pr, n in prec[::2]:
        print(f"  extended z >= {t:5.2f}: precision {pr:6.1%} (n={n})")
    print(f"\n-> Z_EXT_CONFIDENT = {conf}  (extended precision >= core@4, else >= 50%)")
    print(f"-> Z_EXT_SHOW      = {show}  (extended precision >= half of core@4)")
    at = dict((t, pr) for t, pr, n in prec)
    with open(f"{args.serving}/calibration.json", "w") as fh:
        json.dump({"core_precision_at_4": float(target),
                   "z_ext_confident": conf, "z_ext_show": show,
                   "ext_precision_at_confident": float(at[conf]) if conf else None,
                   "ext_precision_at_show": float(at[show]) if show else None,
                   "n_queries": int(len(q)),
                   "note": "precision = top hit is a measured active (lower bound)"},
                  fh, indent=2)


if __name__ == "__main__":
    main()
