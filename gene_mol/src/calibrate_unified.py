"""Can one calibrated ranking replace the two-list panel?

The tiered tool ranks core and extended targets separately because raw
z-scores are not comparable across models: a small model with 30 actives
reaches z = 10 on noise, and a flat cut-off either floods the output or hides
real hits. Two lists avoid the problem instead of solving it, and leave the
molecule tool surfacing a true extended target for only 77 of 150 molecules
(§4.6-4.7 of the report).

This asks whether a *calibrated* score fixes it properly. For each hit we
estimate P(this target is a measured active for this molecule) from the z
score and the model's own support, fit on held-out molecules. Probabilities
are comparable across models by construction, so core and extended can go into
one ranking.

Compared, on the same held-out molecules:
    z-ranked, two lists   what the tool does now
    z-ranked, one list    the naive merge, which dilution ruins
    calibrated, one list  this proposal

Fit and evaluation use disjoint molecule halves, so the calibration is never
scored on what it was fitted to.

Usage:
    python -I src/calibrate_unified.py --n-queries 20000 --n-jobs 32
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
from sklearn.isotonic import IsotonicRegression

import model_store

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _score(files, X):
    return [model_store.LaplacianModel.load(f).zscore(X) for f in files]


def recall_at(rank, ks=(1, 3, 5, 10)):
    return {f"r@{k}": float(np.mean(rank <= k)) for k in ks}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--serving", default="models/serving_v2")
    p.add_argument("--data", default="data/processed_merged3")
    p.add_argument("--n-queries", type=int, default=20000)
    p.add_argument("--n-jobs", type=int, default=32)
    p.add_argument("--z-show", type=float, default=7.75)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--out", default="models/serving_v2/unified.json")
    args = p.parse_args()
    os.chdir(ROOT)

    man = pd.read_csv(f"{args.serving}/manifest.csv")
    man = man[man.list.isin(["core", "extended"])].reset_index(drop=True)
    is_core = (man.list == "core").to_numpy()
    genes = man.gene.to_numpy()
    n_act = man.n_train_active.fillna(30).to_numpy()
    auc = man.roc_auc.fillna(0.9).to_numpy()

    pairs = pd.read_parquet(f"{args.data}/pairs.parquet",
                            columns=["gene", "mol_id", "label", "split"])
    act = pairs[(pairs.label == 1) & (pairs.split == "test")]
    truth = act.groupby("mol_id").gene.apply(set)
    rng = np.random.default_rng(args.seed)
    q = rng.choice(np.array(sorted(truth.index)),
                   size=min(args.n_queries, len(truth)), replace=False)
    ids = np.load(f"{args.data}/fp_r2_ids.npy", allow_pickle=True)
    row = pd.Series(np.arange(len(ids)), index=ids)
    X = sp.load_npz(f"{args.data}/fp_r2.npz")[row.loc[q].to_numpy()]
    print(f"{len(q):,} held-out query molecules", flush=True)

    files = [os.path.join(ROOT, f) for f in man.model_file]
    shards = np.array_split(np.arange(len(files)), args.n_jobs)
    out = Parallel(n_jobs=args.n_jobs)(
        delayed(_score)([files[i] for i in s], X) for s in shards)
    Z = np.column_stack([z for o in out for z in o])
    print(f"scored {Z.shape[0]:,} x {Z.shape[1]:,}", flush=True)

    # Truth matrix, molecules x targets.
    gidx = {g: i for i, g in enumerate(genes)}
    Y = np.zeros_like(Z, dtype=bool)
    for r, m in enumerate(q):
        for g in truth[m]:
            j = gidx.get(g)
            if j is not None:
                Y[r, j] = True

    half = len(q) // 2
    fit, ev = np.arange(half), np.arange(half, len(q))

    # Calibrate within support bands: a z of 8 means something different for a
    # model with 30 actives than for one with 3,000, and binning by support
    # keeps each isotonic fit well populated.
    band = np.digitize(np.log10(np.maximum(n_act, 1)), [1.7, 2.3, 3.0])
    cal = {}
    for b in np.unique(band):
        cols = np.flatnonzero(band == b)
        zz = Z[np.ix_(fit, cols)].ravel()
        yy = Y[np.ix_(fit, cols)].ravel()
        keep = zz > 1.0                      # calibrate where it could matter
        iso = IsotonicRegression(y_min=0, y_max=1, out_of_bounds="clip")
        iso.fit(zz[keep], yy[keep].astype(float))
        cal[int(b)] = iso
        print(f"  band {b}: {len(cols):4d} targets, {int(keep.sum()):,} fit points, "
              f"P at z=4/8/12: {iso.predict([4,8,12]).round(3)}", flush=True)

    P = np.zeros_like(Z)
    for b, iso in cal.items():
        cols = np.flatnonzero(band == b)
        P[:, cols] = iso.predict(Z[:, cols].ravel()).reshape(len(q), len(cols))

    def rank_of(score, rows, mask=None, blocked=None):
        """Best rank of a true target per row, under a scoring matrix."""
        S = score[rows].copy()
        if mask is not None:
            S[:, ~mask] = -np.inf
        if blocked is not None:
            S[blocked[rows]] = -np.inf
        order = np.argsort(-S, axis=1)
        T = Y[rows]
        out = np.full(len(rows), S.shape[1] + 1)
        for i in range(len(rows)):
            hit = np.flatnonzero(T[i][order[i]])
            if len(hit):
                out[i] = hit[0] + 1
        return out

    # Molecules whose only true targets are extended: the molecule_ext shape.
    ext_only = np.array([not any(is_core[gidx[g]] for g in truth[m] if g in gidx)
                         and any(g in gidx for g in truth[m]) for m in q])
    print(f"\nextended-only molecules in eval half: {int(ext_only[ev].sum())}", flush=True)

    # Scheme 1: two lists, as served today. Core ranked by z; an extended hit
    # is visible only above the show bar.
    hidden = (~np.tile(is_core, (len(q), 1))) & (Z < args.z_show)
    res = {}
    res["z two lists (current)"] = rank_of(Z, ev, blocked=hidden)
    res["z one list (naive merge)"] = rank_of(Z, ev)
    res["calibrated one list"] = rank_of(P, ev)
    res["z core only"] = rank_of(Z, ev, mask=is_core)
    # Keep the two-list schema the LLM is trained to read, but choose and rank
    # the extended entries by calibrated probability instead of a flat z bar.
    for pmin in (0.10, 0.20, 0.30):
        hid = (~np.tile(is_core, (len(q), 1))) & (P < pmin)
        S = np.where(np.tile(is_core, (len(q), 1)), Z, P * 0 + Z)   # core by z
        res[f"two lists, P-gate {pmin:.2f}"] = rank_of(S, ev, blocked=hid)

    print("\nall held-out molecules")
    for k, r in res.items():
        print(f"  {k:28s} {recall_at(r)}")

    e = ev[ext_only[ev]]
    if len(e):
        print(f"\nextended-only molecules (n={len(e)}) - the tool's ceiling")
        schemes = [("z two lists (current)", Z, hidden),
                   ("z one list (naive merge)", Z, None),
                   ("calibrated one list", P, None)]
        for pmin in (0.10, 0.20, 0.30):
            schemes.append((f"two lists, P-gate {pmin:.2f}", Z,
                            (~np.tile(is_core, (len(q), 1))) & (P < pmin)))
        for k, sc, blocked in schemes:
            r = rank_of(sc, e, blocked=blocked)
            print(f"  {k:28s} {recall_at(r)}  surfaced@10 {np.mean(r <= 10):.1%}")

    # What the tool actually hands the model is a *displayed set*: the top core
    # hits plus up to EXT_TOP extended ones that clear a bar. Rank is the wrong
    # metric for that; this asks whether a true target appears at all.
    EXT_TOP, CORE_TOP = 5, 10
    core_cols = np.flatnonzero(is_core)
    ext_cols = np.flatnonzero(~is_core)

    def displayed(rows, ext_score, ext_gate, gate_val):
        hit = np.zeros(len(rows), dtype=bool)
        for n, r in enumerate(rows):
            zc = Z[r, core_cols]
            top_c = core_cols[np.argsort(-zc)[:CORE_TOP]]
            es, eg = ext_score[r, ext_cols], ext_gate[r, ext_cols]
            ok = np.flatnonzero(eg >= gate_val)
            top_e = ext_cols[ok[np.argsort(-es[ok])[:EXT_TOP]]] if len(ok) else []
            shown = np.concatenate([top_c, np.asarray(top_e, dtype=int)])
            hit[n] = Y[r, shown].any()
        return hit.mean()

    print("\nwhat the tool DISPLAYS (top 10 core + up to 5 extended)")
    for label, rows in [("all held-out", ev), ("extended-only", ev[ext_only[ev]])]:
        cur = displayed(rows, Z, Z, args.z_show)
        print(f"  {label:16s} current z>={args.z_show:<5g}      {cur:6.1%}")
        for pmin in (0.05, 0.10, 0.20):
            v = displayed(rows, P, P, pmin)
            print(f"  {label:16s} rank+gate by P>={pmin:<4.2f} {v:6.1%}")

    # Persist the isotonic curves as knots so the serving tool can interpolate
    # them with numpy and needs no sklearn at inference.
    bands = {}
    for b, iso in cal.items():
        xs = np.asarray(iso.X_thresholds_, dtype=float)
        ys = np.asarray(iso.y_thresholds_, dtype=float)
        if len(xs) > 256:                      # thin out, keeping the shape
            idx = np.unique(np.linspace(0, len(xs) - 1, 256).astype(int))
            xs, ys = xs[idx], ys[idx]
        bands[str(b)] = {"z": xs.tolist(), "p": ys.tolist()}
    with open(args.out, "w") as fh:
        json.dump({"bands": bands,
                   "band_edges_log10_actives": [1.7, 2.3, 3.0],
                   "n_queries": int(len(q)),
                   "ext_gate_p": 0.05,
                   "note": ("isotonic P(target is a measured active) vs z, per "
                            "support band, fit on a disjoint half of the "
                            "held-out molecules")},
                  fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
