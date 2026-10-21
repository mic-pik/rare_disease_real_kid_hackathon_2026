"""Logistic regression + LightGBM on scalar features and PCA-reduced embeddings; feature-set ablation.

Selection on val only. Test predictions are saved but not printed (see evaluate_test.py).

Usage: .venv-evo2/bin/python scripts/train_classic.py <cache_dir> <out_dir>
"""
import sys
import time
from pathlib import Path

import joblib
import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).parent))
from model_common import consequence_vocab, emb_blocks, evaluate, load_cache, scalar_features  # noqa: E402

PCA_DIMS = {"evo_site_diff": 64, "evo_mean_diff": 64, "evo_ref_site": 32,
            "esm_site_diff": 32, "esm_mean_diff": 32, "esm_wt_site": 16}
PCA_FIT_ROWS = 60_000
FEATURE_SETS = {
    "scalars": [],
    "scalars+esm": ["esm_site_diff", "esm_mean_diff", "esm_wt_site"],
    "scalars+evo": ["evo_site_diff", "evo_mean_diff", "evo_ref_site"],
    "scalars+evo+esm": list(PCA_DIMS),
}
SEED = 0

cache, out = sys.argv[1], Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
t0 = time.time()
meta, evo, esm = load_cache(cache)
vocab = consequence_vocab(meta)
X_scalar = scalar_features(meta, vocab)
split = meta.split.to_numpy()
y = meta.label.to_numpy()
has_esm = meta.has_esm.to_numpy().astype(bool)
tr, va, te = (split == "train"), (split == "val"), (split == "test")
print(f"loaded in {time.time() - t0:.0f}s; {X_scalar.shape[1]} scalar features", flush=True)

# PCA per embedding block, fit on a train subsample (ESM blocks only on missense rows)
rng = np.random.default_rng(SEED)
pcs, pcas = {}, {}
blocks_all = emb_blocks(evo, esm, np.arange(len(meta)))
for name, k in PCA_DIMS.items():
    fit_pool = np.where(tr & (has_esm if name.startswith("esm") else True))[0]
    fit_rows = rng.choice(fit_pool, min(PCA_FIT_ROWS, len(fit_pool)), replace=False)
    pca = PCA(k, svd_solver="randomized", random_state=SEED).fit(blocks_all[name][fit_rows])
    z = pca.transform(blocks_all[name]).astype(np.float32)
    if name.startswith("esm"):
        z[~has_esm] = 0
    pcs[name], pcas[name] = z, pca
    print(f"PCA {name}: {k} comps, {pca.explained_variance_ratio_.sum():.2f} var explained", flush=True)
del blocks_all

results, preds = [], {"var_key": meta.var_key, "split": split, "label": y}
for fs_name, blocks in FEATURE_SETS.items():
    cols = list(X_scalar.columns) + [f"{b}_pc{i}" for b in blocks for i in range(PCA_DIMS[b])]
    X = np.hstack([X_scalar.to_numpy()] + [pcs[b] for b in blocks]).astype(np.float32)

    scaler = StandardScaler().fit(X[tr])
    Xs = scaler.transform(X)
    best = None
    for C in [0.01, 0.1, 1.0]:
        lr = LogisticRegression(C=C, max_iter=3000).fit(Xs[tr], y[tr])
        p = lr.predict_proba(Xs[va])[:, 1]
        auc = evaluate(meta[va].reset_index(drop=True), p, "")[0]["auroc"]
        if best is None or auc > best[0]:
            best = (auc, C, lr)
    name = f"LR[{fs_name}]"
    p_all = best[2].predict_proba(Xs)[:, 1]
    preds[name] = p_all
    results += [dict(r, C=best[1]) for r in evaluate(meta[va].reset_index(drop=True), p_all[va], name)]
    joblib.dump({"model": best[2], "scaler": scaler, "cols": cols}, out / f"lr_{fs_name}.joblib")
    print(f"{name}: val AUROC {best[0]:.4f} (C={best[1]})  [{time.time() - t0:.0f}s]", flush=True)

    gbm = lgb.LGBMClassifier(n_estimators=3000, learning_rate=0.03, num_leaves=63, min_child_samples=50,
                             subsample=0.8, subsample_freq=1, colsample_bytree=0.5, reg_lambda=1.0,
                             n_jobs=32, random_state=SEED, verbose=-1)
    gbm.fit(X[tr], y[tr], eval_set=[(X[va], y[va])], eval_metric="auc",
            callbacks=[lgb.early_stopping(100, verbose=False)])
    name = f"LGBM[{fs_name}]"
    p_all = gbm.predict_proba(X)[:, 1]
    preds[name] = p_all
    r = evaluate(meta[va].reset_index(drop=True), p_all[va], name)
    results += [dict(x, n_trees=gbm.best_iteration_) for x in r]
    joblib.dump({"model": gbm, "cols": cols}, out / f"lgbm_{fs_name}.joblib")
    print(f"{name}: val AUROC {r[0]['auroc']:.4f} ({gbm.best_iteration_} trees)  [{time.time() - t0:.0f}s]", flush=True)

joblib.dump({"pcas": pcas, "pca_dims": PCA_DIMS, "vocab": vocab}, out / "preprocess.joblib")
pd.DataFrame(preds).to_parquet(out / "preds_classic.parquet", index=False)
res = pd.DataFrame(results)
res.to_csv(out / "val_results_classic.tsv", sep="\t", index=False)

# Zero-shot references on val
zs = []
v = meta[va].reset_index(drop=True)
zs += evaluate(v, -v.delta_ll.to_numpy(), "zero-shot Evo2 delta_ll")
m = v.has_esm.to_numpy().astype(bool)
zs += evaluate(v[m].reset_index(drop=True), -v.esm_ens.to_numpy()[m], "zero-shot ESM-1v (missense)")
res = pd.concat([pd.DataFrame(zs), res])
print("\n=== Validation AUROC by ClinVar consequence group ===")
print(res.pivot_table(index="model", columns="group", values="auroc", sort=False).round(3).to_string())
