"""Validate Evo 2 feature chunks: coverage vs subset, uniqueness, NaNs, zero-shot delta_ll AUROC.
Usage: .venv-evo2/bin/python scripts/validate_evo2_features.py <features_dir> <clinvar_subset.tsv.gz>
"""
import glob
import sys

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

feat_dir, subset = sys.argv[1:3]
keys, dl, nan_emb = [], [], 0
for f in sorted(glob.glob(f"{feat_dir}/*.npz")):
    d = np.load(f, allow_pickle=True)
    keys.append(d["var_key"])
    dl.append(d["delta_ll"])
    nan_emb += int(np.isnan(d["emb"].astype(np.float32)).sum())  # streamed; embeddings not kept
keys, dl = np.concatenate(keys), np.concatenate(dl)
s = pd.read_csv(subset, sep="\t", dtype={"chrom": str}).set_index("var_key")
print(f"variants {len(keys):,} / subset {len(s):,}; unique {len(set(keys)):,}; missing from features "
      f"{len(set(s.index) - set(keys)):,}; NaN emb {nan_emb}; NaN delta_ll {int(np.isnan(dl).sum())}")
s = s.loc[keys].assign(d=dl)
print("\nZero-shot Evo 2 delta_ll AUROC (lower = pathogenic):")
for sp in ["val", "test", "train"]:
    g = s[s.split == sp]
    per = {c: round(roc_auc_score(x.label, -x.d), 3) for c, x in g.groupby("csq_group") if x.label.nunique() == 2}
    print(f"  {sp:5s} n={len(g):>7,} overall {roc_auc_score(g.label, -g.d):.3f}  {per}")
mva = s[s.genes.str.contains("BUB1B|CEP57|TRIP13")]
print(f"  MVA genes n={len(mva)} (P/LP {int(mva.label.sum())}) AUROC {roc_auc_score(mva.label, -mva.d):.3f}")
