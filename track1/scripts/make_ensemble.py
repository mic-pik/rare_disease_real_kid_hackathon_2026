"""Build the final ensemble (mean of LightGBM and MLP probabilities) and compare it with its parts on val.

Model choice is made on val only; the test split is evaluated afterwards with evaluate_test.py.

Usage: .venv-evo2/bin/python scripts/make_ensemble.py <cache_dir> <model_dir>
Writes <model_dir>/preds_ensemble.parquet (var_key, split, label, LGBM, MLP, ENS columns).
"""
import sys
from pathlib import Path

import pandas as pd
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).parent))
from model_common import evaluate  # noqa: E402

LGBM, MLP, ENS = "LGBM[scalars+evo+esm]", "MLP[scalars+all]", "ENS[LGBM+MLP]"

cache, model_dir = Path(sys.argv[1]), Path(sys.argv[2])
meta = pd.read_parquet(cache / "meta.parquet")
c = pd.read_parquet(model_dir / "classic/preds_classic.parquet").set_index("var_key")
m = pd.read_parquet(model_dir / "mlp/preds_mlp_all.parquet").set_index("var_key")

preds = c[["split", "label", LGBM]].join(m[MLP])
# Probability average (not rank average) so the ensemble applies to new variants without re-ranking
preds[ENS] = (preds[LGBM] + preds[MLP]) / 2
preds.reset_index().to_parquet(model_dir / "preds_ensemble.parquet", index=False)

va = meta[meta.split == "val"].reset_index(drop=True)
v = preds.loc[va.var_key]
rows = (evaluate(va, v[LGBM].to_numpy(), LGBM) + evaluate(va, v[MLP].to_numpy(), MLP)
        + evaluate(va, (rankdata(v[LGBM]) + rankdata(v[MLP])) / 2, "rank-avg LGBM+MLP")
        + evaluate(va, v[ENS].to_numpy(), ENS))
print("=== Validation AUROC ===")
print(pd.DataFrame(rows).pivot_table(index="model", columns="group", values="auroc", sort=False).round(3).to_string())
