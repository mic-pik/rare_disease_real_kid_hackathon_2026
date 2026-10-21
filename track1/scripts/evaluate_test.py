"""Final held-out test evaluation: chosen model(s) vs zero-shot baselines, per consequence group and MVA genes.

Run once, after model selection on val. Predictions come from the saved preds_*.parquet files.

Usage: .venv-evo2/bin/python scripts/evaluate_test.py <cache_dir> <out_tsv> <preds.parquet>:<model column> [...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from model_common import evaluate  # noqa: E402

cache, out_tsv, specs = sys.argv[1], sys.argv[2], sys.argv[3:]
meta = pd.read_parquet(Path(cache) / "meta.parquet")
te = (meta.split == "test").to_numpy()
tm = meta[te].reset_index(drop=True)

rows = evaluate(tm, -tm.delta_ll.to_numpy(), "zero-shot Evo2 delta_ll")
m = tm.has_esm.to_numpy().astype(bool)
rows += evaluate(tm[m].reset_index(drop=True), -tm.esm_ens.to_numpy()[m], "zero-shot ESM-1v (missense only)")
for spec in specs:
    path, col = spec.rsplit(":", 1)
    p = pd.read_parquet(path, columns=["var_key", col]).set_index("var_key")[col]
    rows += evaluate(tm, p.loc[tm.var_key].to_numpy(), col)

res = pd.DataFrame(rows)
res.to_csv(out_tsv, sep="\t", index=False)
print("=== TEST AUROC ===")
print(res.pivot_table(index="model", columns="group", values="auroc", sort=False).round(3).to_string())
print("\n=== TEST AUPRC ===")
print(res.pivot_table(index="model", columns="group", values="auprc", sort=False).round(3).to_string())
print("\nn / n_pos per group:", res.drop_duplicates("group").set_index("group")[["n", "n_pos"]].to_dict("index"))
