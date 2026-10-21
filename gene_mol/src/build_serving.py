"""Build the serving manifest for the two-list predict_targets tool.

    core      targets of the v1 ChEMBL-only build (913 ChEMBL targets, 911
              genes). Ranked on their own, exactly as the LLM was trained to read.
    extended  every other target of the merged ChEMBL+ExCAPE build. Ranked
              separately, so a noisy small model can never push a core hit down.
    excluded  extended targets whose model cannot separate its own held-out
              test set (ROC-AUC < --min-auc). These are mostly PubChem
              screening readouts (GLP1R, GNAS, GMNN, ATAD5, ...): thousands of
              "actives" from reporter assays, not a ligand-binding signal.

Each row names the model file that serves it, so the core can come from either
build (--core v1 or --core merged) without the tool knowing the difference.

Usage:
    python -I src/build_serving.py --core merged --merged models/lnb_merged3 \
        --out models/serving_v2
"""

from pathlib import Path
import argparse
import glob
import json
import os

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
KG = str(Path(__file__).resolve().parents[2] / "mva_llm/data/kg/genes.parquet")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--core", choices=["v1", "merged"], required=True)
    p.add_argument("--v1", default="models/lnb_r2")
    p.add_argument("--merged", default="models/lnb_merged3")
    p.add_argument("--min-auc", type=float, default=0.75)
    p.add_argument("--out", default="models/serving_v2")
    args = p.parse_args()
    os.chdir(ROOT)
    os.makedirs(args.out, exist_ok=True)

    kg = pd.read_parquet(KG)[["uniprot", "symbol", "name"]].dropna(subset=["symbol"])
    kg["symbol"] = kg.symbol.str.upper()
    sym = dict(zip(kg.uniprot, kg.symbol))
    kg_by_sym = kg.drop_duplicates("symbol").set_index("symbol")

    # ChEMBL identity per gene: the v1 target if there is one, else the ChEMBL
    # single-protein target with the most measurements for that UniProt.
    v1t = pd.read_parquet("data/processed/targets.parquet")
    v1_files = {os.path.basename(f)[:-4] for f in glob.glob(f"{args.v1}/*.npz")}
    v1t = v1t[v1t.target_chembl_id.isin(v1_files)].copy()
    v1t["gene"] = v1t.uniprot.map(sym)
    acts = pd.read_parquet("data/processed/activities.parquet",
                           columns=["target_chembl_id", "target_name", "uniprot"])
    chembl_t = (acts.groupby(["target_chembl_id", "target_name", "uniprot"]).size()
                .rename("n").reset_index())
    chembl_t["gene"] = chembl_t.uniprot.map(sym)
    chembl_t = chembl_t.sort_values("n", ascending=False).drop_duplicates("gene")
    chembl_by_gene = chembl_t.set_index("gene")

    mm = pd.read_csv(f"{args.merged}/per_target_metrics.csv").rename(
        columns={"target_chembl_id": "gene"})
    v1m = pd.read_csv(f"{args.v1}/per_target_metrics.csv")
    core_genes = set(v1t.gene.dropna())

    rows = []
    if args.core == "v1":
        v1m = v1m.set_index("target_chembl_id")
        for t in v1t.itertuples():
            r = v1m.loc[t.target_chembl_id] if t.target_chembl_id in v1m.index else None
            rows.append(dict(
                list="core", gene=t.gene, model_file=f"{args.v1}/{t.target_chembl_id}.npz",
                target_chembl_id=t.target_chembl_id, target_name=t.target_name,
                uniprot=t.uniprot, n_train_active=None if r is None else r.n_train_active,
                support_tier=None, roc_auc=None if r is None else r.roc_auc, reason=None))
    else:
        mmi = mm.set_index("gene")
        v1_by_gene = v1t.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
        for g in sorted(core_genes & set(mm.gene)):
            t, r = v1_by_gene.loc[g], mmi.loc[g]
            rows.append(dict(
                list="core", gene=g, model_file=f"{args.merged}/{g}.npz",
                target_chembl_id=t.target_chembl_id, target_name=t.target_name,
                uniprot=t.uniprot, n_train_active=r.n_train_active,
                support_tier=r.support_tier, roc_auc=r.roc_auc, reason=None))

    for r in mm[~mm.gene.isin(core_genes)].itertuples():
        c = chembl_by_gene.loc[r.gene] if r.gene in chembl_by_gene.index else None
        k = kg_by_sym.loc[r.gene] if r.gene in kg_by_sym.index else None
        ok = pd.notna(r.roc_auc) and r.roc_auc >= args.min_auc
        rows.append(dict(
            list="extended" if ok else "excluded", gene=r.gene,
            model_file=f"{args.merged}/{r.gene}.npz",
            target_chembl_id=None if c is None else c.target_chembl_id,
            target_name=(c.target_name if c is not None else
                         (k["name"] if k is not None else r.gene)),
            uniprot=(c.uniprot if c is not None else
                     (k.uniprot if k is not None else None)),
            n_train_active=r.n_train_active, support_tier=r.support_tier,
            roc_auc=r.roc_auc,
            reason=None if ok else
            (f"held-out ROC-AUC {r.roc_auc:.2f} < {args.min_auc}" if pd.notna(r.roc_auc)
             else "no held-out evaluation")))

    man = pd.DataFrame(rows)
    assert man.model_file.map(os.path.exists).all(), "manifest names a missing model"
    man.to_csv(f"{args.out}/manifest.csv", index=False)
    cfg = vars(args) | {"counts": man.list.value_counts().to_dict()}
    with open(f"{args.out}/config.json", "w") as fh:
        json.dump(cfg, fh, indent=2)
    print(man.list.value_counts().to_string())
    ex = man[man.list == "excluded"].sort_values("n_train_active", ascending=False)
    print("\nexcluded (largest first):")
    print(ex[["gene", "n_train_active", "roc_auc", "support_tier"]].head(25).to_string(index=False))
    print(f"\nwrote {args.out}/manifest.csv")


if __name__ == "__main__":
    main()
