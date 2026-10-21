"""Apply the trained ensemble (LGBM + MLP) to new variants, and rank the patient's candidates.

  --check : re-score ClinVar val variants from the training cache and compare with the saved
            predictions (must match), to prove this code path reproduces training preprocessing.
  default : score patient candidates (candidates.tsv.gz + ESM-1v / Evo 2 feature chunks), write
            ranked variant and gene tables.

Thresholds come from the ClinVar val set (specificity 95% / 99%), never from the patient.

Usage:
  .venv-evo2/bin/python scripts/score_patient.py --check
  .venv-evo2/bin/python scripts/score_patient.py --candidates work/patient/candidates.tsv.gz \
      --esm work/features/patient_esm1v --evo work/features/patient_evo2 --out work/patient/scored
"""
import argparse
import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from torch import nn

sys.path.insert(0, str(Path(__file__).parent))
from model_common import emb_blocks, load_cache, scalar_features  # noqa: E402
from model_data import load_chunks  # noqa: E402

MODEL_DIR = Path("work/model")
LGBM_NAME = "scalars+evo+esm"
MLP_SEEDS = [0, 1, 2]
dev = "cuda" if torch.cuda.is_available() else "cpu"


def predict(meta, evo, esm):
    """meta needs consequence, impact, var_type, ref, alt, delta_ll, esm_score_1..5, esm_ens, has_esm."""
    pre = joblib.load(MODEL_DIR / "classic/preprocess.joblib")
    X_scalar = scalar_features(meta, pre["vocab"]).to_numpy()
    blocks = emb_blocks(evo, esm, np.arange(len(meta)))
    has_esm = meta.has_esm.to_numpy().astype(bool)

    pcs = []
    for name in pre["pca_dims"]:  # same order as training (FEATURE_SETS["scalars+evo+esm"])
        z = pre["pcas"][name].transform(blocks[name]).astype(np.float32)
        if name.startswith("esm"):
            z[~has_esm] = 0
        pcs.append(z)
    lgbm = joblib.load(MODEL_DIR / f"classic/lgbm_{LGBM_NAME}.joblib")
    X = np.hstack([X_scalar] + pcs).astype(np.float32)
    assert X.shape[1] == len(lgbm["cols"]), (X.shape, len(lgbm["cols"]))
    p_lgbm = lgbm["model"].predict_proba(X)[:, 1]

    p_mlp = np.zeros(len(meta))
    for seed in MLP_SEEDS:
        ck = torch.load(MODEL_DIR / f"mlp/mlp_all_seed{seed}.pt", weights_only=False, map_location=dev)
        Xm = np.hstack([X_scalar] + [blocks[b] for b in ck["blocks"]]).astype(np.float32)
        Xm = torch.from_numpy((Xm - ck["mu"]) / ck["sd"]).to(dev, torch.float16)  # as stored in training
        layers, d = [], Xm.shape[1]
        for h in ck["hidden"]:
            layers += [nn.Linear(d, h), nn.GELU(), nn.Dropout(ck["dropout"])]
            d = h
        model = nn.Sequential(*layers, nn.Linear(d, 1)).to(dev)
        model.load_state_dict(ck["state"])
        model.eval()
        with torch.no_grad():
            p = torch.cat([torch.sigmoid(model(Xm[i:i + 8192].float()).squeeze(1)) for i in range(0, len(Xm), 8192)])
        p_mlp += p.cpu().numpy() / len(MLP_SEEDS)
    return pd.DataFrame({"p_lgbm": p_lgbm, "p_mlp": p_mlp, "score": (p_lgbm + p_mlp) / 2}, index=meta.index)


def val_thresholds():
    pr = pd.read_parquet(MODEL_DIR / "preds_ensemble.parquet")
    benign = pr[(pr.split == "val") & (pr.label == 0)]["ENS[LGBM+MLP]"]
    return {"spec95": float(benign.quantile(0.95)), "spec99": float(benign.quantile(0.99))}


def check():
    meta, evo, esm = load_cache(MODEL_DIR / "cache")
    rows = np.where(meta.split == "val")[0][::10]
    p = predict(meta.iloc[rows].reset_index(drop=True), evo[rows], esm[rows])
    saved = pd.read_parquet(MODEL_DIR / "preds_ensemble.parquet").set_index("var_key").loc[meta.var_key.iloc[rows]]
    for col, ref in [("p_lgbm", "LGBM[scalars+evo+esm]"), ("p_mlp", "MLP[scalars+all]"), ("score", "ENS[LGBM+MLP]")]:
        diff = np.abs(p[col].to_numpy() - saved[ref].to_numpy())
        print(f"{col:7s} vs saved {ref:22s}: max |diff| {diff.max():.2e} over {len(rows):,} val variants")
    print("thresholds from val benign:", {k: round(v, 4) for k, v in val_thresholds().items()})


def score(args):
    cand = pd.read_csv(args.candidates, sep="\t", dtype={"chrom": str})
    evo = load_chunks(args.evo, ["delta_ll", "emb"])
    ei = pd.Series(np.arange(len(evo["var_key"])), index=evo["var_key"]).reindex(cand.var_key)
    assert ei.notna().all(), f"{ei.isna().sum()} candidates without Evo 2 features"
    ei = ei.to_numpy().astype(int)
    cand["delta_ll"] = evo["delta_ll"][ei]
    evo_emb = evo["emb"][ei]
    del evo

    esm = load_chunks(args.esm, ["scores", "emb"])
    si = pd.Series(np.arange(len(esm["var_key"])), index=esm["var_key"]).reindex(cand.var_key)
    has = si.notna().to_numpy()
    cand["has_esm"] = has.astype(np.int8)
    sc = np.full((len(cand), 5), np.nan, dtype=np.float32)
    sc[has] = esm["scores"][si[has].astype(int)]
    for i in range(5):
        cand[f"esm_score_{i + 1}"] = sc[:, i]
    cand["esm_ens"] = sc.mean(1)
    esm_emb = np.zeros((len(cand), 5120), dtype=np.float16)
    esm_emb[has] = esm["emb"][si[has].astype(int)]

    cand = cand.join(predict(cand, evo_emb, esm_emb))
    thr = val_thresholds()
    cand["tier"] = np.select([cand.score >= thr["spec99"], cand.score >= thr["spec95"]], ["high", "moderate"], "low")
    cand = cand.sort_values("score", ascending=False)
    cand["rank"] = np.arange(1, len(cand) + 1)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    cand.to_csv(out / "variants_scored.tsv.gz", sep="\t", index=False)

    # Gene level, recessive model: homozygous damaging variant, or >= 2 damaging hets (phase unknown)
    dmg = cand[cand.tier != "low"]
    g = cand.groupby("genes").agg(n_candidates=("var_key", "size"), max_score=("score", "max"),
                                  in_mva_gene=("in_mva_gene", "any"))
    g["n_damaging"] = dmg.groupby("genes").size().reindex(g.index, fill_value=0)
    g["second_score"] = cand.groupby("genes").score.apply(lambda s: s.nlargest(2).iloc[-1] if len(s) > 1 else 0)
    g["hom_damaging"] = dmg[dmg.zygosity == "hom"].groupby("genes").size().reindex(g.index, fill_value=0)
    g["biallelic_candidate"] = (g.hom_damaging > 0) | (g.n_damaging >= 2)
    # Rank recessive candidates by the weaker of the two alleles (or the hom allele)
    g["recessive_score"] = np.where(g.hom_damaging > 0, g.max_score, g.second_score)
    g = g.sort_values(["biallelic_candidate", "recessive_score"], ascending=False)
    g.to_csv(out / "genes_scored.tsv", sep="\t")

    cols = ["rank", "genes", "consequence", "hgvsp", "hgvsc", "zygosity", "score", "tier", "delta_ll", "esm_ens", "max_af", "dp", "gq"]
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 45)
    print(f"Scored {len(cand):,} candidates; thresholds (val specificity 95% / 99%): "
          f"{thr['spec95']:.3f} / {thr['spec99']:.3f}; tiers {cand.tier.value_counts().to_dict()}")
    print("\n=== Top 25 variants ===")
    print(cand[cols].head(25).round(3).to_string(index=False))
    print(f"\n=== Biallelic candidate genes (hom damaging or >= 2 damaging hets): {int(g.biallelic_candidate.sum())} ===")
    print(g[g.biallelic_candidate].head(25).round(3).to_string())
    print("\n=== MVA genes (all candidates) ===")
    print(cand[cand.in_mva_gene][cols].round(3).to_string(index=False))


ap = argparse.ArgumentParser()
ap.add_argument("--check", action="store_true")
ap.add_argument("--candidates")
ap.add_argument("--esm")
ap.add_argument("--evo")
ap.add_argument("--out")
args = ap.parse_args()
check() if args.check else score(args)
