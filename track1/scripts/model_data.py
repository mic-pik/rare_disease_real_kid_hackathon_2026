"""Assemble a training cache from ClinVar subset + VEP + ESM-1v / Evo 2 feature chunks.

Writes to <cache_dir>:
  meta.parquet   var_key, label, split, csq_group (ClinVar, reporting only), var_type, genes,
                 sample_weight, consequence/impact (VEP, used as features), delta_ll, esm scores
  evo_emb.npy    (N, 16384) fp16 [ref_site, alt_site, ref_mean, alt_mean], row-aligned with meta
  esm_emb.npy    (N, 5120)  fp16 [wt_site, mt_site, wt_mean, mt_mean], zeros where has_esm == 0

No gnomAD allele frequency is loaded: ClinVar benign labels are often frequency-based (label leakage).

Usage: .venv-evo2/bin/python scripts/model_data.py <clinvar_subset.tsv.gz> <vep.tsv> <esm_dir> <evo_dir> <cache_dir>
"""
import glob
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def read_vep(path, cols):
    txt = "".join(line for line in Path(path).read_text().splitlines(keepends=True) if not line.startswith("##")).lstrip("#")
    return pd.read_csv(io.StringIO(txt), sep="\t", dtype=str, na_values=["-"], usecols=cols)


def load_chunks(feat_dir, fields):
    parts = {f: [] for f in ["var_key", *fields]}
    for path in sorted(glob.glob(f"{feat_dir}/*.npz")):
        d = np.load(path, allow_pickle=True)
        for f in parts:
            parts[f].append(d[f])
    return {f: np.concatenate(v) for f, v in parts.items()}


def main():
    subset, vep_path, esm_dir, evo_dir, cache = sys.argv[1:6]
    cache = Path(cache)
    cache.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(subset, sep="\t", dtype={"chrom": str},
                       usecols=["var_key", "chrom", "pos", "ref", "alt", "label", "split", "csq_group",
                                "var_type", "genes", "sample_weight"])
    vep = read_vep(vep_path, ["Uploaded_variation", "Consequence", "IMPACT", "SYMBOL", "Gene"])
    vep = vep.rename(columns={"Uploaded_variation": "var_key", "Consequence": "consequence", "IMPACT": "impact",
                              "SYMBOL": "vep_symbol", "Gene": "vep_gene"})
    assert vep.var_key.is_unique
    meta = meta.merge(vep, on="var_key", how="left", validate="1:1")
    assert meta.consequence.notna().all(), "variants without VEP annotation"

    evo = load_chunks(evo_dir, ["delta_ll", "emb"])
    pos = pd.Series(np.arange(len(evo["var_key"])), index=evo["var_key"])
    idx = pos.reindex(meta.var_key).to_numpy()
    assert not np.isnan(idx.astype(float)).any(), "variants missing Evo 2 features"
    meta["delta_ll"] = evo["delta_ll"][idx]
    np.save(cache / "evo_emb.npy", evo["emb"][idx])
    del evo

    esm = load_chunks(esm_dir, ["scores", "emb"])
    pos = pd.Series(np.arange(len(esm["var_key"])), index=esm["var_key"])
    idx = pos.reindex(meta.var_key)
    has = idx.notna().to_numpy()
    meta["has_esm"] = has.astype(np.int8)
    scores = np.full((len(meta), esm["scores"].shape[1]), np.nan, dtype=np.float32)
    scores[has] = esm["scores"][idx[has].astype(int)]
    for i in range(scores.shape[1]):
        meta[f"esm_score_{i + 1}"] = scores[:, i]
    meta["esm_ens"] = scores.mean(1)
    emb = np.zeros((len(meta), esm["emb"].shape[1]), dtype=np.float16)
    emb[has] = esm["emb"][idx[has].astype(int)]
    np.save(cache / "esm_emb.npy", emb)

    meta.to_parquet(cache / "meta.parquet", index=False)
    print(f"{len(meta):,} variants; has_esm {has.sum():,}; splits {meta.split.value_counts().to_dict()}")
    print("cache:", sorted(p.name for p in cache.iterdir()))


if __name__ == "__main__":
    main()
