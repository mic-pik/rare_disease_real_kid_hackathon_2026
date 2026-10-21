"""Shared feature construction and evaluation for the pathogenicity models.

Features (never gnomAD AF, never gene identity):
  scalars  Evo 2 delta_ll, ESM-1v scores (5 + ensemble, 0 if not missense) + has_esm,
           VEP consequence terms (multi-hot), VEP impact, variant type, indel length
  evo      Evo 2 blocks: alt_site - ref_site, alt_mean - ref_mean, ref_site (3 x 4096)
  esm      ESM-1v blocks: mt_site - wt_site, mt_mean - wt_mean, wt_site (3 x 1280), 0 if not missense
"""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

MIN_TERM_COUNT = 20  # consequence terms rarer than this in train are dropped from the vocabulary
REPORT_GROUPS = ["missense", "intronic", "splice_site", "synonymous", "utr", "inframe_indel", "other",
                 "frameshift", "stop_start"]
MVA_GENES = ("BUB1B", "CEP57", "TRIP13")


def load_cache(cache):
    cache = Path(cache)
    meta = pd.read_parquet(cache / "meta.parquet")
    evo = np.load(cache / "evo_emb.npy", mmap_mode="r")
    esm = np.load(cache / "esm_emb.npy", mmap_mode="r")
    return meta, evo, esm


def consequence_vocab(meta):
    terms = meta.loc[meta.split == "train", "consequence"].str.split(",").explode().value_counts()
    return sorted(terms[terms >= MIN_TERM_COUNT].index)


def scalar_features(meta, vocab):
    """Returns DataFrame of scalar features (unscaled). Works for ClinVar and patient meta alike."""
    terms = meta.consequence.str.split(",")
    f = pd.DataFrame({"delta_ll": meta.delta_ll.astype(np.float32)}, index=meta.index)
    for i in range(1, 6):
        f[f"esm_score_{i}"] = meta[f"esm_score_{i}"].fillna(0).astype(np.float32)
    f["esm_ens"] = meta.esm_ens.fillna(0).astype(np.float32)
    f["has_esm"] = meta.has_esm.astype(np.float32)
    for t in vocab:
        f[f"csq_{t}"] = terms.map(lambda ts, t=t: t in ts).astype(np.float32)
    for imp in ["HIGH", "MODERATE", "LOW", "MODIFIER"]:
        f[f"impact_{imp}"] = (meta.impact == imp).astype(np.float32)
    for vt in ["SNV", "deletion", "insertion", "delins"]:
        f[f"vtype_{vt}"] = (meta.var_type == vt).astype(np.float32)
    f["indel_len"] = (meta.alt.str.len() - meta.ref.str.len()).astype(np.float32)
    f["frame_shift_len"] = (f.indel_len.abs() % 3 != 0).astype(np.float32)
    return f


def emb_blocks(evo, esm, rows):
    """Embedding blocks for the given row indices, float32: dict name -> (n, d)."""
    e = np.asarray(evo[rows], dtype=np.float32)
    s = np.asarray(esm[rows], dtype=np.float32)
    d, k = 4096, 1280
    return {
        "evo_site_diff": e[:, d:2 * d] - e[:, :d],
        "evo_mean_diff": e[:, 3 * d:] - e[:, 2 * d:3 * d],
        "evo_ref_site": e[:, :d],
        "esm_site_diff": s[:, k:2 * k] - s[:, :k],
        "esm_mean_diff": s[:, 3 * k:] - s[:, 2 * k:3 * k],
        "esm_wt_site": s[:, :k],
    }


def is_mva(meta):
    return meta.genes.str.split("|").map(lambda gs: any(g in MVA_GENES for g in gs))


def evaluate(meta, pred, name):
    """AUROC/AUPRC overall, per ClinVar consequence group and for MVA genes. Returns list of dict rows."""
    rows = []

    def add(group, m):
        y = meta.label[m]
        if y.nunique() == 2:
            rows.append({"model": name, "group": group, "n": int(m.sum()), "n_pos": int(y.sum()),
                         "auroc": roc_auc_score(y, pred[m]), "auprc": average_precision_score(y, pred[m])})

    add("overall", np.ones(len(meta), bool))
    for g in REPORT_GROUPS:
        add(g, (meta.csq_group == g).to_numpy())
    add("MVA_genes", is_mva(meta).to_numpy())
    return rows
