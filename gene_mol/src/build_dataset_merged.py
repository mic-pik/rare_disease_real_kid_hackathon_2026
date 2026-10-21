"""Rebuild the target-prediction dataset with measured inactives.

The existing dataset labels most negatives by assumption: a molecule active
somewhere else and never tested here is called inactive. ExCAPE-DB supplies
~42M *measured* inactives from ChEMBL20 and PubChem BioAssay, so for many
targets the presumed decoys can be replaced by real ones.

Targets are keyed by gene symbol rather than ChEMBL target id, because that is
the only identifier both sources share.

Molecules are keyed by InChIKey (also used as mol_id). ChEMBL and ExCAPE write
SMILES differently, so keying on the raw string (the merged2 build) stored 365k
compounds twice: truth split across copies, ~5k compounds on both sides of the
split, and a compound active under one id could be sampled as a decoy under the
other. Each compound keeps ChEMBL's SMILES when it has one, so its scaffold
(and split) matches the v1 ChEMBL-only build.

Label policy
    active    ChEMBL pChEMBL >= --active-thr, or ExCAPE Activity_Flag = A
    inactive  ChEMBL pChEMBL <= --inactive-thr, or ExCAPE Activity_Flag = N
    decoy     only used to top a target up to --min-negatives, and recorded as
              a separate source so its effect can be measured

Usage:
    python -I src/build_dataset_merged.py --activities data/processed/activities.parquet \
        --excape data/excape_raw/excape_v2.tsv.xz --outdir data/processed_merged
"""

from pathlib import Path
import argparse
import hashlib
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

KG = str(Path(__file__).resolve().parents[2] / "mva_llm/data/kg/genes.parquet")


def _scaffold(smiles: str, max_atoms: int) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or not (6 <= mol.GetNumHeavyAtoms() <= max_atoms):
        return None
    try:
        s = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return None
    return s or Chem.MolToSmiles(mol)


def _inchikeys(smiles: list[str]) -> list[str | None]:
    RDLogger.DisableLog("rdApp.*")
    out = []
    for s in smiles:
        mol = Chem.MolFromSmiles(s)
        try:
            out.append(Chem.MolToInchiKey(mol) or None if mol is not None else None)
        except Exception:
            out.append(None)
    return out


def _test_split(scaffold: str, frac: float, seed: int) -> bool:
    h = hashlib.blake2b(f"{seed}:{scaffold}".encode(), digest_size=8).digest()
    return (int.from_bytes(h, "big") % 10_000) < frac * 10_000


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--activities", required=True)
    p.add_argument("--excape", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--active-thr", type=float, default=6.0)
    p.add_argument("--inactive-thr", type=float, default=5.0)
    p.add_argument("--min-actives", type=int, default=30)
    p.add_argument("--min-negatives", type=int, default=300,
                   help="floor for negatives per target")
    p.add_argument("--negative-ratio", type=float, default=25.0,
                   help="target negatives = ratio x actives. A flat floor made "
                        "the median target ~50%% active, which caps enrichment "
                        "factor near 2 and does not reflect deployment, where "
                        "actives are rare.")
    p.add_argument("--max-negatives", type=int, default=25000)
    p.add_argument("--max-heavy-atoms", type=int, default=70)
    p.add_argument("--test-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-jobs", type=int, default=32)
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # --- ChEMBL 37 side ------------------------------------------------------
    ch = pd.read_parquet(args.activities)
    g = pd.read_parquet(KG)[["uniprot", "symbol"]].dropna()
    sym = dict(zip(g.uniprot, g.symbol.str.upper()))
    ch["gene"] = ch.uniprot.map(sym)
    ch = ch[ch.gene.notna()].copy()
    ch["label"] = np.where(
        (ch.pchembl >= args.active_thr) & (~ch.any_censored_weak), 1,
        np.where((ch.pchembl <= args.inactive_thr)
                 | (ch.any_censored_weak & (ch.pchembl < args.active_thr)), 0, -1))
    ch = ch[ch.label >= 0]
    ch = ch[["gene", "smiles", "label"]].assign(source="chembl37")
    print(f"ChEMBL37: {len(ch):,} labelled pairs "
          f"({int((ch.label == 1).sum()):,} active)", flush=True)

    # --- ExCAPE side ---------------------------------------------------------
    print("reading ExCAPE ...", flush=True)
    ex = pd.read_csv(args.excape, sep="\t", low_memory=False,
                     usecols=["Gene_Symbol", "Activity_Flag", "Tax_ID", "SMILES"])
    ex = ex[(ex.Tax_ID == 9606) & ex.Gene_Symbol.notna() & ex.SMILES.notna()]
    ex["gene"] = ex.Gene_Symbol.str.upper()
    ex["label"] = (ex.Activity_Flag == "A").astype(int)
    ex = ex[["gene", "SMILES", "label"]].rename(columns={"SMILES": "smiles"})
    ex["source"] = "excape"
    print(f"ExCAPE human: {len(ex):,} pairs "
          f"({int((ex.label == 1).sum()):,} active, "
          f"{int((ex.label == 0).sum()):,} measured inactive)", flush=True)

    both = pd.concat([ch, ex], ignore_index=True)

    # --- molecule identity: InChIKey, ChEMBL spelling preferred --------------
    uniq = pd.DataFrame({"smiles": both.smiles.unique()})
    uniq["from_chembl"] = uniq.smiles.isin(set(ch.smiles))
    print(f"computing InChIKeys for {len(uniq):,} distinct SMILES ...", flush=True)
    chunks = np.array_split(uniq.smiles.to_numpy(), args.n_jobs * 8)
    uniq["ik"] = sum(Parallel(n_jobs=args.n_jobs)(
        delayed(_inchikeys)(list(c)) for c in chunks), [])
    uniq = uniq[uniq.ik.notna()]
    canon = (uniq.sort_values(["from_chembl", "smiles"], ascending=[False, True])
             .drop_duplicates("ik").set_index("ik").smiles)
    rep_smiles = uniq.ik.map(canon)
    print(f"{len(uniq):,} parseable SMILES -> {len(canon):,} compounds "
          f"({len(uniq) - len(canon):,} alternative spellings merged)", flush=True)
    both["smiles"] = both.smiles.map(dict(zip(uniq.smiles, rep_smiles)))
    both = both[both.smiles.notna()]
    del uniq, rep_smiles
    # A pair measured by both: trust an active call over an inactive one, since
    # inactivity is assay-specific and a confirmed active is harder to produce.
    # Among equal labels prefer the ChEMBL record (stable sort, ChEMBL rows come
    # first in `both`), so `source` is meaningful for source-filtered evals.
    both = (both.sort_values("label", ascending=False, kind="stable")
            .drop_duplicates(["gene", "smiles"], keep="first"))
    print(f"merged: {len(both):,} unique gene-molecule pairs", flush=True)

    # --- molecules: scaffold + split ----------------------------------------
    mols = both.drop_duplicates("smiles")[["smiles"]].reset_index(drop=True)
    print(f"computing scaffolds for {len(mols):,} molecules ...", flush=True)
    scafs = Parallel(n_jobs=args.n_jobs, batch_size=4096)(
        delayed(_scaffold)(s, args.max_heavy_atoms) for s in mols.smiles)
    mols["scaffold"] = scafs
    mols = mols[mols.scaffold.notna()].reset_index(drop=True)
    mols["split"] = ["test" if _test_split(s, args.test_frac, args.seed)
                     else "train" for s in mols.scaffold]
    mols["mol_id"] = mols.smiles.map(canon.reset_index().set_index("smiles").ik)
    assert mols.mol_id.notna().all() and mols.mol_id.is_unique
    print(f"usable molecules: {len(mols):,} "
          f"({mols.scaffold.nunique():,} scaffolds)", flush=True)

    both = both.merge(mols[["smiles", "mol_id", "split"]], on="smiles")

    n_act = both[both.label == 1].groupby("gene").mol_id.nunique()
    keep = set(n_act[n_act >= args.min_actives].index)
    both = both[both.gene.isin(keep)]
    print(f"targets with >={args.min_actives} actives: {len(keep):,}", flush=True)

    # --- decoys, only where measured negatives are thin ----------------------
    pool = mols.mol_id.to_numpy()
    measured_neg = both[both.label == 0].groupby("gene").mol_id.nunique()
    seen = both.groupby("gene").mol_id.apply(set)
    rows = []
    for gene in sorted(keep):
        have = int(measured_neg.get(gene, 0))
        n_a = int(n_act.get(gene, 0))
        target_n = min(args.max_negatives,
                       max(args.min_negatives, int(args.negative_ratio * n_a)))
        want = max(0, target_n - have)
        if want <= 0:
            continue
        cand = rng.choice(pool, size=min(len(pool), int(want * 1.3) + 50),
                          replace=False)
        s = seen.get(gene, set())
        cand = [c for c in cand if c not in s][:want]
        rows.append(pd.DataFrame({"gene": gene, "mol_id": cand, "label": 0,
                                  "source": "decoy"}))
    decoys = (pd.concat(rows, ignore_index=True) if rows
              else pd.DataFrame(columns=["gene", "mol_id", "label", "source"]))
    decoys = decoys.merge(mols[["mol_id", "split"]], on="mol_id")
    print(f"decoys added: {len(decoys):,} "
          f"(targets needing top-up: {decoys.gene.nunique() if len(decoys) else 0:,})",
          flush=True)

    pairs = pd.concat([both[["gene", "mol_id", "label", "source", "split"]],
                       decoys], ignore_index=True)
    pairs = pairs.drop_duplicates(["gene", "mol_id"])

    targets = (pairs[pairs.label == 1].groupby("gene").mol_id.nunique()
               .rename("n_actives").reset_index())
    targets["n_measured_inactive"] = targets.gene.map(
        pairs[(pairs.label == 0) & (pairs.source != "decoy")]
        .groupby("gene").mol_id.nunique()).fillna(0).astype(int)
    targets["n_decoy"] = targets.gene.map(
        pairs[pairs.source == "decoy"].groupby("gene").mol_id.nunique()
    ).fillna(0).astype(int)

    mols[mols.mol_id.isin(set(pairs.mol_id))].to_parquet(
        f"{args.outdir}/molecules.parquet", index=False)
    pairs.to_parquet(f"{args.outdir}/pairs.parquet", index=False)
    targets.sort_values("n_actives", ascending=False).to_parquet(
        f"{args.outdir}/targets.parquet", index=False)

    print(f"\npairs {len(pairs):,} over {pairs.gene.nunique():,} targets",
          flush=True)
    print(pairs.source.value_counts().to_string(), flush=True)
    frac = 1 - (targets.n_decoy.sum() / max(1, targets.n_decoy.sum()
                                            + targets.n_measured_inactive.sum()))
    print(f"\nnegatives that are measured rather than presumed: {frac:.1%}",
          flush=True)


if __name__ == "__main__":
    main()
