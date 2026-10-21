"""Turn the extracted ChEMBL activities into a per-target classification set.

Labelling
    active    pChEMBL >= --active-thr  (and not a '>' censored measurement)
    inactive  pChEMBL <= --inactive-thr, or a '>' measurement below the
              active threshold (genuinely measured as weak)
    dropped   everything in the grey zone between the two thresholds

Negatives are the measured inactives topped up with "decoys": molecules
active at some other target but never measured against this one. Targets
need --min-actives actives to be kept.

Splitting is by Bemis-Murcko scaffold at the molecule level, so a scaffold
never appears in both train and test for any target.

Usage:
    python -I src/build_dataset.py --activities data/processed/activities.parquet \
        --outdir data/processed
"""

import argparse
import hashlib
import os

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")


def _scaffold(smiles: str, max_atoms: int) -> str | None:
    """Murcko scaffold SMILES, or None if the molecule is unusable."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumHeavyAtoms() > max_atoms or mol.GetNumHeavyAtoms() < 6:
        return None
    try:
        # Acyclic molecules give an empty scaffold; bucket them by themselves.
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
    except Exception:
        return None
    return scaf if scaf else Chem.MolToSmiles(mol)


def _test_split(scaffold: str, test_frac: float, seed: int) -> bool:
    """Deterministic scaffold -> train/test assignment (no global state)."""
    h = hashlib.blake2b(f"{seed}:{scaffold}".encode(), digest_size=8).digest()
    return (int.from_bytes(h, "big") % 10_000) < test_frac * 10_000


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--activities", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--active-thr", type=float, default=6.0)
    p.add_argument("--inactive-thr", type=float, default=5.0)
    p.add_argument("--min-actives", type=int, default=100)
    p.add_argument("--decoy-ratio", type=float, default=25.0,
                   help="decoys per active, on top of measured inactives")
    p.add_argument("--max-decoys", type=int, default=25000,
                   help="cap per target; keep it well above the typical active "
                        "count or the big targets end up majority-active, which "
                        "makes AUC/enrichment look far better than they are")
    p.add_argument("--max-heavy-atoms", type=int, default=70)
    p.add_argument("--test-frac", type=float, default=0.2)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--n-jobs", type=int, default=32)
    args = p.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = pd.read_parquet(args.activities)
    print(f"loaded {len(df):,} pairs", flush=True)

    # ---- labels -----------------------------------------------------------
    is_active = (df.pchembl >= args.active_thr) & (~df.any_censored_weak)
    is_inactive = (df.pchembl <= args.inactive_thr) | (
        df.any_censored_weak & (df.pchembl < args.active_thr)
    )
    df = df[is_active | is_inactive].copy()
    df["label"] = np.where(is_active[is_active | is_inactive], 1, 0)
    print(f"after thresholding: {len(df):,} pairs "
          f"({df.label.sum():,} active / {(1 - df.label).sum():,} inactive)", flush=True)

    # ---- molecules: scaffold + split --------------------------------------
    mols = df.drop_duplicates("mol_chembl_id")[["mol_chembl_id", "smiles"]].reset_index(drop=True)
    print(f"computing scaffolds for {len(mols):,} molecules ...", flush=True)
    scafs = Parallel(n_jobs=args.n_jobs, batch_size=2048)(
        delayed(_scaffold)(s, args.max_heavy_atoms) for s in mols.smiles
    )
    mols["scaffold"] = scafs
    mols = mols[mols.scaffold.notna()].reset_index(drop=True)
    mols["split"] = [
        "test" if _test_split(s, args.test_frac, args.seed) else "train"
        for s in mols.scaffold
    ]
    print(f"usable molecules: {len(mols):,} "
          f"({(mols.split == 'test').sum():,} in test, "
          f"{mols.scaffold.nunique():,} scaffolds)", flush=True)

    df = df[df.mol_chembl_id.isin(set(mols.mol_chembl_id))]

    # ---- target selection -------------------------------------------------
    n_act = df[df.label == 1].groupby("target_chembl_id").size()
    keep = set(n_act[n_act >= args.min_actives].index)
    df = df[df.target_chembl_id.isin(keep)].copy()
    print(f"targets with >={args.min_actives} actives: {len(keep):,}", flush=True)

    # ---- decoys -----------------------------------------------------------
    # Pool = anything active somewhere; sampling is weighted by nothing, i.e.
    # a plain random draw, which is the usual assumption for target prediction.
    pool = mols.mol_chembl_id.to_numpy()
    measured = df.groupby("target_chembl_id").mol_chembl_id.apply(set)
    rng = np.random.default_rng(args.seed)

    decoy_rows = []
    for tid, grp in df.groupby("target_chembl_id"):
        n_a = int((grp.label == 1).sum())
        n_want = min(int(n_a * args.decoy_ratio), args.max_decoys)
        seen = measured[tid]
        # Oversample then filter, cheaper than set-differencing the whole pool.
        cand = rng.choice(pool, size=min(len(pool), int(n_want * 1.3) + 100), replace=False)
        cand = [c for c in cand if c not in seen][:n_want]
        decoy_rows.append(pd.DataFrame({"target_chembl_id": tid,
                                        "mol_chembl_id": cand,
                                        "label": 0,
                                        "source": "decoy"}))
    df["source"] = np.where(df.label == 1, "measured_active", "measured_inactive")
    pairs = pd.concat([df[["target_chembl_id", "mol_chembl_id", "label", "source"]],
                       *decoy_rows], ignore_index=True)
    pairs = pairs.drop_duplicates(["target_chembl_id", "mol_chembl_id"])
    pairs = pairs.merge(mols[["mol_chembl_id", "split"]], on="mol_chembl_id")

    targets = (
        df[df.label == 1].drop_duplicates("target_chembl_id")[
            ["target_chembl_id", "target_name", "uniprot"]
        ]
        .merge(n_act.rename("n_actives").reset_index(), on="target_chembl_id")
        .sort_values("n_actives", ascending=False)
    )

    mols[mols.mol_chembl_id.isin(set(pairs.mol_chembl_id))].to_parquet(
        f"{args.outdir}/molecules.parquet", index=False)
    pairs.to_parquet(f"{args.outdir}/pairs.parquet", index=False)
    targets.to_parquet(f"{args.outdir}/targets.parquet", index=False)

    tr = pairs[pairs.split == "train"]
    print(f"\npairs: {len(pairs):,} over {pairs.target_chembl_id.nunique():,} targets",
          flush=True)
    print(f"  train {len(tr):,} ({tr.label.sum():,} active)", flush=True)
    print(f"  test  {len(pairs) - len(tr):,} "
          f"({pairs[pairs.split == 'test'].label.sum():,} active)", flush=True)
    print(pairs.source.value_counts().to_string(), flush=True)


if __name__ == "__main__":
    main()
