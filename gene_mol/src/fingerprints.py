"""Unfolded ECFP fingerprints as a CSR matrix over the full 2**32 bit space.

laplaciannb.rdkit_to_csr does exactly this, but single-threaded and via Python
lists, which is slow for ~10^6 molecules. This is the same bit convention
(rdFingerprintGenerator sparse Morgan on-bits reinterpreted as uint32) chunked
across processes, so the output is interchangeable with the library's.

Usage:
    python -I src/fingerprints.py --molecules data/processed/molecules.parquet \
        --out data/processed/fp_r2.npz --radius 2
"""

import argparse

import numpy as np
import pandas as pd
import scipy.sparse as sp
from joblib import Parallel, delayed
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")

N_BITS = 2 ** 32


def _chunk(smiles: list[str], radius: int, use_counts: bool):
    """Return (row_offsets_per_molecule, column_indices) for one chunk."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius)
    counts, cols = [], []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            counts.append(0)
            continue
        bits = gen.GetSparseFingerprint(mol).GetOnBits()
        # Reinterpret the signed 32-bit hash as unsigned, losslessly.
        c = np.fromiter(bits, dtype=np.int64, count=len(bits)) & 0xFFFFFFFF
        c = np.unique(c)
        cols.append(c)
        counts.append(len(c))
    return np.asarray(counts, dtype=np.int64), (
        np.concatenate(cols) if cols else np.empty(0, dtype=np.int64)
    )


def smiles_to_csr(smiles: list[str], radius: int = 2, n_jobs: int = 32,
                  chunk_size: int = 5000, verbose: bool = True) -> sp.csr_matrix:
    chunks = [smiles[i:i + chunk_size] for i in range(0, len(smiles), chunk_size)]
    out = Parallel(n_jobs=n_jobs, verbose=5 if verbose else 0)(
        delayed(_chunk)(c, radius, False) for c in chunks
    )
    counts = np.concatenate([c for c, _ in out])
    indices = np.concatenate([i for _, i in out])
    indptr = np.zeros(len(counts) + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    data = np.ones(len(indices), dtype=np.bool_)
    return sp.csr_matrix((data, indices, indptr), shape=(len(smiles), N_BITS),
                         dtype=np.bool_)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--molecules", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--radius", type=int, default=2, help="2 => ECFP4")
    p.add_argument("--n-jobs", type=int, default=64)
    args = p.parse_args()

    mols = pd.read_parquet(args.molecules)
    X = smiles_to_csr(mols.smiles.tolist(), radius=args.radius, n_jobs=args.n_jobs)
    print(f"{X.shape[0]:,} molecules, {X.nnz:,} set bits, "
          f"{X.nnz / X.shape[0]:.1f} per molecule, "
          f"{len(np.unique(X.indices)):,} distinct features", flush=True)

    sp.save_npz(args.out, X)
    # The merged ChEMBL+ExCAPE build keys molecules by mol_id, the original by
    # mol_chembl_id; accept whichever is present.
    id_col = "mol_chembl_id" if "mol_chembl_id" in mols.columns else "mol_id"
    np.save(args.out.replace(".npz", "_ids.npy"), mols[id_col].to_numpy())
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
