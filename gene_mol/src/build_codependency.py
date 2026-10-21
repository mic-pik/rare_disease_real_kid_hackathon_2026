"""Gene-gene co-dependency from DepMap CRISPR, for "what else does this cell need?".

The pipeline's existing dependency tool (`dependency_targets`) correlates a
*drug's* killing profile against gene dependency profiles. There is no
gene-vs-gene equivalent, so the only way to go from a broken gene to a
druggable vulnerability has been to hand-curate its pathway neighbours - which
for BUB1B returns the mitotic core (CDK1, PLK1, AURKB, the APC/C), every
member of which is pan-essential. Hand curation cannot help there, because the
thing it retrieves is the textbook pathway and the textbook pathway is what
every dividing cell runs.

Co-dependency is the data-driven replacement. Two genes whose Chronos
profiles rise and fall together across 1,178 cell lines are needed under the
same conditions, whether or not anyone has drawn them in the same pathway.

Method (fixed - these numbers are a regression test):
  * drop genes with >10% missing values across the panel
  * fill the remaining gaps with that gene's panel mean
  * Pearson correlation of every gene's profile against every other's
Top BUB1B partners come out CCNB1 0.250, MAD2L1 0.228, TP53 0.204,
FBXO5 0.200, CUL1 0.199, PSMA6 0.198, RACGAP1 0.194.

Storage choice: the full matrix is 17,205 x 17,205 = 296M correlations, ~1.2 GB
as float32 and far more as a long parquet, for a tool that is only ever asked
for one gene's top partners. We therefore compute the whole matrix in row
blocks but persist only the top `--top-n` positive partners per gene
(default 50, ~860k rows, a few MB). Negative co-dependency - genes whose
dependencies anti-correlate - is a different question (synthetic rescue /
paralog buffering) and is deliberately not stored; nothing downstream asks for
it, and keeping it would double the file for an answer no tool reads. The
gene universe is every gene that survives the NaN filter, i.e. it is not
restricted to a curated druggable set: restricting it would reintroduce
exactly the curation bias this tool exists to remove.

Usage:
    python -I src/build_codependency.py \
        --out ../mva_llm/data/kg/codependency.parquet
"""
import argparse
import os
import re

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--crispr", default="data/depmap_raw/crispr_gene_effect.csv")
    p.add_argument("--out", default="../mva_llm/data/kg/codependency.parquet")
    p.add_argument("--max-nan-frac", type=float, default=0.10,
                   help="drop genes missing in more than this fraction of lines")
    p.add_argument("--top-n", type=int, default=50,
                   help="positive partners stored per gene")
    p.add_argument("--block", type=int, default=2000,
                   help="genes per correlation block (memory knob)")
    args = p.parse_args()
    os.chdir(ROOT)

    d = pd.read_csv(args.crispr, index_col=0)
    print(f"{d.shape[0]:,} cell lines x {d.shape[1]:,} genes", flush=True)
    d.columns = [re.sub(r"\s*\(\d+\)$", "", c).strip() for c in d.columns]
    d = d.loc[:, ~d.columns.duplicated()]

    keep = d.columns[d.isna().mean() <= args.max_nan_frac]
    d = d[keep]
    n_lines, n_genes = d.shape
    print(f"after <={args.max_nan_frac:.0%} NaN filter: {n_genes:,} genes", flush=True)

    X = d.to_numpy(dtype=np.float32)
    col_mean = np.nanmean(X, axis=0)
    miss = np.where(np.isnan(X))
    X[miss] = np.take(col_mean, miss[1])

    # z-score per gene so a correlation is a dot product / n
    Z = (X - X.mean(axis=0)) / (X.std(axis=0) + 1e-12)
    genes = np.asarray(d.columns)
    del d, X

    k = args.top_n
    src, dst, cor, rnk = [], [], [], []
    for s in range(0, n_genes, args.block):
        e = min(s + args.block, n_genes)
        R = (Z[:, s:e].T @ Z) / n_lines                   # block x n_genes
        for bi in range(e - s):
            i = s + bi
            r = R[bi]
            r[i] = -np.inf                                # never self
            idx = np.argpartition(-r, k)[:k]
            idx = idx[np.argsort(-r[idx])]
            idx = idx[r[idx] > 0]                         # positive partners only
            src.append(np.full(len(idx), genes[i]))
            dst.append(genes[idx])
            cor.append(r[idx].astype(np.float32))
            rnk.append(np.arange(1, len(idx) + 1, dtype=np.int16))
        print(f"  {e:,}/{n_genes:,}", flush=True)

    out = pd.DataFrame({"gene": np.concatenate(src),
                        "partner": np.concatenate(dst),
                        "correlation": np.round(np.concatenate(cor), 4),
                        "rank": np.concatenate(rnk)})
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    out.to_parquet(args.out, index=False)

    b = out[out.gene == "BUB1B"].head(9)
    print("\nBUB1B regression check:")
    print(b[["partner", "correlation"]].to_string(index=False))
    print(f"\nwrote {args.out} ({len(out):,} rows, "
          f"{out.gene.nunique():,} genes, top {k} each)")


if __name__ == "__main__":
    main()
