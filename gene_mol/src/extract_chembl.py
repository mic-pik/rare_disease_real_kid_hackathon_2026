"""Extract a clean bioactivity table from the ChEMBL SQLite dump.

Keeps only high-confidence binding assays against single proteins with a
pChEMBL value, then collapses replicate measurements to one value per
(target, molecule) pair.

Usage:
    python -I src/extract_chembl.py --db data/chembl37_db/chembl_37.db \
        --out data/processed/activities.parquet
"""

import argparse
import sqlite3

import pandas as pd

QUERY = """
SELECT
    td.chembl_id              AS target_chembl_id,
    td.pref_name              AS target_name,
    cseq.accession            AS uniprot,
    md.chembl_id              AS mol_chembl_id,
    cs.canonical_smiles       AS smiles,
    act.pchembl_value         AS pchembl,
    act.standard_type         AS standard_type,
    act.standard_relation     AS relation
FROM activities act
JOIN assays a               ON act.assay_id = a.assay_id
JOIN target_dictionary td   ON a.tid = td.tid
JOIN target_components tc   ON td.tid = tc.tid
JOIN component_sequences cseq ON tc.component_id = cseq.component_id
JOIN molecule_dictionary md ON act.molregno = md.molregno
JOIN compound_structures cs ON md.molregno = cs.molregno
WHERE td.target_type = 'SINGLE PROTEIN'
  AND td.organism = ?
  AND a.assay_type IN ('B', 'F')
  AND a.confidence_score >= 8
  AND act.pchembl_value IS NOT NULL
  AND act.potential_duplicate = 0
  AND act.data_validity_comment IS NULL
  AND act.standard_type IN ('IC50', 'EC50', 'Ki', 'Kd', 'AC50', 'Potency')
  AND act.standard_relation IN ('=', '<', '>', '<=', '>=')
"""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", required=True, help="path to chembl_XX.db")
    p.add_argument("--out", required=True, help="output parquet path")
    p.add_argument("--organism", default="Homo sapiens")
    args = p.parse_args()

    con = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    print(f"querying {args.db} ...", flush=True)
    df = pd.read_sql_query(QUERY, con, params=(args.organism,))
    con.close()
    print(f"raw rows: {len(df):,}  targets: {df.target_chembl_id.nunique():,}  "
          f"molecules: {df.mol_chembl_id.nunique():,}", flush=True)

    # A '>' or '>=' relation means the true potency is weaker than reported;
    # such rows are only trustworthy as evidence of inactivity.
    df["censored_weak"] = df.relation.isin([">", ">="])

    # Collapse replicates: median pChEMBL per (target, molecule).
    agg = (
        df.groupby(["target_chembl_id", "mol_chembl_id"], as_index=False)
        .agg(
            pchembl=("pchembl", "median"),
            n_meas=("pchembl", "size"),
            pchembl_std=("pchembl", "std"),
            any_censored_weak=("censored_weak", "any"),
        )
    )
    meta_t = df.drop_duplicates("target_chembl_id")[
        ["target_chembl_id", "target_name", "uniprot"]
    ]
    meta_m = df.drop_duplicates("mol_chembl_id")[["mol_chembl_id", "smiles"]]
    out = agg.merge(meta_t, on="target_chembl_id").merge(meta_m, on="mol_chembl_id")

    out.to_parquet(args.out, index=False)
    print(f"wrote {args.out}: {len(out):,} target-molecule pairs, "
          f"{out.target_chembl_id.nunique():,} targets, "
          f"{out.mol_chembl_id.nunique():,} molecules", flush=True)


if __name__ == "__main__":
    main()
