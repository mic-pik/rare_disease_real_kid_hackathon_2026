"""Drug-level safety and development history, from the local ChEMBL dump.

Answers the questions a repurposing nomination has to survive: has this been
withdrawn anywhere, does it carry a black-box warning, what has it actually
been trialled for and how far did it get, and what pharmacological class is it.

ChEMBL 37 tables used: drug_warning (withdrawal / black box, with class and
year), drug_indication (indication + max phase reached for it),
molecule_atc_classification + atc_classification, molecule_dictionary.

Usage:
    python -I src/build_drug_safety.py --out ../mva_llm/data/kg/drug_safety.parquet
"""
import argparse
import os
import sqlite3

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB = ("data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", default=DB)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    os.chdir(ROOT)
    c = sqlite3.connect(args.db)

    warn = pd.read_sql("""
        SELECT md.chembl_id, md.pref_name, w.warning_type, w.warning_class,
               w.warning_description, w.warning_country, w.warning_year
        FROM drug_warning w JOIN molecule_dictionary md USING (molregno)""", c)
    print(f"warnings: {len(warn):,} rows, {warn.chembl_id.nunique():,} drugs", flush=True)

    ind = pd.read_sql("""
        SELECT md.chembl_id, i.mesh_heading, i.efo_term, i.max_phase_for_ind
        FROM drug_indication i JOIN molecule_dictionary md USING (molregno)""", c)
    ind = ind[ind.max_phase_for_ind.notna()]
    print(f"indications: {len(ind):,} rows, {ind.chembl_id.nunique():,} drugs", flush=True)

    atc = pd.read_sql("""
        SELECT md.chembl_id, a.level5, a.who_name, a.level1_description,
               a.level2_description, a.level4_description
        FROM molecule_atc_classification m
        JOIN atc_classification a USING (level5)
        JOIN molecule_dictionary md USING (molregno)""", c)

    mol = pd.read_sql("""
        SELECT chembl_id, pref_name, max_phase, withdrawn_flag, oral, parenteral,
               topical, black_box_warning, first_approval, prodrug, therapeutic_flag
        FROM molecule_dictionary""", c)
    print(f"molecules: {len(mol):,}; withdrawn {int(mol.withdrawn_flag.fillna(0).sum()):,}; "
          f"black box {int(mol.black_box_warning.fillna(0).sum()):,}", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    base = os.path.splitext(args.out)[0]
    mol.to_parquet(f"{base}.parquet", index=False)
    warn.to_parquet(f"{base}_warnings.parquet", index=False)
    ind.to_parquet(f"{base}_indications.parquet", index=False)
    atc.to_parquet(f"{base}_atc.parquet", index=False)
    print(f"wrote {base}{{,_warnings,_indications,_atc}}.parquet")


if __name__ == "__main__":
    main()
