"""Build a complete ClinVar table, including VUS and conflicting records.

rare_real's clinvar_set.tsv.gz is a *training* set for a pathogenicity
classifier: it keeps only P/LP and B/LB with >=1 review star and drops
Uncertain_significance and Conflicting_classifications. That is right for
training a binary classifier and wrong for a lookup tool - asked about a VUS,
the tool answered "not in ClinVar", which is false and is exactly the kind of
confident absence a clinician must not be given.

This parses the raw ClinVar VCF so every submitted record is retrievable, with
its review status.

Usage:
    python -I src/build_clinvar_full.py --vcf <clinvar.vcf.gz> --out data/kg/clinvar_full.parquet
"""

import argparse
import gzip
import re

import pandas as pd

# ClinVar review status -> star rating, as ClinVar itself defines it.
STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_multiple_submitters": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_classification_provided": 0,
    "no_classifications_from_unflagged_records": 0,
    "no_assertion_provided": 0,
}

INFO = re.compile(r"([^=;]+)=([^;]*)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    rows = []
    with gzip.open(args.vcf, "rt") as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            info = dict(INFO.findall(f[7]))
            sig = info.get("CLNSIG")
            if not sig:
                continue
            rev = info.get("CLNREVSTAT", "")
            gene = info.get("GENEINFO", "")
            gene = "|".join(g.split(":")[0] for g in gene.split("|") if g)
            rows.append((f"{chrom}:{pos}:{ref}:{alt}", chrom, int(pos), ref, alt,
                         vid, sig, rev, STARS.get(rev, 0), gene,
                         info.get("MC", ""), info.get("CLNVC", ""),
                         info.get("CLNDN", "")))

    df = pd.DataFrame(rows, columns=[
        "var_key", "chrom", "pos", "ref", "alt", "clinvar_id", "clnsig",
        "review_status", "stars", "genes", "mc", "clnvc", "conditions"])
    df = df.drop_duplicates("var_key")
    df.to_parquet(args.out, index=False)

    print(f"{len(df):,} variants", flush=True)
    print(df.clnsig.value_counts().head(12).to_string(), flush=True)
    print("\nby review stars:", flush=True)
    print(df.stars.value_counts().sort_index().to_string(), flush=True)
    print(f"\nwrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
