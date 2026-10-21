"""Build the ClinVar training set: labels, filters, gene-grouped train/val/test split, sites VCF for VEP.

Usage: .venv/bin/python scripts/build_clinvar_set.py <clinvar.vcf.gz> <genome.fa.fai> <out_dir>
Outputs: <out_dir>/clinvar_set.tsv.gz, <out_dir>/clinvar_set.vcf (unsorted-safe, sorted by genome order)
"""
import subprocess
import sys
from io import StringIO
from pathlib import Path

import numpy as np
import pandas as pd

CLINVAR, FAI, OUT = sys.argv[1], sys.argv[2], Path(sys.argv[3])
OUT.mkdir(parents=True, exist_ok=True)

LABELS = {
    "Pathogenic": 1, "Likely_pathogenic": 1, "Pathogenic/Likely_pathogenic": 1,
    "Benign": 0, "Likely_benign": 0, "Benign/Likely_benign": 0,
}
STARS = {
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "reviewed_by_expert_panel": 3,
    "practice_guideline": 4,
}
GERMLINE_BITS = 1 | 4 | 8 | 16 | 32 | 64 | 128
CHROMS = [str(c) for c in range(1, 23)] + ["X", "Y"]
MAX_LEN = 51  # REF/ALT up to 51 bp, i.e. indels up to 50 bp
FORCE_TEST_GENES = {"BUB1B", "CEP57", "TRIP13"}  # MVA genes: keep out of train/val
SPLIT_FRAC = {"train": 0.8, "val": 0.1}  # remainder -> test
SEED = 0

fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/CLNSIG\t%INFO/CLNREVSTAT\t%INFO/GENEINFO\t%INFO/MC\t%INFO/CLNVC\t%INFO/ORIGIN\n"
raw = subprocess.run(["bcftools", "query", "-f", fmt, CLINVAR], check=True, capture_output=True, text=True).stdout
cols = ["chrom", "pos", "clinvar_id", "ref", "alt", "clnsig", "revstat", "geneinfo", "mc", "clnvc", "origin"]
df = pd.read_csv(StringIO(raw), sep="\t", names=cols, dtype=str, na_values=["."], keep_default_na=False)
df["pos"] = df.pos.astype(int)
n0 = len(df)


def report(step):
    print(f"{step:45s} {len(df):>9,}")


report("ClinVar records")
df = df[df.clnsig.isin(LABELS)]
report("P/LP or B/LB")
df = df[df.revstat.isin(STARS)]
report(">= 1 review star")
origin = df.origin.fillna("0").str.split("|").map(lambda v: int(v[0]) if v[0].isdigit() else 0)
df = df[(origin == 0) | (origin & GERMLINE_BITS > 0)]
report("germline or unknown origin")
df = df[df.chrom.isin(CHROMS)]
report("chr1-22, X, Y")
acgt = df.ref.str.fullmatch("[ACGT]+") & df.alt.str.fullmatch("[ACGT]+")
df = df[acgt & (df.ref.str.len() <= MAX_LEN) & (df.alt.str.len() <= MAX_LEN)]
report("ACGT alleles, indels <= 50 bp")
df = df[df.geneinfo.notna()]
report("has gene")

df["label"] = df.clnsig.map(LABELS)
df["stars"] = df.revstat.map(STARS)
df["genes"] = df.geneinfo.str.split("|").map(lambda gs: sorted({g.split(":")[0] for g in gs}))
df["var_type"] = np.select(
    [(df.ref.str.len() == 1) & (df.alt.str.len() == 1), df.ref.str.len() > df.alt.str.len(), df.ref.str.len() < df.alt.str.len()],
    ["SNV", "deletion", "insertion"], "delins")

# Union-find over genes that co-occur on a variant, so overlapping genes land in the same split
parent = {}


def find(g):
    parent.setdefault(g, g)
    while parent[g] != g:
        parent[g] = parent[parent[g]]
        g = parent[g]
    return g


for gs in df.genes:
    for g in gs[1:]:
        parent[find(g)] = find(gs[0])
df["gene_group"] = df.genes.map(lambda gs: find(gs[0]))

groups = df.gene_group.unique()
rng = np.random.default_rng(SEED)
u = dict(zip(groups, rng.random(len(groups))))
forced = {find(g) for g in FORCE_TEST_GENES if g in parent}
split = {g: "test" if g in forced else "train" if u[g] < SPLIT_FRAC["train"]
         else "val" if u[g] < SPLIT_FRAC["train"] + SPLIT_FRAC["val"] else "test" for g in groups}
df["split"] = df.gene_group.map(split)
df["genes"] = df.genes.str.join("|")

# Stable variant key used across VEP / feature tables
df["var_key"] = df.chrom + ":" + df.pos.astype(str) + ":" + df.ref + ":" + df.alt
assert df.var_key.is_unique, "duplicate var_key"

order = {c: i for i, c in enumerate(pd.read_csv(FAI, sep="\t", header=None, dtype={0: str})[0])}
df = df.sort_values(["chrom", "pos"], key=lambda s: s.map(order) if s.name == "chrom" else s).reset_index(drop=True)

keep = ["var_key", "chrom", "pos", "ref", "alt", "clinvar_id", "label", "clnsig", "stars", "var_type",
        "genes", "gene_group", "split", "mc", "clnvc"]
df[keep].to_csv(OUT / "clinvar_set.tsv.gz", sep="\t", index=False)

with open(OUT / "clinvar_set.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    for line in open(FAI):
        c, length = line.split("\t")[:2]
        if c in CHROMS:
            f.write(f"##contig=<ID={c},length={length}>\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    df[["chrom", "pos", "var_key", "ref", "alt"]].assign(q=".", f=".", i=".").to_csv(f, sep="\t", header=False, index=False)

print(f"\nKept {len(df):,} / {n0:,}; {df.gene_group.nunique():,} gene groups; "
      f"forced-test groups: {sorted(forced)}")
print("\nVariants by split x label:")
print(pd.crosstab(df.split, df.label, margins=True))
print("\nGene groups per split:")
print(df.groupby("split").gene_group.nunique())
print("\nVariants by type x label:")
print(pd.crosstab(df.var_type, df.label, margins=True))
print("\nMVA genes:")
mva = df[df.genes.str.split("|").map(lambda gs: bool(FORCE_TEST_GENES & set(gs)))]
print(pd.crosstab(mva.genes, [mva.split, mva.label]))
