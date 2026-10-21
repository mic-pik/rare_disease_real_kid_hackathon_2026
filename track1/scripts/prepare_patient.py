"""Select patient candidate variants from the VEP-annotated PASS VCF and attach genotypes.

Candidates: PASS, in a transcript (not intergenic), ACGT alleles, indels <= 50 bp, alt allele in GT, and
  (gnomAD max AF < MAX_AF or absent)  OR  in an MVA gene (any frequency: hypomorphic alleles can be
  more common). gnomAD AF is used only here, for filtering — never as a model feature.

Writes:
  <out_dir>/candidates.tsv.gz      var_key, chrom, pos, ref, alt, split ("patient"), label (NaN), var_type,
                                    genes (VEP SYMBOL), consequence, impact, HGVS, genotype, depth, max_af
  <out_dir>/candidates.vep.tsv     the VEP rows of the candidates (input for esm1v_features.py)

Usage: .venv-evo2/bin/python scripts/prepare_patient.py <patient.pass.ids.vcf.gz> <patient_pass.vep.tsv> <out_dir>
"""
import io
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

MAX_AF = 0.01
MAX_LEN = 51
MVA_GENES = {"BUB1B", "CEP57", "TRIP13"}

vcf, vep_path, out = sys.argv[1], sys.argv[2], Path(sys.argv[3])
out.mkdir(parents=True, exist_ok=True)

vep_txt = "".join(line for line in Path(vep_path).read_text().splitlines(keepends=True) if not line.startswith("##")).lstrip("#")
vep = pd.read_csv(io.StringIO(vep_txt), sep="\t", dtype=str, na_values=["-"])
vep = vep.rename(columns={"Uploaded_variation": "var_key"})
assert vep.var_key.is_unique, "duplicate var_key in VEP output"
print(f"VEP rows (PASS): {len(vep):,}")

gt_txt = subprocess.run(["bcftools", "query", "-f", "%ID\t%REF\t%ALT\t[%GT\t%AD\t%DP\t%GQ]\n", vcf],
                        check=True, capture_output=True, text=True).stdout
gt = pd.read_csv(io.StringIO(gt_txt), sep="\t", names=["var_key", "ref", "alt", "gt", "ad", "dp", "gq"],
                 dtype=str, na_values=["."])
assert gt.var_key.is_unique

df = vep.merge(gt, on="var_key", how="inner", validate="1:1")
df["chrom"] = df.var_key.str.split(":").str[0]
df["pos"] = df.var_key.str.split(":").str[1].astype(int)
df["max_af"] = pd.to_numeric(df.MAX_AF, errors="coerce")
df["first_csq"] = df.Consequence.str.split(",").str[0]
df["in_mva_gene"] = df.SYMBOL.isin(MVA_GENES)


def report(step, d):
    print(f"{step:45s} {len(d):>9,}   (MVA-gene variants {int(d.in_mva_gene.sum())})")


report("PASS, VEP-annotated", df)
df = df[df.first_csq != "intergenic_variant"]
report("in a transcript", df)
df = df[df.ref.str.fullmatch("[ACGT]+") & df.alt.str.fullmatch("[ACGT]+")
        & (df.ref.str.len() <= MAX_LEN) & (df.alt.str.len() <= MAX_LEN)]
report("ACGT, indels <= 50 bp", df)
df = df[df["gt"].fillna("").str.contains("1")]  # df.gt is DataFrame.gt (>), use []
report("alt allele in genotype", df)
df = df[(df.max_af.isna() | (df.max_af < MAX_AF)) | df.in_mva_gene]
report(f"gnomAD max AF < {MAX_AF} or absent, or MVA gene", df)

df["var_type"] = np.select(
    [(df.ref.str.len() == 1) & (df.alt.str.len() == 1), df.ref.str.len() > df.alt.str.len(),
     df.ref.str.len() < df.alt.str.len()], ["SNV", "deletion", "insertion"], "delins")
df["zygosity"] = np.where(df["gt"].str.replace("|", "/").isin(["1/1"]), "hom", "het")
cand = pd.DataFrame({
    "var_key": df.var_key, "chrom": df.chrom, "pos": df.pos, "ref": df.ref, "alt": df.alt,
    "split": "patient", "label": np.nan, "var_type": df.var_type, "genes": df.SYMBOL.fillna(df.Gene),
    "gene_id": df.Gene, "consequence": df.Consequence, "impact": df.IMPACT, "hgvsc": df.HGVSc,
    "hgvsp": df.HGVSp, "max_af": df.max_af, "gt": df["gt"], "zygosity": df["zygosity"], "ad": df.ad,
    "dp": df.dp, "gq": df.gq, "in_mva_gene": df.in_mva_gene,
})
cand.to_csv(out / "candidates.tsv.gz", sep="\t", index=False)

# VEP rows of candidates, with the original header, for esm1v_features.py
with open(vep_path) as f_in, open(out / "candidates.vep.tsv", "w") as f_out:
    keep = set(cand.var_key)
    for line in f_in:
        if line.startswith("#") or line.split("\t", 1)[0] in keep:
            f_out.write(line)

print(f"\nCandidates: {len(cand):,}; missense {int((df.first_csq == 'missense_variant').sum()):,}; "
      f"hom {int((cand.zygosity == 'hom').sum()):,}; genes {cand.genes.nunique():,}")
print(cand.impact.value_counts().to_dict())
