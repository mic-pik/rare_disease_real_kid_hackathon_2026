"""Subsample the ClinVar set for embedding extraction: keep all P/LP, up to RATIO x benign per stratum.

Strata = split x consequence group (most severe ClinVar MC term), so benign variants are matched to
pathogenic ones by consequence rather than dominated by synonymous/intronic SNVs.
Truncating/canonical-splice P/LP are capped (PATHOGENIC_CAPS, split proportionally across splits):
they are nearly all pathogenic, add little beyond a consequence rule, and dominate GPU time.
All variants in FORCE_KEEP_GENES are kept. sample_weight = 1 / sampling fraction (restores ClinVar mix).

Usage: .venv/bin/python scripts/subsample_clinvar.py <clinvar_set.tsv.gz> <out_dir> [ratio]
Outputs: <out_dir>/clinvar_subset.tsv.gz, <out_dir>/clinvar_subset.vcf.gz (+ .tbi)
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd

IN, OUT = sys.argv[1], Path(sys.argv[2])
RATIO = float(sys.argv[3]) if len(sys.argv) > 3 else 2.0
SEED = 0
FORCE_KEEP_GENES = {"BUB1B", "CEP57", "TRIP13"}
PATHOGENIC_CAPS = {"frameshift": 20_000, "stop_start": 20_000, "splice_site": 20_000}  # total over splits
SECONDS_PER_VARIANT, N_GPUS = 0.72, 4  # Evo 2 7B, 8 kb ref+alt windows, measured in smoke test (unoptimised)

# Most severe first; first matching group wins
CONSEQUENCE_GROUPS = [
    ("frameshift", {"frameshift_variant"}),
    ("stop_start", {"nonsense", "stop_lost", "initiator_codon_variant"}),
    ("splice_site", {"splice_donor_variant", "splice_acceptor_variant"}),
    ("inframe_indel", {"inframe_deletion", "inframe_insertion", "inframe_indel"}),
    ("missense", {"missense_variant"}),
    ("synonymous", {"synonymous_variant"}),
    ("utr", {"5_prime_UTR_variant", "3_prime_UTR_variant"}),
    ("intronic", {"intron_variant"}),
]


def consequence_group(mc):
    terms = set() if pd.isna(mc) else {t.split("|")[1] for t in mc.split(",")}
    return next((g for g, members in CONSEQUENCE_GROUPS if terms & members), "other")


df = pd.read_csv(IN, sep="\t", dtype={"chrom": str})
df["csq_group"] = df.mc.map(consequence_group)
forced = df.genes.str.split("|").map(lambda gs: bool(FORCE_KEEP_GENES & set(gs)))

n_path_group = df[~forced & (df.label == 1)].groupby("csq_group").size()

parts = []
for (split, csq), g in df[~forced].groupby(["split", "csq_group"]):
    path, benign = g[g.label == 1], g[g.label == 0]
    if csq in PATHOGENIC_CAPS:
        n_path = min(len(path), round(PATHOGENIC_CAPS[csq] * len(path) / n_path_group[csq]))
        path = path.sample(n_path, random_state=SEED).assign(sample_weight=(g.label == 1).sum() / max(n_path, 1))
    else:
        path = path.assign(sample_weight=1.0)
    n_benign = min(len(benign), int(RATIO * len(path)))
    kept = benign.sample(n_benign, random_state=SEED)
    parts += [path, kept.assign(sample_weight=len(benign) / n_benign if n_benign else 1.0)]
sub = pd.concat(parts + [df[forced].assign(sample_weight=1.0)])
sub = df[["var_key"]].merge(sub, on="var_key")  # restore genome order

OUT.mkdir(parents=True, exist_ok=True)
sub.to_csv(OUT / "clinvar_subset.tsv.gz", sep="\t", index=False)

# Sites VCF for downstream tools: reuse the full-set header
full_vcf = Path(IN).with_suffix("").with_suffix(".vcf.gz")
keys = OUT / "clinvar_subset.keys"
sub.var_key.to_csv(keys, index=False, header=False)
subprocess.run(["bcftools", "view", "-i", f"ID=@{keys}", "-Oz", "-o", OUT / "clinvar_subset.vcf.gz", full_vcf], check=True)
subprocess.run(["tabix", "-f", "-p", "vcf", OUT / "clinvar_subset.vcf.gz"], check=True)
keys.unlink()

print(f"Kept {len(sub):,} / {len(df):,} variants (benign:pathogenic <= {RATIO}:1 per split x consequence)")
print("\nSubset by consequence group x label (full set in brackets):")
full = pd.crosstab(df.csq_group, df.label)
subt = pd.crosstab(sub.csq_group, sub.label)
print(subt.astype(str) + " [" + full.loc[subt.index].astype(str) + "]")
print("\nSubset by split x label:")
print(pd.crosstab(sub.split, sub.label, margins=True))
hours = len(sub) * SECONDS_PER_VARIANT / 3600
print(f"\nEstimated Evo 2 7B time: {hours:.0f} GPU-h ({hours / N_GPUS:.1f} h on {N_GPUS} GPUs)")
