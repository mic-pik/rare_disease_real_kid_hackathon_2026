#!/usr/bin/env bash
# Annotate a VCF with Ensembl VEP 115 (offline cache), one picked consequence per variant, tab output.
# Usage: scripts/run_vep.sh <project_dir> <in.vcf.gz> <out.tsv> [forks]
# --distance 0: only transcripts the variant overlaps, so --pick (MANE first) can't choose a neighbouring
# gene's MANE transcript (upstream/downstream) over the overlapped gene's transcript.
# ClinVar/CLIN_SIG columns from the cache are deliberately not output, to avoid label leakage into features.
set -euo pipefail
PROJ=$(realpath "$1"); IN=$(realpath "$2"); OUT=$(realpath -m "$3"); FORKS=${4:-64}

FIELDS="Uploaded_variation,Location,Allele,VARIANT_CLASS,Gene,SYMBOL,Feature,BIOTYPE,STRAND,Consequence,IMPACT"
FIELDS+=",MANE_SELECT,CANONICAL,ENSP,EXON,INTRON,DISTANCE,HGVSc,HGVSp,cDNA_position,CDS_position"
FIELDS+=",Protein_position,Amino_acids,Codons,gnomADe_AF,gnomADg_AF,MAX_AF"

BINDS=$(printf "%s\n" "$PROJ" "$(dirname "$IN")" "$(dirname "$OUT")" | sort -u | paste -sd,)

apptainer exec --bind "$BINDS" "$PROJ/tools/vep_115.0.sif" vep \
    --input_file "$IN" --format vcf \
    --output_file "$OUT" --tab --fields "$FIELDS" --force_overwrite \
    --stats_file "${OUT%.tsv}_summary.html" \
    --offline --cache --dir_cache "$PROJ/ref/vep_cache" --cache_version 115 --assembly GRCh38 \
    --fasta "$PROJ/ref/grch38/GRCh38.primary_assembly.fa" \
    --distance 0 \
    --pick --pick_order mane_select,canonical,appris,tsl,biotype,ccds,rank,length \
    --symbol --canonical --mane --biotype --protein --hgvs --numbers --variant_class \
    --check_existing --af_gnomade --af_gnomadg --max_af \
    --fork "$FORKS" --buffer_size 5000
