#!/usr/bin/env bash
# Fetch VEP 115 container + GRCh38 offline cache; verify Ensembl checksums before extracting.
# Usage: scripts/setup_vep.sh <project_dir>
set -euo pipefail
cd "$1"

ENS=https://ftp.ensembl.org/pub/release-115

# Ensembl CHECKSUMS lines are "<sum> <blocks> <file>" as produced by `sum`
verify() {  # verify <dir> <file>
    local expected actual
    expected=$(awk -v f="$2" '$3 == f {print $1, $2}' "$1/CHECKSUMS")
    actual=$(cd "$1" && sum "$2" | awk '{print $1, $2}')
    [[ -n "$expected" && "$expected" == "$actual" ]] || { echo "CHECKSUM MISMATCH $2: expected '$expected' got '$actual'"; exit 1; }
    echo "checksum OK: $2"
}

verify ref/pep Homo_sapiens.GRCh38.pep.all.fa.gz

[[ -f tools/vep_115.0.sif ]] || apptainer pull tools/vep_115.0.sif docker://ensemblorg/ensembl-vep:release_115.0
apptainer exec tools/vep_115.0.sif vep --help | grep -m1 -i 'ensembl-vep'

mkdir -p ref/vep_cache
cd ref/vep_cache
wget -q -O CHECKSUMS "$ENS/variation/indexed_vep_cache/CHECKSUMS"
# Single-connection FTP is ~0.8 MB/s; parallel ranged connections are much faster
aria2c -q -c -x 16 -s 16 -k 64M "$ENS/variation/indexed_vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz"
cd ../..
verify ref/vep_cache homo_sapiens_vep_115_GRCh38.tar.gz
tar xzf ref/vep_cache/homo_sapiens_vep_115_GRCh38.tar.gz -C ref/vep_cache

ls ref/vep_cache/homo_sapiens
echo DONE
