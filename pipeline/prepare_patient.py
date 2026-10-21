#!/usr/bin/env python3
"""Normalise an authorised local GRCh38 VCF and run the frozen VEP policy."""
import argparse
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]/'track1'


def run(*argv):
    subprocess.run([str(a) for a in argv], check=True, cwd=ROOT)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--vcf', type=Path, required=True)
    ap.add_argument('--cpus', type=int, default=8)
    args = ap.parse_args()
    vcf = args.vcf.resolve()
    if not vcf.is_file(): raise FileNotFoundError(vcf)
    out = ROOT/'work/patient'; out.mkdir(parents=True, exist_ok=True)
    (ROOT/'work/vep').mkdir(parents=True, exist_ok=True)
    rename = out/'contigs.tsv'
    # Support both UCSC and Ensembl naming without changing coordinates.
    names = [str(i) for i in range(1,23)]+['X','Y','M','MT']
    rename.write_text(''.join(f'chr{x}\t{"MT" if x in ("M","MT") else x}\n' for x in names)+'M\tMT\n')
    run('bcftools','annotate','--rename-chrs',rename,vcf,'-Oz','-o',out/'renamed.vcf.gz')
    run('tabix','-f','-p','vcf',out/'renamed.vcf.gz')
    contigs = ','.join([str(i) for i in range(1,23)]+['X','Y','MT'])
    run('bcftools','view','-r',contigs,out/'renamed.vcf.gz','-Oz','-o',out/'primary.vcf.gz')
    run('bcftools','norm','-m','-any','-f',ROOT/'ref/grch38/GRCh38.primary_assembly.fa','-c','e',
        out/'primary.vcf.gz','-Oz','-o',out/'normalised_with_stars.vcf.gz')
    run('bcftools','view','-e','ALT="*"',out/'normalised_with_stars.vcf.gz','-Oz','-o',out/'patient.norm.vcf.gz')
    run('tabix','-f','-p','vcf',out/'patient.norm.vcf.gz')
    run('bcftools','view','-f','PASS','-t',','.join(contigs.split(',')[:-1]),out/'patient.norm.vcf.gz',
        '-Oz','-o',out/'pass.vcf.gz')
    run('bcftools','annotate','--set-id','%CHROM:%POS:%REF:%FIRST_ALT',out/'pass.vcf.gz',
        '-Oz','-o',out/'patient.pass.ids.vcf.gz')
    run('tabix','-f','-p','vcf',out/'patient.pass.ids.vcf.gz')
    run('bash',ROOT/'scripts/run_vep.sh',ROOT,out/'patient.pass.ids.vcf.gz',ROOT/'work/vep/patient_pass.vep.tsv',args.cpus)


if __name__ == '__main__': main()
