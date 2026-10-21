#!/usr/bin/env python3
"""Portable sequential pipeline; --dry-run shows every command without executing it.

Run heavy stages inside an allocation. There are no embedded cluster accounts,
credentials or machine paths. Relative configuration paths resolve from the repo.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]
STAGES = ('download', 'references', 'training-data', 'training-features',
          'train', 'knowledge', 'chemical-models', 'patient', 'patient-features',
          'score', 'reason')


def resolve(value):
    p = Path(value).expanduser()
    return p.resolve() if p.is_absolute() else (ROOT / p).resolve()


def build_plan(config, stages):
    """Return command records; pure planning allows tests without data or GPUs."""
    t = ROOT / 'track1'
    g = ROOT / 'gene_mol'
    m = ROOT / 'mva_llm'
    py = {key: str(resolve(value)) for key, value in config['python'].items()}
    n = str(config.get('cpus', 8))
    plan = []
    def add(stage, cwd, *argv):
        if stage in stages:
            plan.append({'stage': stage, 'cwd': str(cwd), 'argv': [str(x) for x in argv]})
    def script(stage, cwd, env, name, *args):
        add(stage, cwd, py[env], cwd/'scripts'/name if cwd == t else cwd/'src'/name, *args)
    add('download', ROOT, py['reasoning'], ROOT/'pipeline/prepare_public.py', 'download')
    add('references', ROOT, sys.executable, ROOT/'pipeline/prepare_public.py', 'references')
    add('references', t, 'bash', t/'scripts/setup_vep.sh', t)
    script('training-data', t, 'evo', 'build_clinvar_set.py', 'ref/clinvar/clinvar.vcf.gz',
           'ref/grch38/GRCh38.primary_assembly.fa.fai', 'work/clinvar')
    add('training-data', t, 'bgzip', '-f', 'work/clinvar/clinvar_set.vcf')
    add('training-data', t, 'tabix', '-f', '-p', 'vcf', 'work/clinvar/clinvar_set.vcf.gz')
    add('training-data', t, 'bcftools', 'norm', '-f', 'ref/grch38/GRCh38.primary_assembly.fa',
        '-c', 'e', 'work/clinvar/clinvar_set.vcf.gz', '-Ob', '-o', os.devnull)
    script('training-data', t, 'evo', 'subsample_clinvar.py', 'work/clinvar/clinvar_set.tsv.gz', 'work/clinvar', '2.0')
    add('training-data', t, 'bash', t/'scripts/run_vep.sh', t, 'work/clinvar/clinvar_subset.vcf.gz', 'work/vep/clinvar_subset.vep.tsv', n)
    for stage, subset, vep, prefix in [
        ('training-features', 'work/clinvar/clinvar_subset.tsv.gz', 'work/vep/clinvar_subset.vep.tsv', ''),
        ('patient-features', 'work/patient/candidates.tsv.gz', 'work/patient/candidates.vep.tsv', 'patient_')]:
        script(stage, t, 'esm', 'esm1v_features.py', '--subset', subset, '--vep', vep,
               '--pep', 'ref/pep/Homo_sapiens.GRCh38.pep.all.fa.gz', '--models', 'models',
               '--out', 'work/features/'+prefix+'esm1v', '--device', 'cuda', '--batch', '16', '--shard', '0/1')
        script(stage, t, 'evo', 'evo2_features.py', '--subset', subset,
               '--fasta', 'ref/grch38/GRCh38.primary_assembly.fa', '--ckpt', 'models/evo2_7b/evo2_7b.pt',
               '--out', 'work/features/'+prefix+'evo2', '--window', '4096', '--batch', '4', '--shard', '0/1')
    script('train', t, 'evo', 'model_data.py', 'work/clinvar/clinvar_subset.tsv.gz', 'work/vep/clinvar_subset.vep.tsv',
           'work/features/esm1v', 'work/features/evo2', 'work/model/cache')
    script('train', t, 'evo', 'train_classic.py', 'work/model/cache', 'work/model/classic')
    script('train', t, 'evo', 'train_mlp.py', 'work/model/cache', 'work/model/mlp', 'all')
    script('train', t, 'evo', 'make_ensemble.py', 'work/model/cache', 'work/model')
    script('train', t, 'evo', 'evaluate_test.py', 'work/model/cache', 'work/model/test_results.tsv',
           'work/model/classic/preds_classic.parquet:LGBM[scalars]',
           'work/model/classic/preds_classic.parquet:LGBM[scalars+evo+esm]',
           'work/model/mlp/preds_mlp_all.parquet:MLP[scalars+all]',
           'work/model/preds_ensemble.parquet:ENS[LGBM+MLP]')
    add('knowledge', ROOT, py['reasoning'], ROOT/'pipeline/prepare_public.py', 'knowledge-inputs')
    script('knowledge', m, 'reasoning', 'build_kg.py', '--raw', 'data/raw', '--outdir', 'data/kg')
    script('knowledge', m, 'reasoning', 'build_clinvar_full.py', '--vcf', t/'ref/clinvar/clinvar.vcf.gz', '--out', 'data/kg/clinvar_full.parquet')
    add('knowledge', ROOT, py['reasoning'], ROOT/'pipeline/prepare_public.py', 'literature')
    add('knowledge', ROOT, py['reasoning'], ROOT/'pipeline/prepare_public.py', 'supplements')
    db = 'data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db'
    script('chemical-models', g, 'chem', 'extract_chembl.py', '--db', db, '--out', 'data/processed/activities.parquet')
    script('chemical-models', g, 'chem', 'build_dataset.py', '--activities', 'data/processed/activities.parquet', '--outdir', 'data/processed', '--n-jobs', n)
    script('chemical-models', g, 'chem', 'build_dataset_merged.py', '--activities', 'data/processed/activities.parquet',
           '--excape', 'data/excape_raw/excape_v2.tsv.xz', '--outdir', 'data/processed_merged3', '--n-jobs', n)
    for data, model in [('processed', 'lnb_r2'), ('processed_merged3', 'lnb_merged3')]:
        script('chemical-models', g, 'chem', 'fingerprints.py', '--molecules', f'data/{data}/molecules.parquet',
               '--out', f'data/{data}/fp_r2.npz', '--radius', '2', '--n-jobs', n)
        script('chemical-models', g, 'chem', 'train_targets.py', '--data', f'data/{data}', '--fp', f'data/{data}/fp_r2.npz',
               '--outdir', f'models/{model}', '--n-jobs', n)
    script('chemical-models', g, 'chem', 'build_serving.py', '--core', 'v1', '--merged', 'models/lnb_merged3', '--out', 'models/serving_v2')
    for name in ['calibrate_extended.py', 'calibrate_unified.py']:
        script('chemical-models', g, 'chem', name, '--n-jobs', n)
    script('chemical-models', g, 'chem', 'build_essentiality.py', '--out', 'work/essentiality.parquet')
    script('chemical-models', g, 'chem', 'build_codependency.py', '--out', m/'data/kg/codependency.parquet')
    script('chemical-models', g, 'chem', 'build_drug_safety.py', '--db', db, '--out', m/'data/kg/drug_safety.parquet')
    script('chemical-models', g, 'chem', 'build_sibling_diseases.py', '--out', m/'data/kg/sibling_diseases.parquet')
    if 'patient' in stages:
        add('patient', ROOT, sys.executable, ROOT/'pipeline/prepare_patient.py', '--vcf', resolve(config['patient_vcf']), '--cpus', n)
        script('patient', t, 'evo', 'prepare_patient.py', 'work/patient/patient.pass.ids.vcf.gz', 'work/vep/patient_pass.vep.tsv', 'work/patient')
    script('score', t, 'evo', 'score_patient.py', '--check')
    script('score', t, 'evo', 'score_patient.py', '--candidates', 'work/patient/candidates.tsv.gz',
           '--esm', 'work/features/patient_esm1v', '--evo', 'work/features/patient_evo2', '--out', 'work/patient/scored')
    drift = ['--allow-reference-drift'] if config.get('allow_reference_drift', False) else []
    add('reason', ROOT, py['reasoning'], ROOT/'pipeline/check_resources.py',
        '--out', resolve(config['output'])/'resource_check.json', *drift)
    add('reason', ROOT, py['reasoning'], m/'unified/workflow.py', '--config', resolve(config['output'])/'workflow_config.json',
        '--out', resolve(config['output'])/'analysis')
    # Preserve stage order even when commands are declared in shared loops.
    return sorted(plan, key=lambda item: STAGES.index(item['stage']))


def workflow_config(config):
    return {'case': str(resolve(config['case'])),
            'variants': str(ROOT/'track1/work/patient/scored/variants_scored.tsv.gz'),
            'hpo_annotations': str(ROOT/'mva_llm/data/raw/hpo/genes_to_phenotype.txt'),
            'hpo_ontology': str(ROOT/'mva_llm/data/raw/hpo/hp.obo'),
            'model': str(ROOT/'mva_llm/models/reasoning'), 'tp': config.get('tp', 4), 'seed': config.get('seed', 17),
            'fulltext_dir': str(ROOT/'mva_llm/data/unified_public_v1/fulltext'),
            'identity_file': str(ROOT/'mva_llm/data/unified_public_v1/identities.json')}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', required=True)
    ap.add_argument('--stage', choices=['all', *STAGES], required=True)
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    config = json.loads(Path(args.config).read_text())
    stages = STAGES if args.stage == 'all' else [args.stage]
    plan = build_plan(config, stages)
    if args.dry_run:
        for item in plan: print(f"[{item['stage']}] cwd={item['cwd']}\n  {shlex.join(item['argv'])}")
        return
    out = resolve(config['output'])
    out.mkdir(parents=True, exist_ok=True)
    if 'reason' in stages and (out/'analysis').exists() and any((out/'analysis').iterdir()):
        raise ValueError('Analysis output already exists. Choose a new output directory.')
    for rel in ['track1/work/vep', 'track1/work/patient', 'track1/work/clinvar',
                'track1/tools', 'track1/ref/vep_cache', 'gene_mol/data/processed', 'gene_mol/work', 'mva_llm/data/kg']:
        (ROOT/rel).mkdir(parents=True, exist_ok=True)
    if 'reason' in stages:
        (out/'workflow_config.json').write_text(json.dumps(workflow_config(config), indent=2)+'\n')
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(config.get('cpus', 8)), VLLM_WORKER_MULTIPROC_METHOD='spawn')
    for i, item in enumerate(plan):
        entry = {**item, 'started': datetime.now(timezone.utc).isoformat()}
        print(f"[{item['stage']}] {shlex.join(item['argv'])}", flush=True)
        proc = subprocess.run(item['argv'], cwd=item['cwd'], env=env)
        entry['returncode'] = proc.returncode
        with (out/'execution.jsonl').open('a') as f: f.write(json.dumps(entry)+'\n')
        if proc.returncode: raise SystemExit(proc.returncode)


if __name__ == '__main__':
    main()
