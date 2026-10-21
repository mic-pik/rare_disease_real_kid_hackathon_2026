#!/usr/bin/env python3
"""Fail before inference on missing resources or changed historical evidence."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def semantic_digest(path):
    path = Path(path)
    if path.suffix == '.parquet':
        import pandas as pd
        # Preserve row and column order: ordering can affect retrieval tie breaks.
        text = pd.read_parquet(path).to_json(orient='split', index=False, double_precision=15)
    elif path.suffix == '.json':
        data = json.loads(path.read_text())
        def clean(value):
            if isinstance(value, list): return [clean(x) for x in value]
            if isinstance(value, dict):
                return {k: clean(v) for k,v in value.items() if k not in ('original_sha256', 'provenance')}
            return value
        text = json.dumps(clean(data), sort_keys=True, ensure_ascii=False, separators=(',',':'))
    else:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    return hashlib.sha256(text.encode()).hexdigest()


def check(root=ROOT, allow_drift=False):
    expected = json.loads((root/'config/expected_public_semantic.json').read_text())
    missing, changed = [], []
    for rel, digest in expected.items():
        path = root/rel
        if not path.is_file(): missing.append(rel)
        elif semantic_digest(path) != digest: changed.append(rel)
    for rel in ['gene_mol/work/essentiality.parquet','gene_mol/models/serving_v2/manifest.csv',
                'gene_mol/models/serving_v2/unified.json','mva_llm/models/reasoning/config.json']:
        if not (root/rel).is_file(): missing.append(rel)
    if (root/'gene_mol/models/serving_v2/manifest.csv').is_file():
        import csv
        with (root/'gene_mol/models/serving_v2/manifest.csv').open() as f:
            for row in csv.DictReader(f):
                if row['list'] in ('core','extended') and not (root/'gene_mol'/row['model_file']).is_file():
                    missing.append('gene_mol/'+row['model_file'])
    if not list((root/'mva_llm/models/reasoning').glob('*.safetensors')):
        missing.append('mva_llm/models/reasoning/*.safetensors')
    result = {'missing': missing, 'changed_public_resources': changed,
              'historical_reference_match': not missing and not changed,
              'reference_drift_explicitly_allowed': allow_drift}
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--allow-reference-drift', action='store_true')
    ap.add_argument('--out', type=Path)
    args = ap.parse_args()
    result = check(allow_drift=args.allow_reference_drift)
    text = json.dumps(result, indent=2)+'\n'
    print(text)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True); args.out.write_text(text)
    if result['missing'] or (result['changed_public_resources'] and not args.allow_reference_drift):
        raise SystemExit('Resources incomplete or changed. Review the report; do not claim an exact historical replay.')


if __name__ == '__main__': main()
