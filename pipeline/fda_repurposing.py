#!/usr/bin/env python3
"""Public-data FDA candidate retrieval, separate from patient inference.

Requires pandas/pyarrow (the chemical environment). No model calls, network
requests, inferred approval, fuzzy identity joins, or efficacy nominations.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import sqlite3
import subprocess
from datetime import datetime, timezone
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FDA_URL = 'https://www.fda.gov/media/89850/download?attachment='
SCHEMA_VERSION = '1.0'


def norm(value):
    return ' '.join(str(value or '').upper().split())


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def read_fda(archive):
    tables = {}
    with zipfile.ZipFile(archive) as z:
        for name in ('Products', 'Applications', 'Submissions', 'MarketingStatus',
                     'MarketingStatus_Lookup'):
            text = z.read(name + '.txt').decode('cp1252')
            tables[name] = list(csv.DictReader(io.StringIO(text), delimiter='\t'))
    return tables


def fda_products(tables):
    """Retain all products; eligibility requires original approval AND status.

    Application-level ORIG/AP evidence is not sufficient alone: a product may
    be tentative even when another product in its application was approved.
    """
    apps = {r['ApplNo']: r for r in tables['Applications']}
    approved = {}
    for r in tables['Submissions']:
        if r['SubmissionType'].strip() == 'ORIG' and r['SubmissionStatus'].strip() == 'AP':
            approved.setdefault(r['ApplNo'], []).append(r['SubmissionStatusDate'])
    labels = {r['MarketingStatusID']: r['MarketingStatusDescription']
              for r in tables['MarketingStatus_Lookup']}
    status = {}
    for r in tables['MarketingStatus']:
        status.setdefault((r['ApplNo'], r['ProductNo']), set()).add(r['MarketingStatusID'])
    out = []
    for r in tables['Products']:
        app, product = r['ApplNo'], r['ProductNo']
        ids = status.get((app, product), set())
        # Conflicting or unknown marketing statuses require review.
        recognised = bool(ids) and ids <= {'1', '2', '3'}
        confirmed = app in approved and recognised and apps.get(app, {}).get('ApplType', '').strip() in ('NDA', 'ANDA', 'BLA')
        parts = [norm(x) for x in r['ActiveIngredient'].split(';') if norm(x)]
        out.append({
            'product_id': app + ':' + product, 'application_number': app,
            'product_number': product, 'application_type': apps.get(app, {}).get('ApplType'),
            'brand_name': r['DrugName'], 'active_ingredient_text': r['ActiveIngredient'],
            'ingredient_names': parts, 'is_combination': len(parts) > 1,
            'form_route': r['Form'], 'strength': r['Strength'],
            'marketing_status': sorted(labels.get(x, 'Unknown: ' + x) for x in ids),
            'original_approval_dates': sorted(set(approved.get(app, []))),
            'approval_evidence_confirmed': confirmed,
            'eligible_current_product': confirmed and bool(ids & {'1', '2'}) and '3' not in ids,
            'safety_withdrawal_status': 'not_assessed',
            'fda_url': 'https://www.accessdata.fda.gov/scripts/cder/daf/index.cfm?event=overview.process&ApplNo=' + app,
        })
    return out


def resolve_names(con, names):
    """Exact preferred-name/synonym mapping; collapse only ChEMBL parent links.

    No arbitrary largest-fragment or salt-name stripping. Distinct parents
    matching the same name are quarantined even if one is more common.
    """
    wanted = set(names)
    matches = {n: set() for n in wanted}
    for mol, name in con.execute('SELECT molregno, pref_name FROM molecule_dictionary WHERE pref_name IS NOT NULL'):
        key = norm(name)
        if key in wanted:
            matches[key].add(mol)
    for mol, name in con.execute('SELECT molregno, synonyms FROM molecule_synonyms'):
        key = norm(name)
        if key in wanted:
            matches[key].add(mol)
    hierarchy = dict(con.execute('SELECT molregno, parent_molregno FROM molecule_hierarchy WHERE parent_molregno IS NOT NULL'))
    resolved, issues = {}, []
    for name in sorted(wanted):
        parents = {hierarchy.get(m, m) for m in matches[name]}
        if len(parents) != 1:
            issues.append({'ingredient_name': name, 'issue': 'unmatched' if not parents else 'ambiguous_parent',
                           'matching_molregnos': sorted(matches[name]), 'parent_molregnos': sorted(parents)})
            continue
        parent = next(iter(parents))
        row = con.execute('''SELECT m.chembl_id, m.pref_name, s.canonical_smiles, s.standard_inchi_key
            FROM molecule_dictionary m LEFT JOIN compound_structures s USING (molregno)
            WHERE m.molregno=?''', (parent,)).fetchone()
        if row is None or not row[2] or not row[3]:
            issues.append({'ingredient_name': name, 'issue': 'missing_parent_structure',
                           'matching_molregnos': sorted(matches[name]), 'parent_molregnos': [parent]})
            continue
        resolved[name] = {'ingredient_id': row[0], 'canonical_name': row[1],
                          'smiles': row[2], 'inchikey': row[3], 'parent_molregno': parent}
    return resolved, issues


def build_catalog(products, resolved):
    grouped = {}
    for p in products:
        # Keep unresolved components on the product; never convert a partially
        # resolved combination into a single-ingredient product.
        p['mapped_ingredient_ids'] = [resolved[n]['ingredient_id'] if n in resolved else None
                                      for n in p['ingredient_names']]
        for name in p['ingredient_names']:
            if name not in resolved:
                continue
            identity = resolved[name]
            entry = grouped.setdefault(identity['ingredient_id'], {
                **identity, 'fda_names': set(), 'product_ids': set(),
                'approved_product_ids': set(), 'eligible_product_ids': set(),
                'eligible_single_ingredient_product_ids': set(),
                'mapping_method': 'exact_name_or_synonym_then_chembl_parent',
                'safety_withdrawal_status': 'not_assessed'})
            entry['fda_names'].add(name)
            entry['product_ids'].add(p['product_id'])
            if p['approval_evidence_confirmed']:
                entry['approved_product_ids'].add(p['product_id'])
            if p['eligible_current_product']:
                entry['eligible_product_ids'].add(p['product_id'])
                if not p['is_combination']:
                    entry['eligible_single_ingredient_product_ids'].add(p['product_id'])
    return [{k: sorted(v) if isinstance(v, set) else v for k, v in row.items()}
            for _, row in sorted(grouped.items())]


def measured_interactions(con, ingredients, genes):
    """Retrieve original assay rows; do not average away censoring/conflicts."""
    con.execute('CREATE TEMP TABLE selected_parents (parent INTEGER PRIMARY KEY, ingredient_id TEXT)')
    con.executemany('INSERT INTO selected_parents VALUES (?,?)',
                    [(r['parent_molregno'], r['ingredient_id']) for r in ingredients if r['eligible_product_ids']])
    con.execute('CREATE TEMP TABLE selected_genes (gene TEXT PRIMARY KEY)')
    con.executemany('INSERT INTO selected_genes VALUES (?)', [(g,) for g in sorted(set(genes))])
    # Materialise the small target set first, avoiding a scan of all activities.
    con.execute('''CREATE TEMP TABLE selected_targets AS
        SELECT DISTINCT t.tid, t.chembl_id target_chembl_id, t.pref_name target_name,
          t.organism target_organism, s.component_synonym gene, c.accession uniprot
        FROM target_dictionary t JOIN target_components tc USING (tid)
        JOIN component_sequences c USING (component_id)
        JOIN component_synonyms s USING (component_id)
        JOIN selected_genes g ON s.component_synonym=g.gene
        WHERE s.syn_type='GENE_SYMBOL' AND t.organism='Homo sapiens'
          AND t.target_type='SINGLE PROTEIN' ''')
    coverage = pd.read_sql_query('SELECT * FROM selected_targets', con)
    query = '''SELECT DISTINCT sp.ingredient_id, m.chembl_id tested_molecule_chembl_id,
       t.gene, t.target_chembl_id, t.target_name, t.target_organism, t.uniprot,
       a.activity_id, ass.chembl_id assay_chembl_id, ass.assay_type,
       ass.description assay_description, ass.assay_organism, ass.assay_cell_type,
       ass.assay_tissue, ass.confidence_score,
       a.standard_type, a.standard_relation, a.standard_value, a.standard_units,
       a.pchembl_value, a.activity_comment, a.data_validity_comment,
       a.potential_duplicate, a.action_type reported_action_type,
       d.chembl_id document_chembl_id, d.pubmed_id, d.doi, d.year,
       d.title document_title
       FROM selected_targets t JOIN assays ass ON ass.tid=t.tid
       JOIN activities a ON a.assay_id=ass.assay_id
       JOIN molecule_dictionary m ON m.molregno=a.molregno
       LEFT JOIN molecule_hierarchy h ON h.molregno=m.molregno
       JOIN selected_parents sp ON sp.parent=COALESCE(h.parent_molregno,m.molregno)
       LEFT JOIN docs d ON d.doc_id=a.doc_id
       WHERE ass.assay_type IN ('B','F') AND ass.confidence_score >= 8'''
    frame = pd.read_sql_query(query, con)
    frame['source'] = 'ChEMBL original activity record'
    frame['direction_status'] = frame.reported_action_type.apply(
        lambda x: 'reported_in_activity_record' if pd.notna(x) and str(x).strip()
        and str(x).upper() not in ('UNKNOWN', 'OTHER') else 'unknown')
    frame['interpretation'] = 'Measured assay result; not evidence of disease rescue or achievable exposure'
    return frame, coverage


def sibling_candidates(con, ingredients, source, diseases):
    """Retrieve curated indication metadata, not invented trial outcomes."""
    frame = pd.read_parquet(source)
    frame = frame[frame.disease.str.upper().isin({norm(d) for d in diseases})].copy()
    mapping = dict(con.execute('''SELECT m.chembl_id, p.chembl_id FROM molecule_dictionary m
        LEFT JOIN molecule_hierarchy h ON m.molregno=h.molregno
        JOIN molecule_dictionary p ON p.molregno=COALESCE(h.parent_molregno,m.molregno)'''))
    allowed = {r['ingredient_id'] for r in ingredients if r['eligible_product_ids']}
    frame['ingredient_id'] = frame.chembl_id.map(mapping)
    frame = frame[frame.ingredient_id.isin(allowed)].copy()
    frame['source'] = 'Local curated sibling-disease table / ChEMBL indications'
    frame['evidence_level'] = 'indication_metadata_only'
    frame['measured_endpoint'] = 'not_retrieved'
    frame['transfer_status'] = 'unestablished; different lesion and population'
    frame['exposure_transfer'] = 'Exposure history requires primary-source confirmation; tolerability does not transfer automatically'
    return frame


def candidate_record(identity, routes, evidence):
    """Contract also represents historical candidates without FDA confirmation."""
    return {'schema_version': SCHEMA_VERSION, 'identity': identity,
            'retrieval_routes': sorted(routes), 'evidence': evidence,
            'mechanism': {'hypothesis': None, 'required_direction': 'unknown'},
            'assessments': {key: 'unestablished' for key in
                            ('mechanistic_fit', 'functional_rescue', 'exposure_feasibility', 'safety')},
            'opposing_evidence': [], 'unresolved_assumptions': [],
            'decision': 'insufficient_evidence',
            'decision_reason': 'Retrieval only; biological evidence review has not been performed',
            'clinical_suitability': 'unestablished'}


def export_candidates(ingredients, measured, siblings):
    rows = []
    for identity in ingredients:
        cid = identity['ingredient_id']
        m = measured[measured.ingredient_id == cid]
        s = siblings[siblings.ingredient_id == cid]
        if m.empty and s.empty:
            continue
        evidence, routes = [], []
        if not m.empty:
            routes.append('measured_pharmacology')
            evidence.append({'artifact': 'measured_interactions.parquet',
                             'ingredient_id': cid, 'n_records': len(m),
                             'targets': sorted(set(m.gene)), 'role': 'context'})
        if not s.empty:
            routes.append('related_disease')
            evidence.append({'artifact': 'related_disease_candidates.parquet',
                             'ingredient_id': cid, 'n_records': len(s),
                             'diseases': sorted(set(s.disease)),
                             'intents': sorted(set(s.intent.dropna())), 'role': 'context'})
        r = candidate_record(identity, routes, evidence)
        r['unresolved_assumptions'] = ['Direction and disease relevance need review',
             'FDA ingredient eligibility does not establish suitability for this disease',
             'Safety-related withdrawals and relevant exposure not yet assessed']
        rows.append(r)
    return rows


def baseline(out, source):
    # Recorded run artifacts carry case-derived content and are retained privately, so the
    # historical comparison is optional: record its absence rather than an empty candidate set.
    if not source.is_dir():
        write_json(out/'baseline_candidates.json', [])
        write_json(out/'baseline.json', {
            'historical_baseline': 'unavailable',
            'reason': 'Recorded run artifacts are not distributed with this repository',
            'expected_source': str(source),
            'git_revision': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_status': subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True)})
        return
    paths = sorted(source.glob('04_candidate_*.json'))
    result = []
    for p in paths:
        original = json.loads(p.read_text())
        c = original['candidate']
        r = candidate_record({'canonical_name': c['name'], 'original_identity': c.get('identity'),
                              'fda_status': 'not_assessed'}, ['historical_discovery'], c.get('evidence', []))
        r['original_result'] = original
        r['decision_reason'] = 'Historical result preserved; no FDA eligibility inferred'
        result.append(r)
    write_json(out/'baseline_candidates.json', result)
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).split(b'\0')
    hashes = {p.decode(): sha256(ROOT/p.decode()) for p in tracked if p and (ROOT/p.decode()).is_file()}
    write_json(out/'baseline.json', {
        'historical_baseline': 'available',
        'git_revision': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'git_status': subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True),
        'tracked_file_hashes': hashes, 'historical_candidate_hashes': {p.name:sha256(p) for p in paths},
        'original_verdicts': [x['original_result']['gate'] for x in result]})


def run(config, out):
    out = Path(out)
    if out.exists():
        raise ValueError('Output directory must be fresh; previous results are never overwritten')
    required = ('fda_zip', 'chembl_db', 'sibling_table')
    for key in required:
        if not Path(config[key]).is_file():
            raise ValueError('Missing input: ' + key)
    if not config.get('genes') or not config.get('sibling_diseases'):
        raise ValueError('Explicit genes and sibling_diseases are required; no inferred case inputs')
    out.mkdir(parents=True)
    write_json(out/'status.json', {'status': 'running'})
    try:
        baseline(out, ROOT/'submission/audit/final_automatic_run')
        inputs = {}
        for key in required:
            print('Hashing input: ' + key, flush=True)
            inputs[key] = {'path': str(Path(config[key]).resolve()), 'sha256': sha256(config[key])}
        manifest = {'created_utc': datetime.now(timezone.utc).isoformat(), 'schema_version': SCHEMA_VERSION,
                    'inputs': inputs, 'config': config, 'fda_source_url': FDA_URL,
                    'implementation_sha256': sha256(__file__), 'pandas_version': pd.__version__,
                    'limitations': ['FDA safety withdrawal review not implemented',
                       'ChEMBL-only measured retrieval; ExCAPE assay provenance not yet integrated',
                       'Related-disease rows are metadata, not verified clinical endpoints',
                       'No therapeutic ranking or efficacy nomination']}
        write_json(out/'manifest.json', manifest)
        print('Reading FDA products and resolving exact ingredient identities', flush=True)
        products = fda_products(read_fda(config['fda_zip']))
        con = sqlite3.connect(Path(config['chembl_db']).resolve().as_uri() + '?mode=ro', uri=True)
        try:
            names = {n for p in products for n in p['ingredient_names']}
            resolved, issues = resolve_names(con, names)
            ingredients = build_catalog(products, resolved)
            pd.DataFrame(products).to_parquet(out/'fda_products.parquet', index=False)
            pd.DataFrame(ingredients).to_parquet(out/'fda_ingredients.parquet', index=False)
            pd.DataFrame(issues, columns=['ingredient_name','issue','matching_molregnos','parent_molregnos']).to_csv(out/'identity_mapping_issues.csv', index=False)
            print('Retrieving original ChEMBL measurements', flush=True)
            measured, coverage = measured_interactions(con, ingredients, config['genes'])
            measured.to_parquet(out/'measured_interactions.parquet', index=False)
            coverage.to_csv(out/'target_coverage.csv', index=False)
            siblings = sibling_candidates(con, ingredients, config['sibling_table'], config['sibling_diseases'])
            siblings.to_parquet(out/'related_disease_candidates.parquet', index=False)
        finally:
            con.close()
        candidates = export_candidates(ingredients, measured, siblings)
        write_json(out/'candidates.json', candidates)
        pd.DataFrame([{'ingredient_id': r['identity']['ingredient_id'],
                       'name': r['identity']['canonical_name'],
                       'retrieval_routes': ';'.join(r['retrieval_routes']),
                       'has_current_single_ingredient_product': bool(r['identity']['eligible_single_ingredient_product_ids']),
                       'decision': r['decision'], 'decision_reason': r['decision_reason']}
                      for r in candidates], columns=['ingredient_id','name','retrieval_routes',
                       'has_current_single_ingredient_product','decision','decision_reason']).to_csv(out/'candidates.csv', index=False)
        summary = {'products': len(products), 'distinct_fda_names': len(names),
                   'mapped_names': len(resolved), 'quarantined_names': len(issues),
                   'canonical_ingredients': len(ingredients),
                   'eligible_ingredients': sum(bool(r['eligible_product_ids']) for r in ingredients),
                   'measured_records': len(measured), 'related_disease_records': len(siblings),
                   'candidate_ingredients': len(candidates),
                   'candidate_ingredients_with_current_single_product': sum(
                       bool(r['identity']['eligible_single_ingredient_product_ids']) for r in candidates),
                   'measured_records_without_action_direction': int((measured.direction_status == 'unknown').sum()),
                   'genes_without_target_mapping': sorted(set(config['genes']) - set(coverage.gene)),
                   'genes_without_measurements': sorted(set(config['genes']) - set(measured.gene)),
                   'requested_diseases_without_eligible_records': sorted(set(config['sibling_diseases']) - set(siblings.disease))}
        write_json(out/'summary.json', summary)
        manifest['output_hashes'] = {p.name:sha256(p) for p in sorted(out.iterdir())
                                      if p.is_file() and p.name not in ('manifest.json','status.json')}
        write_json(out/'manifest.json', manifest)
        write_json(out/'status.json', {'status':'completed_retrieval_only'})
        print(json.dumps(summary, indent=2))
    except Exception as exc:
        write_json(out/'status.json', {'status':'failed', 'error': str(exc)})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    config_path = Path(args.config).resolve()
    config = json.loads(config_path.read_text())
    # All resource paths are relative to the config, not the caller's cwd.
    for key in ('fda_zip','chembl_db','sibling_table'):
        config[key] = str((config_path.parent/Path(config[key])).resolve())
    run(config, args.out)


if __name__ == '__main__':
    main()
