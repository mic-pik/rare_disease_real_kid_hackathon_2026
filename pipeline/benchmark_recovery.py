#!/usr/bin/env python3
"""Prepare blinded features, rank without labels, and evaluate FDA recovery.

Three separate commands. No LLM, network, learned weights or indication-based
candidate filtering. Nonpositive labels are unlabelled, never known negatives.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3

import numpy as np
import pandas as pd

from fda_repurposing import measured_interactions, sha256, write_json
from review_fda_candidates import assay_triage

FEATURES = {'catalog.parquet', 'genetics.parquet', 'mechanisms.parquet',
            'active_assays.parquet', 'queries.json'}
METHODS = ('target_overlap', 'assay_overlap', 'chemical_similarity', 'fixed_mean',
           'target_degree', 'random_expectation')


def fresh(path):
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    return path


def seal(path, **metadata):
    write_json(path/'manifest.json', {**metadata,
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'output_hashes': {str(p.relative_to(path)): sha256(p)
                          for p in sorted(path.rglob('*')) if p.is_file() and p != path/'manifest.json'}})


def verify(path):
    m = json.loads((path/'manifest.json').read_text())
    for name, digest in m['output_hashes'].items():
        if sha256(path/name) != digest:
            raise ValueError('Hash mismatch: '+name)
    return m


def select_genetics(frame, disease_ids, minimum, maximum):
    # Never use the overall score: it includes known-drug evidence.
    g = frame[frame.datatype.eq('genetic_association') &
              frame.disease_id.isin(disease_ids) & frame.datatype_score.ge(minimum)]
    g = g.groupby(['disease_id', 'gene_id'], as_index=False).datatype_score.max()
    g = g.sort_values(['disease_id', 'datatype_score', 'gene_id'],
                      ascending=[True, False, True])
    return g if maximum is None else g.groupby('disease_id').head(maximum)


def prepare(args):
    cfg = json.loads(Path(args.config).read_text())
    kg, fda = Path(args.kg), Path(args.fda)
    verify(fda)
    queries = cfg['diseases']; ids = [q['disease_id'] for q in queries]
    if len(set(ids)) != len(ids):
        raise ValueError('Duplicate disease query')
    catalog = pd.read_parquet(fda/'fda_ingredients.parquet')
    # Retain every currently eligible mapped ingredient, even with no target data.
    catalog = catalog[catalog.eligible_product_ids.map(len).gt(0)].copy()
    allowed = set(catalog.ingredient_id)
    evidence = pd.read_parquet(kg/'gene_disease_evidence.parquet')
    genetics = select_genetics(evidence, ids, cfg['genetics_min_score'], cfg['max_genes_per_disease'])
    genes = pd.read_parquet(kg/'genes.parquet', columns=['gene_id', 'symbol'])
    genetics = genetics.merge(genes, on='gene_id', how='left', validate='many_to_one')
    if genetics.symbol.isna().any():
        raise ValueError('Missing gene symbol')
    dt = pd.read_parquet(kg/'drug_target.parquet', columns=['drug_id', 'gene_id', 'action_type'])
    dd = pd.read_parquet(kg/'drug_disease.parquet')
    dd = dd[dd.disease_id.isin(ids) & dd.max_phase.eq('APPROVAL')].copy()
    con = sqlite3.connect('file:'+str(Path(args.chembl).resolve())+'?mode=ro', uri=True)
    con.execute('CREATE TEMP TABLE wanted_ids (id TEXT PRIMARY KEY)')
    con.executemany('INSERT INTO wanted_ids VALUES (?)', [(x,) for x in sorted(set(dt.drug_id) | set(dd.drug_id))])
    mapping = dict(con.execute('''SELECT m.chembl_id, p.chembl_id
        FROM wanted_ids w JOIN molecule_dictionary m ON m.chembl_id=w.id
        LEFT JOIN molecule_hierarchy h ON h.molregno=m.molregno
        JOIN molecule_dictionary p ON p.molregno=COALESCE(h.parent_molregno,m.molregno)'''))
    dt['ingredient_id'] = dt.drug_id.map(mapping)
    dt = dt[dt.ingredient_id.isin(allowed)][['ingredient_id', 'gene_id', 'action_type']].drop_duplicates()
    dd['ingredient_id'] = dd.drug_id.map(mapping)
    dd['in_candidate_universe'] = dd.ingredient_id.isin(allowed)
    labels = dd[dd.in_candidate_universe][['disease_id', 'ingredient_id']].drop_duplicates()
    missing = set(ids)-set(labels.disease_id)
    if missing:
        raise ValueError('No eligible approval labels: '+str(sorted(missing)))
    ingredient_records = catalog.to_dict('records')
    for record in ingredient_records:
        record['eligible_product_ids'] = list(record['eligible_product_ids'])
    assays, coverage = measured_interactions(con, ingredient_records, genetics.symbol.unique())
    con.close()
    print(f'Extracted {len(assays)} original assay rows', flush=True)
    triage = pd.DataFrame([assay_triage(r) for r in assays.to_dict('records')])
    assays = pd.concat([assays.reset_index(drop=True), triage], axis=1)
    active = assays[assays.activity_class.eq('active') & assays.eligible_for_support_count]
    active = active[['ingredient_id', 'gene', 'activity_id']].merge(
        genes.rename(columns={'symbol':'gene'}), on='gene', how='inner')
    active = active[['ingredient_id', 'gene_id', 'activity_id']].drop_duplicates()
    out = fresh(args.out); features = fresh(out/'features')
    catalog[['ingredient_id', 'canonical_name', 'smiles', 'inchikey']].sort_values('ingredient_id').to_parquet(features/'catalog.parquet', index=False)
    genetics.to_parquet(features/'genetics.parquet', index=False)
    dt.to_parquet(features/'mechanisms.parquet', index=False)
    active.to_parquet(features/'active_assays.parquet', index=False)
    # Query names only; curated sentinels/indication labels never enter features.
    write_json(features/'queries.json', queries)
    seal(features, allowed_files=sorted(FEATURES), evidence_datatype='genetic_association')
    truth = fresh(out/'truth')
    labels.to_csv(truth/'labels.csv', index=False)
    dd.to_csv(truth/'label_mapping_audit.csv', index=False)
    write_json(truth/'sentinels.json', cfg.get('sentinels', []))
    seal(truth, label_level='approval indication metadata intersected with current FDA ingredient catalogue')
    assays.to_parquet(out/'original_assays_triaged.parquet', index=False)
    coverage.to_csv(out/'target_coverage.csv', index=False)
    write_json(out/'config.json', cfg)
    seal(out, config_sha256=sha256(args.config), implementation_sha256=sha256(__file__),
         fda_manifest_sha256=sha256(fda/'manifest.json'),
         kg_hashes={n:sha256(kg/n) for n in ['gene_disease_evidence.parquet','genes.parquet','drug_target.parquet','drug_disease.parquet']},
         chembl={'path':str(Path(args.chembl).resolve()),'size':Path(args.chembl).stat().st_size,
                 'mtime_ns':Path(args.chembl).stat().st_mtime_ns,
                 'note':'DB not fully hashed here; extracted original assay rows and all feature tables are hashed.'})
    print(json.dumps({'candidates':len(catalog),'queries':len(queries),'positives':len(labels),
                      'selected_genes':genetics.gene_id.nunique(),'active_assay_links':len(active)}))


def overlap(catalog_ids, links, genetics):
    joined = links[['ingredient_id', 'gene_id']].drop_duplicates().merge(
        genetics[['gene_id', 'datatype_score']], on='gene_id', how='inner')
    scores = joined.groupby('ingredient_id').datatype_score.max()
    return scores.reindex(catalog_ids).fillna(0).to_numpy(float)


def score(args):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    src = Path(args.features)
    # Fail closed if a truth table or an unexpected field is supplied.
    actual = {p.name for p in src.iterdir() if p.is_file()}-{'manifest.json'}
    if actual != FEATURES:
        raise ValueError('Feature allowlist mismatch')
    verify(src)
    c = pd.read_parquet(src/'catalog.parquet')
    g = pd.read_parquet(src/'genetics.parquet')
    m = pd.read_parquet(src/'mechanisms.parquet')
    a = pd.read_parquet(src/'active_assays.parquet')
    schemas = [(c, {'ingredient_id','canonical_name','smiles','inchikey'}),
               (g, {'disease_id','gene_id','datatype_score','symbol'}),
               (m, {'ingredient_id','gene_id','action_type'}),
               (a, {'ingredient_id','gene_id','activity_id'})]
    if any(set(frame.columns) != columns for frame,columns in schemas):
        raise ValueError('Unexpected feature column')
    queries = json.loads((src/'queries.json').read_text())
    if any(set(q) != {'disease_id','name'} for q in queries):
        raise ValueError('Unexpected query field')
    ids = c.ingredient_id.tolist()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = []
    for smi in c.smiles:
        mol = Chem.MolFromSmiles(smi)
        fps.append(gen.GetFingerprint(mol) if mol is not None else None)
    sim = np.zeros((len(c),len(c)),dtype=np.float32)
    valid = [i for i,fp in enumerate(fps) if fp is not None]
    for i in valid:
        sim[i,valid] = DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in valid])
    # Similarity excludes self and identical connectivity keys (salt/stereo aliases).
    keys = c.inchikey.str.split('-').str[0].to_numpy()
    sim[keys[:,None] == keys[None,:]] = 0
    degree = m.groupby('ingredient_id').gene_id.nunique().reindex(ids).fillna(0).to_numpy(float)
    rows = []
    for query in queries:
        did = query['disease_id']; weights = g[g.disease_id.eq(did)]
        target = overlap(ids,m,weights); assay = overlap(ids,a,weights)
        chemical = (sim*target[None,:]).max(axis=1)
        values = {'target_overlap':target, 'assay_overlap':assay,
                  'chemical_similarity':chemical, 'fixed_mean':(target+assay+chemical)/3,
                  'target_degree':degree, 'random_expectation':np.zeros(len(c))}
        for method, scores in values.items():
            rows.extend({'disease_id':did,'ingredient_id':cid,'method':method,'score':float(s)}
                        for cid,s in zip(ids,scores))
    out = fresh(args.out)
    pd.DataFrame(rows).to_parquet(out/'scores.parquet',index=False)
    write_json(out/'scoring_policy.json', {
        'methods':list(METHODS), 'chemical_similarity':'max Morgan r2/2048 Tanimoto times target-anchor genetic score; excludes same connectivity',
        'overlap':'maximum selected disease genetic score among linked genes',
        'fixed_mean':'equal mean of target, assay and chemical scores; no fitting',
        'limitations':'direction, exposure and clinical efficacy not modelled',
        'invalid_smiles':int(len(c)-len(valid))})
    seal(out, feature_manifest_sha256=sha256(src/'manifest.json'), implementation_sha256=sha256(__file__))
    print('Ranked all candidates without loading truth or indication tables.',flush=True)


def tie_metrics(scores, positive, ks=(1,5,10,50)):
    """Exact expected recall/hit/MRR under uniform ordering within score ties."""
    scores = np.asarray(scores,float); positive = np.asarray(positive,bool)
    if not np.isfinite(scores).all() or not positive.any():
        raise ValueError('Finite scores and at least one positive required')
    total = int(positive.sum()); result = {}; groups = []
    for s in sorted(set(scores), reverse=True):
        mask = scores == s; groups.append((int(mask.sum()),int(positive[mask].sum())))
    for k in ks:
        remaining = min(k,len(scores)); hits = 0.; no_hit = 1.
        for n,p in groups:
            take = min(remaining,n)
            hits += take*p/n
            for j in range(take):
                no_hit *= max(0.,(n-p-j)/(n-j))
            remaining -= take
            if remaining == 0: break
        result[f'recall@{k}'] = hits/total
        result[f'hit@{k}'] = 1-no_hit
    before = 0; mrr = 0.
    for n,p in groups:
        if p:
            survive = 1.
            for j in range(1,n-p+2):
                prob = survive*p/(n-j+1)
                mrr += prob/(before+j)
                survive *= (n-p-j+1)/(n-j+1)
            break
        before += n
    result['mrr'] = mrr
    result['positive_score_coverage'] = float(np.mean(scores[positive]>0))
    result['n_positives'] = total
    return result


def evaluate(args):
    prepared, scored = Path(args.prepared), Path(args.scores)
    verify(prepared); sm = verify(scored)
    if sm['feature_manifest_sha256'] != sha256(prepared/'features/manifest.json'):
        raise ValueError('Scores belong to different features')
    c = pd.read_parquet(prepared/'features/catalog.parquet')
    labels = pd.read_csv(prepared/'truth/labels.csv')
    scores = pd.read_parquet(scored/'scores.parquet')
    queries = json.loads((prepared/'features/queries.json').read_text())
    if set(scores.disease_id) != {q['disease_id'] for q in queries} or set(scores.method) != set(METHODS):
        raise ValueError('Missing or extra queries/methods')
    rows, failures, tops, sentinel_rows = [], [], [], []
    sentinels = json.loads((prepared/'truth/sentinels.json').read_text())
    for query in queries:
        did = query['disease_id']; truth = set(labels.loc[labels.disease_id.eq(did),'ingredient_id'])
        for method in METHODS:
            f = scores[scores.disease_id.eq(did) & scores.method.eq(method)].copy()
            if f.ingredient_id.duplicated().any() or set(f.ingredient_id) != set(c.ingredient_id):
                raise ValueError('Candidate universe mismatch')
            rows.append({'disease_id':did,'disease':query['name'],'method':method,
                         **tie_metrics(f.score,f.ingredient_id.isin(truth))})
            f['rank_min'] = f.score.rank(ascending=False,method='min').astype(int)
            f['rank_max'] = f.score.rank(ascending=False,method='max').astype(int)
            f = f.merge(c[['ingredient_id','canonical_name']],on='ingredient_id')
            f['known_positive'] = f.ingredient_id.isin(truth)
            tops.extend(f.sort_values(['score','ingredient_id'],ascending=[False,True]).head(10).to_dict('records'))
            for r in f[f.known_positive].to_dict('records'):
                r['diagnosis'] = 'no_positive_score' if r['score'] <= 0 else ('outside_top10' if r['rank_min']>10 else ('tie_crosses_top10' if r['rank_max']>10 else 'recovered_top10'))
                failures.append(r)
            for sentinel in sentinels:
                if sentinel['disease_id'] == did:
                    match = f[f.ingredient_id.eq(sentinel['ingredient_id'])]
                    if match.empty: raise ValueError('Sentinel absent from FDA universe')
                    r = match.iloc[0].to_dict()
                    if r['canonical_name'] != sentinel['name']: raise ValueError('Sentinel identity mismatch')
                    sentinel_rows.append({**sentinel,**r})
    out = fresh(args.out); metrics = pd.DataFrame(rows)
    metrics.to_csv(out/'per_disease_metrics.csv',index=False)
    macro = metrics.groupby('method')[[x for x in metrics if '@' in x or x in ('mrr','positive_score_coverage')]].mean()
    macro.to_csv(out/'macro_metrics.csv')
    pd.DataFrame(failures).to_csv(out/'positive_diagnostics.csv',index=False)
    pd.DataFrame(tops).to_csv(out/'top10.csv',index=False)
    pd.DataFrame(sentinel_rows).to_csv(out/'sentinels.csv',index=False)
    lines = ['# FDA drug recovery: fixed development pilot','',
        f'{len(queries)} diseases; {len(c)} candidates per disease; {len(labels)} unique metadata-positive drug–disease pairs. All other candidates are unlabelled, not known negatives.', '',
        'Indications are held out of scoring. This is not a temporal/publication split or validation of the LLM, clinical efficacy, direction or exposure. Broad labels can include symptomatic and combination uses.', '',
        'Metrics average diseases equally and use exact expected ordering within ties. Random expectation is an all-tied ranking. Top-10 CSV uses ingredient ID only for display; inspect rank intervals.', '',
        '| Method | Recall@10 | Hit@10 | MRR | Positive-score coverage |', '|---|---:|---:|---:|---:|']
    for method,r in macro.iterrows():
        lines.append(f"| {method} | {r['recall@10']:.3f} | {r['hit@10']:.3f} | {r['mrr']:.3f} | {r['positive_score_coverage']:.3f} |")
    lines += ['', 'See per_disease_metrics.csv, positive_diagnostics.csv, sentinels.csv and the prepared label_mapping_audit.csv for coverage and failures. No method was selected or tuned using these outcomes.', '']
    (out/'report.md').write_text('\n'.join(lines))
    seal(out, prepared_manifest_sha256=sha256(prepared/'manifest.json'), score_manifest_sha256=sha256(scored/'manifest.json'), implementation_sha256=sha256(__file__))
    print(macro.to_string())


def diagnose(args):
    """Post-evaluation coverage audit; never feeds scores or retunes policy."""
    prepared, kg = Path(args.prepared), Path(args.kg)
    verify(prepared)
    cfg = json.loads((prepared/'config.json').read_text())
    labels = pd.read_csv(prepared/'truth/labels.csv')
    mechanisms = pd.read_parquet(prepared/'features/mechanisms.parquet')
    selected = pd.read_parquet(prepared/'features/genetics.parquet')
    source = kg/'gene_disease_evidence.parquet'
    expected = json.loads((prepared/'manifest.json').read_text())['kg_hashes'][source.name]
    if sha256(source) != expected:
        raise ValueError('Diagnostic evidence differs from prepared snapshot')
    all_genetics = pd.read_parquet(source, filters=[('datatype','==','genetic_association')])
    all_genetics = all_genetics[all_genetics.disease_id.isin(labels.disease_id.unique())]
    rows = []
    for row in labels.itertuples():
        genes = set(mechanisms.loc[mechanisms.ingredient_id.eq(row.ingredient_id),'gene_id'])
        raw = all_genetics[all_genetics.disease_id.eq(row.disease_id) & all_genetics.gene_id.isin(genes)]
        chosen = selected[selected.disease_id.eq(row.disease_id) & selected.gene_id.isin(genes)]
        best = float(raw.datatype_score.max()) if len(raw) else 0.
        reason = ('selected_target_supported' if len(chosen) else
                  'no_mechanism_annotation' if not genes else
                  'no_direct_genetic_support' if best <= 0 else
                  'below_genetic_threshold' if best < cfg['genetics_min_score'] else
                  'excluded_by_gene_limit')
        rows.append({'disease_id':row.disease_id,'ingredient_id':row.ingredient_id,
                     'mechanism_target_count':len(genes),'best_unfiltered_genetic_score':best,
                     'target_route_failure':reason})
    out = fresh(args.out); frame = pd.DataFrame(rows)
    frame.to_csv(out/'target_route_diagnostics.csv',index=False)
    counts = {str(k):int(v) for k,v in frame.target_route_failure.value_counts().items()}
    write_json(out/'summary.json', counts)
    seal(out, prepared_manifest_sha256=sha256(prepared/'manifest.json'),
         source_sha256=expected, implementation_sha256=sha256(__file__))
    print(json.dumps(counts,indent=2))


def main():
    p = argparse.ArgumentParser(description=__doc__); sub = p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('prepare')
    for name in ('config','kg','fda','chembl','out'): a.add_argument('--'+name,required=True)
    a=sub.add_parser('score')
    for name in ('features','out'): a.add_argument('--'+name,required=True)
    a=sub.add_parser('evaluate')
    for name in ('prepared','scores','out'): a.add_argument('--'+name,required=True)
    a=sub.add_parser('diagnose')
    for name in ('prepared','kg','out'): a.add_argument('--'+name,required=True)
    args=p.parse_args(); {'prepare':prepare,'score':score,'evaluate':evaluate,'diagnose':diagnose}[args.command](args)


if __name__ == '__main__': main()
