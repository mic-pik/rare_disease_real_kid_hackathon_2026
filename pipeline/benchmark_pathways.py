#!/usr/bin/env python3
"""Frozen family-separated test of uncapped genetics and unsigned pathways.

No fitted weights, network access, LLM calls or clinical nominations. Shared
pathway membership is not a signed downstream relationship or rescue evidence.
"""
import argparse
import json
from pathlib import Path
import shutil
import sqlite3

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from benchmark_recovery import fresh, seal, verify, select_genetics, overlap, tie_metrics
from fda_repurposing import sha256, write_json

METHODS = ('target_top50', 'target_uncapped', 'pathway_only', 'pathway_augmented',
           'target_degree', 'random_expectation')
FEATURES = {'catalog.parquet','mechanisms.parquet','genetics.parquet',
            'pathways.parquet','queries.json','policy.json','directions.json'}
UP = {'AGONIST','ACTIVATOR','POSITIVE ALLOSTERIC MODULATOR','POSITIVE MODULATOR','POTENTIATOR'}
DOWN = {'INHIBITOR','ANTAGONIST','BLOCKER','NEGATIVE ALLOSTERIC MODULATOR','NEGATIVE MODULATOR','DEGRADER'}


def validate_families(queries, ontology):
    by_id = {r['id']:r for r in ontology}
    if len({q['disease_id'] for q in queries}) != len(queries):
        raise ValueError('Duplicate disease')
    audit = []
    for q in queries:
        if q['split'] not in {'development','heldout'}:
            raise ValueError('Unknown split')
        own = by_id[q['disease_id']]
        ancestors = set(own['ancestors']) | {q['disease_id']}
        if q['family_root'] not in ancestors:
            raise ValueError('Family root is not an ancestor')
        if q['name'] != own['name']:
            raise ValueError('Disease name mismatch')
        for other in queries:
            if q['split'] == other['split']: continue
            if other['family_root'] in ancestors or other['disease_id'] in ancestors:
                raise ValueError('Cross-split disease family overlap')
        audit.append({**q,'family_name':by_id[q['family_root']]['name'],
                      'ancestors':sorted(ancestors)})
    return audit


def filter_pathways(frame, policy):
    # Exclude whole pathways if any annotation places them in a forbidden class.
    bad = set(frame.loc[frame.top_level.isna() | frame.top_level.isin(policy['excluded_top_levels']),'pathway_id'])
    sizes = frame.groupby('pathway_id').gene_id.nunique()
    good = set(sizes[sizes.between(policy['pathway_min_genes'],policy['pathway_max_genes'])].index)-bad
    clean = frame[frame.pathway_id.isin(good)][['gene_id','pathway_id','pathway']].drop_duplicates()
    clean['pathway_size'] = clean.pathway_id.map(sizes)
    return clean


def validate_constraints(constraints):
    seen = set()
    for row in constraints:
        if set(row) != {'disease_id','gene_id','required_action','source','locator'}:
            raise ValueError('Invalid direction constraint fields')
        key = row['disease_id'],row['gene_id']
        if key in seen or row['required_action'] not in {'increase','decrease'}:
            raise ValueError('Duplicate or invalid direction constraint')
        if not row['source'] or not row['locator']:
            raise ValueError('Direction constraint requires source and locator')
        seen.add(key)


def direction_status(actions, requirement=None):
    """This checks action compatibility only; even compatible does not mean rescue."""
    if requirement is None: return 'unknown_required_direction'
    signs = set()
    for action in actions:
        action = str(action).strip().upper()
        signs.add('increase' if action in UP else 'decrease' if action in DOWN else 'unknown')
    if {'increase','decrease'} <= signs: return 'mixed_drug_actions'
    opposite = 'decrease' if requirement == 'increase' else 'increase'
    if opposite in signs: return 'opposing_action'
    if 'unknown' in signs or not signs: return 'unknown_drug_action'
    return 'compatible_action_only'


def pathway_links(genetics, pathways, target_genes):
    """Best distinct seed-to-target membership connection; no inferred edge sign."""
    weights = dict(zip(genetics.gene_id,genetics.datatype_score))
    rows = []
    for pid, block in pathways.groupby('pathway_id',sort=True):
        genes = set(block.gene_id)
        seeds = sorted(((float(weights[g]),g) for g in genes & weights.keys()),key=lambda x:(-x[0],x[1]))
        if not seeds: continue
        size = int(block.pathway_size.iloc[0])
        for target in sorted(genes & target_genes):
            seed = next(((w,g) for w,g in seeds if g != target),None)
            if seed is None: continue
            weight,gene = seed
            rows.append({'gene_id':target,'seed_gene_id':gene,'pathway_id':pid,
                         'pathway':block.pathway.iloc[0],'pathway_size':size,
                         'pathway_score':weight/np.sqrt(size-1)})
    columns = ['gene_id','seed_gene_id','pathway_id','pathway','pathway_size','pathway_score']
    result = pd.DataFrame(rows,columns=columns)
    return result.sort_values(['pathway_score','pathway_id','seed_gene_id'],ascending=[False,True,True]).drop_duplicates('gene_id')


def prepare(args):
    frozen, base, kg = Path(args.frozen),Path(args.base),Path(args.kg)
    cfg = json.loads((frozen/'config.json').read_text())
    if sha256(frozen/'config.json') != json.loads((frozen/'freeze.json').read_text())['config_sha256']:
        raise ValueError('Frozen configuration changed')
    verify(base)
    # Pin reuse of the exact v1 candidate/mechanism snapshots.
    catalog = pd.read_parquet(base/'features/catalog.parquet')
    mechanisms = pd.read_parquet(base/'features/mechanisms.parquet')
    queries = cfg['queries']; ids = [q['disease_id'] for q in queries]
    ontology = ds.dataset(args.ontology,format='parquet').to_table(columns=['id','name','ancestors']).to_pylist()
    family_audit = validate_families(queries,ontology)
    validate_constraints(cfg['direction_constraints'])
    evidence = pd.read_parquet(kg/'gene_disease_evidence.parquet')
    genetics = select_genetics(evidence,ids,cfg['genetics_min_score'],None)
    pathways = filter_pathways(pd.read_parquet(kg/'gene_pathway.parquet'),cfg)
    labels = pd.read_parquet(kg/'drug_disease.parquet')
    labels = labels[labels.disease_id.isin(ids) & labels.max_phase.eq('APPROVAL')].copy()
    con = sqlite3.connect('file:'+str(Path(args.chembl).resolve())+'?mode=ro',uri=True)
    con.execute('CREATE TEMP TABLE wanted (id TEXT PRIMARY KEY)')
    con.executemany('INSERT INTO wanted VALUES (?)',[(x,) for x in sorted(set(labels.drug_id))])
    mapping = dict(con.execute('''SELECT m.chembl_id,p.chembl_id FROM wanted w
        JOIN molecule_dictionary m ON m.chembl_id=w.id
        LEFT JOIN molecule_hierarchy h ON h.molregno=m.molregno
        JOIN molecule_dictionary p ON p.molregno=COALESCE(h.parent_molregno,m.molregno)'''))
    con.close()
    labels['ingredient_id'] = labels.drug_id.map(mapping)
    labels['in_universe'] = labels.ingredient_id.isin(catalog.ingredient_id)
    truth = labels[labels.in_universe][['disease_id','ingredient_id']].drop_duplicates()
    if set(truth.disease_id) != set(ids): raise ValueError('Missing label coverage; do not silently replace diseases')
    out = fresh(args.out); feature = fresh(out/'features'); labeldir = fresh(out/'truth')
    catalog.to_parquet(feature/'catalog.parquet',index=False)
    mechanisms.to_parquet(feature/'mechanisms.parquet',index=False)
    genetics.to_parquet(feature/'genetics.parquet',index=False)
    pathways.to_parquet(feature/'pathways.parquet',index=False)
    write_json(feature/'queries.json',queries)
    policy = {k:cfg[k] for k in ['genetics_min_score','baseline_gene_limit','pathway_min_genes','pathway_max_genes','pathway_weight','excluded_top_levels']}
    write_json(feature/'policy.json',policy)
    write_json(feature/'directions.json',cfg['direction_constraints'])
    seal(feature)
    truth.to_csv(labeldir/'labels.csv',index=False);labels.to_csv(labeldir/'mapping_audit.csv',index=False)
    seal(labeldir)
    write_json(out/'config.json',cfg);write_json(out/'family_audit.json',family_audit)
    shutil.copy2(frozen/'freeze.json',out/'freeze.json')
    # Capture executable sources, not just hashes, for the uncommitted implementation.
    code = fresh(out/'code')
    for name in ['benchmark_pathways.py','benchmark_recovery.py','fda_repurposing.py','review_fda_candidates.py']:
        shutil.copy2(Path(__file__).with_name(name),code/name)
    seal(out,base_manifest_sha256=sha256(base/'manifest.json'),
         source_hashes={n:sha256(kg/n) for n in ['gene_disease_evidence.parquet','gene_pathway.parquet','drug_disease.parquet']},
         ontology_hashes={p.name:sha256(p) for p in sorted(Path(args.ontology).glob('*.parquet'))})
    print(f'Prepared {len(queries)} queries and {len(catalog)} candidates; policy frozen; no rankings inspected.')


def score(args):
    source = Path(args.features);verify(source)
    if {p.name for p in source.iterdir()} != FEATURES|{'manifest.json'}: raise ValueError('Feature allowlist mismatch')
    c=pd.read_parquet(source/'catalog.parquet');m=pd.read_parquet(source/'mechanisms.parquet')
    g=pd.read_parquet(source/'genetics.parquet');p=pd.read_parquet(source/'pathways.parquet')
    expected=[(c,{'ingredient_id','canonical_name','smiles','inchikey'}),(m,{'ingredient_id','gene_id','action_type'}),
              (g,{'disease_id','gene_id','datatype_score'}),(p,{'gene_id','pathway_id','pathway','pathway_size'})]
    if any(set(frame.columns)!=cols for frame,cols in expected): raise ValueError('Unexpected feature fields')
    queries=json.loads((source/'queries.json').read_text());policy=json.loads((source/'policy.json').read_text())
    if any(set(q)!={'disease_id','name','split','family_root'} for q in queries): raise ValueError('Unexpected query fields')
    directions=json.loads((source/'directions.json').read_text());validate_constraints(directions)
    required={(r['disease_id'],r['gene_id']):r['required_action'] for r in directions}
    ids=c.ingredient_id.tolist();rows=[];route_rows=[]
    mechanisms=m.groupby(['ingredient_id','gene_id']).action_type.agg(lambda x:sorted(set(str(a) for a in x))).reset_index()
    degree=m.groupby('ingredient_id').gene_id.nunique().reindex(ids).fillna(0).to_numpy(float)
    for q in queries:
        did=q['disease_id'];weights=g[g.disease_id.eq(did)].sort_values(['datatype_score','gene_id'],ascending=[False,True])
        capped=overlap(ids,m,weights.head(policy['baseline_gene_limit']))
        uncapped=overlap(ids,m,weights)
        paths=pathway_links(weights,p,set(m.gene_id))
        links=mechanisms.merge(paths,on='gene_id',how='inner')
        links['disease_id']=did
        links['action_type']=links.action_type.map(lambda x:';'.join(x))
        links['direction_status']=[direction_status(r.action_type.split(';'),required.get((did,r.gene_id))) for r in links.itertuples()]
        links['required_action']=[required.get((did,r.gene_id),'unknown') for r in links.itertuples()]
        links['advancement']='requires_functional_rescue_exposure_and_safety_evidence'
        # Retain conflicts in audit, but do not use them to expand retrieval.
        usable=links[~links.direction_status.isin(['opposing_action','mixed_drug_actions'])]
        neighbour=usable.groupby('ingredient_id').pathway_score.max().reindex(ids).fillna(0).to_numpy(float)
        route_rows.extend(links.to_dict('records'))
        values={'target_top50':capped,'target_uncapped':uncapped,'pathway_only':neighbour,
                'pathway_augmented':np.maximum(uncapped,policy['pathway_weight']*neighbour),
                'target_degree':degree,'random_expectation':np.zeros(len(ids))}
        for method,scores in values.items():
            rows.extend({'disease_id':did,'ingredient_id':cid,'method':method,'score':float(v)} for cid,v in zip(ids,scores))
    out=fresh(args.out);pd.DataFrame(rows).to_parquet(out/'scores.parquet',index=False)
    pd.DataFrame(route_rows).to_parquet(out/'pathway_routes.parquet',index=False)
    seal(out,features_manifest_sha256=sha256(source/'manifest.json'),implementation_sha256=sha256(__file__))
    print('Scored frozen methods without reading indication labels.')


def paired_delta(metrics,method,baseline,seed,n):
    pivot=metrics.pivot(index='disease_id',columns='method',values='recall@10')
    delta=(pivot[method]-pivot[baseline]).to_numpy()
    rng=np.random.default_rng(seed)
    samples=delta[rng.integers(0,len(delta),size=(n,len(delta)))].mean(axis=1)
    low,high=np.quantile(samples,[.025,.975])
    return {'method':method,'baseline':baseline,'mean_delta':float(delta.mean()),
            'lower95':float(low),'upper95':float(high),'n_diseases':len(delta),
            'improved':int((delta>0).sum()),'worsened':int((delta<0).sum()),'tied':int((delta==0).sum())}


def evaluate(args):
    prepared,scored=Path(args.prepared),Path(args.scores);verify(prepared);sm=verify(scored)
    if sm['features_manifest_sha256'] != sha256(prepared/'features/manifest.json'):raise ValueError('Score input mismatch')
    # Heldout outcomes must use the same executable frozen before scoring.
    if sha256(__file__) != sha256(prepared/'code/benchmark_pathways.py') or sm['implementation_sha256'] != sha256(__file__):
        raise ValueError('Implementation changed since freeze; use the captured code')
    cfg=json.loads((prepared/'config.json').read_text())
    queries=cfg['queries'];c=pd.read_parquet(prepared/'features/catalog.parquet')
    truth=pd.read_csv(prepared/'truth/labels.csv');scores=pd.read_parquet(scored/'scores.parquet')
    if set(scores.disease_id)!={q['disease_id'] for q in queries} or set(scores.method)!=set(METHODS):raise ValueError('Missing queries or methods')
    rows=[];ranks=[]
    for q in queries:
        positives=set(truth.loc[truth.disease_id.eq(q['disease_id']),'ingredient_id'])
        for method in METHODS:
            f=scores[scores.disease_id.eq(q['disease_id'])&scores.method.eq(method)].copy()
            if f.ingredient_id.duplicated().any() or set(f.ingredient_id)!=set(c.ingredient_id):raise ValueError('Candidate coverage mismatch')
            rows.append({**q,'method':method,**tie_metrics(f.score,f.ingredient_id.isin(positives))})
            f['rank_min']=f.score.rank(ascending=False,method='min').astype(int)
            f['rank_max']=f.score.rank(ascending=False,method='max').astype(int)
            f['known_positive']=f.ingredient_id.isin(positives);f['split']=q['split']
            ranks.extend(f.to_dict('records'))
    out=fresh(args.out);metrics=pd.DataFrame(rows);metrics.to_csv(out/'per_disease.csv',index=False)
    ranks=pd.DataFrame(ranks).merge(c[['ingredient_id','canonical_name']],on='ingredient_id')
    ranks.to_parquet(out/'rankings.parquet',index=False)
    macro=metrics.groupby(['split','method'])[['recall@10','hit@10','mrr','positive_score_coverage']].mean()
    macro.to_csv(out/'macro.csv');comparisons=[]
    for split in ['development','heldout']:
        for method,baseline in [('target_uncapped','target_top50'),('pathway_augmented','target_top50'),('pathway_augmented','target_uncapped')]:
            comparisons.append({'split':split,**paired_delta(metrics[metrics.split.eq(split)],method,baseline,cfg['bootstrap_seed'],cfg['bootstrap_samples'])})
    pd.DataFrame(comparisons).to_csv(out/'paired_comparisons.csv',index=False)
    h=[r for r in comparisons if r['split']=='heldout']
    promote_pathway=all(r['mean_delta']>0 and r['lower95']>0 for r in h if r['method']=='pathway_augmented')
    promote_uncapped=next(r for r in h if r['method']=='target_uncapped')['lower95']>0
    decision={'pathway_retrieval_promotion':promote_pathway,'uncapped_retrieval_promotion':promote_uncapped,
              'clinical_or_MVA_promotion':False,'required_direction_constraints':len(cfg['direction_constraints']),
              'warning':'Five heldout diseases provide limited precision; paired bootstrap is exploratory. Unsigned membership is not compensation.'}
    write_json(out/'decision.json',decision)
    lines=['# Uncapped genetics and pathway benchmark','',
           'Indications withheld from scoring. Fixed methods; five heldout families; no temporal/publication holdout. Unlabelled drugs are not negatives. Direction compatibility is not efficacy.','',
           '| Split | Method | Recall@10 | Hit@10 | MRR |','|---|---|---:|---:|---:|']
    for (split,method),r in macro.iterrows():lines.append(f"| {split} | {method} | {r['recall@10']:.3f} | {r['hit@10']:.3f} | {r['mrr']:.3f} |")
    lines+=['','## Prespecified decision','',json.dumps(decision,indent=2),'',
            'See paired_comparisons.csv for changes and bootstrap intervals, per_disease.csv for individual outcomes, rankings.parquet for complete ranks and the scores directory for auditable pathway connections.','']
    (out/'report.md').write_text('\n'.join(lines))
    seal(out,prepared_manifest_sha256=sha256(prepared/'manifest.json'),scored_manifest_sha256=sha256(scored/'manifest.json'),implementation_sha256=sha256(__file__))
    print(macro.to_string());print(json.dumps(decision))


def main():
    parser=argparse.ArgumentParser(description=__doc__);subs=parser.add_subparsers(dest='command',required=True)
    p=subs.add_parser('prepare')
    for key in ['frozen','base','kg','ontology','chembl','out']:p.add_argument('--'+key,required=True)
    p=subs.add_parser('score')
    for key in ['features','out']:p.add_argument('--'+key,required=True)
    p=subs.add_parser('evaluate')
    for key in ['prepared','scores','out']:p.add_argument('--'+key,required=True)
    args=parser.parse_args();{'prepare':prepare,'score':score,'evaluate':evaluate}[args.command](args)


if __name__=='__main__':main()
