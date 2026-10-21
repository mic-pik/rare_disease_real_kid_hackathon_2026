#!/usr/bin/env python3
"""Exploratory FDA compound reversal of a single-case expression signature."""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


def weights(direction, magnitude=None):
    direction = np.asarray(direction)
    magnitude = np.ones(len(direction)) if magnitude is None else np.asarray(magnitude)
    w = np.zeros(len(direction))
    for sign in [-1, 1]:
        mask = direction == sign
        if not mask.any() or magnitude[mask].sum() <= 0:
            raise ValueError("Both query directions require positive weight")
        w[mask] = -sign * .5 * magnitude[mask] / magnitude[mask].sum()
    return w


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while block := f.read(8 * 1024 * 1024):
            h.update(block)
    return h.hexdigest()


def background_check(conditions, minimum=20):
    """Post-hoc sensitivity: equal-weight other drugs in exactly matched contexts."""
    keys = ['cell_id','pert_time','pert_time_unit','pert_dose','pert_dose_unit']
    rows = []
    for _, group in conditions.groupby(keys,dropna=False):
        for row in group.itertuples(index=False):
            others = group[group.ingredient_id != row.ingredient_id].groupby('ingredient_id').primary.mean()
            d = row._asdict()
            d['background_n_ingredients'] = len(others)
            d['background_median'] = float(others.median()) if len(others)>=minimum else np.nan
            d['background_adjusted'] = row.primary-d['background_median']
            d['background_percentile'] = float(((others<row.primary).sum()+.5*(others==row.primary).sum())/len(others)) if len(others)>=minimum else np.nan
            rows.append(d)
    return pd.DataFrame(rows)


def run(args):
    out = args.output
    out.mkdir(parents=True, exist_ok=True)
    query = pd.read_csv(args.query / 'lincs_measured_landmark_query.csv', dtype={'entrez_id': str})
    info_path = args.lincs / 'GSE92742_Broad_LINCS_gene_info.txt.gz'
    info = pd.read_csv(info_path, sep='\t', dtype={'pr_gene_id': str})
    genes = pd.read_csv(args.query / 'gene_contrasts.csv')
    probes = pd.read_csv(args.query / 'probe_contrasts.csv')
    ids = probes[probes.eligible].groupby('symbol').entrez_id.agg(lambda s: sorted(set(str(int(v)) for v in s.dropna())))
    ids = ids[ids.map(len) == 1].map(lambda s: s[0])
    universe = sorted(set(ids) & set(info.loc[info.pr_is_lm == 1, 'pr_gene_id']))
    assert set(query.entrez_id) <= set(universe) and query.entrez_id.is_unique
    si_path = args.lincs / 'GSE92742_Broad_LINCS_sig_info.txt.gz'
    met_path = args.lincs / 'GSE92742_Broad_LINCS_sig_metrics.txt.gz'
    pert_path = args.lincs / 'GSE92742_Broad_LINCS_pert_info.txt.gz'
    si = pd.read_csv(si_path, sep='\t', low_memory=False)
    met = pd.read_csv(met_path, sep='\t', low_memory=False)
    pert = pd.read_csv(pert_path, sep='\t', low_memory=False)
    fda = pd.read_parquet(args.fda)
    fda = fda[fda.eligible_product_ids.map(len) > 0]
    mapping = pert[pert.pert_type == 'trt_cp'].merge(fda, left_on='inchi_key', right_on='inchikey')
    ambiguous = mapping.groupby('pert_id').ingredient_id.nunique()
    mapping = mapping[~mapping.pert_id.isin(ambiguous[ambiguous > 1].index)]
    mapping = mapping[['pert_id','pert_iname','inchi_key','ingredient_id','canonical_name']].drop_duplicates()
    mapping.to_csv(out / 'identity_mapping.csv', index=False)
    s = si[si.pert_type == 'trt_cp'].merge(mapping.drop(columns='pert_iname'), on='pert_id', validate='many_to_one')
    s = s.merge(met[['sig_id','is_exemplar','distil_cc_q75','distil_nsample','tas']], on='sig_id', validate='one_to_one')
    total = len(s)
    s = s[(s.is_exemplar == 1) & (s.distil_cc_q75 >= .2) & (s.distil_nsample >= 3)].copy()
    quality = len(s)
    # Do not count overlapping constituent experiments as independent profiles.
    used, selected, rejected = set(), [], []
    for row in s.sort_values(['distil_cc_q75','sig_id'], ascending=[False,True]).itertuples():
        distils = set(str(row.distil_id).split('|'))
        if distils & used or not distils or '-666' in distils:
            rejected.append(row.sig_id)
        else:
            used.update(distils)
            selected.append(row.sig_id)
    s = s[s.sig_id.isin(selected)].sort_values('sig_id').reset_index(drop=True)
    pd.DataFrame({'sig_id':rejected}).to_csv(out / 'overlapping_or_missing_constituents.csv', index=False)
    gctx = args.lincs / 'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
    decode = lambda a: [v.decode() if isinstance(v, bytes) else str(v) for v in a]
    with h5py.File(gctx) as f:
        rid = decode(f['0/META/ROW/id'][:]); cid = decode(f['0/META/COL/id'][:])
        ri = {v:i for i,v in enumerate(rid)}; ci = {v:i for i,v in enumerate(cid)}
        assert set(s.sig_id) <= set(ci) and set(universe) <= set(ri)
        columns = [ri[v] for v in universe]
        x = np.empty((len(s),len(universe)), dtype=np.float32)
        for i, sig in enumerate(s.sig_id):
            x[i] = f['0/DATA/0/matrix'][ci[sig],:][columns]
            if i % 500 == 0: print(f'Read {i}/{len(s)} profiles', flush=True)
    assert np.isfinite(x).all()
    # Rank normalization limits outlier z-scores; null queries use this same universe.
    ranks = pd.DataFrame(x).rank(axis=1, method='average').to_numpy()
    ranks = 2 * (ranks - 1) / (len(universe) - 1) - 1
    qidx = [universe.index(v) for v in query.entrez_id]
    variants = {}
    for name, mask, mag in [
        ('primary', np.ones(len(query),dtype=bool), None),
        ('relative_supported', query.relative_direction_supported.to_numpy(), None),
        ('fibro_weighted', np.ones(len(query),dtype=bool), query.fibro_delta.abs().to_numpy()),
        ('lcl_weighted', np.ones(len(query),dtype=bool), query.lcl_delta.abs().to_numpy()),
        ('effect025', (query.minimum_abs_delta >= .25).to_numpy(), None),
    ]:
        w = np.zeros(len(universe))
        w[np.array(qidx)[mask]] = weights(query.robust_direction.to_numpy()[mask], None if mag is None else mag[mask])
        variants[name] = w
        s[name] = ranks @ w
    condition_keys = ['ingredient_id','canonical_name','cell_id','pert_time','pert_time_unit','pert_dose','pert_dose_unit']
    s.to_csv(out / 'signature_scores.csv', index=False)
    conditions = s.groupby(condition_keys, dropna=False)[list(variants)].mean().reset_index()
    conditions.to_csv(out / 'condition_scores.csv', index=False)
    background = background_check(conditions)
    background.to_csv(out / 'matched_drug_background.csv',index=False)
    # Balanced aggregation: each condition within a cell, then each cell, gets equal weight.
    frame = pd.DataFrame(ranks, columns=universe)
    for key in condition_keys: frame[key] = s[key]
    cf = frame.groupby(condition_keys,dropna=False)[universe].mean().reset_index()
    cell = cf.groupby(['ingredient_id','canonical_name','cell_id'])[universe].mean()
    agg = cell.groupby(['ingredient_id','canonical_name']).mean()
    results = pd.DataFrame(index=agg.index)
    for name,w in variants.items(): results[name] = agg.to_numpy() @ w
    results['n_cells'] = cell.groupby(level=[0,1]).size()
    cell_score = cell.to_numpy() @ variants['primary']
    cell_results = cell.reset_index()[['ingredient_id','canonical_name','cell_id']]
    cell_results['primary'] = cell_score
    cell_results.to_csv(out / 'cell_scores.csv',index=False)
    results['positive_cell_fraction'] = cell_results.assign(positive=cell_results.primary>0).groupby(['ingredient_id','canonical_name']).positive.mean()
    results['n_signatures'] = s.groupby(['ingredient_id','canonical_name']).size()
    rng = np.random.default_rng(args.seed)
    nup = int((query.robust_direction==1).sum()); ndown = len(query)-nup
    nullw = np.zeros((len(universe),args.nulls),dtype=np.float32)
    for i in range(args.nulls):
        draw = rng.choice(len(universe),len(query),replace=False)
        nullw[draw[:nup],i] = -.5/nup
        nullw[draw[nup:],i] = .5/ndown
    null = agg.to_numpy() @ nullw
    observed = results.primary.to_numpy()
    results['random_gene_p'] = (1+(null >= observed[:,None]).sum(axis=1))/(args.nulls+1)
    results['random_gene_bh_q'] = bh(results.random_gene_p)
    results['all_sensitivity_positive'] = (results[list(variants)]>0).all(axis=1)
    results['computational_followup'] = (results.random_gene_bh_q<=.1)&(results.n_cells>=2)&(results.positive_cell_fraction>=.75)&results.all_sensitivity_positive
    bgcell = background.dropna(subset=['background_adjusted']).groupby(['ingredient_id','canonical_name','cell_id'])[['background_adjusted','background_percentile']].mean()
    bgsummary = bgcell.groupby(level=[0,1]).mean()
    bgsummary['background_covered_cells'] = bgcell.groupby(level=[0,1]).size()
    results = results.join(bgsummary)
    results['posthoc_background_supported'] = results.computational_followup & (results.background_covered_cells>=2) & (results.background_adjusted>0) & (results.background_percentile>=.75)
    results = results.sort_values(['computational_followup','primary'],ascending=[False,False])
    results.to_csv(out / 'compound_scores.csv')
    results[results.computational_followup].to_csv(out / 'computational_followup.csv')
    np.savez_compressed(out / 'null_scores.npz', scores=null, ingredient_ids=agg.index.get_level_values(0).to_numpy(dtype=str))
    np.savez_compressed(out / 'selected_landmark_matrix.npz', z=x, sig_ids=s.sig_id.to_numpy(dtype=str), gene_ids=np.array(universe))
    summary = {'mapped_ingredients':int(mapping.ingredient_id.nunique()),'mapped_profiles':total,'quality_profiles':quality,
               'nonoverlapping_profiles':len(s),'scored_ingredients':len(results),'null_universe_genes':len(universe),
               'query_up':nup,'query_down':ndown,'nulls':args.nulls,'seed':args.seed,
               'variant_gene_counts':{k:int((v!=0).sum()) for k,v in variants.items()},
               'computational_followup_count':int(results.computational_followup.sum()),
               'posthoc_background_supported_count':int(results.posthoc_background_supported.sum()),
               'decision':'Expression-reversal hypotheses only; no functional rescue or exposure validation',
               'limitations':['single affected donor','unmatched cell types','no growth/viability adjustment','random-gene null does not preserve pathway coexpression','tissue sensitivities reweight the same consensus genes, not independent tissue signatures']}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    inputs = [info_path,si_path,met_path,pert_path,args.fda,args.query/'lincs_measured_landmark_query.csv',args.query/'probe_contrasts.csv',args.query/'gene_contrasts.csv']
    manifest = {'input_sha256':{p.name:digest(p) for p in inputs},'script_sha256':digest(Path(__file__)),
                'gctx':{'name':gctx.name,'size':gctx.stat().st_size,'mtime_ns':gctx.stat().st_mtime_ns,'hash_scope':'selected matrix and identifiers archived and hashed; full 23GB source not rehashed'},
                'output_sha256':{p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'}}
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--query',type=Path,required=True)
    p.add_argument('--lincs',type=Path,required=True)
    p.add_argument('--fda',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--nulls',type=int,default=10000)
    p.add_argument('--seed',type=int,default=20261019)
    run(p.parse_args())
