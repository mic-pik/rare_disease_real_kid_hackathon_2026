#!/usr/bin/env python3
"""Fixed-candidate cycle-gene ablation and reagent-aware BUB1B checks."""
import argparse
import json
import re
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
try:
    from .score_gse22206_lincs import weights, bh, digest
except ImportError:
    from score_gse22206_lincs import weights, bh, digest


def rank_unit(x):
    r = pd.DataFrame(x).rank(axis=1, method='average').to_numpy(copy=True)
    r -= r.mean(axis=1, keepdims=True)
    norm = np.linalg.norm(r, axis=1, keepdims=True)
    return np.divide(r, norm, out=np.zeros_like(r), where=norm > 0)


def reagent_agreement(x):
    if len(x) < 2:
        return np.nan
    u = rank_unit(x)
    return float(np.median((u @ u.T)[np.triu_indices(len(u), 1)]))


def aggregate(frame, values):
    keys=['ingredient_id','canonical_name','cell_id','pert_time','pert_time_unit','pert_dose','pert_dose_unit']
    joined=pd.concat([frame[keys].reset_index(drop=True), values.reset_index(drop=True)],axis=1)
    conditions=joined.groupby(keys,dropna=False)[list(values.columns)].mean().reset_index()
    cell=conditions.groupby(keys[:3])[list(values.columns)].mean()
    return cell.groupby(level=[0,1]).mean(), cell, conditions


def run(a):
    o=a.output; o.mkdir(parents=True,exist_ok=True)
    manifest=json.loads((a.baseline/'manifest.json').read_text())
    for name in ['selected_landmark_matrix.npz','signature_scores.csv','compound_scores.csv']:
        assert digest(a.baseline/name)==manifest['output_sha256'][name]
    baseline=pd.read_csv(a.baseline/'compound_scores.csv')
    s=pd.read_csv(a.baseline/'signature_scores.csv')
    archive=np.load(a.baseline/'selected_landmark_matrix.npz')
    assert np.array_equal(archive['sig_ids'],s.sig_id.to_numpy())
    z=archive['z']; ids=archive['gene_ids'].tolist()
    q=pd.read_csv(a.query/'lincs_measured_landmark_query.csv',dtype={'entrez_id':str})
    gp=pd.read_parquet(a.kg/'gene_pathway.parquet'); genes=pd.read_parquet(a.kg/'genes.parquet')
    banned=gp[gp.top_level.isin(['Cell Cycle','DNA Replication'])].merge(genes[['gene_id','symbol']],on='gene_id')
    banned.to_csv(o/'excluded_pathway_memberships.csv',index=False)
    info=pd.read_csv(a.lincs/'GSE92742_Broad_LINCS_gene_info.txt.gz',sep='\t',dtype={'pr_gene_id':str})
    excluded_ids=set(info.loc[info.pr_gene_symbol.isin(set(banned.symbol)), 'pr_gene_id'])
    retained=np.array([g not in excluded_ids for g in ids])
    q['removed_cycle']=q.entrez_id.isin(excluded_ids)
    q.to_csv(o/'query_ablation.csv',index=False)
    q=q[~q.removed_cycle]
    # Preserve original within-profile rank scale for comparison.
    r=2*(pd.DataFrame(z).rank(axis=1).to_numpy()-1)/(len(ids)-1)-1
    w=np.zeros(len(ids)); w[[ids.index(g) for g in q.entrez_id]]=weights(q.robust_direction)
    agg,cell,conditions=aggregate(s,pd.DataFrame(r,columns=ids))
    result=baseline.set_index(['ingredient_id','canonical_name']).copy()
    result['noncycle_score']=pd.Series(agg.to_numpy()@w,index=agg.index)
    rng=np.random.default_rng(20261019)
    nullw=np.zeros((len(ids),10000)); nup=int((q.robust_direction==1).sum()); ndown=len(q)-nup
    for i in range(10000):
        draw=rng.choice(np.flatnonzero(retained),len(q),replace=False)
        nullw[draw[:nup],i]=-.5/nup;nullw[draw[nup:],i]=.5/ndown
    null=agg.to_numpy()@nullw; observed=agg.to_numpy()@w
    p=(1+(null>=observed[:,None]).sum(axis=1))/10001
    result['noncycle_p']=pd.Series(p,index=agg.index)
    result['noncycle_q']=pd.Series(bh(p),index=agg.index)
    cv=pd.Series(cell.to_numpy()@w,index=cell.index)
    result['noncycle_positive_cell_fraction']=(cv>0).groupby(level=[0,1]).mean()
    result['noncycle_screen_supported']=(result.noncycle_q<=.1)&(result.noncycle_positive_cell_fraction>=.75)&(result.n_cells>=2)
    s['noncycle_score']=r@w
    s.to_csv(o/'drug_profile_scores.csv',index=False)
    # Individual shRNAs, never mixed with consensus or overexpression records.
    si=pd.read_csv(a.lincs/'GSE92742_Broad_LINCS_sig_info.txt.gz',sep='\t',low_memory=False)
    metrics=pd.read_csv(a.lincs/'GSE92742_Broad_LINCS_sig_metrics.txt.gz',sep='\t',low_memory=False)
    kd=si[(si.pert_iname=='BUB1B')&(si.pert_type=='trt_sh')].merge(metrics[['sig_id','is_exemplar','distil_cc_q75','distil_nsample']],on='sig_id',validate='one_to_one')
    kd['quality_pass']=(kd.is_exemplar==1)&(kd.distil_cc_q75>=.2)&(kd.distil_nsample>=3)
    kd.to_csv(o/'kd_inventory.csv',index=False); kd=kd[kd.quality_pass].copy()
    used=set();chosen=[]
    for row in kd.sort_values(['distil_cc_q75','sig_id'],ascending=[False,True]).itertuples():
        constituents=set(str(row.distil_id).split('|'))
        if constituents & used or '-666' in constituents: continue
        used.update(constituents);chosen.append(row.sig_id)
    kd=kd[kd.sig_id.isin(chosen)].sort_values('sig_id').reset_index(drop=True)
    path=a.lincs/'GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx'
    decode=lambda x:[v.decode() if isinstance(v,bytes) else str(v) for v in x]
    with h5py.File(path) as f:
        ri={v:i for i,v in enumerate(decode(f['0/META/ROW/id'][:]))}
        ci={v:i for i,v in enumerate(decode(f['0/META/COL/id'][:]))}
        kz=np.array([f['0/DATA/0/matrix'][ci[g],:][[ri[v] for v in ids]] for g in kd.sig_id])
    assert np.isfinite(kz).all()
    np.savez_compressed(o/'kd_measured_matrix.npz',z=kz,sig_ids=kd.sig_id.to_numpy(dtype=str),gene_ids=np.array(ids))
    kd.to_csv(o/'kd_selected.csv',index=False)
    kdf=pd.DataFrame(kz,columns=ids)
    for key in ['cell_id','pert_time','pert_id']: kdf[key]=kd[key]
    kr=kdf.groupby(['cell_id','pert_time','pert_id'])[ids].mean()
    contexts=[]; comparisons=[]
    drug_u=rank_unit(z[:,retained])
    # Patient agreement uses only the fixed retained query, equal direction weights.
    patient=np.zeros(len(ids)); patient[[ids.index(g) for g in q.entrez_id]]=-weights(q.robust_direction)
    for (cell_id,kd_time),g in kr.groupby(level=[0,1]):
        agreement=reagent_agreement(g.to_numpy()[:,retained])
        usable=len(g)>=2 and agreement>.1
        ku=rank_unit(g.to_numpy()[:,retained])
        context={'cell_id':cell_id,'kd_time_h':kd_time,'n_reagents':len(g),'median_reagent_rank_correlation':agreement if np.isfinite(agreement) else None,'usable':usable}
        k_r=2*(pd.DataFrame(g.to_numpy()).rank(axis=1).to_numpy()-1)/(len(ids)-1)-1
        context['patient_direction_agreement']=float(np.mean(k_r@patient))
        contexts.append(context)
        if not usable: continue
        for i in s.index[s.cell_id==cell_id]:
            scores=-(ku@drug_u[i])
            for reagent,score in zip(g.index.get_level_values('pert_id'),scores):
                comparisons.append({'sig_id':s.loc[i,'sig_id'],'ingredient_id':s.loc[i,'ingredient_id'],'canonical_name':s.loc[i,'canonical_name'],
                                    'cell_id':cell_id,'drug_time_h':s.loc[i,'pert_time'],'drug_dose_uM':s.loc[i,'pert_dose'],
                                    'kd_time_h':kd_time,'reagent':reagent,'kd_reversal':score})
    pd.DataFrame(contexts).to_csv(o/'kd_context_agreement.csv',index=False)
    comp=pd.DataFrame(comparisons);comp.to_csv(o/'drug_vs_individual_shrna.csv',index=False)
    if not comp.empty:
        # Avoid pooling doses, genes or repeated profiles as independent evidence.
        prof=comp.groupby(['ingredient_id','canonical_name','cell_id','kd_time_h','drug_time_h','drug_dose_uM']).agg(kd_reversal=('kd_reversal','mean'),reagent_positive_fraction=('kd_reversal',lambda x:float((x>0).mean())))
        kc=prof.groupby(level=[0,1,2])[['kd_reversal','reagent_positive_fraction']].mean()
        ks=kc.groupby(level=[0,1]).mean();ks['kd_covered_cells']=kc.groupby(level=[0,1]).size()
        ks['kd_positive_cell_fraction']=(kc.kd_reversal>0).groupby(level=[0,1]).mean()
        result=result.join(ks)
        result['kd_direction_supported']=(result.kd_covered_cells>=2)&(result.kd_positive_cell_fraction>=.75)&(result.reagent_positive_fraction>=.75)&(result.kd_reversal>0)
        prof.to_csv(o/'kd_condition_scores.csv')
    result.to_csv(o/'all_compound_checks.csv')
    result[result.posthoc_background_supported].to_csv(o/'fixed20_checks.csv')
    # PRISM is contextual viability evidence, never a matched rescue/safety gate.
    if a.prism:
        cells=pd.read_csv(a.prism/'prism_cells.csv')
        normalize=lambda x: re.sub('[^A-Z0-9]','',x.upper())
        cells['cell_key']=cells.ccle_name.str.split('_').str[0].map(normalize)
        cellmap=cells.groupby('cell_key').depmap_id.agg(lambda x:sorted(set(x.dropna())))
        cellmap=cellmap[cellmap.map(len)==1].map(lambda x:x[0])
        wanted=s[s.ingredient_id.isin(result[result.posthoc_background_supported].index.get_level_values(0))].copy()
        wanted['depmap_id']=wanted.cell_id.map(normalize).map(cellmap)
        cmp=pd.read_csv(a.prism/'prism_compounds.csv').dropna(subset=['IDs'])
        cmp['pert_id']=cmp.IDs.str.extract(r'(BRD-[A-Z0-9]+)',expand=False)
        prism=pd.read_csv(a.prism/'prism_matrix.csv',index_col=0,low_memory=False)
        records=[]
        for row in wanted.itertuples():
            if pd.isna(row.depmap_id) or row.depmap_id not in prism.columns: continue
            for match in cmp[cmp.pert_id==row.pert_id].drop_duplicates(['IDs','dose','screen']).itertuples():
                if match.IDs not in prism.index: continue
                val=prism.at[match.IDs,row.depmap_id]
                if pd.isna(val):continue
                records.append({'ingredient_id':row.ingredient_id,'canonical_name':row.canonical_name,'sig_id':row.sig_id,
                    'cell_id':row.cell_id,'depmap_id':row.depmap_id,'prism_id':match.IDs,'prism_screen':match.screen,
                    'lincs_dose_uM':row.pert_dose,'lincs_time_h':row.pert_time,'prism_dose_metadata':match.dose,
                    'prism_log2_viability':float(val),'below_half_relative_viability':bool(val < -1),
                    'dose_numeric_equal':bool(np.isclose(row.pert_dose,match.dose)),
                    'context':'PRISM aggregated profile; exact Broad compound root, normalized unique cell-name mapping; not a matched-time assay'})
        pv=pd.DataFrame(records);pv.to_csv(o/'prism_context.csv',index=False)
        if not pv.empty:
            pv.drop_duplicates(['ingredient_id','depmap_id','prism_id']).groupby(['ingredient_id','canonical_name']).agg(
                prism_cells=('depmap_id','nunique'),median_prism_log2_viability=('prism_log2_viability','median'),
                fraction_below_half=('below_half_relative_viability','mean')).to_csv(o/'prism_summary.csv')
    summary={'original_query_genes':74,'removed_query_genes':74-len(q),'remaining_up':nup,'remaining_down':ndown,'noncycle_universe_genes':int(retained.sum()),
             'kd_individual_quality_profiles':len(kd),'kd_contexts':contexts,
             'fixed20_noncycle_supported':int(result[result.posthoc_background_supported].noncycle_screen_supported.sum()),
             'fixed20_kd_supported':int(result[result.posthoc_background_supported].kd_direction_supported.sum()),
             'fixed20_both_supported':int((result.posthoc_background_supported&result.noncycle_screen_supported&result.kd_direction_supported).sum()),
             'limits':['No observed growth adjustment','shRNA off-targets unresolved','KD96h versus drugs6/24h are different states','same LINCS platform, not independent disease cohort','cycle exclusion also removes real BUB1B biology']}
    (o/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    inputs=[a.baseline/'manifest.json',a.query/'lincs_measured_landmark_query.csv',a.kg/'gene_pathway.parquet',a.kg/'genes.parquet']+[a.lincs/f'GSE92742_Broad_LINCS_{name}.txt.gz' for name in ['gene_info','sig_info','sig_metrics']]
    if a.prism: inputs += [a.prism/name for name in ['prism_cells.csv','prism_compounds.csv','prism_matrix.csv']]
    (o/'manifest.json').write_text(json.dumps({'inputs_sha256':{p.name:digest(p) for p in inputs},'script_sha256':digest(Path(__file__)),
        'outputs_sha256':{p.name:digest(p) for p in o.iterdir() if p.is_file() and p.name!='manifest.json'}},indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['baseline','query','lincs','kg','output']:p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--prism',type=Path)
    run(p.parse_args())
