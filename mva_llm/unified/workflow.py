#!/usr/bin/env python3
"""Single entry point: scored variants + observed HPO -> research evidence report."""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
import importlib.metadata
from pathlib import Path
from datetime import datetime, timezone

from core import VERSION, digest, write_json, read_rows, load_ontology, rank_variants, make_handoff, audit_candidate, evidence_pmids

from submission import build_submission, write_submission

HERE = Path(__file__).resolve().parent

def now(): return datetime.now(timezone.utc).isoformat()

def input_paths(config):
    paths = [config[k] for k in ('case','variants','hpo_annotations','hpo_ontology')]
    for key in ('identity_file','upstream_manifest'):
        if config.get(key): paths.append(config[key])
    if config.get('fulltext_dir'): paths.extend(str(p) for p in sorted(Path(config['fulltext_dir']).glob('*.json')))
    # Frozen local tools, source data and model config; checkpoint weights are stat-inventoried below.
    root = HERE.parent
    paths.extend(str(p) for p in sorted((root/'src').glob('*.py')))
    paths.extend(str(p) for p in sorted((root/'data/kg').glob('*.parquet')))
    paths.extend(str(p) for p in sorted(HERE.glob('*.py')))
    panel = root.parent/'gene_mol'
    paths.extend(str(p) for p in sorted((panel/'src').glob('*.py')))
    paths.extend(str(p) for p in sorted((panel/'models/serving_v2').glob('*.json')))
    return sorted(set(paths))

def checkpoint(out, name, fn, resume):
    path = out/(name+'.json')
    seal = out/(name+'.sha256')
    if resume and path.exists():
        if not seal.exists() or digest(path) != seal.read_text().strip():
            raise ValueError(f'Checkpoint integrity failure: {name}')
        return json.loads(path.read_text())
    print(f'[{now()}] {name}',flush=True)
    result = fn()
    write_json(path,result)
    seal.write_text(digest(path)+'\n')
    return result

def track2(config, handoff, out, resume):
    from local_model import LocalModel
    engine = LocalModel(config['model'],config.get('tp',4),config.get('seed',17),out,
                        config.get('fulltext_dir'),config.get('identity_file'))
    if handoff['data_class'] == 'synthetic':
        # Artificial coordinates must never be looked up as real patient alleles.
        engine.functions.pop('lookup_variant', None)
        engine.schemas = [s for s in engine.schemas if s['function']['name'] != 'lookup_variant']
    if resume and (out/'sources.json').exists():
        if not (out/'sources.sha256').exists() or digest(out/'sources.json') != (out/'sources.sha256').read_text().strip():
            raise ValueError('Source-ledger integrity failure')
        engine.sources = json.loads((out/'sources.json').read_text())
    gene = handoff['lead_gene']
    context = json.dumps(handoff)
    def discover():
        # Mandatory broad retrieval is dispatched by the driver, not left to model discretion.
        records = {}
        for name,args in [('lookup_gene',{'symbol':gene}),('gene_diseases',{'symbol':gene}),
                          ('protein_regulation',{'gene':gene})]:
            records[name] = engine.call(name,args,'discovery')
        result = engine.ask_json('discovery', '''Generate a diverse set of up to FIVE small-molecule
experimental hypotheses for the constitutional disorder. No candidate or winner is supplied.
Start from the actual lesion uncertainty: distinguish restoration/stabilisation, bypass and
supportive consequences. Exclude cytotoxicity as the constitutional endpoint. Use additional
local tools as needed, including literature and related-disease searches. Do not use a lack
of approved agents as a reason to invent evidence. Return:
{"mechanism": "established links and conditional links", "candidates":
[{"name":"exact molecule name", "rationale":"hypothesis", "target":"molecular target or process",
"source_ids":["PMIDs from retrieved evidence supporting this candidate"]}],
"gaps":["unresolved questions"]}.
A candidate set may be empty. Candidate names must come from retrieved evidence.\nCASE:\n'''+context+
            '\nMANDATORY TOOL OUTPUTS:\n'+json.dumps(records,default=str))
        result['protein_alias'] = records['protein_regulation'].get('protein_alias',gene)
        result['discovery_source_ids'] = list(engine.sources)
        return result
    discovery = checkpoint(out,'03_discovery',discover,resume)
    candidates = discovery.get('candidates')
    if not isinstance(candidates,list) or len(candidates)>5: raise ValueError('Discovery must return <=5 candidates')
    results = []
    for index, proposed in enumerate(candidates):
        if not isinstance(proposed.get('name'),str) or not proposed['name'].strip(): raise ValueError('Candidate requires name')
        name = proposed['name']
        stage = f'04_candidate_{index+1:02}'
        def evaluate():
            identity = engine.call('compound_identity', {'name':name},stage)
            safety = engine.call('drug_safety',{'drug_name':name},stage)
            panel = engine.call('safety_panel',{'smiles':identity['smiles']},stage) if identity.get('resolved') else {'ok':False,'error':'identity unresolved'}
            short = min([name,*identity.get('aliases',[])],key=len)
            protein = discovery.get('protein_alias') or gene
            queries = [f'{short} {gene} {protein} protein abundance rescue', f'{short} {gene} {protein} degradation adverse inhibition',
                       f'{short} {protein} chromosome segregation toxicity',
                       f'{name} {short} human oral trial tolerability', f'{name} {short} children chronic safety']
            searches = [engine.call('literature_evidence',{'query':q,'top':8},stage) for q in queries]
            # Automatically materialise complete abstracts, not only 420-character search snippets.
            anchors = list(proposed.get('source_ids',[])) + discovery.get('discovery_source_ids',[])
            pmids = evidence_pmids(anchors,searches,limit=config.get('source_budget',12))
            sources = [engine.call('paper_text',{'pmid':p},stage) for p in pmids]
            record = engine.ask_json(stage, '''Evaluate ONLY this exact candidate, including evidence against it.
Nominate means prioritise a controlled PRECLINICAL EXPERIMENT, not prescribe or establish
clinical suitability. A lack of human efficacy does not by itself refute a bounded experiment.
Return {"name":"exact proposed name", "nominate":true or false,
"rationale":"reason", "functional_rescue":"unestablished" or "supported_in_stated_model",
"evidence":[{"compound":"actual tested molecule, or none for a genetic manipulation",
"evidence_role":"direct_support" or "direct_opposing" or "context",
"intervention":"actual experimental manipulation",
"system":"organism/tissue/population", "dose":"value or not reported", "route":"value or not reported",
"duration":"value or not reported", "endpoint":"measured endpoint", "direction":"actual measured direction",
"source_id":"PMID", "quote":"verbatim contiguous source passage at least 30 characters",
"interpretation":"what this supports, contradicts, or leaves unknown"}],
"uncertainties":["..."], "experiment":{"design":"comparators and dose range",
"success":"functional and viability criteria", "failure":"what undermines the mechanism",
"dependency_control":"what loss of effect after target disruption would mean",
"mechanism_control":{"perturbation":"block_proposed_mechanism",
"outcome":"loss_of_benefit", "interpretation":"supports_dependency_if_controls_valid"}}}.
Use context for other molecules, genetic experiments, or background mechanisms; they
must not be labelled as direct candidate exposure. At least one actual exposure record
is required for a research nomination. The hypothetical structured control describes
how an interpretable loss of benefit would be interpreted, not an observed result.
Keep drug and genetic interventions distinct. Quote passages containing actual relevant results.
Do not claim clinical suitability. Do not treat model safety scores as safety evidence.
An empty evidence list is acceptable when evidence is absent.\nCASE:\n'''+context+
                '\nCANDIDATE:\n'+json.dumps(proposed)+'\nIDENTITY:\n'+json.dumps(identity)+
                '\nSAFETY METADATA AND PREDICTIONS:\n'+json.dumps([safety,panel],default=str)+
                '\nSEARCHES:\n'+json.dumps(searches,default=str)+'\nSOURCES:\n'+json.dumps(sources,default=str))
            # The driver owns identity and search provenance, not model-generated JSON.
            if record.get('name') != name: raise ValueError('Model changed candidate identity')
            record['identity'] = identity
            record['opposing_searches'] = queries[1:3]
            record['opposing_search_results'] = searches[1:3]
            record['search_results'] = searches
            record['safety_metadata'] = safety
            record['safety_predictions'] = panel
            def review_record(proposal, suffix):
                relevant = {str(e.get('source_id')):engine.sources.get(str(e.get('source_id'))) for e in proposal.get('evidence',[])}
                return engine.ask_json(stage+suffix, '''Independently audit the candidate against the
source text below. Do NOT defer to the proposing model. Verify that cited intervention,
molecule, dose, route, organism, endpoint and direction actually support each interpretation.
Context records may describe other molecules or genetic manipulations if clearly labelled
and not transferred as direct candidate evidence. Direct records must test this molecule.
Check that protein abundance is not equated to function, adult exposure is not paediatric
safety, absent metadata is not no human use, and genetic manipulation is not drug treatment.
Audit the proposed experiment: loss of benefit after blocking a mechanism can support that
mechanism; it does not itself refute the original benefit. Flag all contradictions.
Return {"verdict":"supported" or "unsupported" or "uncertain", "issues":["specific issue"],
"rationale":"source-linked explanation"}. Supported means a bounded experimental hypothesis,
never clinical suitability. If any material assertion is unsupported, do not return supported.
For each issue identify the actual offending claim in the proposal. Do not attribute a
claim the proposal does not make, or treat an explicitly stated limitation as a claim of
efficacy. Distinguish a proposed experimental design from a reported experimental result.
Clinical uncertainty does not by itself invalidate a conditional preclinical hypothesis.
\nCASE:\n'''+context+'\nPROPOSAL:\n'+json.dumps(proposal,default=str)+'\nSOURCES:\n'+json.dumps(relevant,default=str),tools=False)
            review = review_record(record,'_independent_review')
            gate = audit_candidate(record,engine.sources,review)
            initial = {'candidate':record,'review':review,'gate':gate}
            if gate['errors']:
                scientific = ('name','nominate','rationale','functional_rescue','evidence','uncertainties','experiment')
                driver = {k:v for k,v in record.items() if k not in scientific}
                proposal = {k:record[k] for k in scientific if k in record}
                corrected = engine.ask_json(stage+'_correction', '''Correct this research proposal against the exact
source texts. Reviewer comments may themselves be wrong: verify them, do not obey them blindly.
Return the same scientific JSON schema as the proposal below. Keep the exact candidate name.
Copy quotes verbatim from either the supplied full text or its explicitly supplied abstract.
Correct or withdraw unsupported claims. Keep other compounds and genetic interventions in
context records. Protein abundance alone is not functional rescue. Do not extrapolate an
unrelated tumour experiment to established risk in this patient. Make proposed assays clearly
hypothetical, with viability, proliferation and functional endpoints; do not invent clinical
doses from mouse exposures. Nominate means worth a controlled preclinical test, never treatment.
Either retain or withdraw nomination according to the evidence; no preferred answer is required.
If uncertainty cannot be resolved, preserve it. This is the only scientific correction attempt.
\nCASE:\n'''+context+'\nPROPOSAL:\n'+json.dumps(proposal)+
                    '\nAUTOMATIC FLAGS:\n'+json.dumps(gate)+'\nREVIEW:\n'+json.dumps(review)+
                    '\nSOURCE RECORDS:\n'+json.dumps(sources,default=str),tools=False)
                if corrected.get('name') != name: raise ValueError('Correction changed candidate identity')
                record = {**corrected,**driver}
                review = review_record(record,'_correction_review')
                gate = audit_candidate(record,engine.sources,review)
            return {'candidate':record,'review':review,'gate':gate,'initial_assessment':initial,
                    'correction_attempted':bool(initial['gate']['errors'])}
        results.append(checkpoint(out,stage,evaluate,resume))
    return {'discovery':discovery,'candidates':results,
            'research_leads':[r['candidate']['name'] for r in results if r['gate']['status']=='research_lead'],
            'clinical_suitability':'unestablished', 'status':'completed_with_research_leads' if any(r['gate']['status']=='research_lead' for r in results) else 'review_required'}

def report(case, ranked, handoff, result, out):
    lines = [f'# {case["case_id"]}: automatic research workflow', '',
             f'Data class: **{case["data_class"]}**. Workflow {VERSION}.', '',
             '## Track 1', '', '| Rank | Gene | Variant score | Phenotype similarity | Genotype status |',
             '|---|---|---|---|---|']
    for g in ranked[:10]: lines.append(f'| {g["rank"]} | {g["gene"]} | {g["variant_score"]:.3f} | {g["phenotype_score"]:.3f} | {g["genotype_status"]} |')
    lines += ['', 'Ranking is an uncalibrated heuristic. See `01_ranking.json` for source annotations and QC flags.',
              '', '## Handoff', '', f'Lead gene: **{handoff["lead_gene"]}**. Phase unknown; residual function unestablished.',
              '', *['- '+s for s in handoff['limitations']], '', '## Track 2', '']
    if result is None:
        lines += ['Not executed (`--track1-only`). This report is not an end-to-end result.']
    else:
        lines += [f'Status: **{result["status"]}**. Clinical suitability: **unestablished**.',
                  '', 'Research leads: '+(', '.join(result['research_leads']) or 'none passing the automated audit')+'.']
        for r in result['candidates']:
            c = r['candidate']; lines += ['', '### '+c['name'], '', 'Gate: **'+r['gate']['status']+'**.', '', c.get('rationale','')]
            if r.get('correction_attempted'):
                lines += ['', 'One automatic source-based correction and a fresh review were performed; the initial assessment is preserved in the candidate JSON.']
            if r['gate']['errors']: lines += ['', 'Audit flags: '+', '.join(r['gate']['errors'])]
            lines += ['', 'Independent review: '+r['review'].get('rationale','')]
            for e in c.get('evidence',[]):
                sid = str(e.get('source_id',''))
                lines += ['', f'- {e.get("interpretation", "")} ([PMID {sid}](https://pubmed.ncbi.nlm.nih.gov/{sid}/)); '+
                          f'intervention: {e.get("intervention")}; system: {e.get("system")}; endpoint: {e.get("endpoint")}.']
            lines += ['', 'Uncertainties: '+'; '.join(c.get('uncertainties',[])), '', 'Proposed experiment:', '',
                      '```json',json.dumps(c.get('experiment',{}),indent=2),'```']
    lines += ['', '## Submission export', '',
              'Automatically generated: `track1_submission.csv`; calculation audit: `02_submission.json`. '
              'EPCR is the gene ranking score divided by the within-gene hypothesis rank. '
              'These are uncalibrated ranking estimates, not validated causal probabilities. '
              'Up to 10 hypotheses are exported as primary candidates; secondary findings are not inferred.',
              '', '## Limits', '', 'This is a research workflow, not a validated diagnostic or treatment system. '
              'Quote checks establish source membership; model review can still misinterpret science. '
              'The existing literature corpus and several legacy tools are biased toward MVA. '
              'No independent generalisation benchmark has yet been passed.']
    (out/'report.md').write_text('\n'.join(lines)+'\n')

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True); p.add_argument('--out',required=True)
    p.add_argument('--track1-only',action='store_true'); p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    config=json.loads(Path(args.config).read_text())
    case=json.loads(Path(config['case']).read_text())
    if case.get('data_class') not in ('synthetic','public','authorised_patient'):
        raise ValueError('case.data_class must be synthetic, public or authorised_patient')
    if case['data_class']=='authorised_patient' and not case.get('reuse_authorisation'):
        raise ValueError('Record the provenance of patient reuse authorisation in the case manifest')
    out=Path(args.out).resolve()
    if out.exists() and any(out.iterdir()) and not args.resume:
        raise ValueError('Output directory nonempty: choose a new directory or --resume')
    out.mkdir(parents=True,exist_ok=True)
    os.environ.update(HF_HUB_OFFLINE='1',TRANSFORMERS_OFFLINE='1',GENE_MOL_MODELS='tiered',GENE_MOL_EXT_RANK='calibrated')
    fingerprint={'config':config, 'inputs':{str(path):digest(path) for path in input_paths(config)}}
    fingerprint['runtime'] = {'python':sys.version, 'executable':sys.executable,
        'packages':dict(sorted((d.metadata['Name'],d.version) for d in importlib.metadata.distributions() if d.metadata.get('Name')))}
    panel = HERE.parent.parent/'gene_mol/models/serving_v2'
    fingerprint['panel_inventory'] = {str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
        for p in sorted(panel.rglob('*')) if p.is_file()}
    if config.get('model'):
        fingerprint['model_files']={str(p):{'size':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns}
                                    for p in sorted(Path(config['model']).rglob('*')) if p.is_file()}
    manifest=out/'manifest.json'
    if args.resume and manifest.exists():
        if json.loads(manifest.read_text())['fingerprint'] != fingerprint:
            raise ValueError('Inputs, code, configuration or model inventory changed; use a fresh run directory')
    else: write_json(manifest,{'created':now(),'fingerprint':fingerprint,'data_class':case['data_class']})
    snapshot = out/'source_snapshot'
    snapshot.mkdir(exist_ok=True)
    for source in HERE.glob('*.py'):
        target = snapshot/source.name
        if not target.exists(): target.write_bytes(source.read_bytes())
    try:
        ranked=checkpoint(out,'01_ranking',lambda:rank_variants(read_rows(config['variants']),case['observed_hpo'],
            read_rows(config['hpo_annotations']),load_ontology(config['hpo_ontology'])),args.resume)
        handoff=checkpoint(out,'02_handoff',lambda:make_handoff(case,ranked),args.resume)
        submission=checkpoint(out,'02_submission',lambda:build_submission(case,ranked,config.get('submission')),args.resume)
        write_submission(submission,out/'track1_submission.csv')
        (out/'track1_submission.sha256').write_text(digest(out/'track1_submission.csv')+'\n')
        result=None
        if not args.track1_only:
            result=checkpoint(out,'05_track2',lambda:track2(config,handoff,out,args.resume),args.resume)
        report(case,ranked,handoff,result,out)
        write_json(out/'status.json',{'status':'track1_complete' if result is None else result['status'],'finished':now()})
        print(f'Wrote {out}/report.md',flush=True)
    except Exception as exc:
        write_json(out/'status.json',{'status':'failed','error':f'{type(exc).__name__}: {exc}','time':now()})
        raise

if __name__=='__main__': main()
