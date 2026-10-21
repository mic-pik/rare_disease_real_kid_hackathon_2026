#!/usr/bin/env python3
"""Conservative assay triage plus explicit, source-linked curator judgments.

The ordered output is a research review queue, not a drug efficacy prediction.
No language model calls or undocumented clinical-exposure imputations.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re

import pandas as pd

from fda_repurposing import sha256, write_json

POLICY = {
    'version': '1.0', 'active_potency_nM': 1000, 'weak_potency_nM': 10000,
    'screen_active_inhibition_percent': 50, 'screen_weak_inhibition_percent': 20,
    'screen_active_max_test_nM': 1000, 'screen_weak_min_test_nM': 10000,
    'meaning': 'Operational retrieval thresholds, not universal pharmacological or efficacy cutoffs',
}
UNITS = {'pm': .001, 'nm': 1, 'um': 1000, 'mm': 1000000, 'm': 1e9}


def text(x):
    return '' if x is None or pd.isna(x) else str(x).strip()


def number(x):
    try:
        n = float(x)
        return n if math.isfinite(n) else None
    except (ValueError, TypeError):
        return None


def direction(row):
    action = text(row.get('reported_action_type')).upper()
    if action in ('INHIBITOR', 'ANTAGONIST', 'BLOCKER', 'NEGATIVE ALLOSTERIC MODULATOR'):
        return 'inhibition', 'structured_activity_annotation'
    if action in ('ACTIVATOR', 'AGONIST', 'POSITIVE ALLOSTERIC MODULATOR'):
        return 'activation', 'structured_activity_annotation'
    # Binding constants do not establish functional direction, even when a
    # generic assay description says "inhibition".
    if text(row.get('standard_type')) in ('Kd','Ki'):
        return 'unknown', 'binding_measurement'
    description = text(row.get('assay_description')).lower()
    if re.search(r'phosphorylation|acetylation|cell (?:growth|proliferation)|mib binding', description):
        return 'unknown', 'indirect_endpoint_requires_review'
    if text(row.get('standard_type')) == 'IC50' and description.startswith('inhibition of '):
        return 'inhibition', 'assay_description_inference_not_independent_validation'
    return 'unknown', 'not_established'


def potency_class(value, relation):
    if relation in ('=', '<', '<=') and value <= POLICY['active_potency_nM']:
        return 'active', 'potency_bound_at_or_below_screening_threshold'
    if relation in ('=', '>', '>=') and value >= POLICY['weak_potency_nM']:
        return 'inactive_at_threshold', 'weak_potency_bound_not_universal_inactivity'
    return 'ambiguous', 'intermediate_potency_or_bound_crosses_threshold'


def assay_triage(row):
    result = {'activity_class': 'uninterpretable', 'classification_reason': '',
              'value_nM': None, 'test_concentration_nM': None,
              'inhibition_percent': None, 'quality_flags': [],
              'eligible_for_support_count': False}
    result['inferred_direction'], result['direction_basis'] = direction(row)
    if text(row.get('data_validity_comment')):
        result['quality_flags'].append('data_validity:' + text(row['data_validity_comment']))
    if number(row.get('potential_duplicate')) not in (None, 0):
        result['quality_flags'].append('potential_duplicate')
    value = number(row.get('standard_value'))
    relation = text(row.get('standard_relation'))
    typ = text(row.get('standard_type'))
    units = text(row.get('standard_units')).replace('µ','u').replace('μ','u').lower()
    desc = text(row.get('assay_description')).lower().replace('µ','u').replace('μ','u')
    if value is None or relation not in ('=', '<', '<=', '>', '>='):
        result['classification_reason'] = 'missing_numeric_value_or_supported_relation'
    elif typ in ('IC50','EC50','AC50','Ki','Kd','Potency'):
        if value <= 0 or units not in UNITS:
            result['classification_reason'] = 'nonpositive_potency_or_unknown_units'
        else:
            value *= UNITS[units]
            result['value_nM'] = value
            result['activity_class'], result['classification_reason'] = potency_class(value, relation)
    elif units == '%' and typ in ('Inhibition','% of inhibition','Residual Activity','% Control','% of control','Activity'):
        # Do not read phosphorylation or generic "Activity" as inhibition.
        if any(w in desc for w in ('phosphorylation','acetylation','cell viability','cell growth')):
            result['classification_reason'] = 'indirect_cellular_endpoint'
        else:
            residual = typ in ('Residual Activity','% Control','% of control') or bool(re.search(r'(activity remaining|residual activity|residual kinase activity)', desc))
            inhibitory = typ in ('Inhibition','% of inhibition')
            doses = re.findall(r'\bat\s+([0-9]+(?:\.[0-9]+)?)\s*(nm|um|mm)\b', desc)
            if len(doses) != 1 or not (residual or inhibitory) or not 0 <= value <= 100:
                result['classification_reason'] = 'ambiguous_percent_polarity_dose_or_range'
            else:
                dose = float(doses[0][0]) * UNITS[doses[0][1]]
                effect = 100-value if residual else value
                rel = {'>':'<', '>=':'<=', '<':'>', '<=':'>=', '=':'='}[relation] if residual else relation
                result.update(test_concentration_nM=dose, inhibition_percent=effect)
                if dose > 0 and dose <= POLICY['screen_active_max_test_nM'] and effect >= 50 and rel in ('=','>','>='):
                    result.update(activity_class='active', classification_reason='single_point_effect_at_low_test_concentration')
                elif dose >= POLICY['screen_weak_min_test_nM'] and effect <= 20 and rel in ('=','<','<='):
                    result.update(activity_class='inactive_at_threshold', classification_reason='small_effect_at_high_test_concentration')
                else:
                    result.update(activity_class='ambiguous', classification_reason='single_point_result_not_decisive_at_policy_thresholds')
    else:
        result['classification_reason'] = 'endpoint_requires_assay_specific_manual_interpretation'
    comment = text(row.get('activity_comment')).lower()
    if result['activity_class'] == 'active' and comment in ('inactive','not active','inconclusive','uncertain'):
        result['quality_flags'].append('numeric_activity_conflicts_with_source_comment')
    if result['activity_class'] == 'inactive_at_threshold' and comment == 'active':
        result['quality_flags'].append('numeric_inactivity_conflicts_with_source_comment')
    result['eligible_for_support_count'] = not result['quality_flags'] and result['activity_class'] in ('active','inactive_at_threshold')
    result['quality_flags'] = ';'.join(result['quality_flags'])
    return result


def aggregate(assays):
    rows = []
    for (ingredient, gene), group in assays.groupby(['ingredient_id','gene'], sort=True):
        good = group[group.eligible_for_support_count]
        active = good[good.activity_class == 'active']
        weak = good[good.activity_class == 'inactive_at_threshold']
        # Distinct documents are not necessarily independent experiments.
        rows.append({'ingredient_id':ingredient, 'gene':gene, 'n_records':len(group),
                     'n_active_unflagged':len(active), 'n_weak_unflagged':len(weak),
                     'mixed_active_weak_evidence':bool(len(active) and len(weak)),
                     'active_document_count':active.document_chembl_id.dropna().nunique(),
                     'has_inhibitory_active':bool((active.inferred_direction == 'inhibition').any()),
                     'has_binding_active':bool(active.standard_type.isin(['Kd','Ki']).any()),
                     'has_activating_active':bool((active.inferred_direction == 'activation').any())})
    return pd.DataFrame(rows, columns=['ingredient_id','gene','n_records','n_active_unflagged',
         'n_weak_unflagged','mixed_active_weak_evidence','active_document_count',
         'has_inhibitory_active','has_binding_active','has_activating_active'])


def default_decision(summary, intents):
    if 'instability_directed' in intents:
        return 2, 'defer_primary_evidence_review', 'Related-disease intent is metadata; efficacy and transfer need primary-source review.'
    if (summary.n_active_unflagged > 0).any():
        if summary.has_inhibitory_active.any():
            return 3, 'defer_direction_conflict', 'Inhibitory activity in the selected checkpoint/regulation panel does not establish restoration; no functional rescue evidence supplied.'
        return 3, 'defer_direction_unknown', 'Assay activity/binding found; beneficial action direction and disease rescue remain unestablished.'
    if intents:
        return 4, 'reject_retrieval_rationale', 'Retrieved sibling uses address consequences or symptoms; they do not support the constitutional segregation-rescue objective.'
    if len(summary) and (summary.n_weak_unflagged > 0).any():
        return 5, 'deprioritize_no_active_support', 'No unflagged active evidence at the operational thresholds; weak measurements are not positive repurposing evidence.'
    return 5, 'defer_uninterpretable_assays', 'No decisive unflagged assay evidence; generic fold-change, missing values or ambiguous endpoints require assay-specific interpretation.'


def validate_review(review, ids):
    sources = review['sources']
    for cid, item in review['candidate_reviews'].items():
        if cid not in ids:
            raise ValueError('Reviewed identity absent from retrieval: ' + cid)
        if not item['source_ids'] or not set(item['source_ids']) <= set(sources):
            raise ValueError('Missing evidence source for ' + cid)
        if item['decision'] not in ('conditional_experiment_priority','defer_transfer','reject_retrieval_rationale'):
            raise ValueError('Unsupported curator decision')
        if item['decision'] == 'conditional_experiment_priority' and not item.get('decisive_experiment'):
            raise ValueError('Conditional priority requires a falsifiable experiment')


def run(input_dir, review_path, out):
    input_dir, review_path, out = Path(input_dir), Path(review_path), Path(out)
    if out.exists():
        raise ValueError('Output directory must be fresh')
    manifest = json.loads((input_dir/'manifest.json').read_text())
    # Verify before interpreting any upstream result.
    for name, expected in manifest['output_hashes'].items():
        if sha256(input_dir/name) != expected:
            raise ValueError('Upstream hash mismatch: '+name)
    candidates = json.loads((input_dir/'candidates.json').read_text())
    review = json.loads(review_path.read_text())
    validate_review(review, {c['identity']['ingredient_id'] for c in candidates})
    measured = pd.read_parquet(input_dir/'measured_interactions.parquet')
    siblings = pd.read_parquet(input_dir/'related_disease_candidates.parquet')
    assays = pd.concat([measured.reset_index(drop=True), pd.DataFrame(
        [assay_triage(row) for row in measured.to_dict('records')])], axis=1)
    pairs = aggregate(assays)
    rows, detailed = [], []
    for candidate in candidates:
        identity = candidate['identity']; cid = identity['ingredient_id']
        summary = pairs[pairs.ingredient_id == cid]
        s = siblings[siblings.ingredient_id == cid]
        intents = set(s.intent.dropna())
        tier, decision, reason = default_decision(summary, intents)
        manual = review['candidate_reviews'].get(cid)
        if manual:
            tier, decision, reason = manual['tier'], manual['decision'], manual['reason']
        row = {'ingredient_id':cid, 'name':identity['canonical_name'], 'review_tier':tier,
               'decision':decision, 'reason':reason,
               'review_depth':'source_level_curated' if manual else 'assay_and_retrieval_triage',
               'single_ingredient_product':bool(identity['eligible_single_ingredient_product_ids']),
               'n_active_records':int(summary.n_active_unflagged.sum()),
               'n_weak_records':int(summary.n_weak_unflagged.sum()),
               'mixed_target_evidence':';'.join(summary.loc[summary.mixed_active_weak_evidence,'gene']),
               'sibling_intents':';'.join(sorted(intents)),
               'exposure_status':manual['exposure_status'] if manual else 'not_assessed',
               'safety_status':manual['safety_status'] if manual else 'not_assessed',
               'clinical_suitability':'unestablished'}
        rows.append(row)
        detailed.append({'review':row, 'original_candidate':candidate,
                         'curator_dossier':manual, 'functional_rescue_in_MVA':'unestablished'})
    frame = pd.DataFrame(rows).sort_values(['review_tier','name']).reset_index(drop=True)
    # Ordinal is presentation order; ties within tiers are alphabetical.
    frame.insert(0,'review_order',range(1,len(frame)+1))
    out.mkdir(parents=True)
    assays.to_parquet(out/'assay_triage.parquet',index=False)
    pairs.to_csv(out/'target_evidence.csv',index=False)
    frame.to_csv(out/'review_queue.csv',index=False)
    frame[frame.decision=='conditional_experiment_priority'].to_csv(out/'experimental_shortlist.csv',index=False)
    write_json(out/'candidate_reviews.json',detailed)
    write_json(out/'review_evidence.json',review)
    write_json(out/'policy.json',POLICY)
    dossiers = out/'dossiers'; dossiers.mkdir()
    for cid, item in review['candidate_reviews'].items():
        name = frame.loc[frame.ingredient_id==cid,'name'].iloc[0]
        lines = [f'# {name}', '', f"Decision: **{item['decision']}**. Clinical suitability: unestablished.",
                 '', item['reason'], '']
        for label,key in [('Supporting evidence','support'),('Counter-evidence and limitations','counter'),
                          ('Exposure feasibility','exposure'),('Safety','safety'),
                          ('Decisive experiment','decisive_experiment')]:
            lines += ['## '+label, '', item[key], '']
        lines += ['## Sources', '']
        for source_id in item['source_ids']:
            source = review['sources'][source_id]
            link = source.get('url') or os.path.relpath(
                Path(__file__).resolve().parents[1]/source['path'], dossiers.resolve())
            lines.append(f"- [{source['title']}]({link})")
        (dossiers/(cid+'.md')).write_text('\n'.join(lines)+'\n')
    summary = {'candidate_ingredients':len(frame), 'assay_records':len(assays),
               'assay_classes':dict(Counter(assays.activity_class)),
               'quality_flagged_records':int(assays.quality_flags.ne('').sum()),
               'decisions':dict(Counter(frame.decision)),
               'source_level_dossiers':len(review['candidate_reviews']),
               'clinical_leads':0}
    write_json(out/'summary.json',summary)
    lines = ['# FDA candidate evidence review', '',
        'The queue prioritises further research; it is not an efficacy ranking or treatment recommendation.', '',
        f"Reviewed {len(frame)} ingredients and {len(assays)} assay rows. {len(review['candidate_reviews'])} candidates received source-level dossiers; others received explicit assay/route triage.", '',
        '## Decisions', '', '| Decision | Ingredients |','|---|---:|']
    lines += [f'| {k} | {v} |' for k,v in summary['decisions'].items()]
    lines += ['', '## Curated review queue', '', '| Ingredient | Decision | Reason |','|---|---|---|']
    for r in frame[frame.review_depth=='source_level_curated'].itertuples():
        lines += [f'| [{r.name}](dossiers/{r.ingredient_id}.md) | {r.decision} | {r.reason} |']
    lines += ['', '## Interpretation limits', '',
        '- Active means an operational assay threshold was met. It does not mean rescue, safety, or clinical exposure is established.',
        '- Weak/inactive-at-threshold results apply to the measured target and conditions, not every possible mechanism.',
        '- Mixed active/weak evidence is retained; different assay contexts can explain disagreement.',
        '- No unreviewed compound receives an exposure margin or safety endorsement.',
        '- The shortlist is deliberately not padded to five compounds when the evidence does not justify five.',
        '- Curator judgments are explicit versioned inputs, not independently validated model predictions.',
        '- FDA approval comes from the upstream frozen catalogue; this review does not independently recheck all withdrawals.', '']
    (out/'report.md').write_text('\n'.join(lines))
    write_json(out/'manifest.json',{'created_utc':datetime.now(timezone.utc).isoformat(),
        'upstream_manifest_sha256':sha256(input_dir/'manifest.json'),
        'review_evidence_sha256':sha256(review_path), 'implementation_sha256':sha256(__file__),
        'output_hashes':{str(p.relative_to(out)):sha256(p) for p in sorted(out.rglob('*')) if p.is_file()}})
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True); p.add_argument('--review',required=True); p.add_argument('--out',required=True)
    a = p.parse_args(); run(a.input,a.review,a.out)
