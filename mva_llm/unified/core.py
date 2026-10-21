"""Pure-Python contracts, phenotype ranking and conservative research gates."""
from __future__ import annotations
import csv
import gzip
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path

VERSION = '0.5.0'

def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()

def write_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)

def read_rows(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def number(value, default=None):
    if value in (None, '', '.', '-', 'nan', 'NA'): return default
    x = float(value)
    if not math.isfinite(x): raise ValueError('Non-finite numeric input')
    return x

def load_ontology(path):
    parents, names, aliases = {}, {}, {}
    current = None
    for line in Path(path).read_text().splitlines():
        if line == '[Term]': current = None
        elif line.startswith('id: HP:'):
            current = line[4:]; parents[current] = set()
        elif current and line.startswith('name: '): names[current] = line[6:]
        elif current and line.startswith('is_a: '): parents[current].add(line.split()[1])
        elif current and line.startswith('alt_id: '): aliases[line[8:]] = current
        elif line.startswith('['): current = None
    cache = {}
    def ancestors(term, seen=None):
        term = aliases.get(term, term)
        if term in cache: return cache[term]
        seen = set() if seen is None else seen
        if term in seen: raise ValueError('Cycle in HPO ontology')
        result = {term}
        for p in parents.get(term, ()):
            result |= ancestors(p, seen | {term})
        cache[term] = result
        return result
    for term in parents: ancestors(term)
    return cache, names, aliases

def rank_variants(rows, phenotype_ids, annotations, ontology):
    """Transparent heuristic, not a calibrated disease probability.

    Resnik-style information-content similarity, normalised per query term.
    No diagnostic labels, disease name, expected gene or compound enter ranking.
    """
    ancestors, names, aliases = ontology
    query = list(dict.fromkeys(aliases.get(t, t) for t in phenotype_ids))
    if not query or any(t not in ancestors for t in query):
        raise ValueError('Supply non-empty observed HPO IDs present in the ontology')
    gene_terms = defaultdict(set)
    inheritance = defaultdict(set)
    for a in annotations:
        term = aliases.get(a['hpo_id'], a['hpo_id'])
        if term in {'HP:0000006', 'HP:0000007', 'HP:0001417', 'HP:0001419'}:
            inheritance[a['gene_symbol']].add(term)
        # Exclude inheritance/clinical modifier terms; only phenotype root descendants.
        if 'HP:0000118' in ancestors.get(term, set()):
            gene_terms[a['gene_symbol']].update(ancestors[term])
    counts = defaultdict(int)
    for terms in gene_terms.values():
        for term in terms: counts[term] += 1
    n = max(1, len(gene_terms))
    ic = {t: -math.log((v + 1) / (n + 1)) for t, v in counts.items()}
    groups = defaultdict(list)
    seen = set()
    for row in rows:
        vid = row.get('var_key') or row.get('variant')
        gene = row.get('genes') or row.get('gene')
        score = number(row.get('score'))
        if not vid or not gene or score is None or not 0 <= score <= 1:
            raise ValueError('Every variant requires var_key, genes, score in [0,1]')
        if vid in seen: raise ValueError(f'Duplicate variant: {vid}')
        seen.add(vid)
        v = dict(row, variant=vid, score=score)
        v['dp'] = number(row.get('dp')); v['gq'] = number(row.get('gq'))
        v['max_af'] = number(row.get('max_af'))
        v['qc_flags'] = []
        if v['dp'] is None or v['gq'] is None: v['qc_flags'].append('missing_depth_or_genotype_quality')
        if v['dp'] is not None and v['dp'] < 10: v['qc_flags'].append('low_depth')
        if v['gq'] is not None and v['gq'] < 20: v['qc_flags'].append('low_genotype_quality')
        v['eligible'] = not any(x.startswith('low_') for x in v['qc_flags'])
        if v['max_af'] is not None and v['max_af'] >= .01: v['eligible'] = False
        for symbol in re.split('[;,]', gene):
            if symbol.strip(): groups[symbol.strip()].append(dict(v, gene=symbol.strip()))
    ranked = []
    for gene, variants in groups.items():
        variants.sort(key=lambda v: (-v['score'], v['variant']))
        good = [v for v in variants if v['eligible']]
        if not good: continue
        matches = []
        for term in query:
            common = ancestors[term] & gene_terms.get(gene, set())
            best = max(common, key=lambda t: (ic.get(t, 0), t), default=None)
            denom = ic.get(term, math.log(n + 1))
            similarity = ic.get(best, 0) / denom if denom else 0
            matches.append({'query': term, 'query_name': names.get(term), 'common_ancestor': best,
                            'similarity': round(similarity, 6)})
        pheno = sum(m['similarity'] for m in matches) / len(query)
        high = [v for v in good if v.get('tier') in ('high', 'moderate')]
        hom = [v for v in high if v.get('zygosity') in ('hom', 'homozygous', '1/1', '1|1')]
        # Multiple nearby frameshifts can represent one event; flag, never count as proven alleles.
        cluster = False
        for i, a in enumerate(high):
            for b in high[i+1:]:
                pa, pb = a['variant'].split(':'), b['variant'].split(':')
                if len(pa) == len(pb) == 4 and pa[0] == pb[0] and abs(int(pa[1])-int(pb[1])) <= 20:
                    cluster = True
        recessive = 'HP:0000007' in inheritance.get(gene, set())
        genotype = ('homozygous_damaging_candidate' if hom else
                    'multiple_damaging_variants_phase_unknown' if len(high) >= 2 and not cluster else
                    'one_damaging_allele_second_unresolved' if high and recessive else
                    'inheritance_unresolved')
        # Neither count nor a high score establishes phase or functional rescue substrate.
        score = .55 * good[0]['score'] + .45 * pheno
        partners = [v for v in good if v['variant'] != good[0]['variant']]
        def partner_key(v):
            coding = bool(re.search('missense|stop_gained|frameshift|splice_(acceptor|donor)', v.get('consequence', '')))
            return (-int(coding), -(v['max_af'] is None), -v['score'], v['variant'])
        partners.sort(key=partner_key)
        ranked.append({'gene': gene, 'rank_score': round(score, 8), 'variant_score': good[0]['score'],
                       'phenotype_score': round(pheno, 8), 'phenotype_matches': matches,
                       'inheritance_annotations': sorted(inheritance.get(gene, ())),
                       'genotype_status': genotype, 'clustered_calls': cluster,
                       'lead_variant': good[0], 'second_allele_hypotheses': partners[:5],
                       'n_variants': len(variants), 'n_damaging': len(high),
                       'follow_up': ['phase_and_residual_function_unresolved'] +
                       (['depth_SV_MEI_follow_up_needed'] if recessive and not hom else [])})
    ranked.sort(key=lambda g: (-g['rank_score'], g['gene']))
    if not ranked: raise ValueError('No eligible variants remain')
    for i, g in enumerate(ranked, 1): g['rank'] = i
    return ranked

def make_handoff(case, ranked):
    lead = ranked[0]
    fields = ('variant','gene','score','tier','consequence','hgvsc','hgvsp',
              'zygosity','max_af','dp','gq','qc_flags','eligible')
    variants = [{k: v[k] for k in fields if k in v}
                for v in [lead['lead_variant'], *lead['second_allele_hypotheses']]]
    return {'schema_version': VERSION, 'case_id': case['case_id'],
            'data_class': case['data_class'], 'phenotype': case.get('phenotype', ''),
            'input_note': case.get('note', ''),
            'observed_hpo': case['observed_hpo'], 'scope': 'constitutional_research',
            'family_history_hpo': case.get('family_history_hpo', []),
            'lead_gene': lead['gene'], 'variants': variants,
            'alternative_genes': [{k: g[k] for k in ('gene','rank_score','genotype_status')} for g in ranked[1:6]],
            'genotype_status': lead['genotype_status'], 'phase': 'unknown',
            'residual_function': 'unestablished', 'follow_up': lead['follow_up'],
            'lesion_hypotheses': [
                {'interpretation': 'lead allele may contribute to phenotype', 'status': 'hypothesis'},
                {'interpretation': 'a second allele could leave rescuable residual product',
                 'status': 'conditional_on_phase_and_function'}],
            'limitations': ['Gene ranking is heuristic, not a calibrated causal probability.',
                           'Second-allele candidates retain their original pathogenicity scores.',
                           'Two variants do not establish trans phase or biallelic loss of function.',
                           'No raw-read artefact exclusion or SV/MEI execution in this milestone.']}

def normalise_text(s): return ' '.join(s.split()).casefold()

def evidence_pmids(anchors, searches, limit=25):
    """Keep discovery citations; interleave fresh searches to avoid truncation bias."""
    ordered = list(dict.fromkeys(str(p) for p in anchors if p))
    if len(ordered) > limit:
        raise ValueError('Too many discovery sources for the evidence budget; do not silently discard them')
    hits = [s.get('hits',[]) for s in searches]
    for i in range(max((len(h) for h in hits),default=0)):
        for group in hits:
            if i < len(group) and group[i].get('pmid'):
                pmid = str(group[i]['pmid'])
                if pmid not in ordered and len(ordered) < limit: ordered.append(pmid)
    return ordered

def audit_candidate(candidate, sources, review):
    """Fail closed: provenance checks plus independent model interpretation review.

    Quote validation proves source membership, NOT scientific entailment.
    """
    errors = []
    if type(candidate.get('nominate')) is not bool: errors.append('nomination_must_be_boolean')
    records = candidate.get('evidence', [])
    required = ('compound', 'intervention', 'system', 'dose', 'route', 'duration', 'endpoint',
                'direction', 'source_id', 'quote', 'interpretation', 'evidence_role')
    aliases = {normalise_text(str(s)) for s in [candidate.get('name',''),
        candidate.get('identity',{}).get('name',''),*candidate.get('identity',{}).get('aliases',[])] if s}
    direct_support = False
    if not records: errors.append('no_evidence_records')
    for i, record in enumerate(records):
        for key in required:
            if not isinstance(record.get(key), str) or not record[key].strip(): errors.append(f'evidence_{i}_missing_{key}')
        source = sources.get(str(record.get('source_id')))
        quote = record.get('quote', '')
        if source is None: errors.append(f'evidence_{i}_unknown_source')
        elif len(quote) < 30 or not any(normalise_text(quote) in normalise_text(t)
                for t in [source['text'],source.get('abstract_text','')]):
            errors.append(f'evidence_{i}_quote_not_in_source')
        role = record.get('evidence_role')
        if role not in ('direct_support','direct_opposing','context'):
            errors.append(f'evidence_{i}_invalid_role')
        if role in ('direct_support','direct_opposing'):
            if normalise_text(record.get('compound','')) not in aliases:
                errors.append(f'evidence_{i}_compound_mismatch')
            elif role == 'direct_support': direct_support = True
    if not direct_support and candidate.get('nominate') is True: errors.append('no_direct_compound_support')
    control = candidate.get('experiment',{}).get('mechanism_control',{})
    if not all(control.get(k) for k in ('perturbation','outcome','interpretation')):
        errors.append('missing_structured_mechanism_control')
    if (control.get('outcome') == 'loss_of_benefit'
        and control.get('interpretation') != 'supports_dependency_if_controls_valid'):
        errors.append('inverted_dependency_control_logic')
    if review.get('verdict') != 'supported': errors.append('independent_review_not_supported')
    if review.get('issues'): errors.append('independent_review_issues')
    if not candidate.get('identity', {}).get('resolved'): errors.append('unresolved_compound_identity')
    if not candidate.get('opposing_searches'): errors.append('opposing_search_missing')
    completed = candidate.get('opposing_search_results', [])
    if len(completed) != 2 or not all(r.get('ok') is True for r in completed):
        errors.append('opposing_search_execution_failed')
    if candidate.get('functional_rescue') not in ('unestablished', 'supported_in_stated_model'):
        errors.append('invalid_functional_rescue_status')
    return {'status': ('research_lead' if candidate.get('nominate') is True else 'not_nominated') if not errors else 'review_required',
            'errors': errors, 'clinical_suitability': 'unestablished',
            'caveat': 'Automatic checks and model review do not establish experimental validity or clinical suitability.'}
