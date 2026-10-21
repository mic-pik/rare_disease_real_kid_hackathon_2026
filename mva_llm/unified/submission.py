"""Deterministic Track 1 submission export; EPCR values are uncalibrated priorities."""
from __future__ import annotations
import csv
import math
import re
from pathlib import Path

FIELDS = ['proband_id', 'chrom_1', 'pos_1', 'ref_1', 'alt_1',
          'chrom_2', 'pos_2', 'ref_2', 'alt_2', 'epcr', 'finding_type', 'notes']
POLICY = 'gene_score_divided_by_local_hypothesis_rank'
HET = {'het', '0/1', '1/0', '0|1', '1|0'}


def coordinates(key):
    parts = key.split(':')
    if len(parts) != 4:
        raise ValueError(f'Invalid variant key: {key}')
    chrom, pos, ref, alt = parts
    chrom = chrom.removeprefix('chr')
    if chrom == 'MT': chrom = 'M'
    if chrom not in {*map(str, range(1, 23)), 'X', 'Y', 'M'} or not pos.isdigit() or int(pos) < 1:
        raise ValueError(f'Invalid coordinates: {key}')
    if not re.fullmatch('[ACGT]+', ref) or not re.fullmatch('[ACGT]+', alt) or ref == alt:
        raise ValueError(f'Expected normalised sequence alleles: {key}')
    return [f'chr{chrom}', str(int(pos)), ref, alt]


def build_submission(case, ranked, options=None):
    options = options or {}
    if set(options) - {'genome_build', 'max_rows'}:
        raise ValueError('Unknown submission option')
    if options.get('genome_build', 'GRCh38') != 'GRCh38':
        raise ValueError('Submission requires GRCh38 coordinates; no implicit liftover')
    count = options.get('max_rows', 10)
    if type(count) is not int or not 1 <= count <= 10:
        raise ValueError('Submission max_rows must be an integer in [1, 10]')
    case_id = case.get('case_id')
    if not isinstance(case_id, str) or not case_id.strip():
        raise ValueError('Nonempty case_id required')
    hypotheses = []
    for gene in ranked:
        score = float(gene['rank_score'])
        if not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError('Gene score must be finite and within [0, 1]')
        if score == 0: continue
        lead = gene['lead_variant']
        if not lead.get('eligible', False): continue
        first = coordinates(lead['variant'])
        pairs = []
        # Gene-level AR annotation permits a hypothesis, never establishes mode or phase.
        if ('HP:0000007' in gene.get('inheritance_annotations', [])
                and lead.get('zygosity') in HET and not gene.get('clustered_calls')):
            seen = {tuple(first)}
            for partner in gene.get('second_allele_hypotheses', []):
                if not partner.get('eligible', False) or partner.get('zygosity') not in HET:
                    continue
                second = coordinates(partner['variant'])
                if tuple(second) in seen: continue
                seen.add(tuple(second))
                if first[0] != second[0] or abs(int(first[1])-int(second[1])) <= 20:
                    continue
                pairs.append(partner)
        # Preserve the existing first-partner preference, then the single-variant alternative.
        local = [pairs[0], None, *pairs[1:]] if pairs else [None]
        for local_rank, partner in enumerate(local, 1):
            second = coordinates(partner['variant']) if partner else ['', '', '', '']
            value = score / local_rank
            if value < 1e-9: continue  # explicit CSV precision floor, never round to zero
            description = ('candidate compound-het; phase unknown; functional effect unresolved'
                           if partner else 'single-variant hypothesis; inheritance interpretation unresolved')
            note = (f"{gene['gene']}: {description}. Lead score={lead['score']:.6f}"
                    + (f"; partner score={partner['score']:.6f}" if partner else '')
                    + f". EPCR=gene score {score:.8f}/local hypothesis rank {local_rank}; "
                      'uncalibrated priority, not a validated causal probability.')
            if gene.get('clustered_calls'): note += ' Nearby damaging calls flagged; artefact review required.'
            hypotheses.append({'gene': gene['gene'], 'gene_score': score,
                'local_hypothesis_rank': local_rank, 'raw_epcr': value,
                'lead_score': lead['score'], 'partner_score': partner['score'] if partner else None,
                'row': dict(zip(FIELDS, [case_id, *first, *second, f'{value:.9f}', 'primary', note]))})
    hypotheses.sort(key=lambda h: (-h['raw_epcr'], h['gene'],
                                    tuple(h['row'][f] for f in FIELDS[1:9])))
    unique, seen = [], set()
    for h in hypotheses:
        row = h['row']
        a = tuple(row[f] for f in FIELDS[1:5])
        b = tuple(row[f] for f in FIELDS[5:9])
        key = tuple(sorted([a, b])) if b[0] else (a,)
        if key not in seen:
            unique.append(h); seen.add(key)
    if not unique: raise ValueError('No eligible submission hypotheses')
    selected = unique[:count]
    return {'policy': POLICY, 'genome_build': 'GRCh38', 'max_rows': count,
            'epcr_interpretation': 'Uncalibrated ranking estimate; not a posterior probability. Values need not sum to one.',
            'formula': 'EPCR = gene rank_score / within-gene hypothesis rank (1-based)',
            'selection': 'First eligible AR heterozygous pair, single lead, remaining eligible pairs per gene; global descending EPCR; deduplicate; keep max_rows.',
            'secondary_findings': 'Not inferred. All exported hypotheses compete as primary explanations.',
            'candidate_count': len(unique), 'omitted_count': len(unique)-len(selected),
            'hypotheses': selected, 'rows': [h['row'] for h in selected]}


def write_submission(result, path):
    """Render the sealed JSON stage deterministically, including on resume."""
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, lineterminator='\n')
        writer.writeheader()
        writer.writerows(result['rows'])
    temporary.replace(path)
