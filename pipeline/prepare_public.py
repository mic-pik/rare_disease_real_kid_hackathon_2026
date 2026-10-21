#!/usr/bin/env python3
"""Fetch public resources and rebuild local evidence. Never reads case inputs.

Downloads are atomic. config/sources.json contains versioned URLs and available
historical checksums. Credentials for gated model access come from HF_TOKEN or
the user's Hugging Face credential store, never from a configuration file.
"""
import argparse
import gzip
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]


def download(url, path, checksum=None, algorithm='sha256'):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def valid():
        if not path.is_file() or path.stat().st_size == 0: return False
        if not checksum: return True
        h = hashlib.new(algorithm)
        with path.open('rb') as f:
            for block in iter(lambda: f.read(8 << 20), b''): h.update(block)
        return h.hexdigest() == checksum
    if valid(): return path
    part = path.with_name(path.name+'.part')
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'WorkingVariant-reproduction/1.0'})
            with urllib.request.urlopen(req, timeout=180) as source, part.open('wb') as target:
                shutil.copyfileobj(source, target, 8 << 20)
            part.replace(path)
            if not valid():
                raise ValueError(f'Checksum mismatch: {path.name}. Do not substitute a newer release silently.')
            return path
        except (OSError, ValueError) as exc:
            if attempt == 3: raise RuntimeError(f'Download failed for {path.name}; check URL/access in config/sources.json') from exc
            time.sleep(2 ** attempt)


def extract(archive, dest):
    dest = Path(dest).resolve()
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            target = (dest/member.name).resolve()
            if not target.is_relative_to(dest) or member.issym() or member.islnk() or member.isdev():
                raise ValueError('Unsafe archive member')
        tar.extractall(dest)


def resources(group):
    records = json.loads((ROOT/'config/sources.json').read_text())
    for r in records:
        if r['group'] == group:
            download(r['url'], ROOT/r['destination'], r.get('checksum'), r.get('algorithm', 'sha256'))


def models():
    from huggingface_hub import snapshot_download
    for record in json.loads((ROOT/'config/models.json').read_text()):
        snapshot_download(repo_id=record['repo_id'], revision=record['revision'],
                          local_dir=ROOT/record['destination'], ignore_patterns=['tf_model.h5', '*.msgpack'])


def references():
    resources('references')
    genome = ROOT/'track1/ref/grch38/GRCh38.primary_assembly.fa'
    if not genome.exists():
        with gzip.open(genome.with_suffix('.fa.gz'), 'rb') as source, genome.open('wb') as target:
            shutil.copyfileobj(source, target)
    subprocess.run(['samtools', 'faidx', str(genome)], check=True)
    subprocess.run(['tabix', '-f', '-p', 'vcf', str(ROOT/'track1/ref/clinvar/clinvar.vcf.gz')], check=True)


def knowledge_inputs():
    resources('knowledge')
    cfg = json.loads((ROOT/'config/opentargets.json').read_text())
    for table in cfg['tables']:
        url = cfg['base_url'].rstrip('/')+'/'+table+'/'
        with urllib.request.urlopen(url, timeout=180) as response: listing = response.read().decode()
        names = sorted(set(re.findall(r'href="([^"/]+\.parquet)"', listing)))
        if not names: raise ValueError(f'No parquet partitions in {table}; inspect the pinned release URL')
        for name in names: download(url+name, ROOT/'mva_llm/data/raw/opentargets'/table/name)
    archive = ROOT/'gene_mol/data/chembl37_raw/chembl_37_sqlite.tar.gz'
    db = ROOT/'gene_mol/data/chembl37_db/chembl_37/chembl_37_sqlite/chembl_37.db'
    if not db.exists(): extract(archive, ROOT/'gene_mol/data/chembl37_db')
    # Only the assessed germline gene associations consumed by sibling_diseases.
    import pandas as pd
    xml = ET.parse(ROOT/'gene_mol/data/orphanet/en_product6.xml').getroot()
    if xml.findtext('Availability/Licence/ShortIdentifier') != 'CC-BY-4.0':
        raise ValueError('Unexpected Orphadata licence; review the resource before reuse')
    types = {'Disease-causing germline mutation(s) in', 'Disease-causing germline mutation(s) (loss of function) in',
             'Disease-causing germline mutation(s) (gain of function) in', 'Major susceptibility factor in'}
    rows = []
    for disorder in xml.iter('Disorder'):
        for association in disorder.iter('DisorderGeneAssociation'):
            if association.findtext('DisorderGeneAssociationStatus/Name') != 'Assessed': continue
            if association.findtext('DisorderGeneAssociationType/Name') not in types: continue
            gene = association.find('Gene')
            refs = {x.findtext('Source'): x.findtext('Reference') for x in gene.iter('ExternalReference')}
            rows.append({'orpha': disorder.findtext('OrphaCode'), 'disorder': disorder.findtext('Name'),
                         'symbol': gene.findtext('Symbol'), 'ensembl': refs.get('Ensembl'),
                         'assoc_type': association.findtext('DisorderGeneAssociationType/Name'),
                         'assoc_status': 'Assessed'})
    pd.DataFrame(rows).to_csv(ROOT/'gene_mol/data/orphanet/genes.tsv', sep='\t', index=False)


def literature():
    """Fetch fixed public IDs, preserving historical corpus order and query labels."""
    import pandas as pd
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location('lit_builder', ROOT/'mva_llm/src/build_literature_index.py')
    builder = module_from_spec(spec); spec.loader.exec_module(builder)
    records = json.loads((ROOT/'config/literature_ids.json').read_text())
    raw = ROOT/'mva_llm/data/raw/literature_frozen'
    rows = []
    for start in range(0, len(records), 100):
        batch = records[start:start+100]
        query = 'EXT_ID:('+ ' OR '.join(r['pmid'] for r in batch)+') AND SRC:MED'
        url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode(
            {'query': query, 'resultType': 'core', 'pageSize': 100, 'format': 'json'})
        data = json.loads(download(url, raw/f'{start:06}.json').read_text())
        found = {str(h.get('pmid')): h for h in data.get('resultList', {}).get('result', [])}
        for record in batch:
            if record['pmid'] not in found or not found[record['pmid']].get('abstractText'):
                raise ValueError(f"Missing historical abstract {record['pmid']}; corpus reproduction incomplete")
            h = found[record['pmid']]
            rows.append({'pmid': record['pmid'], 'doi': h.get('doi'), 'title': builder.clean_text(h.get('title') or ''),
                         'abstract': builder.clean_text(h['abstractText']), 'year': h.get('pubYear'),
                         'journal': h.get('journalTitle'), 'licence': h.get('license'), 'corpus_query': record['corpus_query']})
        print(f'Literature {len(rows)}/{len(records)}', flush=True)
        time.sleep(.2)
    target = ROOT/'mva_llm/data/kg/literature.parquet'
    target.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_parquet(target, index=False)


class Paragraphs(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True); self.rows = []; self.stack = []
    def handle_starttag(self, tag, attrs):
        if tag == 'p':
            self.rows.append([]); self.stack.append(len(self.rows)-1)
    def handle_data(self, data):
        for index in self.stack: self.rows[index].append(data)
    def handle_endtag(self, tag):
        if tag == 'p' and self.stack: self.stack.pop()
    @property
    def paragraphs(self):
        # Preserve inline superscripts and empty paragraphs as in the original
        # extraction; inserting spaces around tags changes source quotations.
        return [re.sub(r'\s+', ' ', ''.join(parts)).strip() for parts in self.rows]


def supplements():
    import pandas as pd
    resources('supplements')
    out = ROOT/'mva_llm/data/unified_public_v1'
    (out/'fulltext').mkdir(parents=True, exist_ok=True)
    html = ROOT/'mva_llm/data/raw/supplements/PMC4194088.html'
    parser = Paragraphs(); parser.feed(html.read_text())
    if len(parser.paragraphs) < 30 or not any('BubR1' in p and 'NMN' in p for p in parser.paragraphs):
        raise ValueError('PMC response lacks the expected paper; it may be a bot challenge. Supply the public HTML snapshot.')
    (out/'fulltext/24825348.json').write_text(json.dumps({'source_id': '24825348',
        'url': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC4194088/', 'level': 'fulltext',
        'text': '\n\n'.join(parser.paragraphs), 'provenance': 'Rebuilt from public PMC HTML; inspect snapshot differences.'}, indent=2)+'\n')
    pubchem = json.loads((ROOT/'mva_llm/data/raw/supplements/nmn_pubchem.json').read_text())['PropertyTable']['Properties'][0]
    smiles = pubchem.get('SMILES') or pubchem.get('IsomericSMILES')
    if not smiles: raise ValueError('PubChem response lacks an isomeric structure')
    identities = [{'name':'beta-nicotinamide mononucleotide','id':'PubChem:14180','smiles':smiles,
        'aliases':['NMN','nicotinamide mononucleotide','beta-NMN','β-nicotinamide mononucleotide',
                   'nicotinamide mononucleotide (NMN)','NMN (nicotinamide mononucleotide)','β-nicotinamide mononucleotide (NMN)'],
        'source':'https://pubchem.ncbi.nlm.nih.gov/compound/14180'}]
    drugs = pd.read_parquet(ROOT/'mva_llm/data/kg/drugs.parquet')
    match = drugs[drugs.drug_id == 'CHEMBL1140']
    if len(match) != 1: raise ValueError('Nicotinamide identity absent or ambiguous')
    identities.append({'name':'nicotinamide','id':'CHEMBL1140','smiles':match.iloc[0].smiles,
        'aliases':['NAM','niacinamide','nicotinamide (niacinamide)','nicotinamide (NAM)'],
        'source':'local ChEMBL / Open Targets drugs.parquet CHEMBL1140'})
    (out/'identities.json').write_text(json.dumps(identities, indent=2)+'\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('stage', choices=['download','references','knowledge-inputs','literature','supplements'])
    args = ap.parse_args()
    {'download': models, 'references': references, 'knowledge-inputs': knowledge_inputs,
     'literature': literature, 'supplements': supplements}[args.stage]()


if __name__ == '__main__': main()
