#!/usr/bin/env python3
"""Post-run evaluation only. Expected answers are NEVER read by workflow.py."""
import argparse
import json
from pathlib import Path
from core import write_json

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',required=True); ap.add_argument('--expected',required=True); ap.add_argument('--out',required=True)
    args=ap.parse_args(); run=Path(args.run)
    expected=json.loads(Path(args.expected).read_text())
    handoff=json.loads((run/'02_handoff.json').read_text())
    ranks=json.loads((run/'01_ranking.json').read_text())
    result={'evaluation_type':'development replay, not independent validation',
            'expected_gene_rank':next((r['rank'] for r in ranks if r['gene']==expected['gene']),None),
            'lead_variant_match':handoff['variants'][0]['variant']==expected['lead_variant'],
            'first_partner_match':len(handoff['variants'])>1 and handoff['variants'][1]['variant']==expected['partner_variant'],
            'phase_preserved_unknown':handoff['phase']=='unknown'}
    if (run/'05_track2.json').exists():
        t2=json.loads((run/'05_track2.json').read_text())
        names=[r['candidate']['name'] for r in t2['candidates']]
        aliases={x.casefold() for x in expected.get('compound_aliases',[])}
        result['hypothesis_retrieved']=any(n.casefold() in aliases for n in names)
        result['hypothesis_passed_audit']=any(n.casefold() in aliases for n in t2['research_leads'])
        result['clinical_suitability']=t2['clinical_suitability']
    write_json(args.out,result); print(json.dumps(result,indent=2))

if __name__=='__main__': main()
