import copy
import csv
import io
import math
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from submission import build_submission, write_submission, FIELDS


def variant(key, score=.1, zygosity='het', eligible=True):
    return dict(variant=key, score=score, zygosity=zygosity, eligible=eligible)


def gene(symbol='GENE_A', score=.8):
    return dict(gene=symbol, rank_score=score, inheritance_annotations=['HP:0000007'],
                lead_variant=variant('1:100:A:G', .9),
                second_allele_hypotheses=[variant('1:1000:C:T'),variant('1:2000:G:A')])


class SubmissionTests(unittest.TestCase):
    def build(self, genes=None, options=None):
        return build_submission({'case_id':'SYNTHETIC'},genes or [gene()],options)

    def test_pair_single_alternative_and_formula(self):
        r=self.build()
        self.assertEqual([float(x['epcr']) for x in r['rows']], [.8,.4,.266666667])
        self.assertEqual([x['pos_2'] for x in r['rows']], ['1000','','2000'])
        self.assertEqual(r['hypotheses'][0]['partner_score'], .1)
        self.assertIn('phase unknown', r['rows'][0]['notes'])

    def test_no_pairs_for_homozygous_lead_or_without_ar(self):
        for update in ('hom', 'inheritance'):
            g=gene()
            if update=='hom': g['lead_variant']['zygosity']='hom'
            else: g['inheritance_annotations']=[]
            r=self.build([g])
            self.assertEqual(len(r['rows']), 1)
            self.assertEqual(r['rows'][0]['pos_2'], '')

    def test_ineligible_homozygous_nearby_and_cross_chrom_partners_excluded(self):
        g=gene()
        g['second_allele_hypotheses']=[variant('1:100:A:G'), variant('1:110:T:C'),
            variant('2:1000:C:T'), variant('1:1000:C:T', zygosity='hom'),
            variant('1:2000:G:A', eligible=False)]
        self.assertEqual(len(self.build([g])['rows']),1)
        g=gene(); g['clustered_calls']=True
        self.assertEqual(len(self.build([g])['rows']),1)

    def test_deduplication_limits_and_determinism(self):
        a=gene(); b=copy.deepcopy(a); b['gene']='GENE_B'
        self.assertEqual(len(self.build([a,b])['rows']),3)
        self.assertEqual(self.build([a,b]),self.build([b,a]))
        genes=[]
        for i in range(15):
            g=gene(f'G{i}',.8-i*.01); g['lead_variant']=variant(f'1:{100+i}:A:C')
            g['second_allele_hypotheses']=[]; genes.append(g)
        r=self.build(genes)
        self.assertEqual(len(r['rows']),10)
        self.assertEqual(r['omitted_count'],5)
        self.assertTrue(all(x['finding_type']=='primary' for x in r['rows']))

    def test_no_probability_claim_or_legacy_answer_dependency(self):
        g=gene(); original=copy.deepcopy(g)
        a=self.build([g]); g['rank_score']=.6; b=self.build([g])
        self.assertEqual(float(b['rows'][0]['epcr']),.6)
        self.assertIn('not a posterior', a['epcr_interpretation'])
        self.assertEqual(original['lead_variant'],g['lead_variant'])

    def test_validation(self):
        for options in ({'max_rows':11},{'max_rows':0},{'max_rows':True}, {'genome_build':'GRCh37'}, {'unknown':1}):
            with self.assertRaises(ValueError): self.build(options=options)
        for score in (float('nan'),float('inf'),-1,1.1):
            with self.assertRaises(ValueError): self.build([gene(score=score)])
        for key in ('1:0:A:G','1:100:A:<DEL>','1:100:A:A','99:100:A:G'):
            g=gene();g['lead_variant']['variant']=key
            with self.assertRaises(ValueError): self.build([g])
        g=gene(); g['lead_variant']['eligible']=False
        with self.assertRaises(ValueError): self.build([g])

    def test_csv_contract_and_repeatable_render(self):
        with tempfile.TemporaryDirectory() as directory:
            p=Path(directory)/'out.csv'; r=self.build(); write_submission(r,p)
            first=p.read_bytes(); write_submission(r,p)
            self.assertEqual(first,p.read_bytes())
            reader=csv.DictReader(io.StringIO(p.read_text()))
            self.assertEqual(reader.fieldnames,FIELDS)
            rows=list(reader); self.assertEqual(rows,r['rows'])
            self.assertTrue(all(0<float(row['epcr'])<=1 for row in rows))
            self.assertTrue(all(row['chrom_1']=='chr1' for row in rows))
