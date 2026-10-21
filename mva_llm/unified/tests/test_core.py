import json
import sys
import tempfile
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from core import rank_variants,make_handoff,audit_candidate,evidence_pmids
from workflow import checkpoint

class Contracts(unittest.TestCase):
    def setUp(self):
        self.ont=({'HP:0000118':{'HP:0000118'},'HP:0010001':{'HP:0010001','HP:0000118'},
                   'HP:0010002':{'HP:0010002','HP:0000118'}},{}, {})
        self.ann=[{'gene_symbol':'MATCH','hpo_id':'HP:0010001'},
                  {'gene_symbol':'MATCH','hpo_id':'HP:0000007'},
                  {'gene_symbol':'OTHER','hpo_id':'HP:0010002'}]
        self.rows=[dict(var_key='1:100:A:G',genes='MATCH',score='.95',tier='high',dp='40',gq='99',
                        consequence='stop_gained',zygosity='het',max_af='.0001'),
                   dict(var_key='1:10000:C:T',genes='MATCH',score='.05',tier='low',dp='40',gq='99',
                        consequence='missense_variant',zygosity='het',max_af=''),
                   dict(var_key='2:100:A:G',genes='OTHER',score='.99',tier='high',dp='40',gq='99',
                        consequence='stop_gained',zygosity='het',max_af='.0001')]
    def rank(self): return rank_variants(self.rows,['HP:0010001'],self.ann,self.ont)
    def test_phenotype_changes_lead_without_gene_bonus(self):
        self.assertEqual(self.rank()[0]['gene'],'MATCH')
        result=rank_variants(self.rows,['HP:0010002'],self.ann,self.ont)
        self.assertEqual(result[0]['gene'],'OTHER')
    def test_second_hit_not_promoted(self):
        r=self.rank()[0]
        self.assertEqual(r['genotype_status'],'one_damaging_allele_second_unresolved')
        self.assertEqual(r['second_allele_hypotheses'][0]['score'],.05)
        h=make_handoff({'case_id':'synthetic','data_class':'synthetic','observed_hpo':['HP:0010001']},self.rank())
        self.assertEqual(h['phase'],'unknown'); self.assertEqual(h['residual_function'],'unestablished')
    def test_order_independent(self):
        first=self.rank(); self.rows.reverse(); self.assertEqual(first,self.rank())
    def test_low_quality_excludes_variant(self):
        self.rows[0]['dp']='2'
        self.assertEqual(self.rank()[0]['gene'],'OTHER')
    def test_duplicate_and_unknown_term_fail(self):
        with self.assertRaises(ValueError): rank_variants(self.rows,['HP:9999999'],self.ann,self.ont)
        self.rows.append(self.rows[0])
        with self.assertRaises(ValueError): self.rank()
    def test_cluster_not_proven_biallelic(self):
        self.rows[1].update(var_key='1:105:C:T',score='.98',tier='high')
        r=self.rank()[0]
        self.assertTrue(r['clustered_calls'])
        self.assertNotEqual(r['genotype_status'],'multiple_damaging_variants_phase_unknown')
    def test_resume_integrity(self):
        with tempfile.TemporaryDirectory() as t:
            out=Path(t)
            checkpoint(out,'stage',lambda:{'x':1},False)
            self.assertEqual(checkpoint(out,'stage',lambda:1/0,True),{'x':1})
            (out/'stage.json').write_text('{}')
            with self.assertRaises(ValueError): checkpoint(out,'stage',lambda:{},True)

class Evidence(unittest.TestCase):
    def test_discovery_sources_survive_fresh_search(self):
        searches=[{'hits':[{'pmid':str(i)} for i in range(10)]}]
        self.assertEqual(evidence_pmids(['primary'],searches,3),['primary','0','1'])
    def test_searches_share_budget_and_duplicates_do_not_consume_it(self):
        searches=[{'hits':[{'pmid':'a'},{'pmid':'b'}]}, {'hits':[{'pmid':'c'},{'pmid':'d'}]}]
        self.assertEqual(evidence_pmids(['a'],searches,3),['a','c','b'])
    def setUp(self):
        self.quote='Compound A increased protein abundance in cultured mouse cells.'
        self.sources={'123':{'text':self.quote}}
        self.c={'name':'Compound A','nominate':True,'identity':{'resolved':True},
                'opposing_searches':['degradation'], 'functional_rescue':'unestablished',
                'opposing_search_results':[{'ok':True},{'ok':True}],
                'experiment':{'mechanism_control':{'perturbation':'block_proposed_mechanism',
                    'outcome':'loss_of_benefit','interpretation':'supports_dependency_if_controls_valid'}},
                'evidence':[dict(compound='Compound A',evidence_role='direct_support',intervention='compound A exposure',system='mouse cells',
                    dose='not reported',route='in vitro',duration='not reported',endpoint='protein abundance',
                    direction='increased',source_id='123',quote=self.quote,interpretation='abundance hypothesis')]}
        self.review={'verdict':'supported','issues':[]}
    def audit(self): return audit_candidate(self.c,self.sources,self.review)
    def test_bounded_lead_never_clinical(self):
        self.assertEqual(self.audit()['status'],'research_lead')
        self.assertEqual(self.audit()['clinical_suitability'],'unestablished')
    def test_missing_paper_abstains(self):
        self.sources.clear(); self.assertEqual(self.audit()['status'],'review_required')
    def test_invented_quote_abstains(self):
        self.c['evidence'][0]['quote']='This molecule rescued affected children permanently.'
        self.assertIn('evidence_0_quote_not_in_source',self.audit()['errors'])
    def test_abstract_quote_remains_valid_with_fulltext(self):
        self.sources['123']['abstract_text']=self.quote
        self.sources['123']['text']='The full text contains different notation and additional methods.'
        self.assertNotIn('evidence_0_quote_not_in_source',self.audit()['errors'])
    def test_quote_cannot_span_source_versions(self):
        self.sources['123']['abstract_text']='a sentence ending here'
        self.sources['123']['text']='and starting somewhere else in the full text'
        self.c['evidence'][0]['quote']='ending here and starting somewhere else'
        self.assertIn('evidence_0_quote_not_in_source',self.audit()['errors'])
    def test_cross_compound_transfer_abstains(self):
        self.c['evidence'][0]['compound']='Compound B'
        self.assertIn('evidence_0_compound_mismatch',self.audit()['errors'])
    def test_ambiguous_identity_abstains(self):
        self.c['identity']['resolved']=False
        self.assertEqual(self.audit()['status'],'review_required')
    def test_review_disagreement_abstains(self):
        self.review['issues']=['Genetic-intervention result assigned to drug']
        self.assertEqual(self.audit()['status'],'review_required')
    def test_negative_nomination_not_overridden(self):
        self.c['nominate']=False
        self.assertEqual(self.audit()['status'],'not_nominated')
    def test_context_not_transferred_to_direct_evidence(self):
        self.c['evidence'][0].update(compound='none',evidence_role='context')
        self.assertIn('no_direct_compound_support',self.audit()['errors'])
        self.assertNotIn('evidence_0_compound_mismatch',self.audit()['errors'])
    def test_inverted_dependency_control_rejected(self):
        self.c['experiment']['mechanism_control']['interpretation']='falsifies_mechanism'
        self.assertIn('inverted_dependency_control_logic',self.audit()['errors'])
    def test_failed_opposing_search_is_not_a_completed_search(self):
        self.c['opposing_search_results'][0]={'ok':False,'error':'corpus unavailable'}
        self.assertIn('opposing_search_execution_failed',self.audit()['errors'])

if __name__=='__main__': unittest.main()
