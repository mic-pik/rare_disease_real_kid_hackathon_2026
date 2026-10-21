import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import pandas as pd
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline'))
import benchmark_pathways as b
from benchmark_recovery import select_genetics, seal


class PathwayTests(unittest.TestCase):
    def test_uncapped_retains_all_above_threshold(self):
        f=pd.DataFrame({'gene_id':['G'+str(i) for i in range(60)],'disease_id':['D']*60,
                        'datatype':['genetic_association']*60,'datatype_score':[.2]*60})
        self.assertEqual(len(select_genetics(f,['D'],.1,None)),60)
        self.assertEqual(len(select_genetics(f,['D'],.1,50)),50)

    def test_sibling_family_leakage_rejected(self):
        q=[{'disease_id':'A','name':'a','split':'development','family_root':'F'},
           {'disease_id':'B','name':'b','split':'heldout','family_root':'F'}]
        ontology=[{'id':'A','name':'a','ancestors':['F']},{'id':'B','name':'b','ancestors':['F']},
                  {'id':'F','name':'family','ancestors':[]}]
        with self.assertRaisesRegex(ValueError,'overlap'):b.validate_families(q,ontology)

    def test_disease_pathways_excluded_wholesale(self):
        f=pd.DataFrame([['G1','P','p','Signal Transduction'],['G2','P','p','Disease'],
                        ['G3','Q','q','Metabolism'],['G4','Q','q','Metabolism']],
                       columns=['gene_id','pathway_id','pathway','top_level'])
        out=b.filter_pathways(f,{'excluded_top_levels':['Disease','Drug ADME'],
                                'pathway_min_genes':2,'pathway_max_genes':100})
        self.assertEqual(set(out.pathway_id),{'Q'})

    def test_distinct_seed_required_and_large_paths_penalized(self):
        g=pd.DataFrame({'gene_id':['G1'],'datatype_score':[1.]})
        p=pd.DataFrame([['G1','P','p',2],['G2','P','p',2],['G1','Q','q',3],
                        ['G3','Q','q',3],['G4','Q','q',3]],
                       columns=['gene_id','pathway_id','pathway','pathway_size'])
        links=b.pathway_links(g,p,{'G1','G2','G3'}).set_index('gene_id')
        self.assertNotIn('G1',links.index)
        self.assertEqual(links.loc['G2','pathway_score'],1)
        self.assertLess(links.loc['G3','pathway_score'],1)

    def test_no_inferred_direction_from_membership(self):
        self.assertEqual(b.direction_status(['INHIBITOR']),'unknown_required_direction')
        self.assertEqual(b.direction_status(['INHIBITOR'],'increase'),'opposing_action')
        self.assertEqual(b.direction_status(['ACTIVATOR'],'increase'),'compatible_action_only')
        self.assertEqual(b.direction_status(['INHIBITOR','ACTIVATOR'],'increase'),'mixed_drug_actions')
        self.assertEqual(b.direction_status(['MODULATOR'],'increase'),'unknown_drug_action')
        self.assertEqual(b.direction_status(['PARTIAL AGONIST'],'increase'),'unknown_drug_action')

    def test_direction_requires_source(self):
        with self.assertRaisesRegex(ValueError,'source'):
            b.validate_constraints([{'disease_id':'D','gene_id':'G','required_action':'increase','source':'','locator':''}])

    def test_bootstrap_pairs_same_diseases(self):
        m=pd.DataFrame([['A','new',.4],['B','new',.2],['B','old',.1],['A','old',.3]],
                       columns=['disease_id','method','recall@10'])
        result=b.paired_delta(m,'new','old',1,100)
        self.assertAlmostEqual(result['mean_delta'],.1)
        self.assertAlmostEqual(result['lower95'],.1)

    def fixture(self,path):
        pd.DataFrame({'ingredient_id':['C1','C2'],'canonical_name':['one','two'],
                      'smiles':['CCO','CC'],'inchikey':['A','B']}).to_parquet(path/'catalog.parquet')
        pd.DataFrame({'ingredient_id':['C1','C2'],'gene_id':['G2','G3'],
                      'action_type':['INHIBITOR','ACTIVATOR']}).to_parquet(path/'mechanisms.parquet')
        pd.DataFrame({'disease_id':['D'],'gene_id':['G1'],'datatype_score':[1.]}).to_parquet(path/'genetics.parquet')
        pd.DataFrame({'gene_id':['G1','G2','G3'],'pathway_id':['P']*3,'pathway':['p']*3,
                      'pathway_size':[3]*3}).to_parquet(path/'pathways.parquet')
        (path/'queries.json').write_text(json.dumps([{'disease_id':'D','name':'d','split':'development','family_root':'F'}]))
        (path/'policy.json').write_text(json.dumps({'baseline_gene_limit':50,'pathway_weight':.5}))
        (path/'directions.json').write_text(json.dumps([{'disease_id':'D','gene_id':'G2','required_action':'increase',
                                                       'source':'synthetic-test','locator':'synthetic fixture'}]))
        seal(path)

    def test_opposing_route_audited_but_excluded_and_truth_cannot_change_scores(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);p=root/'features';p.mkdir();self.fixture(p)
            (root/'truth.json').write_text('{"D":"C1"}');b.score(SimpleNamespace(features=p,out=root/'one'))
            (root/'truth.json').write_text('{"D":"C2"}');b.score(SimpleNamespace(features=p,out=root/'two'))
            s=pd.read_parquet(root/'one/scores.parquet')
            pd.testing.assert_frame_equal(s,pd.read_parquet(root/'two/scores.parquet'))
            ss=s[s.method.eq('pathway_augmented')].set_index('ingredient_id').score
            self.assertEqual(ss['C1'],0);self.assertGreater(ss['C2'],0)
            routes=pd.read_parquet(root/'one/pathway_routes.parquet').set_index('ingredient_id')
            self.assertEqual(routes.loc['C1','direction_status'],'opposing_action')
            self.assertEqual(routes.loc['C2','direction_status'],'unknown_required_direction')
            (p/'labels.csv').write_text('D,C1')
            with self.assertRaisesRegex(ValueError,'allowlist'):b.score(SimpleNamespace(features=p,out=root/'three'))


if __name__=='__main__':unittest.main()
