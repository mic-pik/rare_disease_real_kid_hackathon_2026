import itertools
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline'))
import benchmark_recovery as b


class RecoveryBenchmarkTests(unittest.TestCase):
    def test_ties_match_exhaustive_permutations(self):
        # Two negatives above a tie containing two positives and two negatives.
        values=[]
        for perm in itertools.permutations([True, True, False, False]):
            labels=[False,False,*perm]
            values.append((sum(labels[:3])/2, float(any(labels[:3])),
                           1/(labels.index(True)+1)))
        result=b.tie_metrics([2,2,1,1,1,1],[False,False,True,True,False,False],ks=(3,))
        for key,expected in zip(['recall@3','hit@3','mrr'],np.mean(values,axis=0)):
            self.assertAlmostEqual(result[key],expected)

    def test_recall_is_fraction_of_all_positives_not_any_hit(self):
        result=b.tie_metrics([4,3,2,1],[True,False,False,True],ks=(1,))
        self.assertEqual(result['recall@1'],.5)
        self.assertEqual(result['hit@1'],1)

    def test_zero_scores_do_not_gain_coverage_through_id_order(self):
        result=b.tie_metrics([0]*100,[True]+[False]*99,ks=(10,))
        self.assertAlmostEqual(result['recall@10'],.1)
        self.assertAlmostEqual(result['hit@10'],.1)
        self.assertEqual(result['positive_score_coverage'],0)

    def test_missing_truth_or_nonfinite_score_rejected(self):
        for scores,labels in [([1],[False]),([np.nan],[True])]:
            with self.assertRaises(ValueError):b.tie_metrics(scores,labels)

    def test_known_drug_and_literature_evidence_never_select_genes(self):
        frame=pd.DataFrame([
            ['D','G1','known_drug',1],['D','G2','literature',1],
            ['D','G3','genetic_association',.2],['D','G4','genetic_association',.3],
            ['OTHER','G5','genetic_association',1]],
            columns=['disease_id','gene_id','datatype','datatype_score'])
        result=b.select_genetics(frame,['D'],.1,1)
        self.assertEqual(result.gene_id.tolist(),['G4'])

    def test_missing_target_does_not_remove_candidate(self):
        links=pd.DataFrame({'ingredient_id':['C1'],'gene_id':['G']})
        genetics=pd.DataFrame({'gene_id':['G'],'datatype_score':[.7]})
        np.testing.assert_equal(b.overlap(['C1','C2'],links,genetics),[.7,0])

    def fixture(self, path):
        pd.DataFrame({'ingredient_id':['C1','C2','C3'],
                      'canonical_name':['one','two','three'],
                      'smiles':['CCO','CCO','CCCO'],
                      'inchikey':['SAME-A','SAME-B','OTHER-A']}).to_parquet(path/'catalog.parquet')
        pd.DataFrame({'disease_id':['D'],'gene_id':['G'],
                      'datatype_score':[1.],'symbol':['GENE']}).to_parquet(path/'genetics.parquet')
        pd.DataFrame({'ingredient_id':['C1'],'gene_id':['G'],
                      'action_type':['INHIBITOR']}).to_parquet(path/'mechanisms.parquet')
        pd.DataFrame({'ingredient_id':['C1'],'gene_id':['G'],
                      'activity_id':[1]}).to_parquet(path/'active_assays.parquet')
        (path/'queries.json').write_text(json.dumps([{'disease_id':'D','name':'disease'}]))
        b.seal(path)

    def test_separate_scoring_and_no_identity_similarity_shortcut(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); features=root/'features';features.mkdir();self.fixture(features)
            # Truth deliberately exists outside the feature directory.
            (root/'truth.json').write_text('{"D":"C3"}')
            b.score(SimpleNamespace(features=features,out=root/'first'))
            (root/'truth.json').write_text('{"D":"C1"}')
            b.score(SimpleNamespace(features=features,out=root/'second'))
            one=pd.read_parquet(root/'first/scores.parquet')
            pd.testing.assert_frame_equal(one,pd.read_parquet(root/'second/scores.parquet'))
            similarity=one[one.method.eq('chemical_similarity')].set_index('ingredient_id').score
            self.assertEqual(similarity['C1'],0)
            self.assertEqual(similarity['C2'],0)
            self.assertGreater(similarity['C3'],0)
            (features/'labels.csv').write_text('D,C1')
            with self.assertRaisesRegex(ValueError,'allowlist'):
                b.score(SimpleNamespace(features=features,out=root/'third'))

    def test_tampered_inputs_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp);(p/'data').write_text('before');b.seal(p)
            (p/'data').write_text('after')
            with self.assertRaisesRegex(ValueError,'Hash mismatch'):b.verify(p)


if __name__=='__main__':unittest.main()
