import importlib.util
from pathlib import Path
import sys
import unittest

try:
    import pandas as pd
except ImportError:
    pd = None
if pd is not None:
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'pipeline'))
    import review_fda_candidates as review


@unittest.skipIf(pd is None,'Requires chemical environment')
class AssayReview(unittest.TestCase):
    def row(self, **overrides):
        r = dict(standard_type='IC50', standard_value=100, standard_units='nM',
                 standard_relation='=', assay_description='Inhibition of human target',
                 reported_action_type=None, data_validity_comment=None, potential_duplicate=0,
                 activity_comment=None, ingredient_id='C1', gene='G1', document_chembl_id='D1')
        r.update(overrides)
        return r

    def test_censored_potencies_do_not_become_exact_values(self):
        for value,rel,expected in [(100,'>','ambiguous'),(30000,'<','ambiguous'),
                (1000,'<=','active'),(10000,'>=','inactive_at_threshold')]:
            self.assertEqual(review.assay_triage(self.row(standard_value=value,standard_relation=rel))['activity_class'],expected)

    def test_units_zero_and_unknown_units(self):
        r=review.assay_triage(self.row(standard_value=.5,standard_units='µM'))
        self.assertEqual(r['value_nM'],500)
        self.assertEqual(r['activity_class'],'active')
        for kw in [{'standard_value':0},{'standard_units':'mg/mL'},{'standard_value':float('nan')}]:
            self.assertEqual(review.assay_triage(self.row(**kw))['activity_class'],'uninterpretable')

    def test_binding_is_not_activation_or_inhibition(self):
        r=review.assay_triage(self.row(standard_type='Kd'))
        self.assertEqual(r['activity_class'],'active')
        self.assertEqual(r['inferred_direction'],'unknown')

    def test_source_comment_conflict_and_duplicates_not_counted(self):
        for kw in [{'activity_comment':'Inactive'},{'potential_duplicate':1},
                   {'data_validity_comment':'Outside typical range'}]:
            r=review.assay_triage(self.row(**kw))
            self.assertEqual(r['activity_class'],'active')
            self.assertFalse(r['eligible_for_support_count'])

    def test_residual_percent_polarity_and_censoring(self):
        r=self.row(standard_type='Residual Activity',standard_units='%',standard_value=90,
                   standard_relation='>',assay_description='Residual activity of kinase at 10 uM')
        out=review.assay_triage(r)
        self.assertEqual(out['inhibition_percent'],10)
        self.assertEqual(out['activity_class'],'inactive_at_threshold')
        r.update(standard_relation='=',standard_value=10,assay_description='Residual activity of kinase at 1 uM')
        self.assertEqual(review.assay_triage(r)['activity_class'],'active')

    def test_percent_requires_unambiguous_dose_and_endpoint(self):
        for desc in ['Inhibition of target', 'Inhibition at 1 uM; control inhibitor at 10 uM',
                     'Increase in target phosphorylation at 1 uM']:
            r=self.row(standard_type='Inhibition',standard_units='%',standard_value=90,assay_description=desc)
            self.assertEqual(review.assay_triage(r)['activity_class'],'uninterpretable')

    def test_fold_change_cannot_be_used_as_potency(self):
        r=review.assay_triage(self.row(standard_type='FC',standard_units=None,standard_value=0))
        self.assertEqual(r['activity_class'],'uninterpretable')

    def test_mixed_evidence_retained(self):
        rows=[self.row(),self.row(standard_value=30000,standard_relation='>',document_chembl_id='D2')]
        frame=pd.DataFrame([{**r,**review.assay_triage(r)} for r in rows])
        summary=review.aggregate(frame)
        self.assertTrue(summary.iloc[0].mixed_active_weak_evidence)
        self.assertEqual(summary.iloc[0].n_active_unflagged,1)
        self.assertEqual(summary.iloc[0].n_weak_unflagged,1)
        tier,decision,_=review.default_decision(summary,set())
        self.assertEqual(decision,'defer_direction_conflict')

    def test_source_free_priority_not_allowed(self):
        r={'sources':{},'candidate_reviews':{'C1':{'source_ids':[]}}}
        with self.assertRaisesRegex(ValueError,'Missing evidence'):
            review.validate_review(r,{'C1'})


if __name__=='__main__':unittest.main()
