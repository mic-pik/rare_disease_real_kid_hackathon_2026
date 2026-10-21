import unittest
import numpy as np
import pandas as pd
from pipeline.check_lincs_followup import rank_unit, reagent_agreement, aggregate


class FollowupTests(unittest.TestCase):
    def test_opposing_reagents_fail_agreement(self):
        self.assertAlmostEqual(reagent_agreement(np.array([[1,2,3],[3,2,1]])),-1)
        self.assertTrue(np.isnan(reagent_agreement(np.array([[1,2,3]]))))

    def test_constant_signature_is_not_reversal(self):
        np.testing.assert_array_equal(rank_unit([[1,1,1]]),[[0,0,0]])

    def test_cells_not_weighted_by_number_of_profiles(self):
        frame=pd.DataFrame({'ingredient_id':['a']*3,'canonical_name':['A']*3,'cell_id':['X','X','Y'],
            'pert_time':[6]*3,'pert_time_unit':['h']*3,'pert_dose':[1,10,1],'pert_dose_unit':['uM']*3})
        a,_,_=aggregate(frame,pd.DataFrame({'score':[1,1,-1]}))
        self.assertEqual(a.score.iloc[0],0)
