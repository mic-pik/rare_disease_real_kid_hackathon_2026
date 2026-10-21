import unittest
import numpy as np
import pandas as pd
from pipeline.score_gse22206_lincs import weights, bh, background_check


class ReversalTests(unittest.TestCase):
    def test_reversal_sign_and_balanced_directions(self):
        w=weights([1,1,1,-1])
        self.assertAlmostEqual(w[:3].sum(),-.5)
        self.assertAlmostEqual(w[3],.5)
        self.assertAlmostEqual(np.array([-1,-1,-1,1])@w,1)
        self.assertAlmostEqual(np.ones(4)@w,0)

    def test_single_direction_rejected(self):
        with self.assertRaises(ValueError): weights([1,1])

    def test_bh_retains_original_order(self):
        np.testing.assert_allclose(bh([.04,.001,.03]),[.04,.003,.04])

    def test_background_excludes_self_and_requires_matched_conditions(self):
        d=pd.DataFrame({'ingredient_id':['a','b','c'], 'primary':[1.,0.,100.],
                        'cell_id':['X','X','Y'], 'pert_time':[6,6,6],
                        'pert_time_unit':['h']*3,'pert_dose':[10]*3,'pert_dose_unit':['uM']*3})
        out=background_check(d,minimum=1).set_index('ingredient_id')
        self.assertEqual(out.loc['a','background_adjusted'],1.)
        self.assertEqual(out.loc['b','background_adjusted'],-1.)
        self.assertTrue(np.isnan(out.loc['c','background_adjusted']))
