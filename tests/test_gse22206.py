import unittest

import pandas as pd

from pipeline.analyze_gse22206 import contrasts, robust_direction


class ContrastTests(unittest.TestCase):
    def test_tissue_baseline_and_relatives_do_not_create_case_effect(self):
        x = pd.DataFrame([[10, 10, 100, 100, 10, 10, 100, 100, 10, 10, 1000, 1000, -100, -100]])
        c = contrasts(x)
        self.assertEqual(c.fibro_delta.iloc[0], 0)
        self.assertEqual(c.lcl_delta.iloc[0], 0)
        self.assertEqual(robust_direction(c).iloc[0], 0)

    def test_conflicting_control_donor_prevents_consensus(self):
        x = pd.DataFrame([[5, 5, 5, 5, 1, 1, 1, 1, 9, 9, 0, 0, 0, 0]])
        self.assertEqual(robust_direction(contrasts(x)).iloc[0], 0)

    def test_replicate_overlap_prevents_consensus(self):
        x = pd.DataFrame([[2, 10, 10, 10, 3, 3, 1, 1, 3, 3, 0, 0, 0, 0]])
        c = contrasts(x)
        self.assertGreater(c.lcl_delta.iloc[0], 0)
        self.assertEqual(robust_direction(c).iloc[0], 0)


if __name__ == "__main__":
    unittest.main()
