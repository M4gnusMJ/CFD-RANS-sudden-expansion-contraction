"""Run with: python3 -m unittest discover -s postprocess -p 'test_*.py'."""
import tempfile
import unittest
from pathlib import Path

import numpy as np

from gci import LEVELS, pressure_profile_table


class PressureGCITest(unittest.TestCase):
    def write_profiles(self, directory, values, coordinates=None):
        x = np.linspace(0, 1, len(values[0]))
        for i, level in enumerate(LEVELS):
            np.savetxt(Path(directory) / f"{level}_axis.csv",
                       np.column_stack((x if coordinates is None else coordinates[i], values[i])),
                       delimiter=",")

    def test_unequal_refinement_second_order_and_pressure_datum(self):
        cells = np.array([100, 400, 2500])
        h = np.sqrt(1 / cells)
        exact = np.array([0., 2., 3.])
        values = exact + h[:, None] ** 2
        # A fine-grid zero must retain a finite absolute error bar.
        values[:, 0] -= h[-1] ** 2
        with tempfile.TemporaryDirectory() as directory:
            self.write_profiles(directory, values)
            table = pressure_profile_table(directory, cells, 1.)
            np.testing.assert_allclose(table.p_local, 2., atol=1e-8)
            np.testing.assert_allclose(table["GCI absolute [m2/s2]"], 1.25 * h[-1] ** 2)
            self.assertTrue(np.isnan(table["GCI relative [%]"].iloc[0]))
            self.write_profiles(directory, values + 10)
            shifted = pressure_profile_table(directory, cells, 1.)
            np.testing.assert_allclose(shifted["GCI absolute [m2/s2]"],
                                       table["GCI absolute [m2/s2]"])

    def test_oscillatory_diverging_and_degenerate_points(self):
        values = np.array([[16., 1., 2.], [-4., 4., 2.], [1., 16., 2.]])
        with tempfile.TemporaryDirectory() as directory:
            self.write_profiles(directory, values)
            table = pressure_profile_table(directory, [1, 4, 16], 1.)
            np.testing.assert_allclose(table.p_local[:2], 2.)
            np.testing.assert_allclose(table.p_average, 2.)
            self.assertEqual(list(table["grid trend"]),
                             ["Oscillatory", "Non-convergent", "Undetermined"])
            np.testing.assert_allclose(table["GCI absolute [m2/s2]"], [1.25 * 5 / 3, 5, 0])

    def test_reject_mismatched_samples_and_unresolved_mean(self):
        with tempfile.TemporaryDirectory() as directory:
            values = np.ones((3, 3))
            self.write_profiles(directory, values, [[0, .5, 1], [0, .6, 1], [0, .5, 1]])
            with self.assertRaisesRegex(ValueError, "same x coordinates"):
                pressure_profile_table(directory, [1, 4, 16], 1.)
            self.write_profiles(directory, values)
            with self.assertRaisesRegex(ValueError, "No resolved local"):
                pressure_profile_table(directory, [1, 4, 16], 1.)


if __name__ == "__main__":
    unittest.main()
