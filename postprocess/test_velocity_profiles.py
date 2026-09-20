import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from sample_wall_profiles import latest, metadata, dictionary
from velocity_profiles import scale_profile, profile_change


class VelocityDiagnosticsTest(unittest.TestCase):
    def test_wall_scaling_from_kinematic_shear(self):
        # Exact linear viscous profile: nu=.01, u_tau=.2, U=4*distance.
        radius = np.array([.2, .5, .9])
        p = pd.DataFrame({"r": radius, "Ux": 4 * (1 - radius), "nut": [0., 0., 0.]})
        wall = pd.DataFrame({"x": [0., 2.], "r": [1., 1.],
                             "tau_x_over_rho": [-.04, -.04], "wall_yplus": [.5, .5]})
        station = {"x": 1., "wall_radius": 1., "bulk_velocity": 1.}
        result, yp, shear = scale_profile(p, wall, station, .01)
        np.testing.assert_allclose(result.u_tau, .2)
        np.testing.assert_allclose(result.uplus, result.yplus)
        self.assertEqual(yp, .5)
        self.assertEqual(shear, -.04)
        with self.assertRaisesRegex(ValueError, "bracket"):
            scale_profile(p, wall, dict(station, x=3.), .01)
        with self.assertRaisesRegex(ValueError, "nonzero wall shear"):
            scale_profile(p, wall.assign(tau_x_over_rho=0.), station, .01)

    def test_area_weighted_development_difference(self):
        eta = np.linspace(0, 1, 1001)[1:-1]
        p = pd.DataFrame({"r_over_R": eta, "U_over_Ubulk": 2 * (1 - eta ** 2)})
        self.assertEqual(profile_change(p, p), (0., 0.))
        q = p.assign(U_over_Ubulk=1.1 * p.U_over_Ubulk)
        rms, maximum = profile_change(p, q)
        self.assertAlmostEqual(rms, 20 / np.sqrt(3), places=3)
        self.assertAlmostEqual(maximum, 20., places=3)

    def test_numeric_latest_and_requested_stations(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("9", "100", "20", "system"):
                (root / name).mkdir()
            self.assertEqual(latest(root).name, "100")
        meta = metadata()
        self.assertAlmostEqual(meta["stations"]["expansion_3D"]["x"], 2.7)
        self.assertAlmostEqual(meta["stations"]["contraction_3D"]["x"], 8.4)
        self.assertEqual(len(meta["stations"]), 6)
        self.assertIn("cellPointFace", dictionary(meta))


if __name__ == "__main__":
    unittest.main()
