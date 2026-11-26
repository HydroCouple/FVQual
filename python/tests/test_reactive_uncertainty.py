import unittest
from fvqual.reactive_uncertainty import gaussian_diagnostics


class GaussianDiagnosticsTests(unittest.TestCase):
    def test_identifiable_and_confounding_directions_with_independent_oracle(self):
        try:import numpy
        except ImportError:self.skipTest("Optional NumPy observation-design diagnostics")
        local={"mean":[2.,3.],"mode_sensitivities":[[2.,0.],[0.,3.]]}
        result=gaussian_diagnostics(local,[1.,1.],[3.,4.])
        self.assertEqual(result["numerical_rank"],2)
        self.assertAlmostEqual(result["posterior_mode_covariance"][0][0],1/5)
        self.assertAlmostEqual(result["posterior_mode_covariance"][1][1],1/10)
        self.assertAlmostEqual(result["posterior_mode_mean"][0],2/5)
        self.assertAlmostEqual(result["posterior_mode_mean"][1],3/10)
        confounded=gaussian_diagnostics({"mean":[0.,0.],"mode_sensitivities":[[1.,2.],[2.,4.]]},[1.,1.])
        self.assertEqual(confounded["numerical_rank"],1);self.assertIsNone(confounded["condition_number"])
        self.assertAlmostEqual(confounded["sensitivity_cosines"][0][1],1)
        with self.assertRaises(ValueError):gaussian_diagnostics(local,[0.,1.])
        with self.assertRaises(ValueError):gaussian_diagnostics(local,[1.,1.],[1.])
