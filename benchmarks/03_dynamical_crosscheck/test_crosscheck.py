"""Fast reproducibility tests. Run: python -m unittest -v test_crosscheck.py."""
import math
import unittest
import numpy as np
from oscillator_crosscheck import (Case, Config, algebra_checks,
    bounded_beta_interval, drift_analysis, hamiltonian, packet, trajectory)

class CrossCheckTests(unittest.TestCase):
    def test_symbolic_identifiability(self):
        result=algebra_checks()
        self.assertEqual(result['rank_full'],3)
        self.assertEqual(result['rank_kappa_known'],3)
        self.assertEqual(result['null_direction'],[1,-1,-1,0])

    def test_deformation_and_complete_mimic_same_generator(self):
        a=Case('signal',beta=.0025)
        b=Case('mimic',eta=.0025,kappa=.0025)
        ha,_,_=hamiltonian(a,48)
        hb,_,_=hamiltonian(b,48)
        np.testing.assert_array_equal(ha,hb)
        self.assertEqual(a.loop,b.loop)

    def test_initial_slope_is_not_spring_response(self):
        times=np.linspace(0,.6,9)
        a=trajectory(Case('spring',lam=.01),48,2.,times)
        self.assertAlmostEqual(a['initial_slope'],2.,places=11)
        b=trajectory(Case('kinetic',beta=.0025,lam=.01),48,2.,times)
        self.assertAlmostEqual(b['initial_slope'],2.+4*.0025/3*(8+3),places=11)

    def test_window_extrapolation(self):
        c=Config(dimension=48)
        for case in [Case('null'),Case('signal',beta=.0025),Case('spring',lam=.01)]:
            result=drift_analysis(case,c,tmax=.4)
            self.assertLess(abs(result['bias']),1e-6)

    def test_no_calibration_means_no_attribution_interval(self):
        r=bounded_beta_interval(.0025,.0025,.0005,.0006,math.inf,math.inf)
        self.assertTrue(r['lower_unbounded'] and r['upper_unbounded'])
        self.assertFalse(r['zero_excluded'])

    def test_certificate_is_conditional_not_automatic(self):
        tight=bounded_beta_interval(.0025,.0025,.0005,.0006,.0001,.0001)
        loose=bounded_beta_interval(.0025,.0025,.0005,.0006,.003,.003)
        self.assertTrue(tight['zero_excluded'])
        self.assertFalse(loose['zero_excluded'])
        mismatch=bounded_beta_interval(.1,0,.0001,.0001,0,0)
        self.assertTrue(mismatch['empty'])
        self.assertFalse(mismatch['zero_excluded'])

if __name__=='__main__':
    unittest.main()
