import parachron
from ctypes import c_double, byref
import unittest
from time import sleep

class TestRandom(unittest.TestCase):
    def test_uniform(self):
        parachron.clibrebound.reb_random_uniform.restype = c_double
        simulation = parachron.Simulation()
        for simp in [0, byref(simulation)]:
            for i in range(10):
                r = parachron.clibrebound.reb_random_uniform(simp, c_double(0.0), c_double(1.0))
                self.assertLess(r, 1.0)
                self.assertGreater(r, 0.0)
    
    def test_powerlaw(self):
        parachron.clibrebound.reb_random_powerlaw.restype = c_double
        simulation = parachron.Simulation()
        for simp in [0, byref(simulation)]:
            for slope in [-1.0, 0.0, 1.0]:
                for i in range(10):
                    r = parachron.clibrebound.reb_random_powerlaw(simp, c_double(1.0), c_double(2.0), c_double(slope))
                    self.assertLess(r, 2.0)
                    self.assertGreater(r, 1.0)
    
    def test_normal(self):
        parachron.clibrebound.reb_random_normal.restype = c_double
        simulation = parachron.Simulation()
        for simp in [0, byref(simulation)]:
            for i in range(10):
                r = parachron.clibrebound.reb_random_normal(simp, c_double(1.0))
                self.assertLess(r, 1e5)
                self.assertGreater(r, -1e5)
    
    def test_rayleigh(self):
        parachron.clibrebound.reb_random_rayleigh.restype = c_double
        simulation = parachron.Simulation()
        for simp in [0, byref(simulation)]:
            for i in range(10):
                r = parachron.clibrebound.reb_random_rayleigh(simp, c_double(1.0))
                self.assertLess(r, 1e5)
                self.assertGreater(r, 0.0)

    def test_reproducible(self):
        parachron.clibrebound.reb_random_uniform.restype = c_double
        sim1 = parachron.Simulation()
        sim2 = sim1.copy()
        sleep(0.05) # Windows implementation is very slow, so we wait to get a different seed
        sim3 = parachron.Simulation()
        self.assertEqual(sim1.rand_seed, sim2.rand_seed)
        self.assertNotEqual(sim2.rand_seed, sim3.rand_seed)
        r1 = parachron.clibrebound.reb_random_uniform(byref(sim1), c_double(0.0), c_double(1.0))
        r2 = parachron.clibrebound.reb_random_uniform(byref(sim2), c_double(0.0), c_double(1.0))
        r3 = parachron.clibrebound.reb_random_uniform(byref(sim3), c_double(0.0), c_double(1.0))
        self.assertEqual(r1, r2)
        self.assertNotEqual(r2, r3)
        r2 = parachron.clibrebound.reb_random_uniform(byref(sim2), c_double(0.0), c_double(1.0))
        r3 = parachron.clibrebound.reb_random_uniform(byref(sim3), c_double(0.0), c_double(1.0))
        self.assertNotEqual(r2, r3)
    
if __name__ == "__main__":
    unittest.main()
