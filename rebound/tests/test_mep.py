import rebound
import unittest
import os
import sys
import tempfile
import math


def solar_system():
    sim = rebound.Simulation()
    sim.add("solar system")
    sim.integrator = "whfast512"
    sim.dt = 6. / 365.25 * 2. * math.pi
    sim.exact_finish_time = False
    sim.integrator.gr_potential = 1
    sim.integrator.corrector = 17
    return sim


@unittest.skipUnless(sys.platform.startswith("linux"), "MEP is only available on Linux")
class TestMEP(unittest.TestCase):
    # short slices keep the tests fast; production runs use the default slice_steps = 60000
    kw = dict(cores=2, slice_steps=600, helpers=1, pin=False)

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.dir.cleanup()

    def test_matches_serial(self):
        sim = solar_system()
        mep = rebound.MEP(sim, **self.kw)
        mep.integrate(8 * mep.window_length)
        self.assertEqual(mep.windows, 8)
        self.assertAlmostEqual(sim.t, 8 * mep.window_length, delta=1e-9 * sim.t)
        ref = solar_system()
        ref.integrator.jacobi_recursive = 1
        ref.integrator.concatenate_steps = 600
        for _ in range(16):
            ref.steps(1)
        for i in range(1, 9):
            d = math.dist(sim.particles[i].xyz, ref.particles[i].xyz)
            self.assertLess(d, 1e-6)

    def test_resume_bitwise(self):
        fa = os.path.join(self.dir.name, "a.mep")
        fb = os.path.join(self.dir.name, "b.mep")
        a = rebound.MEP(solar_system(), checkpoint=fa, **self.kw)
        a.integrate(10 * a.window_length)
        b = rebound.MEP(solar_system(), checkpoint=fb, **self.kw)
        b.integrate(4 * b.window_length)
        del b
        c = rebound.MEP.from_checkpoint(fb)
        self.assertEqual(c.windows, 4)
        c.integrate(10 * c.window_length)
        self.assertEqual(c.windows, 10)
        self.assertEqual(a.state, c.state)
        self.assertEqual(a.t, c.t)

    def test_partial_window_not_integrated(self):
        mep = rebound.MEP(solar_system(), **self.kw)
        mep.integrate(0.5 * mep.window_length)
        self.assertEqual(mep.windows, 0)
        mep.integrate(2.5 * mep.window_length)
        self.assertEqual(mep.windows, 2)

    def test_requirements(self):
        sim = rebound.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.integrator = "whfast512"
        with self.assertRaises(RuntimeError):
            rebound.MEP(sim, **self.kw)

    def test_bad_checkpoint(self):
        fn = os.path.join(self.dir.name, "junk.mep")
        with open(fn, "wb") as f:
            f.write(b"not a checkpoint")
        with self.assertRaises(RuntimeError):
            rebound.MEP.from_checkpoint(fn)


if __name__ == "__main__":
    unittest.main()
