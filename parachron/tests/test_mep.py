# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
import parachron
import unittest
import os
import sys
import tempfile
import math


def solar_system():
    sim = parachron.Simulation()
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
        mep = parachron.MEP(sim, **self.kw)
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
        a = parachron.MEP(solar_system(), checkpoint=fa, **self.kw)
        a.integrate(10 * a.window_length)
        b = parachron.MEP(solar_system(), checkpoint=fb, **self.kw)
        b.integrate(4 * b.window_length)
        del b
        c = parachron.MEP.from_checkpoint(fb)
        self.assertEqual(c.windows, 4)
        c.integrate(10 * c.window_length)
        self.assertEqual(c.windows, 10)
        self.assertEqual(a.state, c.state)
        self.assertEqual(a.t, c.t)

    def test_calibration_resume_bitwise(self):
        # the calibration fit (from the stored per-window defect sums) survives checkpoint and resume bitwise
        # (4 slices: windows accepted at k = 2 have corrected slices n >= 2, whose defects are fitted)
        kw = dict(self.kw, cores=4, k_min=2, tol=3e-4, calibration=3, coarse_giants_steps=30)
        fa = os.path.join(self.dir.name, "ca.mep")
        fb = os.path.join(self.dir.name, "cb.mep")
        a = parachron.MEP(solar_system(), checkpoint=fa, **kw)
        a.integrate(10 * a.window_length)
        self.assertTrue(0. < a.calibration_residual < 1.)
        b = parachron.MEP(solar_system(), checkpoint=fb, **kw)
        b.integrate(5 * b.window_length)
        del b
        c = parachron.MEP.from_checkpoint(fb)
        c.integrate(10 * c.window_length)
        self.assertEqual(a.state, c.state)
        self.assertEqual(a.calibration_residual, c.calibration_residual)
        # the calibration changes the coarse step, so the trajectory differs at the tolerance level, not more
        d = parachron.MEP(solar_system(), **dict(kw, calibration=0))
        d.integrate(10 * d.window_length)
        self.assertNotEqual(a.state, d.state)
        for i in range(1, 9):
            self.assertLess(math.dist(a.sim.particles[i].xyz, d.sim.particles[i].xyz), 1e-3)

    def test_boundaries_in_heartbeat(self):
        m = parachron.MEP(solar_system(), **self.kw)
        rec = []
        m.heartbeat = lambda mm: rec.append((mm.boundaries, mm.state))
        m.integrate(3 * m.window_length)
        self.assertEqual(len(rec), 3)
        for b, st in rec:
            self.assertEqual(len(b), m.cores)
            self.assertEqual(b[-1], st)
        self.assertIsNone(m.boundaries)

    def test_checkpoint_version1(self):
        # a checkpoint of format version 1 (Parachron 1.0: no calibration options, no calibration block) is read
        import struct, ctypes
        f2 = os.path.join(self.dir.name, "v2.mep")
        f1 = os.path.join(self.dir.name, "v1.mep")
        a = parachron.MEP(solar_system(), checkpoint=f2, **self.kw)
        a.integrate(4 * a.window_length)
        raw = open(f2, "rb").read()
        nopt_v2 = ctypes.sizeof(parachron.MEPOptions)
        nopt_v1 = parachron.MEPOptions.calibration.offset
        head, opt, rest = raw[:36], raw[36:36 + nopt_v2], raw[36 + nopt_v2:]
        self.assertEqual(struct.unpack("<i", rest[-4:])[0], 0)       # no calibration block
        open(f1, "wb").write(head[:32] + struct.pack("<i", 1) + opt[:nopt_v1] + rest[:-4])
        c = parachron.MEP.from_checkpoint(f1)
        self.assertEqual(c.windows, 4)
        self.assertEqual(c.state, a.state)
        a.integrate(6 * a.window_length)
        c.integrate(6 * c.window_length)
        self.assertEqual(a.state, c.state)

    def test_heartbeat(self):
        # one call with a heartbeat after every window gives the same states as one call per window
        a = parachron.MEP(solar_system(), **self.kw)
        rec = []
        a.heartbeat = lambda m: rec.append((m.windows, m.t, m.state))
        a.integrate(6 * a.window_length)
        b = parachron.MEP(solar_system(), **self.kw)
        ref = []
        for w in range(1, 7):
            b.integrate(w * b.window_length)
            ref.append((b.windows, b.t, b.state))
        self.assertEqual(rec, ref)
        a.heartbeat = None
        a.integrate(7 * a.window_length)
        self.assertEqual(len(rec), 6)

    def test_partial_window_not_integrated(self):
        mep = parachron.MEP(solar_system(), **self.kw)
        mep.integrate(0.5 * mep.window_length)
        self.assertEqual(mep.windows, 0)
        mep.integrate(2.5 * mep.window_length)
        self.assertEqual(mep.windows, 2)

    def test_requirements(self):
        sim = parachron.Simulation()
        sim.add(m=1.)
        sim.add(m=1e-3, a=1.)
        sim.integrator = "whfast512"
        with self.assertRaises(RuntimeError):
            parachron.MEP(sim, **self.kw)

    def test_bad_checkpoint(self):
        fn = os.path.join(self.dir.name, "junk.mep")
        with open(fn, "wb") as f:
            f.write(b"not a checkpoint")
        with self.assertRaises(RuntimeError):
            parachron.MEP.from_checkpoint(fn)


if __name__ == "__main__":
    unittest.main()
