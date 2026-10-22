# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
import parachron
import unittest
import math
import random


def disc(n=2000, seed=3):
    rnd = random.Random(seed)
    sim = parachron.Simulation()
    sim.G = 1.
    sim.softening = 0.02
    sim.root_size = 20.
    sim.N_root_x = sim.N_root_y = sim.N_root_z = 1
    sim.add(m=1.)
    for i in range(n):
        a = rnd.uniform(1., 5.)
        phi = rnd.uniform(0., 2.*math.pi)
        sim.add(m=0.1/n, x=a*math.cos(phi), y=a*math.sin(phi), z=0.05*a*rnd.uniform(-0.5, 0.5),
                vx=-math.sin(phi)/math.sqrt(a), vy=math.cos(phi)/math.sqrt(a))
    return sim


def update(sim):
    from ctypes import byref
    parachron.clibrebound.reb_simulation_update_acceleration(byref(sim))
    return [(p.ax, p.ay, p.az) for p in sim.particles]


class TestTreeGroup(unittest.TestCase):

    def test_accuracy_vs_direct(self):
        sim = disc()
        sim.gravity = "basic"
        ad = update(sim)
        sim.gravity = "tree"
        sim.opening_angle2 = 0.5**2
        for G, tol_median in ((0, 2e-3), (32, 2e-3), (8, 2e-3)):
            sim.tree_group_size = G
            at = update(sim)
            err = sorted(math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))/math.sqrt(sum(y*y for y in b)) for a, b in zip(at[1:], ad[1:]))
            self.assertLess(err[len(err)//2], tol_median)
            self.assertLess(err[int(0.99*len(err))], 2e-2)

    def test_group_more_accurate_than_default_at_same_theta(self):
        sim = disc()
        sim.gravity = "basic"
        ad = update(sim)
        sim.gravity = "tree"
        sim.opening_angle2 = 0.5**2
        med = {}
        for G in (0, 32):
            sim.tree_group_size = G
            at = update(sim)
            err = sorted(math.sqrt(sum((x-y)**2 for x, y in zip(a, b)))/math.sqrt(sum(y*y for y in b)) for a, b in zip(at[1:], ad[1:]))
            med[G] = err[len(err)//2]
        self.assertLess(med[32], med[0])

    def test_deterministic_and_default_unchanged(self):
        s1, s2, s3 = disc(500), disc(500), disc(500)
        for s in (s1, s2, s3):
            s.gravity = "tree"
            s.integrator = "leapfrog"
            s.dt = 1e-2
        s1.tree_group_size = 16
        s2.tree_group_size = 16
        s1.steps(20)
        s2.steps(20)
        self.assertEqual([p.x for p in s1.particles], [p.x for p in s2.particles])
        s3.steps(20)
        self.assertNotEqual([p.x for p in s1.particles], [p.x for p in s3.particles])

    def test_saved_in_binary(self):
        sim = disc(100)
        sim.gravity = "tree"
        sim.tree_group_size = 24
        sim2 = sim.copy()
        self.assertEqual(sim2.tree_group_size, 24)


if __name__ == "__main__":
    unittest.main()
