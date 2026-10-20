import unittest
import os
import random
import rebound


def disc(block_order, n=300, n_active=None, tptype=0, integrator="whfast"):
    random.seed(9)
    sim = rebound.Simulation()
    sim.integrator = integrator
    sim.dt = 0.01
    sim.gravity_block_order = block_order
    sim.testparticle_hidewarnings = 1
    sim.add(m=1.)
    for i in range(n-1):
        sim.add(m=1e-8, a=1.+2.*random.random(), e=0.02*random.random(), inc=0.01*random.random(), f=6.28*random.random(), primary=sim.particles[0])
    if n_active:
        sim.N_active = n_active
        sim.testparticle_type = tptype
    sim.move_to_com()
    return sim


class TestGravityBlockOrder(unittest.TestCase):
    def compare(self, a, b, tol):
        for pa, pb in zip(a.particles, b.particles):
            self.assertAlmostEqual(pa.x, pb.x, delta=tol)
            self.assertAlmostEqual(pa.vy, pb.vy, delta=tol)

    def test_same_physics_whfast(self):
        a, b = disc(0), disc(1)
        a.integrate(5.); b.integrate(5.)
        self.compare(a, b, 1e-12)

    def test_same_physics_ias15(self):
        a, b = disc(0, n=150, integrator="ias15"), disc(1, n=150, integrator="ias15")
        a.integrate(1.); b.integrate(1.)
        self.assertEqual(a.steps_done, b.steps_done)
        self.compare(a, b, 1e-12)

    def test_testparticles(self):
        for tptype in (0, 1):
            a, b = disc(0, n=400, n_active=40, tptype=tptype), disc(1, n=400, n_active=40, tptype=tptype)
            a.integrate(3.); b.integrate(3.)
            self.compare(a, b, 1e-12)

    def test_saved_in_archive(self):
        sim = disc(1, n=50)
        sim.integrate(1.)
        fn = "test_gravity_block_order.bin"
        sim.save_to_file(fn, delete_file=True)
        sim.integrate(2.)
        sim2 = rebound.Simulation(fn)
        os.remove(fn)
        self.assertEqual(sim2.gravity_block_order, 1)
        sim2.integrate(2.)
        for p1, p2 in zip(sim.particles, sim2.particles):
            self.assertEqual(p1.x, p2.x)


if __name__ == "__main__":
    unittest.main()
