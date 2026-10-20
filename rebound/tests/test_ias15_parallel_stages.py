import unittest
import random
import os
import rebound


def make_sim(parallel_stages, N=40):
    random.seed(5)
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.integrator.parallel_stages = parallel_stages
    sim.add(m=1.)
    for _ in range(N-1):
        sim.add(m=1e-6, a=1.+random.random(), e=0.1*random.random(), inc=0.05*random.random(),
                Omega=6.28*random.random(), f=6.28*random.random(), primary=sim.particles[0])
    sim.move_to_com()
    return sim


class TestIAS15ParallelStages(unittest.TestCase):
    def test_same_solution(self):
        a = make_sim(0)
        b = make_sim(1)
        E0 = a.energy()
        a.integrate(30.)
        b.integrate(30.)
        self.assertEqual(a.steps_done, b.steps_done)
        for pa, pb in zip(a.particles, b.particles):
            self.assertAlmostEqual(pa.x, pb.x, delta=1e-12)
            self.assertAlmostEqual(pa.vy, pb.vy, delta=1e-12)
        self.assertLess(abs((b.energy()-E0)/E0), 1e-13)

    def test_reproducible_after_restart(self):
        sim = make_sim(1)
        sim.integrate(5.)
        fn = "test_ias15_parallel_stages.bin"
        sim.save_to_file(fn, delete_file=True)
        sim.integrate(10.)
        sim2 = rebound.Simulation(fn)
        os.remove(fn)
        self.assertEqual(sim2.integrator.parallel_stages, 1)
        sim2.integrate(10.)
        for p1, p2 in zip(sim.particles, sim2.particles):
            self.assertEqual(p1.x, p2.x)
            self.assertEqual(p1.vz, p2.vz)

    def _compare(self, setup, tmax, tol=1e-11):
        a = rebound.Simulation(); b = rebound.Simulation()
        for sim, ps in ((a, 0), (b, 1)):
            random.seed(11)
            sim.integrator = "ias15"
            sim.integrator.parallel_stages = ps
            setup(sim)
            sim.integrate(tmax)
        self.assertEqual(a.steps_done, b.steps_done)
        for pa, pb in zip(a.particles, b.particles):
            self.assertAlmostEqual(pa.x, pb.x, delta=tol)
            self.assertAlmostEqual(pa.vy, pb.vy, delta=tol)

    def _planets_and_testparticles(self, sim, n_active, n_tp, tptype):
        sim.testparticle_hidewarnings = 1
        sim.add(m=1.)
        for _ in range(n_active-1):
            sim.add(m=1e-6, a=1.+random.random(), e=0.05*random.random(), f=6.28*random.random(), primary=sim.particles[0])
        for _ in range(n_tp):
            sim.add(m=1e-10, a=1.+2.*random.random(), e=0.05*random.random(), f=6.28*random.random(), primary=sim.particles[0])
        sim.N_active = n_active
        sim.testparticle_type = tptype
        sim.move_to_com()

    def test_testparticles_type0(self):
        self._compare(lambda s: self._planets_and_testparticles(s, 20, 150, 0), 3.)

    def test_testparticles_type1(self):
        self._compare(lambda s: self._planets_and_testparticles(s, 20, 150, 1), 3.)

    def test_compensated_gravity(self):
        def setup(s):
            s.gravity = "compensated"
            self._planets_and_testparticles(s, 30, 60, 1)
        self._compare(setup, 3.)

    def test_ghost_boxes(self):
        def setup(s):
            s.root_size = 1.
            s.N_root_x = 1; s.N_root_y = 1; s.N_root_z = 1
            s.boundary = "periodic"
            s.N_ghost_x = 1; s.N_ghost_y = 1; s.N_ghost_z = 0
            s.softening = 0.02
            for _ in range(80):
                s.add(m=1e-4*(0.5+random.random()), x=random.random()-0.5, y=random.random()-0.5, z=0.05*(random.random()-0.5),
                      vx=0.01*(random.random()-0.5), vy=0.01*(random.random()-0.5))
        self._compare(setup, 0.5, tol=1e-10)

    def test_unsupported_configuration_falls_back(self):
        # additional forces: the flag is ignored, results equal the default mode bitwise
        def drag(simp):
            for p in simp.contents.particles[1:simp.contents.N]:
                p.ax -= 1e-6*p.vx
        a = make_sim(0, N=10)
        b = make_sim(1, N=10)
        a.additional_forces = drag
        b.additional_forces = drag
        a.integrate(3.)
        b.integrate(3.)
        for pa, pb in zip(a.particles, b.particles):
            self.assertEqual(pa.x, pb.x)


if __name__ == "__main__":
    unittest.main()
