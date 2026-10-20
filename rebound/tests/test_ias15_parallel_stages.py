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
