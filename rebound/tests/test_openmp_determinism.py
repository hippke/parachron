# Results must not depend on the number of OpenMP threads (bitwise). These tests only run against
# an OpenMP build of the library (REBOUND_OPENMP=1 pip install ...); in a serial build they are skipped.
import unittest
import math
import random
import rebound
from rebound.simulation import clibrebound

OPENMP = hasattr(clibrebound, "omp_get_max_threads")


def state(sim):
    return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]


@unittest.skipUnless(OPENMP, "needs an OpenMP build")
class TestOpenMPDeterminism(unittest.TestCase):
    def run_threads(self, make, steps, threads=(1, 2, 3)):
        results = []
        for t in threads:
            rebound.omp_set_num_threads(t)
            sim = make()
            sim.steps(steps)
            results.append(state(sim))
        rebound.omp_set_num_threads(1)
        for t, r in zip(threads[1:], results[1:]):
            ndiff = sum(1 for a, b in zip(results[0], r) if a != b)
            self.assertTrue(ndiff == 0 and len(r) == len(results[0]),
                            "%d of %d particles differ between 1 and %d threads" % (ndiff, len(r), t))

    def test_tree_collisions_shearing_sheet(self):
        def make():
            random.seed(4)
            sim = rebound.Simulation()
            sim.rand_seed = 7                  # collision resolution order is shuffled with it
            sim.integrator = "sei"
            sim.OMEGA = 0.00013143527
            sim.boundary = "shear"
            sim.gravity = "tree"
            sim.collision = "tree"
            sim.collision_resolve = "hardsphere"
            sim.G = 6.67428e-11
            sim.softening = 0.1
            sim.dt = 1e-3*2.*math.pi/0.00013143527
            sim.opening_angle2 = 0.5
            sim.root_size = 50.
            sim.N_root_x = 2; sim.N_root_y = 2; sim.N_root_z = 1
            sim.N_ghost_x = 1; sim.N_ghost_y = 1; sim.N_ghost_z = 0
            for _ in range(1500):
                x = 100.*(random.random()-0.5); y = 100.*(random.random()-0.5)
                r = 1.+random.random()
                sim.add(m=400.*4./3.*math.pi*r**3, r=r, x=x, y=y, z=random.gauss(0., 1.), vy=-1.5*x*0.00013143527)
            return sim
        self.run_threads(make, 5)

    def test_ias15_stage_mode(self):
        def make():
            random.seed(5)
            sim = rebound.Simulation()
            sim.integrator = "ias15"
            sim.integrator.parallel_stages = 1
            sim.add(m=1.)
            for _ in range(299):
                sim.add(m=1e-7, a=1.+random.random(), e=0.02*random.random(), f=6.28*random.random(), primary=sim.particles[0])
            sim.move_to_com()
            return sim
        self.run_threads(make, 3)

    def test_block_order_gravity(self):
        def make():
            random.seed(8)
            sim = rebound.Simulation()
            sim.integrator = "whfast"
            sim.dt = 0.01
            sim.gravity_block_order = 1
            sim.testparticle_hidewarnings = 1
            sim.add(m=1.)
            for _ in range(599):
                sim.add(m=1e-8, a=1.+2.*random.random(), e=0.02*random.random(), f=6.28*random.random(), primary=sim.particles[0])
            sim.N_active = 300
            sim.testparticle_type = 1
            sim.move_to_com()
            return sim
        self.run_threads(make, 5)

    def test_parallel_basic_gravity(self):
        def make():
            random.seed(6)
            sim = rebound.Simulation()
            sim.integrator = "whfast"
            sim.dt = 0.01
            sim.add(m=1.)
            for _ in range(799):
                sim.add(m=1e-8, a=1.+2.*random.random(), e=0.02*random.random(), f=6.28*random.random(), primary=sim.particles[0])
            sim.move_to_com()
            return sim
        self.run_threads(make, 5)


if __name__ == "__main__":
    unittest.main()
