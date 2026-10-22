import parachron
import unittest
import math
import warnings
    
def gr_potential(sim):
    C2 = 10065.32 * 10065.32
    source = sim.particles[0]
    mu = sim.G*source.m;
    prefac = 3.*mu*mu/C2
    H = 0.0

    for p in sim.particles[1:]:
        dx = p.x - source.x
        dy = p.y - source.y
        dz = p.z - source.z
        r2 = dx*dx + dy*dy + dz*dz
        H -= prefac*p.m/r2
    return H
    
class TestIntegratorWHFast512(unittest.TestCase):
    # Kernel under test: 0 picks the fastest one available. The subclasses at the end of this
    # file repeat all tests for each kernel (1: AVX-512 assembly, 2: C AVX-512, 3: C AVX2, 4: C generic).
    kernel = 0

    def kernel_available(self):
        try:
            sim = parachron.Simulation()
            sim.add("solar system")
            sim.integrator = "whfast512"
            sim.integrator.kernel = self.kernel
            sim.integrator.kernel = self.kernel
            sim.dt = 6/365.25*2*math.pi
            sim.exact_finish_time = 0
            sim.integrator.concatenate_steps = 1
            sim.steps(1)
        except RuntimeError:
            return False
        return True

    def require_kernel(self):
        if not self.kernel_available():
            self.skipTest("WHFast512 kernel %d is not available on this machine" % self.kernel)

    def test_whfast512_avx512_available(self):
        if not parachron.avx512_available:
            print("\n::notice:: AVX512 is not supported by this CPU. WHFast512 runs with its C kernels.")
        else:
            print("\n::notice:: AVX512 is supported by this CPU.")
        # The generic C kernel runs everywhere.
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.kernel = 4
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.integrator.concatenate_steps = 10
        sim.steps(1)
        self.assertEqual(sim.integrator.kernel_used, 4)

    def test_whfast512_basic(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy()
        sim.integrator.concatenate_steps = 1
        sim.steps(1)
        self.assertTrue(sim.is_synchronized)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),5e-11)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy()
        self.assertEqual(sim.t, sim.dt*(1+1e5))
        self.assertLess(math.fabs((e0-e1)/e0),2e-9)

    def test_whfast512_corrector(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.corrector = 17
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy()
        sim.integrator.concatenate_steps = 1
        sim.steps(1)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),6e-15)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy()
        self.assertLess(math.fabs((e0-e1)/e0),2e-12)

    def test_whfast512_gr(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.corrector = 17
        sim.integrator.gr_potential = 1
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy() + gr_potential(sim)
        sim.integrator.concatenate_steps = 1
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        # One step (corrector, step, inverse corrector) is round-off dominated: over 300 ICs
        # perturbed by 1e-9 AU, |dE/E| has median 1.7e-15 and max 6.4e-15 for the assembly kernel
        # (1 % above 6e-15); for the C kernels median 1.6e-15, max 6.2e-15 (0.3 %
        # above 6e-15). This particular IC gives 6.03e-15 with the C kernels, hence 7e-15 (was 6e-15).
        self.assertLess(math.fabs((e0-e1)/e0),7e-15)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),2e-12)

    
    def test_whfast512_fewer_than_8_particles(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("outer solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.corrector = 17
        sim.integrator.gr_potential = 1
        sim.dt = 20.0/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy() + gr_potential(sim)
        sim.integrator.concatenate_steps = 1
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),4e-14)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),1e-11)

    def test_whfast512_independent_kepler_solvers(self):
        self.require_kernel()
        def getSim(innera):
            sim = parachron.Simulation()
            sim.add(m=1)
            sim.add(a=innera, e=0.8)
            sim.add(a=1.2, e=0.01)
            sim.integrator = "whfast512"
            sim.integrator.kernel = self.kernel
            sim.dt = 6.0/365.25*2*math.pi
            sim.exact_finish_time = 0
            # 10^4 steps (not the default 10^6): the inner test particles (e = 0.8, dt up to ~P/2)
            # are far outside the Kepler solver's accurate range (z = beta X^2 >> 1, where all
            # kernels, the assembly included, violate Kepler energy conservation); over 10^6 steps
            # such an orbit can run off to infinity, depending on round-off, and 0*inf in the
            # interaction step then gives NaN in every lane. Lane independence is tested either way.
            sim.integrator.concatenate_steps = 10000
            sim.steps(1)
            return sim
        sim1 = getSim(0.01)
        sim2 = getSim(0.1)
        sim3 = getSim(1.0)
        self.assertEqual(sim1.particles[2].x, sim2.particles[2].x)
        self.assertEqual(sim1.particles[2].x, sim3.particles[2].x)
        self.assertEqual(sim1.particles[2].vx, sim2.particles[2].vx)
        self.assertEqual(sim1.particles[2].vx, sim3.particles[2].vx)
        self.assertNotEqual(sim1.particles[1].x, sim2.particles[1].x)
        self.assertNotEqual(sim1.particles[1].x, sim3.particles[1].x)
        self.assertNotEqual(sim1.particles[1].vx, sim2.particles[1].vx)
        self.assertNotEqual(sim1.particles[1].vx, sim3.particles[1].vx)


    def test_whfast512_high_e_kepler_solver(self):
        self.require_kernel()
        def getSim(e, integrator, f):
            sim = parachron.Simulation()
            sim.add(m=1)
            sim.add(a=1, e=e, f=f, omega=0.53, inc=0.3)
            sim.integrator = integrator
            sim.dt = 6/365.25*2*math.pi
            sim.exact_finish_time = 0
            if sim.integrator=="whfast512":
                sim.integrator.kernel = self.kernel
                sim.integrator.concatenate_steps = 1
            sim.steps(1)
            return sim
        for f in [0,math.pi,0.1234,0.2355]:
            for i in range(100):
                e = i/100
                sim = getSim(e, "whfast", f)
                sim512 = getSim(e, "whfast512", f)
                self.assertAlmostEqual(sim.particles[1].x, sim512.particles[1].x, 15)

    def test_whfast512_Nsystems_2(self):
        self.require_kernel()
        for Nplanets in [1,2,3,4]:
            def getSim(Nplanets):
                sim = parachron.Simulation()
                sim.add(m=1)
                sim.add(m=1e-3, a=1, e=0.1, f=0.5, omega=0.53, inc=0.3)
                if Nplanets>1:
                    sim.add(m=1e-3, a=2.3, e=0.1, f=0.15, omega=1.53, inc=0.12)
                if Nplanets>2:
                    sim.add(m=1e-3, a=3.4, e=0.051, f=0.15, omega=2.53, inc=0.2)
                if Nplanets>3:
                    sim.add(m=1e-3, a=5.4, e=0.01, f=0.35, omega=3.53, inc=0.0)
                return sim
            
            sim = getSim(Nplanets)
            sim.integrator = "whfast512"
            sim.integrator.kernel = self.kernel
            sim.integrator.corrector = 1
            sim.integrator.gr_potential = 1
            sim.integrator.concatenate_steps = 1
            sim.dt = 6/365.25*2*math.pi
            sim.exact_finish_time = 0

            sim2 = sim.copy()
            sim2.integrator.N_systems = 2
            for p in sim.particles:
                sim2.add(p)

            sim.steps(1)
            sim2.steps(1)

            self.assertEqual(sim.particles[0].x, sim2.particles[0].x)
            self.assertEqual(sim.particles[0].x, sim2.particles[1+Nplanets].x)
            for i in range(Nplanets):
                self.assertEqual(sim.particles[1+i].x, sim2.particles[1+i].x)
                self.assertEqual(sim.particles[1+i].x, sim2.particles[2+Nplanets+i].x)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[1+i].vx)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[2+Nplanets+i].vx)

    def test_whfast512_Nsystems_4(self):
        self.require_kernel()
        for Nplanets in [1,2]:
            def getSim(Nplanets):
                sim = parachron.Simulation()
                sim.add(m=1)
                sim.add(m=1e-3, a=1, e=0.1, f=0.5, omega=0.53, inc=0.3)
                if Nplanets>1:
                    sim.add(m=1e-3, a=2.3, e=0.1, f=0.15, omega=1.53, inc=0.12)
                return sim
            
            sim = getSim(Nplanets)
            sim.integrator = "whfast512"
            sim.integrator.kernel = self.kernel
            sim.integrator.corrector = 1
            sim.integrator.gr_potential = 1
            sim.integrator.concatenate_steps = 1
            sim.dt = 6/365.25*2*math.pi
            sim.exact_finish_time = 0

            sim2 = sim.copy()
            sim2.integrator.N_systems = 4
            for i in range(3):
                for p in sim.particles:
                    sim2.add(p)

            sim.steps(1)
            sim2.steps(1)

            self.assertEqual(sim.particles[0].x, sim2.particles[0].x)
            self.assertEqual(sim.particles[0].x, sim2.particles[1+Nplanets].x)
            self.assertEqual(sim.particles[0].x, sim2.particles[2+2*Nplanets].x)
            self.assertEqual(sim.particles[0].x, sim2.particles[3+3*Nplanets].x)
            for i in range(Nplanets):
                self.assertEqual(sim.particles[1+i].x, sim2.particles[1+i].x)
                self.assertEqual(sim.particles[1+i].x, sim2.particles[2+Nplanets+i].x)
                self.assertEqual(sim.particles[1+i].x, sim2.particles[3+2*Nplanets+i].x)
                self.assertEqual(sim.particles[1+i].x, sim2.particles[4+3*Nplanets+i].x)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[1+i].vx)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[2+Nplanets+i].vx)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[3+2*Nplanets+i].vx)
                self.assertEqual(sim.particles[1+i].vx, sim2.particles[4+3*Nplanets+i].vx)

    def test_whfast512_com(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("solar system")
        for p in sim.particles:
            p.vx += 0.1
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.corrector = 17
        sim.integrator.gr_potential = 1
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        e0 = sim.energy() + gr_potential(sim)
        sim.integrator.concatenate_steps = 100000
        sim.steps(1)
        e1 = sim.energy() + gr_potential(sim)
        self.assertLess(math.fabs((e0-e1)/e0),4e-14)
        com = sim.com()
        self.assertAlmostEqual(com.vx, 0.1, delta=1e-16) # round-off level (a few ulp of 0.1)
        self.assertAlmostEqual(com.x, com.vx*sim.t, 12)


    def test_whfast512_restart(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.concatenate_steps = 5
        sim.integrator.gr_potential = True
        sim.integrator.corrector = 17
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.steps(3)
        sim.save_to_file("test.bin",delete_file=True)
        sim.steps(4)
        sim2 = parachron.Simulation("test.bin")
        sim2.steps(4)
      
        self.assertEqual(sim.integrator, sim2.integrator)
        self.assertEqual(sim.t, sim2.t)
        for i in range(sim.N):
            self.assertEqual(sim.particles[i].m, sim2.particles[i].m)
            self.assertEqual(sim.particles[i].x, sim2.particles[i].x)
            self.assertEqual(sim.particles[i].vx, sim2.particles[i].vx)
    
    def test_whfast512_exit_min_distance(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add(m=1)
        sim.add(m=1e-8,a=1)
        sim.add(m=1e-8,a=1.2,f=math.pi)
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.concatenate_steps = 1000
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.exit_min_distance = 0.5
        with self.assertRaises(parachron.Encounter):
            sim.steps(1)
        self.assertAlmostEqual(sim.t,11.456, 2)
    
    def test_whfast512_exit_min_distance_star(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add(m=1)
        sim.add(m=1e-8,a=1, e=0.9, f=math.pi)
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.concatenate_steps = 1000
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.exit_min_distance = 0.5
        with self.assertRaises(parachron.Encounter):
            sim.steps(1)
        self.assertAlmostEqual(sim.t,2.99322, 2)

    def test_whfast512_exit_max_distance(self):
        self.require_kernel()
        sim = parachron.Simulation()
        sim.add(m=1)
        sim.add(m=1e-8,a=1, e=0.9, f=0)
        sim.integrator = "whfast512"
        sim.integrator.kernel = self.kernel
        sim.integrator.concatenate_steps = 1000
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.exit_max_distance = 1
        with self.assertRaises(parachron.Escape):
            sim.steps(1)
        self.assertAlmostEqual(sim.t,0.722502, 2)
       
class TestIntegratorWHFast512Kernels(unittest.TestCase):
    def run_kernel(self, kernel, setup, N_systems):
        sim = setup()
        sim.integrator = "whfast512"
        sim.integrator.N_systems = N_systems
        sim.integrator.kernel = kernel
        sim.integrator.gr_potential = 1
        sim.integrator.corrector = 17
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        sim.integrator.concatenate_steps = 1000
        try:
            sim.steps(3)
        except RuntimeError:
            return None
        return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]

    def test_c_kernels_bitwise_identical(self):
        # The C kernels use the same operations (with explicit fused multiply-adds) on every
        # instruction set, so they agree bit for bit. The assembly agrees to round-off.
        def ss():
            sim = parachron.Simulation()
            sim.add("solar system")
            return sim
        def two_systems():
            sim = parachron.Simulation()
            for s in range(2):
                sim.add(m=1)
                sim.add(m=1e-3, a=1, e=0.1+0.05*s, f=0.5, omega=0.53, inc=0.3, primary=sim.particles[-1])
                sim.add(m=1e-4, a=2.3, e=0.1, f=0.15, omega=1.53, inc=0.12, primary=sim.particles[-2])
            return sim
        def three_planets():
            # One system with at most four planets: the AVX2 and generic kernels use their
            # four-lane variants, which must agree bit for bit with the eight-lane kernels.
            sim = parachron.Simulation()
            sim.add(m=1)
            sim.add(m=1e-5, a=1, e=0.05, inc=0.01)
            sim.add(m=2e-5, a=1.3, e=0.02, f=1.1, inc=0.02)
            sim.add(m=1e-5, a=1.7, e=0.08, f=2.2, inc=0.03)
            sim.move_to_com()
            return sim
        for setup, N_systems in [(ss, 1), (two_systems, 2), (three_planets, 1)]:
            results = {k: self.run_kernel(k, setup, N_systems) for k in [1, 2, 3, 4]}
            self.assertTrue(all(math.isfinite(u) for p in results[4] for u in p))
            ref = results[4]
            self.assertIsNotNone(ref)
            for k in [2, 3]:
                if results[k] is not None:
                    self.assertEqual(results[k], ref)
            if results[1] is not None:
                for a, b in zip(results[1], ref):
                    for u, v in zip(a, b):
                        self.assertAlmostEqual(u, v, delta=1e-9*max(1., abs(v)))

    def test_unavailable_kernel_raises(self):
        sim = parachron.Simulation()
        sim.add("solar system")
        sim.integrator = "whfast512"
        sim.integrator.kernel = 99
        sim.dt = 6/365.25*2*math.pi
        sim.exact_finish_time = 0
        with self.assertRaises(RuntimeError):
            sim.steps(1)

class TestIntegratorWHFast512StepsPair(unittest.TestCase):
    def make(self, kind, kernel, variant=0, **kw):
        sim = parachron.Simulation()
        if kind == "ss":
            sim.add("solar system")
            sim.particles[1].x += 1e-6*variant
        elif kind == "3p":
            sim.add(m=1)
            for k in range(3):
                sim.add(m=3e-6, a=1.12**k, e=0.01, inc=0.01*k, l=1.3*k+0.1*variant)
            sim.move_to_com()
        elif kind == "4x2":
            for s in range(4):
                sim.add(m=1)
                for k in range(2):
                    sim.add(m=3e-6, a=1.12**k*(1+0.01*s), e=0.01, l=1.3*k+0.2*s+0.1*variant, primary=sim.particles[-1-k])
        sim.integrator = "whfast512"
        sim.integrator.kernel = kernel
        if kind == "4x2":
            sim.integrator.N_systems = 4
        sim.integrator.gr_potential = kw.get("gr", 1)
        sim.integrator.corrector = kw.get("corrector", 17)
        sim.dt = kw.get("dt", 6/365.25*2*math.pi)
        sim.exact_finish_time = 0
        if "exit_min_distance" in kw:
            sim.exit_min_distance = kw["exit_min_distance"]
        return sim

    def separate(self, sim, N):
        sim.integrator.concatenate_steps = N
        try:
            sim.steps(1)
        except parachron.Encounter:
            return "encounter"
        return "ok"

    def state(self, sim):
        return [sim.t] + [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles]

    def check(self, kind, kernel, N=3000, **kw):
        a, b = self.make(kind, kernel, 0, **kw), self.make(kind, kernel, 1, **kw)
        a2, b2 = self.make(kind, kernel, 0, **kw), self.make(kind, kernel, 1, **kw)
        try:
            a.steps_pair(b, N)
        except parachron.Encounter:
            pass
        self.separate(a2, N)
        self.separate(b2, N)
        self.assertEqual(self.state(a), self.state(a2))
        self.assertEqual(self.state(b), self.state(b2))
        self.assertEqual(a._status, a2._status)
        self.assertEqual(b._status, b2._status)

    def kernel_ok(self, kernel):
        try:
            s = self.make("ss", kernel)
            s.integrator.concatenate_steps = 1
            s.steps(1)
        except RuntimeError:
            return False
        return True

    def test_pair_identical_to_separate(self):
        for kernel in [2, 3, 4, 1, 0]:   # 1 (assembly) and 0 may not be pairable: same results via the fallback
            if not self.kernel_ok(kernel):
                continue
            for kind in ["ss", "3p", "4x2"]:
                N = 300 if kernel == 4 else 3000
                self.check(kind, kernel, N=N)
                self.check(kind, kernel, N=N, gr=0, corrector=0)

    def test_pair_different_dt(self):
        a, b = self.make("ss", 4, 0), self.make("ss", 4, 1, dt=3/365.25*2*math.pi)
        a2, b2 = self.make("ss", 4, 0), self.make("ss", 4, 1, dt=3/365.25*2*math.pi)
        a.steps_pair(b, 200)
        self.separate(a2, 200)
        self.separate(b2, 200)
        self.assertEqual(self.state(a), self.state(a2))
        self.assertEqual(self.state(b), self.state(b2))

    def test_pair_one_encounter(self):
        # Only the first simulation has a close encounter; the second one continues to N_steps.
        def mk(variant, kernel):
            sim = parachron.Simulation()
            sim.add(m=1)
            sim.add(m=1e-8, a=1)
            sim.add(m=1e-8, a=1.2 if variant == 0 else 2.5, f=math.pi)
            sim.integrator = "whfast512"
            sim.integrator.kernel = kernel
            sim.dt = 6/365.25*2*math.pi
            sim.exact_finish_time = 0
            sim.exit_min_distance = 0.5
            return sim
        for kernel in [3, 4]:
            if not self.kernel_ok(kernel):
                continue
            a, b, a2, b2 = mk(0, kernel), mk(1, kernel), mk(0, kernel), mk(1, kernel)
            with self.assertRaises(parachron.Encounter):
                a.steps_pair(b, 1000)
            self.assertAlmostEqual(a.t, 11.456, 2)
            self.assertEqual(b.t, 1000*b.dt)
            self.assertEqual(self.separate(a2, 1000), "encounter")
            self.separate(b2, 1000)
            self.assertEqual(self.state(a), self.state(a2))
            self.assertEqual(self.state(b), self.state(b2))

class TestIntegratorWHFast512KernelASM(TestIntegratorWHFast512):
    kernel = 1
class TestIntegratorWHFast512KernelCAVX512(TestIntegratorWHFast512):
    kernel = 2
class TestIntegratorWHFast512KernelCAVX2(TestIntegratorWHFast512):
    kernel = 3
class TestIntegratorWHFast512KernelCGeneric(TestIntegratorWHFast512):
    kernel = 4


if __name__ == "__main__":
    unittest.main()
