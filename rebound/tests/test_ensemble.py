# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
import rebound
import unittest
import struct
import hashlib
import threading


def state_hash(sims):
    m = hashlib.md5()
    for s in sims:
        m.update(struct.pack("d", s.t))
        for p in s.particles:
            m.update(struct.pack("7d", p.x, p.y, p.z, p.vx, p.vy, p.vz, p.m))
    return m.hexdigest()


def make(kind, i):
    s = rebound.Simulation()
    s.rand_seed = 1000 + i  # deterministic random numbers (MEGNO initial variations)
    s.add(m=1.)
    if kind == "megno":
        s.add(m=1e-3, a=1.)
        s.add(m=1e-3, a=1.5 + 0.01*i, e=0.05)
        s.integrator = "whfast"
        s.dt = 0.02
        s.init_megno()
    elif kind == "ias15":
        for k in range(4):
            s.add(m=1e-4, a=1. + 0.3*k + 1e-3*i, e=0.1, inc=0.05*k)
        s.integrator = "ias15"
    elif kind == "mercurius":
        s.add(m=1e-3, a=1.)
        for k in range(30):
            s.add(m=1e-7, r=1e-4, a=1.2 + 0.02*k + 1e-4*i, e=0.05, f=0.7*k)
        s.integrator = "mercurius"
        s.dt = 0.02
        s.collision = "direct"
        s.collision_resolve = "merge"
    s.move_to_com()
    return s


class TestEnsemble(unittest.TestCase):

    def check_kind(self, kind, tmax):
        a = [make(kind, i) for i in range(6)]
        for s in a:
            s.integrate(tmax)
        b = [make(kind, i) for i in range(6)]
        st = rebound.integrate_simulations(b, tmax, n_threads=3)
        self.assertEqual(st, [s._status for s in a])
        self.assertEqual(state_hash(a), state_hash(b))
        if kind == "megno":
            self.assertEqual([s.megno() for s in a], [s.megno() for s in b])

    def test_whfast_megno(self):
        # Variational equations in threads (cf. upstream issue #476)
        self.check_kind("megno", 50.)

    def test_ias15(self):
        self.check_kind("ias15", 20.)

    def test_mercurius_collisions(self):
        self.check_kind("mercurius", 20.)

    def test_python_threads(self):
        # sim.integrate releases the interpreter lock; plain Python threads give the same results.
        a = [make("ias15", i) for i in range(4)]
        for s in a:
            s.integrate(10.)
        b = [make("ias15", i) for i in range(4)]
        th = [threading.Thread(target=s.integrate, args=(10.,)) for s in b]
        for t in th:
            t.start()
        for t in th:
            t.join()
        self.assertEqual(state_hash(a), state_hash(b))

    def test_steps(self):
        a = [make("megno", i) for i in range(5)]
        for s in a:
            s.steps(500)
        b = [make("megno", i) for i in range(5)]
        rebound.steps_simulations(b, 500, n_threads=2)
        self.assertEqual(state_hash(a), state_hash(b))

    def test_tmax_list_and_status(self):
        b = [make("ias15", i) for i in range(3)]
        b[1].exit_max_distance = 1.2
        st = rebound.integrate_simulations(b, [1., 50., 3.])
        self.assertEqual(st[0], 0)
        self.assertEqual(st[1], 4)  # escape
        self.assertEqual(b[0].t, 1.)
        self.assertEqual(b[2].t, 3.)
        with self.assertRaises(ValueError):
            rebound.integrate_simulations(b, [1., 2.])
        self.assertEqual(rebound.integrate_simulations([], 1.), [])


if __name__ == "__main__":
    unittest.main()
