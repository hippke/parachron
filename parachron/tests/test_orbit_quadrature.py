import parachron
import math
import unittest

# Cartesian states (G = 1, primary of mass 1 at rest at the origin, planet mass 1e-7) whose eccentric anomaly is
# within ~1e-7 rad of pi/2. Before the fix, sin(E) computed from them exceeded 1 by round-off and
# reb_orbit_from_particle returned NaN for E, M and the mean longitude l.
STATES = [
    ("-0x1.26c3da9fa0f1p-1", "0x1.5a97347a89a57p-1", "0x1.3fa90f74e89a5p-4", "-0x1.dbe379427ec8ep-1", "-0x1.031b9dd5d43bp-1", "-0x1.ffc844e7c4589p-6"),
    ("-0x1.423f1a63d1b7fp-2", "0x1.37fbb62848691p-2", "0x1.28d16417a6b1bp-5", "-0x1.52fdc07e1c8d5p+0", "-0x1.71246bb4a7a5p-1", "-0x1.6c8f3b80012f5p-5"),
    ("-0x1.552148ae00626p+0", "0x1.903a57da03819p+0", "0x1.713fd24f67451p-3", "-0x1.39003f310fc31p-1", "-0x1.54d720d0e3209p-2", "-0x1.509be2f432f6cp-6"),
]


class TestOrbitQuadrature(unittest.TestCase):
    def test_known_states(self):
        for st in STATES:
            x, y, z, vx, vy, vz = (float.fromhex(v) for v in st)
            sim = parachron.Simulation()
            sim.add(m=1.)
            sim.add(m=1e-7, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
            o = sim.particles[1].orbit(primary=sim.particles[0])
            for q in (o.E, o.M, o.l, o.f):
                self.assertTrue(math.isfinite(q))
            self.assertAlmostEqual(math.cos(o.E), 0., delta=1e-6)

    def test_scan_across_quadrature(self):
        sim = parachron.Simulation()
        sim.add(m=1.)
        for i in range(4000):
            e = 0.02 + 0.28 * (i % 97) / 97.
            E = math.pi / 2 + (i - 2000) * 5e-11 * (1 if i % 2 else -1)
            f = 2. * math.atan(math.sqrt((1. + e) / (1. - e)) * math.tan(E / 2.))
            p = parachron.Particle(simulation=sim, primary=sim.particles[0], m=1e-7, a=0.5 + (i % 13) / 7., e=e, inc=0.1, Omega=0.2, omega=0.3, f=f)
            o = p.orbit(primary=sim.particles[0], G=sim.G)
            self.assertTrue(math.isfinite(o.l), msg=f"l is NaN for i={i}")


if __name__ == "__main__":
    unittest.main()
