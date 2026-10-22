# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
# The ensemble runner (rebound.steps_simulations / reb_simulations_steps) steps WHFast512
# simulations two at a time with the paired kernel when possible. Results must be bitwise
# identical to stepping each simulation on its own.
import unittest
import rebound

SS = [  # Sun + 8 planets (G = 1, AU, yr/2pi), as in benchmarks/perf
    (1.0, (-0.008816286905115728, -0.0010954624020472503, 0.0002143887550311963), (0.00030366461216799816, -0.0004844036474574087, -6.271127210658849e-06)),
    (1.6601208254808336e-07, (-0.4005015932430153, -0.1717825232600333, 0.02213496003871717), (0.42148131096005145, -1.410098437567698, -0.15349108993869006)),
    (2.447838287784771e-06, (0.6965750451393025, -0.2028719447054862, -0.04298548683058286), (0.32470124542779734, 1.1260969773655014, -0.0040195087604289335)),
    (3.0404326489511185e-06, (-0.9779521627476298, 0.2005298813779153, 0.0001536219574093047), (-0.20841707233386245, -0.9913854283436436, 4.2810116116963296e-05)),
    (3.2271560828978514e-07, (0.6113037476689052, 1.3829508325618787, 0.01377796023041797), (-0.6908541024516895, 0.43030638318698925, 0.018813543022624015)),
    (0.0009547919099366768, (3.590557108006218, -3.485622651312009, -0.06606227034262536), (0.29883609765296565, 0.3335574047451432, -0.008080213946271302)),
    (0.0002858856700231729, (5.836138545391186, -8.10708045469558, -0.09155318103668262), (0.2433082689003069, 0.19125094883087373, -0.012857622616101758)),
    (4.366249613200406e-05, (14.70069810024564, 13.04566034049498, -0.14245405553548557), (-0.15287862212946733, 0.1601553779126656, 0.0026049064826932225)),
    (5.151383772628957e-05, (29.42823543749598, -5.384055011101007, -0.5675216046021788), (0.03184432633916516, 0.18139574116224567, -0.0044542009119127255)),
]


def make(k, kernel, nplanets=8, gr=1, exit_min=0.):
    sim = rebound.Simulation()
    sim.integrator = "whfast512"
    sim.integrator.gr_potential = gr
    sim.integrator.concatenate_steps = 100
    sim.integrator.kernel = kernel
    sim.dt = 5./365.25*6.283185307179586
    sim.exact_finish_time = 0
    for m, x, v in SS[:nplanets+1]:
        sim.add(m=m, x=x[0], y=x[1], z=x[2], vx=v[0], vy=v[1], vz=v[2])
    sim.particles[1].x += 1e-12*k
    if exit_min:
        sim.exit_min_distance = exit_min
    return sim


def state(sim):
    return [(p.x, p.y, p.z, p.vx, p.vy, p.vz) for p in sim.particles] + [sim.t, sim.steps_done, sim._status]


def kernel_ok(kernel):
    sim = make(0, kernel)
    try:
        sim.steps(1)
    except Exception:
        return False
    return sim.integrator.kernel_used == kernel


class TestEnsemblePairing(unittest.TestCase):
    def check(self, kernel, nplanets, n_sims=5, steps=30, exit_min=0.):
        if not kernel_ok(kernel):
            self.skipTest("kernel %d not available" % kernel)
        ref = [make(k, kernel, nplanets, exit_min=exit_min if k == 1 else 0.) for k in range(n_sims)]
        for s in ref:
            try:
                s.steps(steps)
            except rebound.Encounter:
                pass
        ens = [make(k, kernel, nplanets, exit_min=exit_min if k == 1 else 0.) for k in range(n_sims)]
        rebound.steps_simulations(ens, steps, n_threads=2)
        for a, b in zip(ref, ens):
            self.assertEqual(state(a), state(b))

    def test_c_avx512_solar_system(self):
        self.check(2, 8)

    def test_c_avx2_four_lanes(self):
        self.check(3, 3)

    def test_c_generic(self):
        self.check(4, 8, n_sims=3, steps=10)

    def test_assembly_not_paired(self):
        self.check(1, 8)

    def test_encounter_in_one_simulation(self):
        self.check(3, 3, exit_min=0.5)


if __name__ == "__main__":
    unittest.main()
