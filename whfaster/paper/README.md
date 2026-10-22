# Parallel-in-time integration of the Solar System: material for a paper

Self-contained write-up of the multicore work on one long WHFast512 Solar System trajectory (Sun + 8 planets,
dt = 6 d, GR, 17th-order corrector; the setup of arXiv 2609.25339), state of Oct 22, 2026.

| part | content | audience |
|---|---|---|
| [`01_physics.md`](01_physics.md) | why parareal fails for planetary orbits (Kepler shear × short-period terms) and how a mean-element coarse propagator makes it converge; what limits convergence; accuracy yardstick, ensembles, energy | physics (main part) |
| [`02_algorithm.md`](02_algorithm.md) | windowed, pipelined parareal; stopping rule; parallel and staged coarse step; time budget; speed-ups per accuracy class; rejected alternatives; cost analysis of a more accurate coarse step | numerical methods |
| [`03_implementation.md`](03_implementation.md) | code, reproducibility, concurrency details | brief |

**Headline results (16-core AMD Ryzen 9 7950X).**
* Default accuracy (secular frequencies within ≈ 1e-7 ″/yr of serial runs, i.e. far inside the serial integrator's
  own dt = 6 d error of up to 5e-6 ″/yr): **2.9–3.0× faster**; the 5 Gyr run in ≈ 3.2 h instead of 9.1 h.
* Looser, ensemble-validated tolerances: 3.0–3.6×; at the accuracy of a dt = 8–12 d serial run: 3.7–4.7×
  (2.4–2.9× faster than such a serial run).
* An ensemble of 16 × 200 Myr shows serial and parareal dynamics to be statistically indistinguishable.
* Found on the way: a systematic energy offset of +7.5e-17 per synchronisation in WHFast512 (upstream), fixed by an
  opt-in recursive Jacobi transformation.

**Provenance.** Every number is traceable to `../RESEARCH_LOG.md` (iterations cited as “it. N”; rounds 1–3 by the
previous sessions, iterations 28–48 on the 7950X) and `../pint_history.csv`. Earlier rounds measured speed on a
2-vCPU cloud host and gave modelled speed-ups; all speed-ups quoted here are measured on 16 cores unless marked as a
model. Open items: residual parareal energy drift (−1.25e-10 per Gyr), generalisation to other planetary systems,
thorough study of faster tolerances (`../RESEARCH_PLAN_R5.md`).
