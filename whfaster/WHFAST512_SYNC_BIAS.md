# WHFast512: systematic energy change at every synchronization (draft bug report)

*Found Oct 22, 2026 while validating long parallel-in-time runs (`RESEARCH_LOG.md` iteration 36). Draft text for an
upstream issue; nothing has been submitted.*

## Summary
Every WHFast512 integration call (i.e. every synchronization, e.g. every output) changes the total energy of the
Solar System by a fixed **dE/E ≈ +7e-17**, independent of the state. Runs that synchronize often therefore drift
linearly in energy. Example: a 5 Gyr integration with an output every 1000 yr makes 5e6 calls → ≈ +3.5e-10 of drift on
top of the round-off random walk. The cause is the pair of precomputed 8×8 inertial ↔ Jacobi matrices, which are not
exact inverses in floating point.

## Measurement
Same number of steps done in C calls instead of one call, averaged over many start states along a trajectory
(`whfaster/pint/call_bias.c`; Solar System, 8 planets, dt = 6 d):

| kernel | corrector | GR | extra dE/E per call |
|---|---|---|---|
| asm | 17 | on | +7.2e-17 ± 1.5e-17 |
| C AVX-512 | 17 | on | +1.1e-16 ± 1.5e-17 |
| C AVX-512 | 0 | on | +7.3e-17 ± 0.6e-17 |
| C AVX-512 | 17 | off | +7.4e-17 ± 1.5e-17 |
| C AVX-512 | 0 | off | +6.5e-17 ± 0.6e-17 |
| C AVX2 | 17 | on | +7.4e-17 ± 0.5e-17 |

Serial runs with 100 calls per 986 yr drift linearly by +7.1e-17 (asm) / +7.7e-17 (C) per call. A 1 Gyr serial run
with one call per 986 yr (1e6 calls) ends at +6.8e-11, nearly linearly.

## Cause
`inertial_to_jacobi_posvel()` and `jacobi_to_inertial_posvel_and_com()` in `src/integrator_whfast512.c` use
`mat8_inertial_to_jacobi` and `mat8_jacobi_to_inertial` (built in `recalculate_constants_compute()`). Their product,
evaluated in long double for the Solar System masses, is diagonal to 1e-19 with diag − 1 =
[−1.37, −0.22, +0.07, −0.05, +0.45, −1.37, −0.17, −0.77] × 1e-16 (Mercury … Neptune). A round trip therefore scales
each planet's Jacobi position and velocity by (1 + e_i), the same for every state, changing its energy by ≈ 2 e_i |U_i|.
Weighted with the planets' potential energies (Jupiter +0.45, Saturn −1.37 dominate) this predicts +7.6e-17 per call,
as measured.

## Fix (in this fork, opt-in)
New `struct reb_integrator_whfast512_state` field `jacobi_recursive` (default 0 = unchanged, bitwise identical to the
current code): 1 uses REBOUND's recursive Jacobi transformations (`reb_transformations_inertial_to_jacobi_posvel`,
`reb_transformations_jacobi_to_inertial_posvel`, per system in the COM frame) at the start and end of each call. Their
round-off depends on the state and does not accumulate systematically:

| kernel | corrector / GR | extra dE/E per call with `jacobi_recursive = 1` |
|---|---|---|
| asm | 17 / on | −0.3e-17 ± 1.1e-17 |
| C AVX-512 | 17 / on | −1.5e-17 ± 1.1e-17 |
| C AVX2 (30000 states) | 17 / on | −0.9e-17 ± 0.5e-17 |
| C AVX2 (30000 states) | 17 / off | +0.5e-17 ± 0.5e-17 |
| C AVX2 (30000 states) | 0 / on | +0.1e-17 ± 0.2e-17 |

Results agree with the matrix version to ≤ 8e-16 after one call for N_systems = 1, 2, 4 (`whfaster/pint/test_jr.c`);
speed is unchanged (the transformation runs once per call). In 1 Gyr parallel-in-time runs (≈ 2e6 calls) the energy
drift went from +7.5e-10 per Gyr to −6.7e-11 per Gyr, the level of serial runs.

Possible upstream alternatives: make the recursive transformation the default (changes results at round-off level),
or construct `mat8_jacobi_to_inertial` so that the diagonal of the product is exactly 1 for the given masses.
