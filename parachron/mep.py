# Copyright (c) 2026 Michael Hippke. Part of Parachron (https://github.com/hippke/parachron), a modified version
# of REBOUND. Licensed under the GNU General Public License v3 or later; see LICENSE.
"""
Mean-element parareal (MEP): parallel-in-time integration of one WHFast512 simulation on many cores.

The time axis is cut into windows of `cores` slices. All slices of a window are integrated at the same time with
WHFast512 (one thread per slice), while a cheap coarse propagator in mean orbital elements links them (parareal).
After a few iterations every window has converged to the serial WHFast512 trajectory within the tolerance `tol`.
For the Solar System this is about 3x faster than serial WHFast512 on 16 cores at reference accuracy.

Requirements: Linux, a central object and exactly 8 planets ordered by distance (the four inner planets are treated
by perturbation theory in the coarse step, the four outer ones are integrated directly), WHFast512, G = 1 (units AU,
Msun, yr/2pi), N_systems = 1. The method was designed and validated for the Solar System.

Runs can be checkpointed to a small file and resumed later (also in another process); a resumed run is bitwise
identical to an uninterrupted one.
"""
from ctypes import Structure, CFUNCTYPE, c_int, c_int64, c_uint, c_double, c_char, c_char_p, c_void_p, POINTER, byref
from . import clibrebound
import sys

REB_MEP_KHIST = 12
REB_MEP_PATH = 1024
REB_MEP_MSG = 256


class MEPOptions(Structure):
    """Mirror of struct reb_mep_options (see rebound.h)."""
    _fields_ = [("cores", c_int), ("helpers", c_int), ("pin", c_int), ("k_min", c_int), ("k_max", c_int),
                ("tol", c_double), ("slice_steps", c_int64), ("checkpoint_every", c_int64), ("walltime", c_double),
                ("calibration", c_int), ("coarse_giants_steps", c_int), ("second_order", c_int),
                ("partial", c_int)]


class _MEPState(Structure):
    """Mirror of struct reb_mep (see rebound.h)."""


_HEARTBEAT = CFUNCTYPE(None, POINTER(_MEPState))
_MEPState._fields_ = [("opt", MEPOptions), ("t", c_double), ("t_start", c_double), ("dt", c_double),
                      ("windows", c_int64), ("khist", c_int64 * REB_MEP_KHIST), ("masses", c_double * 9),
                      ("E", (c_double * 6) * 8), ("gr_potential", c_uint), ("corrector", c_uint), ("kernel", c_uint),
                      ("walltime_last", c_double), ("checkpoint", c_char * REB_MEP_PATH), ("message", c_char * REB_MEP_MSG),
                      ("heartbeat", _HEARTBEAT), ("ref", c_void_p), ("boundaries", POINTER(c_double)),
                      ("calibration_residual", c_double), ("internal", c_void_p), ("boundaries_count", c_int)]


_create = clibrebound.reb_mep_create
_create.restype = POINTER(_MEPState)
_create_ckpt = clibrebound.reb_mep_create_from_checkpoint
_create_ckpt.restype = POINTER(_MEPState)
_create_ckpt.argtypes = [c_char_p, POINTER(MEPOptions)]
_integrate = clibrebound.reb_mep_integrate
_integrate.restype = c_int
_integrate.argtypes = [POINTER(_MEPState), c_double]
_update = clibrebound.reb_mep_update_simulation
_free = clibrebound.reb_mep_free
_free.argtypes = [POINTER(_MEPState)]


class MEP:
    """
    Mean-element parareal integration of a WHFast512 simulation on many cores.

    Examples
    --------
    >>> sim = ...                                   # Sun + 8 planets, sim.integrator = "whfast512", sim.dt = 6 days
    >>> mep = parachron.MEP(sim, checkpoint="run.mep")
    >>> mep.integrate(1e6 * 2 * np.pi)              # 1 Myr; sim is updated to the end state
    >>> mep = parachron.MEP.from_checkpoint("run.mep")  # later, possibly in another process
    >>> mep.integrate(2e6 * 2 * np.pi)              # continues bitwise identically
    """

    def __init__(self, sim, cores=None, tol=2e-6, k_min=3, k_max=6, slice_steps=60000, helpers=None, pin=True,
                 checkpoint=None, checkpoint_every=100, calibration=0, coarse_giants_steps=0, second_order=0, partial=0):
        """
        Arguments
        ---------
        sim : parachron.Simulation
            Start state. A central object and exactly 8 planets, integrator "whfast512" (its settings gr_potential,
            corrector and kernel are used), G = 1. `sim` is updated in place after every call of `integrate`.
        cores : int, optional
            Slices per window, one thread each. Default: the number of physical cores available to the process.
        tol : float
            Tolerance [AU] of the stopping rule: a window is accepted when no planet position changes by more than
            `tol` between two iterations. Default 2e-6 (reference accuracy: secular frequencies within a few 1e-8
            arcsec/yr of serial WHFast512 over 10 Myr); 3e-6 and 1e-5 are 4% and 10% faster and were validated by
            ensembles. The tolerance is the same for any number of cores, but the accuracy it gives is not: a window
            of more slices iterates more, and a run has fewer window ends where residual enters. At 64 cores the
            default is about three times more accurate than at 16 (10 Myr secular frequencies: rms 6e-9 instead of
            1.9e-8 arcsec/yr), and tol = 5e-6 still matches the 16-core default (1.4e-8) with about 5% fewer
            iterations (~4% faster; tol = 1e-5: ~8% faster, rms 3.4e-8). On many cores a looser tolerance therefore
            buys speed at the accuracy of the 16-core default.
        k_min, k_max : int
            Minimum and maximum number of parareal iterations per window (defaults 3 and 6). A window that reaches
            k_max is accepted even if it has not met `tol`; on 32 or more cores use k_max = 8 (at 64 cores about 3% of
            the windows need 7 iterations at the default tolerance).
        slice_steps : int
            Time steps per slice (even). Default 60000, i.e. 986 yr for dt = 6 days.
        helpers : int, optional
            Helper threads of the coarse step. Default min(8, cores // 2).
        pin : bool
            Pin the threads to physical cores (default True).
        checkpoint : str, optional
            File to which the run state is written (atomically) every `checkpoint_every` windows and at the end of
            every `integrate` call. Resume with `MEP.from_checkpoint`.
        checkpoint_every : int
            Windows between checkpoints (0: only at the end of each `integrate` call). Default 100.
        calibration : int
            Secular calibration of the coarse step: number of previous windows whose coarse defects are fitted by a
            linear map of the eccentricity and inclination vectors (at most 64; 0, default: off). Removes the
            secular-frequency bias of windows that stop after two iterations; use with k_min = 2 (e.g. 40).
        coarse_giants_steps : int
            Fine time steps per step of the outer planets' coarse integration (0, default: 360 days). Use 30 (180
            days) with k_min = 2.
        second_order : int
            1: add the second-order Earth-Mars-Jupiter near-resonance term to the coarse step (fewer iterations,
            mostly on many cores). 0 (default): off.
        partial : int
            Partial windows: a window whose first slices have converged (every boundary update <= `tol`) while at
            most `partial` boundaries at its end have not is accepted up to its last converged boundary; the next
            window starts there. Saves the iterations that only the last slices would need. 0 (default): off.

        Accuracy classes for the Solar System (16 cores, 7950X): tol 2e-6 (default, 3.2x); tol 1e-5 (3.5x);
        k_min = 2, tol 3e-4, calibration = 40, coarse_giants_steps = 30 (4.0x, secular frequencies within a few
        1e-6 arcsec/yr of serial WHFast512, unbiased).
        """
        if not sys.platform.startswith("linux"):
            raise RuntimeError("MEP (mean-element parareal) is only available on Linux.")
        o = MEPOptions()
        clibrebound.reb_mep_options_default(byref(o))
        o.cores = 0 if cores is None else int(cores)
        o.helpers = -1 if helpers is None else int(helpers)
        o.pin = 1 if pin else 0
        o.k_min, o.k_max, o.tol, o.slice_steps, o.checkpoint_every = int(k_min), int(k_max), float(tol), int(slice_steps), int(checkpoint_every)
        o.calibration, o.coarse_giants_steps, o.second_order = int(calibration), int(coarse_giants_steps), int(second_order)
        o.partial = int(partial)
        ptr = _create(byref(sim), byref(o), c_char_p(checkpoint.encode("utf-8")) if checkpoint else None)
        sim.process_messages()
        if not ptr:
            raise RuntimeError("Could not set up the MEP integration.")
        self._p = ptr
        self.sim = sim

    @classmethod
    def from_checkpoint(cls, filename, helpers=None, pin=True, checkpoint_every=None):
        """
        Resume a run from a checkpoint file. The options that change the result (cores, tol, k_min, k_max,
        slice_steps, calibration and its fit state, coarse_giants_steps, second_order, partial, time step, masses, WHFast512 settings) are
        taken from the file; `helpers`, `pin` and
        `checkpoint_every` may be changed. Returns a new MEP object; its `sim` holds the checkpointed state.
        """
        if not sys.platform.startswith("linux"):
            raise RuntimeError("MEP (mean-element parareal) is only available on Linux.")
        o = MEPOptions()
        clibrebound.reb_mep_options_default(byref(o))
        o.helpers = -1 if helpers is None else int(helpers)
        o.pin = 1 if pin else 0
        if checkpoint_every is not None:
            o.checkpoint_every = int(checkpoint_every)
        ptr = _create_ckpt(filename.encode("utf-8"), byref(o))
        if not ptr:
            raise RuntimeError(f"Could not read the MEP checkpoint file '{filename}'.")
        from .simulation import Simulation
        self = cls.__new__(cls)
        self._p = ptr
        self.sim = Simulation()
        _update(self._p, byref(self.sim))
        return self

    def integrate(self, tmax, walltime=None):
        """
        Integrate whole windows up to the last window boundary <= tmax (see `window_length`), then update `sim`.
        With partial windows, windows are started as long as a whole window ends <= tmax, so the run ends at a
        slice boundary less than one window length before tmax (`t` gives the time reached).

        walltime : float, optional
            Stop after the first window that ends after this many seconds of wall clock (the checkpoint is
            written). Use this to split a long run into pieces.

        Returns True if tmax (rounded down to a window boundary) was reached, False if the run stopped because of
        `walltime`. Raises KeyboardInterrupt if interrupted (Ctrl-C); the unfinished window is discarded, and `sim`
        and the checkpoint hold the last completed window.
        """
        self._p.contents.opt.walltime = 0. if walltime is None else float(walltime)
        status = _integrate(self._p, c_double(tmax))
        _update(self._p, byref(self.sim))
        if status == 1:
            raise RuntimeError(self._p.contents.message.decode("utf-8") or "MEP integration failed.")
        if status == 6:
            raise KeyboardInterrupt
        return status == 0

    @property
    def heartbeat(self):
        """
        Function called as heartbeat(mep) after every completed window, or None. `mep.t`, `mep.windows` and
        `mep.state` are current inside it (`mep.sim` is updated only at the end of `integrate`). It runs on a worker
        thread inside the parallel schedule and holds up the next window while it runs, so keep it short (e.g. store
        `mep.state`). Exceptions inside it are printed and ignored.
        """
        return getattr(self, "_heartbeat_py", None)

    @heartbeat.setter
    def heartbeat(self, func):
        if func is None:
            self._heartbeat_c = None
            self._p.contents.heartbeat = _HEARTBEAT()
        else:
            self._heartbeat_c = _HEARTBEAT(lambda _p: func(self))     # keep a reference: C holds only the pointer
            self._p.contents.heartbeat = self._heartbeat_c
        self._heartbeat_py = func

    @property
    def boundaries(self):
        """
        Inside the heartbeat: the end states of the slices of the completed window (`cores` of them, fewer for a
        partial window), as a list of 8 x 6 lists of Jacobi elements (see `state`); None outside the heartbeat.
        """
        b = self._p.contents.boundaries
        if not b:
            return None
        n = self._p.contents.boundaries_count
        return [[[b[(s * 8 + i) * 6 + e] for e in range(6)] for i in range(8)] for s in range(n)]

    @property
    def calibration_residual(self):
        """Fraction of the eccentricity-defect variance left unexplained by the last calibration fit (0: off)."""
        return self._p.contents.calibration_residual

    @property
    def t(self):
        """Current time (end of the last completed window)."""
        return self._p.contents.t

    @property
    def windows(self):
        """Number of completed windows since the start of the run."""
        return self._p.contents.windows

    @property
    def cores(self):
        """Slices per window (threads)."""
        return self._p.contents.opt.cores

    @property
    def window_length(self):
        """Simulated time per window (cores x slice_steps x dt)."""
        c = self._p.contents
        return c.opt.cores * c.opt.slice_steps * c.dt

    @property
    def iterations(self):
        """Histogram {k: number of windows that needed k parareal iterations}."""
        return {k: int(n) for k, n in enumerate(self._p.contents.khist) if n}

    @property
    def mean_iterations(self):
        h = self.iterations
        return sum(k * n for k, n in h.items()) / max(1, sum(h.values()))

    @property
    def walltime_last(self):
        """Wall clock [s] of the last `integrate` call."""
        return self._p.contents.walltime_last

    @property
    def state(self):
        """Exact run state: Jacobi elements (ln a, lambda, k, h, q, p) of the 8 planets."""
        return [list(row) for row in self._p.contents.E]

    def __del__(self):
        if getattr(self, "_p", None):
            _free(self._p)
            self._p = None

    def __repr__(self):
        return f"<parachron.MEP cores={self.cores} t={self.t:.6g} windows={self.windows} mean_iterations={self.mean_iterations:.2f}>"
