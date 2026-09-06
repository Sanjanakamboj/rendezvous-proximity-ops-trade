"""Milestone 5: independent nonlinear (inertial, Cartesian two-body) validation
model, and LVLH <-> inertial frame kinematics.

This module is deliberately independent of the CW closed-form STM (cw.py):
it integrates the RAW two-body equations of motion

    rddot = -mu * r / ||r||^3

in an inertial (ECI-like) frame for the chief and deputy SEPARATELY with
`scipy.integrate.solve_ivp` (DOP853), and reconstructs the LVLH relative
state at each requested time from first principles. No CW equation, STM
entry, or linearization is used anywhere in this module -- it exists to
check the CW model (cw.py, rendezvous.py), not to reuse it.

Frame convention (matches DESIGN.md section 3, restated in inertial terms):
  x_hat: radial, outward from Earth's center through the chief
  y_hat: along-track, direction of chief's orbital velocity
  z_hat: cross-track, along the chief's orbital angular momentum h = r x v

Rotating-frame kinematics: the LVLH frame's inertial angular velocity is
    omega_inertial(t) = h_c(t) / ||r_c(t)||^2
(exact for any orbit, not only circular -- h_c/r_c^2 is the general
osculating orbit-plane rotation rate). For the deputy inertial velocity
given LVLH relative position/velocity (r_rel, v_rel):

    v_dep = v_chief + C @ v_rel + omega_inertial x (C @ r_rel)

where C = [x_hat | y_hat | z_hat] is the LVLH-to-inertial rotation matrix
(columns are the LVLH basis vectors expressed in inertial coordinates).
The inverse transform solves this same equation for v_rel.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .constraints import (
    DEFAULT_SAMPLE_DT_S,
    R_KOZ,
    CHIEF_CROSSING_TOL_M,
    CHIEF_CROSSING_Y_LIMIT,
    TrajectoryConstraintResult,
    classify_point,
    inside_approach_corridor,
)
from .orbit import MU_EARTH

__all__ = [
    "two_body_ode",
    "propagate_two_body",
    "chief_circular_initial_state",
    "lvlh_basis",
    "lvlh_to_inertial",
    "inertial_to_lvlh",
    "NonlinearTrajectory",
    "propagate_relative_nonlinear",
    "propagate_relative_nonlinear_dense",
    "classify_nonlinear_trajectory",
]


def two_body_ode(t: float, state: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Raw inertial two-body Cartesian ODE: d/dt[r,v] = [v, -mu r/|r|^3].

    Independent of, and never calls, any CW equation.
    """
    r = state[:3]
    v = state[3:]
    r_norm = np.linalg.norm(r)
    a = -mu * r / r_norm**3
    return np.concatenate([v, a])


def propagate_two_body(
    state0: np.ndarray, t_span: tuple[float, float], t_eval: np.ndarray | None = None,
    mu: float = MU_EARTH, rtol: float = 1e-12, atol: float = 1e-9, method: str = "DOP853",
):
    """Integrate the raw two-body ODE with solve_ivp; returns the OdeResult.

    Default tolerances (rtol=1e-12, atol=1e-9 -- atol in meters/[m/s], loose
    relative to the km-scale state but tight relative to the meter-scale
    relative motion being resolved) are documented and used consistently;
    an alternate tolerance is exercised in the integrator-convergence test.
    """
    return solve_ivp(
        two_body_ode, t_span, state0, args=(mu,),
        method=method, rtol=rtol, atol=atol,
        t_eval=t_eval, dense_output=(t_eval is None),
    )


def chief_circular_initial_state(a: float, mu: float = MU_EARTH) -> np.ndarray:
    """Inertial state [r,v] (6,) for a circular chief orbit of radius `a`.

    Placed at inertial position (a,0,0) with velocity (0, sqrt(mu/a), 0) --
    i.e. orbiting in the inertial x-y plane, angular momentum along +z.
    This is an analytic circular-orbit initial condition; it is then
    verified (not assumed) to stay circular under numerical two-body
    propagation in the M5 test suite (conservation checks).
    """
    v_circ = np.sqrt(mu / a)
    return np.array([a, 0.0, 0.0, 0.0, v_circ, 0.0])


def lvlh_basis(r_c: np.ndarray, v_c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """LVLH-to-inertial rotation matrix C and inertial angular velocity omega.

    Returns (C, omega_inertial): C is 3x3 with columns [x_hat, y_hat, z_hat]
    (LVLH basis vectors in inertial coordinates); omega_inertial = h/||r||^2
    is the LVLH frame's inertial angular velocity vector (3,).
    """
    r_norm = np.linalg.norm(r_c)
    x_hat = r_c / r_norm
    h = np.cross(r_c, v_c)
    h_norm = np.linalg.norm(h)
    z_hat = h / h_norm
    y_hat = np.cross(z_hat, x_hat)
    C = np.column_stack([x_hat, y_hat, z_hat])
    omega_inertial = h / r_norm**2
    return C, omega_inertial


def lvlh_to_inertial(
    r_c: np.ndarray, v_c: np.ndarray, r_rel: np.ndarray, v_rel: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a deputy LVLH relative state (r_rel, v_rel) to inertial (r_dep, v_dep).

    v_dep = v_chief + C @ v_rel + omega x (C @ r_rel)  -- see module docstring.
    """
    C, omega = lvlh_basis(r_c, v_c)
    r_rel_inertial = C @ r_rel
    r_dep = r_c + r_rel_inertial
    v_dep = v_c + C @ v_rel + np.cross(omega, r_rel_inertial)
    return r_dep, v_dep


def inertial_to_lvlh(
    r_c: np.ndarray, v_c: np.ndarray, r_dep: np.ndarray, v_dep: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert inertial chief/deputy states to a deputy LVLH relative state.

    Inverse of lvlh_to_inertial; solves the same rotating-frame velocity
    relation for v_rel using C^T (C is orthonormal, so C^-1 = C^T).
    """
    C, omega = lvlh_basis(r_c, v_c)
    r_rel_inertial = r_dep - r_c
    r_rel = C.T @ r_rel_inertial
    v_rel_inertial = v_dep - v_c - np.cross(omega, r_rel_inertial)
    v_rel = C.T @ v_rel_inertial
    return r_rel, v_rel


@dataclass(frozen=True)
class NonlinearTrajectory:
    """Result of a nonlinear two-body relative-motion propagation."""

    times: np.ndarray  # (N,)
    lvlh_states: np.ndarray  # (N, 6): [x,y,z,vx,vy,vz] per row
    chief_states: np.ndarray  # (N, 6) inertial
    deputy_states: np.ndarray  # (N, 6) inertial


def propagate_relative_nonlinear(
    r0_rel: np.ndarray, v0_rel_plus: np.ndarray, times: np.ndarray,
    a: float, mu: float = MU_EARTH, rtol: float = 1e-12, atol: float = 1e-9,
) -> NonlinearTrajectory:
    """Propagate a deputy (given its LVLH state right after burn 1) nonlinearly.

    Chief starts on the analytic circular initial condition
    (chief_circular_initial_state); the deputy's inertial initial state is
    derived from (r0_rel, v0_rel_plus) via lvlh_to_inertial at t=0. Both
    are then propagated INDEPENDENTLY with the raw two-body ODE (no CW
    equations), and the deputy-relative-to-chief LVLH state is
    reconstructed at each time in `times` via inertial_to_lvlh using each
    body's own propagated (osculating) state -- not the initial frame.
    """
    times = np.asarray(times, dtype=float)
    chief0 = chief_circular_initial_state(a, mu)
    r_c0, v_c0 = chief0[:3], chief0[3:]
    r_dep0, v_dep0 = lvlh_to_inertial(r_c0, v_c0, r0_rel, v0_rel_plus)
    deputy0 = np.concatenate([r_dep0, v_dep0])

    t_span = (float(times.min()), float(times.max()))
    # solve_ivp requires t_eval sorted within t_span and starting logic
    # handles t0 possibly not being 0; here times always starts at 0.
    sol_chief = propagate_two_body(chief0, (0.0, t_span[1]), t_eval=times, rtol=rtol, atol=atol)
    sol_deputy = propagate_two_body(deputy0, (0.0, t_span[1]), t_eval=times, rtol=rtol, atol=atol)
    if not (sol_chief.success and sol_deputy.success):
        raise RuntimeError(f"nonlinear propagation failed: chief={sol_chief.message}, deputy={sol_deputy.message}")

    chief_states = sol_chief.y.T  # (N,6)
    deputy_states = sol_deputy.y.T

    lvlh_states = np.zeros((len(times), 6))
    for i in range(len(times)):
        r_c, v_c = chief_states[i, :3], chief_states[i, 3:]
        r_d, v_d = deputy_states[i, :3], deputy_states[i, 3:]
        r_rel, v_rel = inertial_to_lvlh(r_c, v_c, r_d, v_d)
        lvlh_states[i, :3] = r_rel
        lvlh_states[i, 3:] = v_rel

    return NonlinearTrajectory(
        times=times, lvlh_states=lvlh_states,
        chief_states=chief_states, deputy_states=deputy_states,
    )


def propagate_relative_nonlinear_dense(
    r0_rel: np.ndarray, v0_rel_plus: np.ndarray, t: float,
    a: float, mu: float = MU_EARTH, rtol: float = 1e-12, atol: float = 1e-9,
):
    """Like propagate_relative_nonlinear, but returns dense-output interpolants.

    Returns (sol_chief, sol_deputy), each a scipy OdeResult with a callable
    `.sol(t)` continuous interpolant, so the LVLH relative state can be
    queried cheaply at arbitrary times in [0, t] without re-integrating --
    used for continuous-time M4-style constraint checking on the nonlinear
    trajectory (classify_nonlinear_trajectory, below).
    """
    chief0 = chief_circular_initial_state(a, mu)
    r_c0, v_c0 = chief0[:3], chief0[3:]
    r_dep0, v_dep0 = lvlh_to_inertial(r_c0, v_c0, r0_rel, v0_rel_plus)
    deputy0 = np.concatenate([r_dep0, v_dep0])

    sol_chief = propagate_two_body(chief0, (0.0, t), t_eval=None, rtol=rtol, atol=atol)
    sol_deputy = propagate_two_body(deputy0, (0.0, t), t_eval=None, rtol=rtol, atol=atol)
    if not (sol_chief.success and sol_deputy.success):
        raise RuntimeError(f"nonlinear propagation failed: chief={sol_chief.message}, deputy={sol_deputy.message}")
    return sol_chief, sol_deputy


def _lvlh_state_at(sol_chief, sol_deputy, tt: float) -> np.ndarray:
    chief_state = sol_chief.sol(tt)
    deputy_state = sol_deputy.sol(tt)
    r_rel, v_rel = inertial_to_lvlh(chief_state[:3], chief_state[3:], deputy_state[:3], deputy_state[3:])
    return np.concatenate([r_rel, v_rel])


def classify_nonlinear_trajectory(
    r0_rel: np.ndarray, v0_rel_plus: np.ndarray, t: float,
    a: float, mu: float = MU_EARTH,
    sample_dt: float = DEFAULT_SAMPLE_DT_S,
    refine: bool = True,
    rtol: float = 1e-12, atol: float = 1e-9,
) -> TrajectoryConstraintResult:
    """M4-identical continuous-time constraint screening, applied to the
    NONLINEAR two-body trajectory instead of the CW closed-form STM.

    Mirrors constraints.classify_trajectory_constraints exactly in method
    (dense sample -> classify each sample -> bounded-scalar refine the
    closest-approach and worst-outside-corridor-clearance candidates) and
    reuses its exact per-point classification rule (classify_point) and
    constants (R_KOZ, corridor bounds, chief-crossing limit/tolerance) --
    it does not redefine or duplicate the constraint definitions, only the
    trajectory SOURCE differs (nonlinear two-body vs. CW STM).
    """
    from scipy.optimize import minimize_scalar

    sol_chief, sol_deputy = propagate_relative_nonlinear_dense(r0_rel, v0_rel_plus, t, a, mu, rtol, atol)

    n_samples = max(2, int(np.ceil(t / sample_dt)) + 1)
    times = np.linspace(0.0, t, n_samples)
    positions = np.array([_lvlh_state_at(sol_chief, sol_deputy, tt)[:3] for tt in times])

    classifications = [classify_point(p) for p in positions]
    distances = np.array([c.distance_m for c in classifications])
    ys = positions[:, 1]

    i_min_dist = int(np.argmin(distances))
    min_distance_m = float(distances[i_min_dist])
    time_of_min_distance_s = float(times[i_min_dist])

    outside_mask = np.array([not c.inside_corridor for c in classifications])
    if outside_mask.any():
        clearance_outside = distances[outside_mask] - R_KOZ
        i_outside = np.where(outside_mask)[0][np.argmin(clearance_outside)]
        min_clearance_outside_corridor_m = float(distances[i_outside] - R_KOZ)
    else:
        min_clearance_outside_corridor_m = float("nan")

    def distance_at(tt: float) -> float:
        return float(np.linalg.norm(_lvlh_state_at(sol_chief, sol_deputy, tt)[:3]))

    _NOT_A_CANDIDATE = 1.0e9

    def clearance_outside_at(tt: float) -> float:
        r = _lvlh_state_at(sol_chief, sol_deputy, tt)[:3]
        if inside_approach_corridor(r):
            return _NOT_A_CANDIDATE
        return float(np.linalg.norm(r)) - R_KOZ

    refined_clearance_violation_time_s = None
    if refine:
        lo = float(times[max(0, i_min_dist - 1)])
        hi = float(times[min(len(times) - 1, i_min_dist + 1)])
        if hi > lo:
            res = minimize_scalar(distance_at, bounds=(lo, hi), method="bounded")
            if res.fun < min_distance_m:
                min_distance_m = float(res.fun)
                time_of_min_distance_s = float(res.x)

        if outside_mask.any() and np.isfinite(min_clearance_outside_corridor_m):
            lo_c = float(times[max(0, i_outside - 1)])
            hi_c = float(times[min(len(times) - 1, i_outside + 1)])
            if hi_c > lo_c:
                res_c = minimize_scalar(clearance_outside_at, bounds=(lo_c, hi_c), method="bounded")
                if res_c.fun < _NOT_A_CANDIDATE and res_c.fun < min_clearance_outside_corridor_m:
                    min_clearance_outside_corridor_m = float(res_c.fun)
                    refined_clearance_violation_time_s = float(res_c.x)

    min_y_m = float(np.min(ys))
    max_y_m = float(np.max(ys))

    koz_violation = any(c.koz_violation for c in classifications)
    corridor_violation = any(c.corridor_violation for c in classifications)
    chief_crossing_violation = bool(max_y_m > (CHIEF_CROSSING_Y_LIMIT + CHIEF_CROSSING_TOL_M))

    refinement_found_koz_violation = (
        np.isfinite(min_clearance_outside_corridor_m) and min_clearance_outside_corridor_m < 0.0
    )
    if refinement_found_koz_violation:
        koz_violation = True

    first_violation_time_s = None
    reasons: list[str] = []

    def _note(reason: str, tt: float) -> None:
        nonlocal first_violation_time_s
        if reason not in reasons:
            reasons.append(reason)
        if first_violation_time_s is None or tt < first_violation_time_s:
            first_violation_time_s = tt

    for c, tt in zip(classifications, times):
        if c.koz_violation:
            _note("koz_violation", float(tt))
        if c.corridor_violation:
            _note("corridor_violation", float(tt))
        if c.chief_crossing_violation:
            _note("chief_crossing_violation", float(tt))
    if refinement_found_koz_violation and refined_clearance_violation_time_s is not None:
        _note("koz_violation", refined_clearance_violation_time_s)
    if chief_crossing_violation and "chief_crossing_violation" not in reasons:
        idx = int(np.argmax(ys > (CHIEF_CROSSING_Y_LIMIT + CHIEF_CROSSING_TOL_M)))
        _note("chief_crossing_violation", float(times[idx]))

    violation_reason = "none" if not reasons else ",".join(reasons)
    safe = not (koz_violation or corridor_violation or chief_crossing_violation)

    corridor_entry_time_s = None
    for c, tt in zip(classifications, times):
        if c.inside_corridor:
            corridor_entry_time_s = float(tt)
            break

    return TrajectoryConstraintResult(
        t=t,
        sample_dt_s=sample_dt,
        safe=safe,
        min_distance_m=min_distance_m,
        time_of_min_distance_s=time_of_min_distance_s,
        min_y_m=min_y_m,
        max_y_m=max_y_m,
        min_clearance_outside_corridor_m=min_clearance_outside_corridor_m,
        koz_violation=koz_violation,
        corridor_violation=corridor_violation,
        chief_crossing_violation=chief_crossing_violation,
        first_violation_time_s=first_violation_time_s,
        violation_reason=violation_reason,
        corridor_entry_time_s=corridor_entry_time_s,
    )
