"""Clohessy-Wiltshire (CW) / Hill-frame closed-form state-transition matrix.

Implements the four 3x3 STM blocks and the full 6x6 STM exactly as derived
in DESIGN.md section 5, for a deputy in linearized relative motion about a
circular chief orbit of mean motion `n`. State convention throughout this
project (DESIGN.md section 3):

    state = [x, y, z, vx, vy, vz]

with x radial, y along-track, z cross-track (right-handed LVLH/Hill frame).

These are closed-form expressions (not derived here by matrix
exponentiation or ODE integration), so they are independently
cross-checked in the test suite against:
  - the CW ODE residual,
  - STM composition Phi(t1+t2) = Phi(t2) @ Phi(t1),
  - the exact cross-track harmonic solution,
  - an independent scipy.integrate.solve_ivp propagation of the raw ODE,
  - the t -> 0 identity limit,
  - Phi(-t) @ Phi(t) = I (time reversal).

No numerical ODE solver is used in the production propagation path; that
is reserved for the verification suite only (per M2 scope).
"""

from __future__ import annotations

import numpy as np


def Phi_rr(t: float, n: float) -> np.ndarray:
    """Position-due-to-initial-position block (3x3), DESIGN.md section 5."""
    c = np.cos(n * t)
    s = np.sin(n * t)
    nt = n * t
    return np.array(
        [
            [4.0 - 3.0 * c, 0.0, 0.0],
            [6.0 * (s - nt), 1.0, 0.0],
            [0.0, 0.0, c],
        ]
    )


def Phi_rv(t: float, n: float) -> np.ndarray:
    """Position-due-to-initial-velocity block (3x3), DESIGN.md section 5.

    Singular (exactly) whenever sin(n*t) = 0, i.e. n*t = k*pi for integer
    k >= 0 -- see conditioning.py for the detection/rejection policy used
    before this block is ever inverted (via linear solve) by the
    two-impulse solver.
    """
    c = np.cos(n * t)
    s = np.sin(n * t)
    nt = n * t
    return np.array(
        [
            [s / n, 2.0 * (1.0 - c) / n, 0.0],
            [-2.0 * (1.0 - c) / n, (4.0 * s - 3.0 * nt) / n, 0.0],
            [0.0, 0.0, s / n],
        ]
    )


def Phi_vr(t: float, n: float) -> np.ndarray:
    """Velocity-due-to-initial-position block (3x3), DESIGN.md section 5."""
    c = np.cos(n * t)
    s = np.sin(n * t)
    return np.array(
        [
            [3.0 * n * s, 0.0, 0.0],
            [6.0 * n * (c - 1.0), 0.0, 0.0],
            [0.0, 0.0, -n * s],
        ]
    )


def Phi_vv(t: float, n: float) -> np.ndarray:
    """Velocity-due-to-initial-velocity block (3x3), DESIGN.md section 5."""
    c = np.cos(n * t)
    s = np.sin(n * t)
    return np.array(
        [
            [c, 2.0 * s, 0.0],
            [-2.0 * s, 4.0 * c - 3.0, 0.0],
            [0.0, 0.0, c],
        ]
    )


def Phi(t: float, n: float) -> np.ndarray:
    """Full 6x6 CW state-transition matrix at time `t` for mean motion `n`.

    state(t) = Phi(t, n) @ state0, with state = [x, y, z, vx, vy, vz].

    At t -> 0, Phi(0, n) = I6 exactly (Phi_rr -> I, Phi_vv -> I,
    Phi_rv -> 0, Phi_vr -> 0), verified in the test suite to machine
    precision.
    """
    top = np.hstack([Phi_rr(t, n), Phi_rv(t, n)])
    bottom = np.hstack([Phi_vr(t, n), Phi_vv(t, n)])
    return np.vstack([top, bottom])


def propagate_cw(state0: np.ndarray, t: float, n: float) -> np.ndarray:
    """Propagate a 6-vector CW state by time `t` using the closed-form STM.

    state0 = [x, y, z, vx, vy, vz]. Returns state(t) = Phi(t, n) @ state0.
    This is the sole production propagation path for M2 (no ODE solver).
    """
    state0 = np.asarray(state0, dtype=float)
    if state0.shape != (6,):
        raise ValueError(f"state0 must be shape (6,), got {state0.shape}")
    return Phi(t, n) @ state0
