"""Two-impulse CW boundary-value rendezvous solver.

Implements the method of DESIGN.md section 6:

    v0_plus  = Phi_rv(T)^-1 @ (rf - Phi_rr(T) @ r0)      [solved, not inverted]
    dv1      = v0_plus - v0_minus
    vT_minus = Phi_vr(T) @ r0 + Phi_vv(T) @ v0_plus
    dv2      = vf_plus - vT_minus
    dv_total = norm(dv1) + norm(dv2)

Phi_rv(T) is never explicitly inverted; `numpy.linalg.solve` is used
instead (better numerically conditioned than forming an explicit
inverse), and is only attempted after the conditioning policy in
conditioning.py has cleared the transfer time as safe.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .conditioning import PhiRvConditioning, phi_rv_conditioning
from .cw import Phi_rr, Phi_rv, Phi_vr, Phi_vv


class SingularTransferError(ValueError):
    """Raised when the requested transfer time is singular/unsafe for Phi_rv."""


@dataclass(frozen=True)
class TwoImpulseResult:
    """Structured result of a two-impulse CW rendezvous solve."""

    t: float
    n: float
    r0: np.ndarray
    v0_minus: np.ndarray
    rf: np.ndarray
    vf_plus: np.ndarray

    v0_plus: np.ndarray
    dv1: np.ndarray
    dv1_mag: float

    vT_minus: np.ndarray
    dv2: np.ndarray
    dv2_mag: float

    dv_total: float

    final_state: np.ndarray  # propagated [r(T); v(T)] using v0_plus
    terminal_position_residual: np.ndarray  # final_state[:3] - rf
    terminal_position_residual_norm: float

    conditioning: PhiRvConditioning


def solve_two_impulse(
    r0: np.ndarray,
    v0_minus: np.ndarray,
    rf: np.ndarray,
    vf_plus: np.ndarray,
    t: float,
    n: float,
) -> TwoImpulseResult:
    """Solve the two-impulse CW rendezvous boundary-value problem.

    Parameters
    ----------
    r0 : initial relative position (3,), meters.
    v0_minus : initial (pre-burn) relative velocity (3,), m/s.
    rf : target relative position at time t (3,), meters.
    vf_plus : desired (post-second-burn) relative velocity at time t (3,), m/s.
    t : transfer time, seconds (t > 0).
    n : chief mean motion, rad/s.

    Returns
    -------
    TwoImpulseResult

    Raises
    ------
    SingularTransferError
        If Phi_rv(t, n) is singular or exceeds the condition-number policy
        threshold (conditioning.COND_THRESHOLD) -- see conditioning.py.
    ValueError
        If any input vector is not shape (3,), or t <= 0.
    """
    r0 = np.asarray(r0, dtype=float)
    v0_minus = np.asarray(v0_minus, dtype=float)
    rf = np.asarray(rf, dtype=float)
    vf_plus = np.asarray(vf_plus, dtype=float)
    for name, vec in (("r0", r0), ("v0_minus", v0_minus), ("rf", rf), ("vf_plus", vf_plus)):
        if vec.shape != (3,):
            raise ValueError(f"{name} must be shape (3,), got {vec.shape}")
    if t <= 0:
        raise ValueError(f"transfer time t must be > 0, got {t}")

    cond = phi_rv_conditioning(t, n)
    if cond.unsafe:
        raise SingularTransferError(
            f"Transfer time t={t!r} s is singular/unsafe for Phi_rv "
            f"(n*t={n * t!r} rad, det={cond.det!r}, cond={cond.cond!r}, "
            f"threshold={cond.cond_threshold!r}). Refusing to solve; choose "
            f"a transfer time away from n*t = k*pi (DESIGN.md section 6.1)."
        )

    prr = Phi_rr(t, n)
    prv = Phi_rv(t, n)
    pvr = Phi_vr(t, n)
    pvv = Phi_vv(t, n)

    rhs = rf - prr @ r0
    v0_plus = np.linalg.solve(prv, rhs)

    dv1 = v0_plus - v0_minus
    dv1_mag = float(np.linalg.norm(dv1))

    vT_minus = pvr @ r0 + pvv @ v0_plus
    dv2 = vf_plus - vT_minus
    dv2_mag = float(np.linalg.norm(dv2))

    dv_total = dv1_mag + dv2_mag

    r_final = prr @ r0 + prv @ v0_plus
    v_final = pvr @ r0 + pvv @ v0_plus  # == vT_minus, before the 2nd impulse
    final_state = np.concatenate([r_final, v_final])

    residual = r_final - rf
    residual_norm = float(np.linalg.norm(residual))

    return TwoImpulseResult(
        t=t,
        n=n,
        r0=r0,
        v0_minus=v0_minus,
        rf=rf,
        vf_plus=vf_plus,
        v0_plus=v0_plus,
        dv1=dv1,
        dv1_mag=dv1_mag,
        vT_minus=vT_minus,
        dv2=dv2,
        dv2_mag=dv2_mag,
        dv_total=dv_total,
        final_state=final_state,
        terminal_position_residual=residual,
        terminal_position_residual_norm=residual_norm,
        conditioning=cond,
    )
