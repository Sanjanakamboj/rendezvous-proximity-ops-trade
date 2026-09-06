"""Transfer-time / Delta-v trade-study utilities (Milestone 3).

This module does NOT re-derive or duplicate any CW/rendezvous equations.
Every Delta-v number here comes from calling the existing, verified M2
solver (`rendezvous.solve_two_impulse`) and M2 conditioning diagnostics
(`conditioning.phi_rv_conditioning`) directly. This module only:

  - wraps a single-transfer-time evaluation in a structured record that
    also carries the safe/unsafe flag instead of raising for unsafe
    times (so a sweep can carry unsafe points through as explicit
    excluded/NaN entries rather than crashing or being silently
    skipped),
  - sweeps a set of transfer times and collects those records,
  - searches the *already-evaluated safe records* for a minimum-Delta-v
    time, then locally refines that single time with a bounded scalar
    minimizer restricted to a safe sub-interval (never spanning the
    singularity -- see find_best_safe_transfer).

Scope: fixed M1 boundary conditions (r0, v0_minus, rf, vf_plus) and the
fixed 400 km circular chief orbit are assumed throughout; nothing here
changes the M1/M2 model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence

import numpy as np
from scipy.optimize import minimize_scalar

from .conditioning import COND_THRESHOLD, phi_rv_conditioning
from .rendezvous import SingularTransferError, TwoImpulseResult, solve_two_impulse


@dataclass(frozen=True)
class TradeRecord:
    """One transfer-time evaluation for the Delta-v/time trade."""

    t_s: float
    t_min: float  # t_s / 60
    t_over_p: float  # t_s / chief_period
    safe: bool
    cond_phi_rv: float
    det_phi_rv: float

    # None for unsafe points (solver was not invoked / result is not
    # meaningful) -- keeps unsafe points explicit rather than fabricating
    # a number.
    dv1: Optional[np.ndarray] = None
    dv2: Optional[np.ndarray] = None
    dv1_mag: float = float("nan")
    dv2_mag: float = float("nan")
    dv_total: float = float("nan")
    closure_residual_m: float = float("nan")

    notes: str = ""


def evaluate_transfer_time(
    r0: np.ndarray,
    v0_minus: np.ndarray,
    rf: np.ndarray,
    vf_plus: np.ndarray,
    t: float,
    n: float,
    chief_period: float,
) -> TradeRecord:
    """Evaluate one transfer time via the M2 solver; never raises.

    Returns a TradeRecord with `safe=False` and NaN Delta-v fields (but
    real `cond_phi_rv`/`det_phi_rv` diagnostics) when the M2 conditioning
    policy (conditioning.COND_THRESHOLD) marks the transfer unsafe,
    instead of propagating SingularTransferError -- so a sweep can record
    every requested time explicitly.
    """
    diag = phi_rv_conditioning(t, n)
    t_over_p = t / chief_period

    if diag.unsafe:
        return TradeRecord(
            t_s=t,
            t_min=t / 60.0,
            t_over_p=t_over_p,
            safe=False,
            cond_phi_rv=diag.cond,
            det_phi_rv=diag.det,
            notes=f"unsafe: cond(Phi_rv)={diag.cond:.3e} > {COND_THRESHOLD:.1e}",
        )

    try:
        result: TwoImpulseResult = solve_two_impulse(r0, v0_minus, rf, vf_plus, t, n)
    except SingularTransferError as exc:
        # Defensive: should not happen given the diag.unsafe check above
        # (both paths use the same conditioning.py policy), but recorded
        # explicitly rather than silently swallowed if it ever does.
        return TradeRecord(
            t_s=t,
            t_min=t / 60.0,
            t_over_p=t_over_p,
            safe=False,
            cond_phi_rv=diag.cond,
            det_phi_rv=diag.det,
            notes=f"solver raised SingularTransferError: {exc}",
        )

    return TradeRecord(
        t_s=t,
        t_min=t / 60.0,
        t_over_p=t_over_p,
        safe=True,
        cond_phi_rv=result.conditioning.cond,
        det_phi_rv=result.conditioning.det,
        dv1=result.dv1,
        dv2=result.dv2,
        dv1_mag=result.dv1_mag,
        dv2_mag=result.dv2_mag,
        dv_total=result.dv_total,
        closure_residual_m=result.terminal_position_residual_norm,
    )


def sweep_transfer_times(
    r0: np.ndarray,
    v0_minus: np.ndarray,
    rf: np.ndarray,
    vf_plus: np.ndarray,
    times: Iterable[float],
    n: float,
    chief_period: float,
) -> list[TradeRecord]:
    """Evaluate `evaluate_transfer_time` for each t in `times`.

    Each time is evaluated independently (no shared mutable state), so
    the returned list does not depend on the order `times` is given in.
    """
    return [
        evaluate_transfer_time(r0, v0_minus, rf, vf_plus, t, n, chief_period)
        for t in times
    ]


def find_best_safe_transfer(
    records: Sequence[TradeRecord],
    r0: np.ndarray,
    v0_minus: np.ndarray,
    rf: np.ndarray,
    vf_plus: np.ndarray,
    n: float,
    chief_period: float,
    refine: bool = True,
) -> TradeRecord:
    """Find the minimum-Delta-v transfer among `records`, optionally refined.

    Only considers records with `safe=True`. The refinement step brackets
    the coarse minimum using its two immediate (safe) sweep neighbors and
    runs a bounded scalar minimizer (`scipy.optimize.minimize_scalar`,
    method="bounded") strictly within that bracket -- so refinement never
    crosses into an unsafe/singular sub-interval, and never searches
    outside the interval actually swept by `records`.

    This is a search over the explicitly scanned interval only; see
    DESIGN.md section 13 for the precise, deliberately narrow meaning of
    "minimum" used throughout this module (not a claim of global or
    mission-level optimality).
    """
    safe_records = [r for r in records if r.safe]
    if not safe_records:
        raise ValueError("no safe records to search")

    coarse_best = min(safe_records, key=lambda r: r.dv_total)

    if not refine:
        return coarse_best

    sorted_safe = sorted(safe_records, key=lambda r: r.t_s)
    idx = sorted_safe.index(coarse_best)
    lo = sorted_safe[idx - 1].t_s if idx > 0 else coarse_best.t_s
    hi = sorted_safe[idx + 1].t_s if idx < len(sorted_safe) - 1 else coarse_best.t_s

    if hi <= lo:
        return coarse_best

    def objective(t: float) -> float:
        rec = evaluate_transfer_time(r0, v0_minus, rf, vf_plus, t, n, chief_period)
        return rec.dv_total if rec.safe else np.inf

    res = minimize_scalar(objective, bounds=(lo, hi), method="bounded")
    refined_t = float(res.x)

    refined_record = evaluate_transfer_time(
        r0, v0_minus, rf, vf_plus, refined_t, n, chief_period
    )
    if not refined_record.safe or refined_record.dv_total > coarse_best.dv_total:
        # Refinement failed to improve (or landed unsafe); fall back to
        # the coarse grid best rather than reporting a worse "refined" answer.
        return coarse_best
    return refined_record
