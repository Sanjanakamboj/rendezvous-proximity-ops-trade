"""Phi_rv conditioning diagnostics and singularity policy.

The two-impulse solver (rendezvous.py) must invert (via linear solve)
Phi_rv(T, n), which is exactly singular whenever sin(n*T) = 0, i.e.
n*T = k*pi for integer k >= 0 (DESIGN.md section 6.1). This module
implements the detection policy so the solver can refuse transfer times
that are singular or numerically unsafe, rather than silently returning a
numerically produced but physically meaningless Delta-v.

Policy chosen (documented, not merely an absolute determinant cutoff):

    A transfer time T is UNSAFE if the 2-norm condition number of
    Phi_rv(T, n) exceeds COND_THRESHOLD = 1e6 (equivalently, the
    reciprocal condition number 1/cond falls below 1e-6).

Rationale:
  - The determinant alone is a poor safety signal: Phi_rv's entries scale
    with 1/n (~1e3 for a 400 km LEO chief), so det(Phi_rv) is already
    large (~1e9) at well-conditioned transfer times (see DESIGN.md
    section 7.2, T = 1800 s). An absolute determinant threshold would
    need to be re-tuned per orbit regime; the condition number is
    scale-invariant and directly measures how much a small error in the
    right-hand side (rf - Phi_rr @ r0) is amplified into v0_plus, which
    is exactly the failure mode we care about.
  - cond(Phi_rv) = 1e6 means a solve can lose up to ~6 decimal digits of
    accuracy relative to double precision's ~16 digits -- comfortably
    still meaningful for this project's mm/s-level Delta-v targets, while
    still being a conservative (not overly permissive) cutoff. This is a
    standard engineering rule of thumb for "safe to invert" (Golub & Van
    Loan-style loss-of-precision reasoning), not a value tuned to make
    any particular test pass.
  - Numerically (verified for the M1 400 km chief scenario, half-period
    singularity at T = 2776.812 s): cond(Phi_rv) is exactly ~8.9e16 at
    the singular point itself, falls to ~9.6e6 only 1 ms away, ~9.6e5 at
    10 ms, and ~9.6e3 at 1 s away. The COND_THRESHOLD = 1e6 policy above
    therefore rejects a guard band of roughly +/-10 ms around each
    singular time for this scenario -- a deliberately tight band, because
    the conditioning degrades extremely sharply (not gradually) near
    n*T = k*pi, so a tight guard band still excludes every numerically
    dangerous point while excluding a negligible sliver of the traded
    interval.

This threshold is a project-level policy constant, reused unchanged by
the M3 trade sweep (not yet implemented).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .cw import Phi_rv

#: Condition-number policy threshold (see module docstring for rationale).
COND_THRESHOLD = 1.0e6


@dataclass(frozen=True)
class PhiRvConditioning:
    """Conditioning diagnostics for Phi_rv(T, n) at a given transfer time."""

    t: float
    n: float
    det: float
    cond: float
    unsafe: bool
    cond_threshold: float = COND_THRESHOLD


def phi_rv_conditioning(t: float, n: float) -> PhiRvConditioning:
    """Compute determinant + condition number of Phi_rv(t, n) and flag safety.

    Returns a PhiRvConditioning record. `unsafe` is True when Phi_rv is
    exactly singular (det == 0) or when its 2-norm condition number
    exceeds COND_THRESHOLD (see module docstring for the policy and its
    justification).
    """
    m = Phi_rv(t, n)
    det = float(np.linalg.det(m))
    if det == 0.0:
        cond = float("inf")
    else:
        cond = float(np.linalg.cond(m))
    unsafe = (not np.isfinite(cond)) or (cond > COND_THRESHOLD)
    return PhiRvConditioning(t=t, n=n, det=det, cond=cond, unsafe=unsafe)
