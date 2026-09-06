"""Shared fixtures / reference data for the M2 verification suite."""

from __future__ import annotations

import numpy as np
import pytest

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rendezvous_cw.orbit import chief_mean_motion_and_period  # noqa: E402


@pytest.fixture(scope="session")
def chief_n() -> float:
    """M1/M2 baseline chief mean motion, rad/s (400 km circular LEO)."""
    n, _ = chief_mean_motion_and_period()
    return n


@pytest.fixture(scope="session")
def chief_period() -> float:
    """M1/M2 baseline chief orbital period, s (400 km circular LEO)."""
    _, p = chief_mean_motion_and_period()
    return p


@pytest.fixture(scope="session")
def m1_scenario(chief_n):
    """Exact M1 scenario state vectors (DESIGN.md section 2), meters/m/s."""
    return {
        "n": chief_n,
        "r0": np.array([0.0, -1000.0, 30.0]),
        "v0_minus": np.array([0.0, 0.0, 0.0]),
        "rf": np.array([0.0, -50.0, 0.0]),
        "vf_plus": np.array([0.0, 0.0, 0.0]),
        "T": 1800.0,
    }


def cw_ode_rhs(t, state, n):
    """Raw first-order CW ODE right-hand side (reference implementation).

    Independent of the closed-form STM in cw.py -- used only in the
    verification suite (see DESIGN.md section 9, checks B and J), never
    in the production propagation path.

        xdot, ydot, zdot, 2n*ydot + 3n^2*x, -2n*xdot, -n^2*z
    """
    x, y, z, vx, vy, vz = state
    ax = 2.0 * n * vy + 3.0 * n**2 * x
    ay = -2.0 * n * vx
    az = -(n**2) * z
    return np.array([vx, vy, vz, ax, ay, az])
