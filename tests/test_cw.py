"""M2 verification suite for the closed-form CW STM (cw.py).

Covers DESIGN.md section 9 checks A (t=0 identity), B (ODE residual),
C (STM composition), D (planar invariance), E (exact cross-track
oscillator), and K (time-reversal / inverse STM).
"""

from __future__ import annotations

import numpy as np
import pytest

from rendezvous_cw.cw import Phi, Phi_rr, Phi_rv, Phi_vr, Phi_vv, propagate_cw

from conftest import cw_ode_rhs

ATOL = 1e-9


# ---------------------------------------------------------------------------
# A. STM identity at t=0
# ---------------------------------------------------------------------------


def test_phi_rr_zero_is_identity(chief_n):
    np.testing.assert_allclose(Phi_rr(0.0, chief_n), np.eye(3), atol=ATOL)


def test_phi_vv_zero_is_identity(chief_n):
    np.testing.assert_allclose(Phi_vv(0.0, chief_n), np.eye(3), atol=ATOL)


def test_phi_rv_zero_is_zero(chief_n):
    np.testing.assert_allclose(Phi_rv(0.0, chief_n), np.zeros((3, 3)), atol=ATOL)


def test_phi_vr_zero_is_zero(chief_n):
    np.testing.assert_allclose(Phi_vr(0.0, chief_n), np.zeros((3, 3)), atol=ATOL)


def test_full_stm_zero_is_identity6(chief_n):
    np.testing.assert_allclose(Phi(0.0, chief_n), np.eye(6), atol=ATOL)


# ---------------------------------------------------------------------------
# B. CW differential-equation residual
#
# Independent check: numerically (finite-difference) differentiate the
# STM-propagated trajectory's velocity components to obtain accelerations,
# then compare against the CW ODE right-hand side (a separately written
# formula, not the STM closed form) evaluated at the same state. This does
# not compare the STM against itself: the accelerations come from a
# central-difference derivative of numerically evaluated states, and the
# "expected" accelerations come from the independent cw_ode_rhs formula.
# ---------------------------------------------------------------------------

_TEST_STATES = [
    np.array([100.0, -500.0, 20.0, 0.05, -0.02, 0.01]),
    np.array([0.0, -1000.0, 30.0, 0.0, 0.0, 0.0]),
    np.array([-250.0, 300.0, -40.0, -0.03, 0.04, -0.02]),
]
_TEST_TIMES = [0.0, 250.0, 900.0, 1800.0, 4000.0]


@pytest.mark.parametrize("state0", _TEST_STATES)
@pytest.mark.parametrize("t", _TEST_TIMES)
def test_cw_ode_residual(chief_n, state0, t):
    n = chief_n
    h = 1e-2  # seconds, central-difference step

    s_minus = propagate_cw(state0, t - h, n)
    s_plus = propagate_cw(state0, t + h, n)
    s_mid = propagate_cw(state0, t, n)

    # Numerically differentiate the propagated state -> [xdot..vzdot]
    d_state_dt_numeric = (s_plus - s_minus) / (2.0 * h)

    # Independent formula: d/dt [x,y,z,vx,vy,vz] = [vx,vy,vz, ax,ay,az]
    d_state_dt_expected = cw_ode_rhs(t, s_mid, n)

    np.testing.assert_allclose(
        d_state_dt_numeric, d_state_dt_expected, atol=1e-6, rtol=1e-6
    )


# ---------------------------------------------------------------------------
# C. STM composition: Phi(t1+t2) == Phi(t2) @ Phi(t1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "t1,t2",
    [
        (100.0, 200.0),
        (900.0, 900.0),
        (1800.0, 500.0),
        (2000.0, 3000.0),
        (50.0, 5000.0),
    ],
)
def test_stm_composition(chief_n, t1, t2):
    n = chief_n
    lhs = Phi(t1 + t2, n)
    rhs = Phi(t2, n) @ Phi(t1, n)
    np.testing.assert_allclose(lhs, rhs, atol=1e-6, rtol=1e-9)


# ---------------------------------------------------------------------------
# D. Planar invariance: z0 = vz0 = 0 => z(t) = vz(t) = 0 for all t
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t", [10.0, 500.0, 1800.0, 4000.0, 5500.0])
def test_planar_invariance(chief_n, t):
    state0 = np.array([120.0, -800.0, 0.0, 0.01, -0.05, 0.0])
    state_t = propagate_cw(state0, t, chief_n)
    assert abs(state_t[2]) < 1e-9  # z
    assert abs(state_t[5]) < 1e-9  # vz


# ---------------------------------------------------------------------------
# E. Exact cross-track harmonic oscillator
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("z0,vz0", [(30.0, 0.0), (0.0, 0.5), (-15.0, 0.2), (30.0, -0.1)])
@pytest.mark.parametrize("t", [0.0, 300.0, 1800.0, 2776.812135626114, 5000.0])
def test_exact_cross_track_oscillator(chief_n, z0, vz0, t):
    n = chief_n
    state0 = np.array([0.0, 0.0, z0, 0.0, 0.0, vz0])
    state_t = propagate_cw(state0, t, n)

    z_expected = z0 * np.cos(n * t) + (vz0 / n) * np.sin(n * t)
    vz_expected = -z0 * n * np.sin(n * t) + vz0 * np.cos(n * t)

    assert state_t[2] == pytest.approx(z_expected, abs=1e-9)
    assert state_t[5] == pytest.approx(vz_expected, abs=1e-9)


# ---------------------------------------------------------------------------
# K. Time-reversal / inverse STM: Phi(-t) @ Phi(t) == I
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t", [50.0, 300.0, 1800.0, 4000.0, 5000.0])
def test_time_reversal_inverse(chief_n, t):
    n = chief_n
    lhs = Phi(-t, n) @ Phi(t, n)
    np.testing.assert_allclose(lhs, np.eye(6), atol=1e-6, rtol=1e-9)
