"""M2 verification: two-impulse solver against the M1 scenario.

Covers DESIGN.md section 9 checks F (terminal position closure),
G (burn reconstruction), H (trivial zero-transfer-state case), plus the
M1 hand-calculation regression comparison (task section 7).
"""

from __future__ import annotations

import numpy as np
import pytest

from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.rendezvous import solve_two_impulse

# M1 DESIGN.md section 7.3 rounded hand-calculation targets (mm/s).
M1_DV1_MM_S = 531.7
M1_DV2_MM_S = 532.8
M1_DVTOTAL_MM_S = 1064.5


def test_m1_scenario_matches_hand_calculation_target(m1_scenario):
    result = solve_two_impulse(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        m1_scenario["T"],
        m1_scenario["n"],
    )

    dv1_mm_s = result.dv1_mag * 1000.0
    dv2_mm_s = result.dv2_mag * 1000.0
    dvtotal_mm_s = result.dv_total * 1000.0

    # M1 values were hand-rounded to 4 significant figures; agreement to
    # 0.1 mm/s is well within that rounding.
    assert dv1_mm_s == pytest.approx(M1_DV1_MM_S, abs=0.1)
    assert dv2_mm_s == pytest.approx(M1_DV2_MM_S, abs=0.1)
    assert dvtotal_mm_s == pytest.approx(M1_DVTOTAL_MM_S, abs=0.1)


# ---------------------------------------------------------------------------
# F. Terminal position closure
# ---------------------------------------------------------------------------


def test_terminal_position_closure(m1_scenario):
    result = solve_two_impulse(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        m1_scenario["T"],
        m1_scenario["n"],
    )
    assert result.terminal_position_residual_norm < 1e-8  # meters
    np.testing.assert_allclose(
        result.final_state[:3], m1_scenario["rf"], atol=1e-8
    )


# ---------------------------------------------------------------------------
# G. Burn reconstruction
# ---------------------------------------------------------------------------


def test_burn_reconstruction_consistency(m1_scenario):
    result = solve_two_impulse(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        m1_scenario["T"],
        m1_scenario["n"],
    )

    # dv1 = v0_plus - v0_minus
    np.testing.assert_allclose(
        result.dv1, result.v0_plus - m1_scenario["v0_minus"], atol=1e-12
    )

    # dv2 = vf_plus - vT_minus
    np.testing.assert_allclose(
        result.dv2, m1_scenario["vf_plus"] - result.vT_minus, atol=1e-12
    )

    # Independently propagate the post-first-burn state [r0, v0_plus]
    # forward by T using the STM directly, and confirm applying both
    # burns reproduces exactly the specified boundary states (rf, vf_plus).
    post_burn1_state = np.concatenate([m1_scenario["r0"], result.v0_plus])
    propagated = propagate_cw(post_burn1_state, m1_scenario["T"], m1_scenario["n"])
    r_at_T, v_at_T_minus = propagated[:3], propagated[3:]

    np.testing.assert_allclose(r_at_T, m1_scenario["rf"], atol=1e-8)
    np.testing.assert_allclose(v_at_T_minus, result.vT_minus, atol=1e-8)

    v_after_second_burn = v_at_T_minus + result.dv2
    np.testing.assert_allclose(
        v_after_second_burn, m1_scenario["vf_plus"], atol=1e-12
    )


# ---------------------------------------------------------------------------
# H. Trivial zero-transfer-state case
#
# A dynamically consistent no-burn case: choose r0 such that propagating
# it (with v0 = 0) under the *unforced* CW dynamics for time T already
# lands exactly on rf with vf = vT_minus, so the required v0_plus equals
# v0_minus (both zero) and dv1 = dv2 = 0. We construct rf/vf directly from
# the free propagation of r0=0, v0=0 (the trivial equilibrium: r stays at
# the origin for all time, since [0,0,0,0,0,0] is a fixed point of the CW
# ODE), which avoids inventing an unverified stationary CW state.
# ---------------------------------------------------------------------------


def test_trivial_zero_separation_requires_zero_delta_v(chief_n):
    r0 = np.array([0.0, 0.0, 0.0])
    v0 = np.array([0.0, 0.0, 0.0])
    T = 1800.0

    # The origin with zero velocity is a fixed point of the CW ODE
    # (verified separately by test_cw.py's ODE-residual check at the
    # zero state), so propagating it forward analytically reproduces
    # (0,0,0,0,0,0) again -- used here only to state rf/vf, not asserted.
    rf = np.array([0.0, 0.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])

    result = solve_two_impulse(r0, v0, rf, vf, T, chief_n)

    assert result.dv1_mag < 1e-9
    assert result.dv2_mag < 1e-9
    assert result.dv_total < 1e-9
    assert result.terminal_position_residual_norm < 1e-9
