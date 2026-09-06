"""M5 verification suite for the nonlinear two-body validation model
(nonlinear.py). Covers task checks A-G.

nonlinear.py never calls cw.py; comparisons here are against the
already-verified M2 CW solver/STM only to quantify agreement, not to
derive the nonlinear model from it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rendezvous_cw.constraints import classify_trajectory_constraints
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.nonlinear import (
    chief_circular_initial_state,
    classify_nonlinear_trajectory,
    inertial_to_lvlh,
    lvlh_basis,
    lvlh_to_inertial,
    propagate_relative_nonlinear,
    propagate_two_body,
)
from rendezvous_cw.orbit import CHIEF_ALTITUDE, MU_EARTH, orbital_radius
from rendezvous_cw.rendezvous import solve_two_impulse

R0 = np.array([0.0, -1000.0, 30.0])
V0_MINUS = np.array([0.0, 0.0, 0.0])
RF = np.array([0.0, -50.0, 0.0])
VF_PLUS = np.array([0.0, 0.0, 0.0])

M4_SELECTED_T = 386.56385551988984
M4_RECOMMENDED_T = 380.0  # final M5-recommended local-robustness candidate


@pytest.fixture(scope="module")
def a_chief() -> float:
    return orbital_radius(CHIEF_ALTITUDE)


# ---------------------------------------------------------------------------
# A. Nonlinear circular-chief propagation (conservation checks)
# ---------------------------------------------------------------------------


def test_chief_orbit_conserves_radius_energy_angular_momentum_and_period(a_chief, chief_period):
    chief0 = chief_circular_initial_state(a_chief)
    times = np.linspace(0.0, chief_period, 400)
    sol = propagate_two_body(chief0, (0.0, chief_period), t_eval=times)
    assert sol.success

    r = sol.y[:3, :]
    v = sol.y[3:, :]
    radii = np.linalg.norm(r, axis=0)
    speeds = np.linalg.norm(v, axis=0)

    # Radius stays constant (circular orbit) to a tight absolute tolerance.
    assert (radii.max() - radii.min()) < 1e-3  # meters, over a 6778 km orbit

    # Specific orbital energy eps = v^2/2 - mu/r is conserved.
    energy = 0.5 * speeds**2 - MU_EARTH / radii
    assert (energy.max() - energy.min()) < 1e-3  # relative to |energy| ~ 3e7

    # Specific angular momentum vector magnitude conserved.
    h = np.cross(r.T, v.T)
    h_mag = np.linalg.norm(h, axis=1)
    assert (h_mag.max() - h_mag.min()) < 1.0  # relative to |h| ~ 5.2e10

    # After one full chief period, position returns close to the start.
    np.testing.assert_allclose(r[:, -1], chief0[:3], atol=1e-2)


# ---------------------------------------------------------------------------
# B. Frame transformation round trips
# ---------------------------------------------------------------------------


def test_lvlh_basis_orthonormal_and_right_handed(a_chief):
    chief0 = chief_circular_initial_state(a_chief)
    C, omega = lvlh_basis(chief0[:3], chief0[3:])
    x_hat, y_hat, z_hat = C[:, 0], C[:, 1], C[:, 2]

    for v in (x_hat, y_hat, z_hat):
        assert np.linalg.norm(v) == pytest.approx(1.0, abs=1e-12)
    assert np.dot(x_hat, y_hat) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(y_hat, z_hat) == pytest.approx(0.0, abs=1e-12)
    assert np.dot(z_hat, x_hat) == pytest.approx(0.0, abs=1e-12)

    # Right-handedness: x_hat cross y_hat == z_hat.
    np.testing.assert_allclose(np.cross(x_hat, y_hat), z_hat, atol=1e-12)
    assert np.linalg.det(C) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize(
    "r_rel,v_rel",
    [
        (np.array([10.0, -500.0, 5.0]), np.array([0.01, -0.02, 0.005])),
        (np.array([-30.0, 800.0, -12.0]), np.array([0.0, 0.0, 0.0])),
        (np.array([0.0, -1000.0, 30.0]), np.array([0.5, -0.2, 0.03])),
    ],
)
def test_round_trip_lvlh_to_eci_to_lvlh(a_chief, r_rel, v_rel):
    chief0 = chief_circular_initial_state(a_chief)
    r_c, v_c = chief0[:3], chief0[3:]
    r_dep, v_dep = lvlh_to_inertial(r_c, v_c, r_rel, v_rel)
    r_rel2, v_rel2 = inertial_to_lvlh(r_c, v_c, r_dep, v_dep)
    np.testing.assert_allclose(r_rel2, r_rel, atol=1e-9)
    np.testing.assert_allclose(v_rel2, v_rel, atol=1e-9)  # includes the omega x r term


def test_zero_relative_state_maps_to_chief_state(a_chief):
    chief0 = chief_circular_initial_state(a_chief)
    r_c, v_c = chief0[:3], chief0[3:]
    r_dep, v_dep = lvlh_to_inertial(r_c, v_c, np.zeros(3), np.zeros(3))
    np.testing.assert_allclose(r_dep, r_c, atol=1e-9)
    np.testing.assert_allclose(v_dep, v_c, atol=1e-9)


def test_plus_x_is_radially_outward(a_chief):
    chief0 = chief_circular_initial_state(a_chief)
    r_c, v_c = chief0[:3], chief0[3:]
    r_dep, _ = lvlh_to_inertial(r_c, v_c, np.array([100.0, 0.0, 0.0]), np.zeros(3))
    delta = r_dep - r_c
    r_hat = r_c / np.linalg.norm(r_c)
    np.testing.assert_allclose(delta, 100.0 * r_hat, atol=1e-9)


def test_plus_y_is_along_track(a_chief):
    chief0 = chief_circular_initial_state(a_chief)
    r_c, v_c = chief0[:3], chief0[3:]
    r_dep, _ = lvlh_to_inertial(r_c, v_c, np.array([0.0, 100.0, 0.0]), np.zeros(3))
    delta = r_dep - r_c
    v_hat = v_c / np.linalg.norm(v_c)
    np.testing.assert_allclose(delta, 100.0 * v_hat, atol=1e-9)


def test_plus_z_is_orbit_normal(a_chief):
    chief0 = chief_circular_initial_state(a_chief)
    r_c, v_c = chief0[:3], chief0[3:]
    r_dep, _ = lvlh_to_inertial(r_c, v_c, np.array([0.0, 0.0, 100.0]), np.zeros(3))
    delta = r_dep - r_c
    h_hat = np.cross(r_c, v_c) / np.linalg.norm(np.cross(r_c, v_c))
    np.testing.assert_allclose(delta, 100.0 * h_hat, atol=1e-9)


# ---------------------------------------------------------------------------
# C. CW small-time consistency
# ---------------------------------------------------------------------------


def test_nonlinear_matches_cw_closely_for_short_transfer_small_separation(chief_n, a_chief):
    r0_small = np.array([0.0, -50.0, 1.5])  # 1/20th of the M1 separation
    result = solve_two_impulse(r0_small, V0_MINUS, np.array([0.0, -10.0, 0.0]), VF_PLUS, 60.0, chief_n)
    times = np.linspace(0.0, 60.0, 20)
    nl = propagate_relative_nonlinear(r0_small, result.v0_plus, times, a_chief)
    state0 = np.concatenate([r0_small, result.v0_plus])
    cw_states = np.array([propagate_cw(state0, t, chief_n) for t in times])
    pos_err = np.linalg.norm(nl.lvlh_states[:, :3] - cw_states[:, :3], axis=1)
    assert pos_err.max() < 0.01  # meters -- very close agreement expected


# ---------------------------------------------------------------------------
# D. Separation-scaling trend
# ---------------------------------------------------------------------------


def test_smaller_separation_does_not_produce_larger_cw_error(chief_n, a_chief):
    r0_unit_mag = np.linalg.norm(R0)
    T = 1800.0
    errs = {}
    for mag in (250.0, 1000.0):
        r0_scaled = R0 * (mag / r0_unit_mag)
        result = solve_two_impulse(r0_scaled, V0_MINUS, RF, VF_PLUS, T, chief_n)
        times = np.linspace(0.0, T, 40)
        nl = propagate_relative_nonlinear(r0_scaled, result.v0_plus, times, a_chief)
        state0 = np.concatenate([r0_scaled, result.v0_plus])
        cw_states = np.array([propagate_cw(state0, t, chief_n) for t in times])
        pos_err = np.linalg.norm(nl.lvlh_states[:, :3] - cw_states[:, :3], axis=1)
        errs[mag] = pos_err[-1]
    assert errs[250.0] < errs[1000.0]


# ---------------------------------------------------------------------------
# E. Integrator tolerance convergence
# ---------------------------------------------------------------------------


def test_nonlinear_vs_cw_error_not_dominated_by_integrator_tolerance(chief_n, a_chief):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M4_RECOMMENDED_T, chief_n)
    times = np.linspace(0.0, M4_RECOMMENDED_T, 30)
    state0 = np.concatenate([R0, result.v0_plus])
    cw_states = np.array([propagate_cw(state0, t, chief_n) for t in times])

    errs = []
    for rtol, atol in [(1e-12, 1e-9), (1e-13, 1e-10)]:
        nl = propagate_relative_nonlinear(R0, result.v0_plus, times, a_chief, rtol=rtol, atol=atol)
        pos_err = np.linalg.norm(nl.lvlh_states[:, :3] - cw_states[:, :3], axis=1)
        errs.append(pos_err[-1])

    # The terminal CW-vs-nonlinear discrepancy should be essentially
    # identical across integrator tolerances (agreeing to a small fraction
    # of its own magnitude), confirming it reflects real linearization
    # error, not integrator noise.
    assert abs(errs[0] - errs[1]) < 0.01 * max(errs[0], errs[1]) + 1e-6


# ---------------------------------------------------------------------------
# F. M4 regression preservation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t,expected_mm_s", [(1800.0, 1064.5203), (M4_SELECTED_T, 5052.7305)])
def test_m4_cw_values_unchanged_by_m5_module(chief_n, t, expected_mm_s):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, t, chief_n)
    cc = classify_trajectory_constraints(R0, result.v0_plus, t, chief_n, sample_dt=1.0)
    assert result.dv_total * 1000.0 == pytest.approx(expected_mm_s, abs=0.05)
    if t == M4_SELECTED_T:
        assert cc.safe  # M4's own regression: still feasible under CW


# ---------------------------------------------------------------------------
# G. Final selected-case geometry (M5-recommended transfer)
# ---------------------------------------------------------------------------


def test_m5_recommended_transfer_cw_metrics(chief_n):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M4_RECOMMENDED_T, chief_n)
    cc = classify_trajectory_constraints(R0, result.v0_plus, M4_RECOMMENDED_T, chief_n, sample_dt=1.0)
    assert result.dv_total * 1000.0 == pytest.approx(5136.02, abs=0.5)
    assert cc.safe
    assert cc.min_clearance_outside_corridor_m == pytest.approx(0.8752, abs=1e-3)


def test_m5_recommended_transfer_nonlinear_terminal_error(chief_n, a_chief):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M4_RECOMMENDED_T, chief_n)
    times = np.array([0.0, M4_RECOMMENDED_T])
    nl = propagate_relative_nonlinear(R0, result.v0_plus, times, a_chief)
    pos_err = np.linalg.norm(nl.lvlh_states[-1, :3] - RF)
    assert pos_err == pytest.approx(0.0109, abs=0.002)


def test_m5_recommended_transfer_nonlinear_constraint_classification(chief_n, a_chief):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M4_RECOMMENDED_T, chief_n)
    nl_cc = classify_nonlinear_trajectory(R0, result.v0_plus, M4_RECOMMENDED_T, a_chief, sample_dt=1.0)
    assert nl_cc.safe
    assert nl_cc.min_clearance_outside_corridor_m == pytest.approx(0.9022, abs=1e-3)


def test_m4_boundary_optimum_nonlinear_margin_is_small_but_positive(chief_n, a_chief):
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M4_SELECTED_T, chief_n)
    nl_cc = classify_nonlinear_trajectory(R0, result.v0_plus, M4_SELECTED_T, a_chief, sample_dt=1.0)
    assert nl_cc.safe  # passes, but only just
    assert 0.0 < nl_cc.min_clearance_outside_corridor_m < 0.1  # small margin, meters
