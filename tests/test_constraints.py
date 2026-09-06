"""M4 verification suite for proximity-geometry constraints (constraints.py).

Covers the M4 task checks A-J. All Delta-v/trajectory numbers ultimately
come from the verified M2 solver (rendezvous.solve_two_impulse) and M2
closed-form STM (cw.propagate_cw) -- constraints.py adds only geometric
classification on top, and duplicates no dynamics.
"""

from __future__ import annotations

import random

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from rendezvous_cw.conditioning import phi_rv_conditioning
from rendezvous_cw.constraints import (
    R_KOZ,
    classify_point,
    classify_trajectory_constraints,
    inside_approach_corridor,
)
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.rendezvous import solve_two_impulse

from conftest import cw_ode_rhs

# M3 authoritative reference transfer times (DESIGN.md section 13).
M2_REF_T = 1800.0
M3_SHORT_MIN_T = 2539.592622
M3_LONG_MIN_T = 4868.130154465808


# ---------------------------------------------------------------------------
# A. Constraint primitive tests
# ---------------------------------------------------------------------------


def test_point_at_150m_outside_corridor_is_allowed():
    # 150 m along -y: outside the 100 m KOZ, outside the corridor's y-range,
    # and y <= -50 m so it does not separately trip chief-crossing either.
    r = np.array([0.0, -150.0, 0.0])
    c = classify_point(r)
    assert c.safe
    assert not c.koz_violation


def test_point_at_90m_outside_corridor_is_koz_violation():
    r = np.array([90.0, 0.0, 0.0])
    c = classify_point(r)
    assert not c.safe
    assert c.koz_violation


def test_point_at_0_neg75_0_is_inside_corridor_and_allowed():
    r = np.array([0.0, -75.0, 0.0])
    c = classify_point(r)
    assert inside_approach_corridor(r)
    assert c.safe
    assert not c.koz_violation
    assert not c.corridor_violation


def test_point_25_neg75_0_is_outside_lateral_corridor():
    r = np.array([25.0, -75.0, 0.0])
    assert not inside_approach_corridor(r)
    c = classify_point(r)
    # within 100 m, within y-range, but x=25 > 20 m lateral bound.
    assert not c.safe
    assert c.corridor_violation
    assert not c.koz_violation


def test_point_0_neg75_25_is_outside_lateral_corridor():
    r = np.array([0.0, -75.0, 25.0])
    assert not inside_approach_corridor(r)
    c = classify_point(r)
    assert not c.safe
    assert c.corridor_violation


def test_point_with_y_greater_than_neg50_is_chief_crossing_violation():
    r = np.array([0.0, -40.0, 0.0])  # y = -40 > -50
    c = classify_point(r)
    assert c.chief_crossing_violation
    assert not c.safe


def test_target_point_is_valid_endpoint():
    r = np.array([0.0, -50.0, 0.0])
    c = classify_point(r)
    assert inside_approach_corridor(r)
    assert c.safe
    assert not c.koz_violation
    assert not c.corridor_violation
    assert not c.chief_crossing_violation


# ---------------------------------------------------------------------------
# B. Endpoint tolerance
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "noise",
    [1e-13, -1e-13, 1e-10, -1e-10, 1e-7, -1e-7],
)
def test_target_endpoint_tolerant_to_floating_point_noise(noise):
    r = np.array([0.0 + noise, -50.0 + noise, 0.0 + noise])
    c = classify_point(r)
    assert c.safe, f"noise={noise} incorrectly rejected the valid endpoint"


def test_endpoint_tolerance_does_not_mask_real_violations():
    # A real, non-noise-scale violation must still be caught.
    r = np.array([0.0, -50.0 + 1.0, 0.0])  # 1 m over the line -- real
    c = classify_point(r)
    assert c.chief_crossing_violation
    assert not c.safe


# ---------------------------------------------------------------------------
# C. Continuous minimum vs. coarse sampling
# ---------------------------------------------------------------------------


def test_refinement_recovers_true_minimum_missed_by_coarse_sampling(chief_n):
    # Use the M3 short-branch minimum transfer, which is known (from the
    # M4 exploration) to dip closer to the chief than either 5 s-spaced
    # coarse sample straddling its true closest-approach time. A very
    # coarse sample (50 s) is virtually guaranteed to miss the exact
    # continuous minimum; refinement must find a value <= the coarse one.
    result = solve_two_impulse(
        np.array([0.0, -1000.0, 30.0]), np.array([0.0, 0.0, 0.0]),
        np.array([0.0, -50.0, 0.0]), np.array([0.0, 0.0, 0.0]),
        M3_LONG_MIN_T, chief_n,
    )
    coarse = classify_trajectory_constraints(
        np.array([0.0, -1000.0, 30.0]), result.v0_plus, M3_LONG_MIN_T, chief_n,
        sample_dt=50.0, refine=False,
    )
    fine = classify_trajectory_constraints(
        np.array([0.0, -1000.0, 30.0]), result.v0_plus, M3_LONG_MIN_T, chief_n,
        sample_dt=50.0, refine=True,
    )
    assert fine.min_distance_m <= coarse.min_distance_m + 1e-9
    # The refined continuous minimum should be meaningfully tighter than
    # what the same coarse grid saw without refinement, for this case.
    assert fine.min_distance_m < coarse.min_distance_m


# ---------------------------------------------------------------------------
# D. Determinism / sweep-order independence
# ---------------------------------------------------------------------------


def test_constraint_classification_independent_of_evaluation_order(chief_n):
    r0 = np.array([0.0, -1000.0, 30.0])
    v0 = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    times = list(np.arange(300.0, 5100.0 + 1.0, 100.0))

    def evaluate_all(order):
        out = {}
        for t in order:
            diag = phi_rv_conditioning(t, chief_n)
            if diag.unsafe:
                continue
            result = solve_two_impulse(r0, v0, rf, vf, t, chief_n)
            cc = classify_trajectory_constraints(r0, result.v0_plus, t, chief_n, sample_dt=5.0)
            out[t] = (cc.safe, cc.min_distance_m, cc.violation_reason)
        return out

    ascending = evaluate_all(times)
    descending = evaluate_all(list(reversed(times)))
    shuffled_times = list(times)
    random.Random(7).shuffle(shuffled_times)
    shuffled = evaluate_all(shuffled_times)

    assert set(ascending) == set(descending) == set(shuffled)
    for t in ascending:
        assert ascending[t] == descending[t] == shuffled[t]


# ---------------------------------------------------------------------------
# E. M3 regression preservation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "t,expected_mm_s",
    [(1800.0, 1064.5203), (M3_SHORT_MIN_T, 678.6092), (M3_LONG_MIN_T, 163.4574)],
)
def test_m3_dv_values_unchanged_through_m4_pipeline(chief_n, t, expected_mm_s):
    result = solve_two_impulse(
        np.array([0.0, -1000.0, 30.0]), np.array([0.0, 0.0, 0.0]),
        np.array([0.0, -50.0, 0.0]), np.array([0.0, 0.0, 0.0]),
        t, chief_n,
    )
    # Evaluating constraints must not perturb the solver's own result.
    classify_trajectory_constraints(
        np.array([0.0, -1000.0, 30.0]), result.v0_plus, t, chief_n, sample_dt=5.0
    )
    assert result.dv_total * 1000.0 == pytest.approx(expected_mm_s, abs=0.05)


# ---------------------------------------------------------------------------
# F. M3 unsafe region remains excluded
# ---------------------------------------------------------------------------


def test_m3_unsafe_times_never_rescued_by_geometry(chief_n, chief_period):
    t_sing = chief_period / 2.0
    for dt in (0.0, 0.001, 0.005):
        t = t_sing + dt
        diag = phi_rv_conditioning(t, chief_n)
        assert diag.unsafe
        with pytest.raises(Exception):
            # solve_two_impulse itself refuses; M4 must never attempt to
            # "rescue" this via geometry -- there is no dv/trajectory to
            # classify at all for an M3-unsafe time.
            solve_two_impulse(
                np.array([0.0, -1000.0, 30.0]), np.array([0.0, 0.0, 0.0]),
                np.array([0.0, -50.0, 0.0]), np.array([0.0, 0.0, 0.0]),
                t, chief_n,
            )


# ---------------------------------------------------------------------------
# G. Final closure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t", [305.0, 340.0, 380.0])
def test_terminal_closure_for_geometrically_feasible_points(chief_n, t):
    r0 = np.array([0.0, -1000.0, 30.0])
    v0 = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    result = solve_two_impulse(r0, v0, rf, vf, t, chief_n)
    cc = classify_trajectory_constraints(r0, result.v0_plus, t, chief_n, sample_dt=1.0)
    assert cc.safe  # sanity: these are the known-feasible fast transfers
    assert result.terminal_position_residual_norm < 1e-8


# ---------------------------------------------------------------------------
# H. Constraint-margin convergence for the selected M4 transfer
# ---------------------------------------------------------------------------


def test_constraint_margin_convergence_at_selected_transfer(chief_n):
    r0 = np.array([0.0, -1000.0, 30.0])
    v0 = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    t_selected = 380.0  # a comfortably feasible fast transfer

    result = solve_two_impulse(r0, v0, rf, vf, t_selected, chief_n)

    cc_2s = classify_trajectory_constraints(r0, result.v0_plus, t_selected, chief_n, sample_dt=2.0)
    cc_1s = classify_trajectory_constraints(r0, result.v0_plus, t_selected, chief_n, sample_dt=1.0)
    cc_refined = classify_trajectory_constraints(
        r0, result.v0_plus, t_selected, chief_n, sample_dt=1.0, refine=True
    )

    # All three should agree the point is safe, and the closest-approach
    # distance should converge as resolution increases / refinement is applied.
    assert cc_2s.safe and cc_1s.safe and cc_refined.safe
    assert abs(cc_2s.min_distance_m - cc_1s.min_distance_m) < 1.0
    assert cc_refined.min_distance_m <= cc_1s.min_distance_m + 1e-9


# ---------------------------------------------------------------------------
# I. Independent trajectory check (STM vs. ODE) for the selected transfer
# ---------------------------------------------------------------------------


def test_independent_ode_check_for_selected_transfer(chief_n):
    r0 = np.array([0.0, -1000.0, 30.0])
    v0 = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    t_selected = 380.0

    result = solve_two_impulse(r0, v0, rf, vf, t_selected, chief_n)
    post_burn1 = np.concatenate([r0, result.v0_plus])

    cc_stm = classify_trajectory_constraints(r0, result.v0_plus, t_selected, chief_n, sample_dt=1.0)

    check_times = np.linspace(0.0, t_selected, 6)
    for tt in check_times:
        sol = solve_ivp(
            cw_ode_rhs, (0.0, tt), post_burn1, args=(chief_n,),
            method="DOP853", rtol=1e-12, atol=1e-12,
        )
        assert sol.success
        state_ode = sol.y[:, -1]
        state_stm = propagate_cw(post_burn1, tt, chief_n)
        np.testing.assert_allclose(state_stm, state_ode, atol=1e-6, rtol=1e-8)

    # Constraint classification is unchanged by cross-checking against the ODE.
    assert cc_stm.safe


# ---------------------------------------------------------------------------
# J. Corridor-entry sanity
# ---------------------------------------------------------------------------


def test_corridor_entry_occurs_inside_authorized_corridor_for_feasible_transfer(chief_n):
    r0 = np.array([0.0, -1000.0, 30.0])
    v0 = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf = np.array([0.0, 0.0, 0.0])
    t_selected = 380.0

    result = solve_two_impulse(r0, v0, rf, vf, t_selected, chief_n)
    cc = classify_trajectory_constraints(r0, result.v0_plus, t_selected, chief_n, sample_dt=0.5)

    assert cc.safe
    assert cc.corridor_entry_time_s is not None

    state0 = np.concatenate([r0, result.v0_plus])
    r_at_entry = propagate_cw(state0, cc.corridor_entry_time_s, chief_n)[:3]
    assert inside_approach_corridor(r_at_entry)

    # And confirm the trajectory is never inside the 100 m sphere before
    # the corridor-entry time (KOZ only ever authorized inside the box).
    for tt in np.linspace(0.0, cc.corridor_entry_time_s - 0.01, 20):
        r_t = propagate_cw(state0, tt, chief_n)[:3]
        assert np.linalg.norm(r_t) >= R_KOZ - 1e-6
