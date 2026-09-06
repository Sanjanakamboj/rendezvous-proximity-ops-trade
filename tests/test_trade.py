"""M3 verification suite for the transfer-time / Delta-v trade study (trade.py).

Covers the M3 task checks A-I. All Delta-v numbers are obtained by
calling trade.py, which itself calls only the verified M2 solver
(rendezvous.solve_two_impulse) -- no CW equations are duplicated here.
"""

from __future__ import annotations

import random

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from rendezvous_cw.conditioning import COND_THRESHOLD
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.rendezvous import solve_two_impulse
from rendezvous_cw.trade import evaluate_transfer_time, find_best_safe_transfer, sweep_transfer_times

from conftest import cw_ode_rhs


# ---------------------------------------------------------------------------
# A. M2 regression preservation
# ---------------------------------------------------------------------------


def test_m1_scenario_via_trade_module_matches_m2_authoritative_values(m1_scenario, chief_period):
    rec = evaluate_transfer_time(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        m1_scenario["T"],
        m1_scenario["n"],
        chief_period,
    )
    assert rec.safe
    assert rec.dv1_mag * 1000.0 == pytest.approx(531.7191, abs=1e-3)
    assert rec.dv2_mag * 1000.0 == pytest.approx(532.8013, abs=1e-3)
    assert rec.dv_total * 1000.0 == pytest.approx(1064.5203, abs=1e-3)


# ---------------------------------------------------------------------------
# B. Trade module vs. direct solver
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("t", [350.0, 750.0, 1234.5, 2000.0, 3300.0, 4500.0, 5050.0])
def test_trade_module_matches_direct_solver(chief_n, chief_period, m1_scenario, t):
    rec = evaluate_transfer_time(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        t,
        chief_n,
        chief_period,
    )
    direct = solve_two_impulse(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        t,
        chief_n,
    )
    assert rec.safe
    assert rec.dv1_mag == pytest.approx(direct.dv1_mag, rel=1e-12)
    assert rec.dv2_mag == pytest.approx(direct.dv2_mag, rel=1e-12)
    assert rec.dv_total == pytest.approx(direct.dv_total, rel=1e-12)
    assert rec.closure_residual_m == pytest.approx(
        direct.terminal_position_residual_norm, abs=1e-15
    )


# ---------------------------------------------------------------------------
# C. Unsafe-zone preservation
# ---------------------------------------------------------------------------


def test_unsafe_zone_marked_and_excluded_from_minimum(chief_n, chief_period, m1_scenario):
    t_sing = chief_period / 2.0
    # Dense probe spanning the known ~19 ms guard band.
    probe_times = t_sing + np.arange(-0.05, 0.05 + 1e-9, 0.001)
    records = sweep_transfer_times(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        probe_times,
        chief_n,
        chief_period,
    )

    unsafe = [r for r in records if not r.safe]
    safe = [r for r in records if r.safe]
    assert len(unsafe) > 0, "expected at least one unsafe record in the probe band"
    for r in unsafe:
        assert r.cond_phi_rv > COND_THRESHOLD or not np.isfinite(r.cond_phi_rv)
        assert np.isnan(r.dv_total)
        assert np.isnan(r.dv1_mag)
        assert np.isnan(r.dv2_mag)

    # Unsafe records must never be selected as the minimum.
    best = find_best_safe_transfer(
        records,
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        chief_n,
        chief_period,
        refine=False,
    )
    assert best.safe
    assert best in safe


# ---------------------------------------------------------------------------
# D. Terminal closure across sweep
# ---------------------------------------------------------------------------


def test_terminal_closure_across_full_sweep(chief_n, chief_period, m1_scenario):
    times = np.arange(300.0, 5100.0 + 1.0, 25.0)
    records = sweep_transfer_times(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        times,
        chief_n,
        chief_period,
    )
    safe = [r for r in records if r.safe]
    assert len(safe) > 0
    worst = max(r.closure_residual_m for r in safe)
    assert worst < 1e-6, f"worst closure residual too large: {worst}"


# ---------------------------------------------------------------------------
# E. Burn-magnitude consistency
# ---------------------------------------------------------------------------


def test_burn_magnitude_consistency_across_sweep(chief_n, chief_period, m1_scenario):
    times = np.arange(300.0, 5100.0 + 1.0, 25.0)
    records = sweep_transfer_times(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        times,
        chief_n,
        chief_period,
    )
    for r in records:
        if not r.safe:
            continue
        assert r.dv1_mag >= 0.0
        assert r.dv2_mag >= 0.0
        assert np.isfinite(r.dv_total)
        assert r.dv_total == pytest.approx(r.dv1_mag + r.dv2_mag, rel=1e-9)


# ---------------------------------------------------------------------------
# F. Grid-resolution convergence of the minimum (branch 2: T > P/2)
# ---------------------------------------------------------------------------


def test_grid_resolution_convergence_of_minimum(chief_n, chief_period, m1_scenario):
    t_sing = chief_period / 2.0
    results = {}
    for step in (20.0, 10.0, 5.0):
        times = np.arange(300.0, 5100.0 + 1.0, step)
        records = sweep_transfer_times(
            m1_scenario["r0"],
            m1_scenario["v0_minus"],
            m1_scenario["rf"],
            m1_scenario["vf_plus"],
            times,
            chief_n,
            chief_period,
        )
        branch2 = [r for r in records if r.safe and r.t_s > t_sing]
        best = min(branch2, key=lambda r: r.dv_total)
        results[step] = best

    # Coarser grids must not do better than finer grids by more than a
    # small amount, and the minimum T location should stabilize.
    dv_20, dv_10, dv_5 = (results[s].dv_total for s in (20.0, 10.0, 5.0))
    assert dv_10 <= dv_20 + 1e-6
    assert dv_5 <= dv_10 + 1e-6
    assert abs(results[10.0].t_s - results[5.0].t_s) <= 10.0

    # Now refine from the finest coarse grid and confirm the refined
    # answer is at least as good as every coarse grid.
    times_fine = np.arange(300.0, 5100.0 + 1.0, 5.0)
    records_fine = sweep_transfer_times(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        times_fine,
        chief_n,
        chief_period,
    )
    branch2_fine = [r for r in records_fine if r.safe and r.t_s > t_sing]
    refined = find_best_safe_transfer(
        branch2_fine,
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        chief_n,
        chief_period,
        refine=True,
    )
    assert refined.dv_total <= dv_5 + 1e-9
    assert refined.dv_total <= dv_10 + 1e-9
    assert refined.dv_total <= dv_20 + 1e-9


# ---------------------------------------------------------------------------
# G. Independence from sweep ordering
# ---------------------------------------------------------------------------


def test_sweep_independent_of_time_ordering(chief_n, chief_period, m1_scenario):
    times = list(np.arange(300.0, 5100.0 + 1.0, 50.0))

    ascending = sweep_transfer_times(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], times, chief_n, chief_period,
    )
    descending = sweep_transfer_times(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], list(reversed(times)), chief_n, chief_period,
    )
    shuffled_times = list(times)
    random.Random(42).shuffle(shuffled_times)
    shuffled = sweep_transfer_times(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], shuffled_times, chief_n, chief_period,
    )

    by_t_ascending = {r.t_s: r for r in ascending}
    by_t_descending = {r.t_s: r for r in descending}
    by_t_shuffled = {r.t_s: r for r in shuffled}

    assert set(by_t_ascending) == set(by_t_descending) == set(by_t_shuffled)
    for t in by_t_ascending:
        a, d, s = by_t_ascending[t], by_t_descending[t], by_t_shuffled[t]
        assert a.safe == d.safe == s.safe
        assert a.dv_total == d.dv_total == s.dv_total or (
            np.isnan(a.dv_total) and np.isnan(d.dv_total) and np.isnan(s.dv_total)
        )
        assert a.cond_phi_rv == d.cond_phi_rv == s.cond_phi_rv


# ---------------------------------------------------------------------------
# H. Singularity-side behavior
# ---------------------------------------------------------------------------


def test_singularity_side_behavior(chief_n, chief_period, m1_scenario):
    t_sing = chief_period / 2.0
    offsets = [10.0, 1.0, 0.1, 0.02]  # approaching from both sides, staying safe

    for sign in (-1.0, 1.0):
        conds = []
        for dt in offsets:
            rec = evaluate_transfer_time(
                m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
                m1_scenario["vf_plus"], t_sing + sign * dt, chief_n, chief_period,
            )
            assert rec.safe
            conds.append(rec.cond_phi_rv)
        # Condition number must increase monotonically as offset shrinks
        # (offsets are given largest-to-smallest above).
        for larger, smaller in zip(conds, conds[1:]):
            assert smaller > larger, conds

    # The unsafe region genuinely interrupts a dense probe crossing it.
    probe_times = t_sing + np.arange(-0.03, 0.03 + 1e-9, 0.001)
    probe = sweep_transfer_times(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], probe_times, chief_n, chief_period,
    )
    assert any(not r.safe for r in probe)
    assert any(r.safe for r in probe if r.t_s < t_sing - 0.02)
    assert any(r.safe for r in probe if r.t_s > t_sing + 0.02)


# ---------------------------------------------------------------------------
# I. Representative trajectory closure at the refined minimum
# ---------------------------------------------------------------------------


def test_refined_minimum_trajectory_closure(chief_n, chief_period, m1_scenario):
    t_sing = chief_period / 2.0
    times = np.arange(300.0, 5100.0 + 1.0, 10.0)
    records = sweep_transfer_times(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], times, chief_n, chief_period,
    )
    branch2 = [r for r in records if r.safe and r.t_s > t_sing]
    refined = find_best_safe_transfer(
        branch2, m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], chief_n, chief_period, refine=True,
    )

    direct = solve_two_impulse(
        m1_scenario["r0"], m1_scenario["v0_minus"], m1_scenario["rf"],
        m1_scenario["vf_plus"], refined.t_s, chief_n,
    )

    # STM-based propagation must land exactly on rf.
    assert direct.terminal_position_residual_norm < 1e-8

    # Independent solve_ivp cross-check of the post-burn-1 state.
    post_burn1 = np.concatenate([m1_scenario["r0"], direct.v0_plus])
    sol = solve_ivp(
        cw_ode_rhs, (0.0, refined.t_s), post_burn1, args=(chief_n,),
        method="DOP853", rtol=1e-12, atol=1e-12,
    )
    assert sol.success
    state_ode = sol.y[:, -1]
    state_stm = propagate_cw(post_burn1, refined.t_s, chief_n)

    np.testing.assert_allclose(state_stm, state_ode, atol=1e-6, rtol=1e-8)
    np.testing.assert_allclose(state_stm[:3], m1_scenario["rf"], atol=1e-6)
    np.testing.assert_allclose(state_stm[3:], direct.vT_minus, atol=1e-8)
