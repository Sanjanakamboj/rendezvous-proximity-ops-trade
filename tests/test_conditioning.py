"""M2 verification: Phi_rv conditioning / singularity detection.

Covers DESIGN.md section 9 check I: probe exactly at n*t = pi, very near
it on both sides, and comfortably away from it; verify the solver refuses
unsafe points and works normally away from the singularity, and report
condition-number growth approaching it.
"""

from __future__ import annotations

import numpy as np
import pytest

from rendezvous_cw.conditioning import COND_THRESHOLD, phi_rv_conditioning
from rendezvous_cw.rendezvous import SingularTransferError, solve_two_impulse


def test_exactly_at_half_period_is_unsafe(chief_n, chief_period):
    t_sing = chief_period / 2.0
    diag = phi_rv_conditioning(t_sing, chief_n)
    assert diag.unsafe
    assert not np.isfinite(diag.cond) or diag.cond > COND_THRESHOLD


@pytest.mark.parametrize("dt", [1e-5, 1e-4, 1e-3])
def test_very_near_half_period_is_unsafe(chief_n, chief_period, dt):
    """Within the documented guard band (~10 ms, conditioning.py docstring)."""
    t_sing = chief_period / 2.0
    for sign in (-1.0, 1.0):
        diag = phi_rv_conditioning(t_sing + sign * dt, chief_n)
        assert diag.unsafe, f"expected unsafe at dt={sign * dt}, cond={diag.cond}"


@pytest.mark.parametrize("dt", [0.1, 1.0, 10.0, 100.0, 500.0])
def test_comfortably_away_from_half_period_is_safe(chief_n, chief_period, dt):
    t_sing = chief_period / 2.0
    for sign in (-1.0, 1.0):
        diag = phi_rv_conditioning(t_sing + sign * dt, chief_n)
        assert not diag.unsafe, f"expected safe at dt={sign * dt}, cond={diag.cond}"


def test_condition_number_grows_monotonically_approaching_singularity(chief_n, chief_period):
    """Report/verify condition number growth as t -> half-period singularity."""
    t_sing = chief_period / 2.0
    offsets = [500.0, 100.0, 10.0, 1.0, 0.1, 0.01, 0.001]
    conds = [phi_rv_conditioning(t_sing - dt, chief_n).cond for dt in offsets]
    # As dt shrinks (approaching the singularity), condition number must
    # strictly increase.
    for smaller, larger in zip(conds, conds[1:]):
        assert larger > smaller, conds


def test_solver_refuses_singular_transfer_time(chief_n, chief_period, m1_scenario):
    t_sing = chief_period / 2.0
    with pytest.raises(SingularTransferError):
        solve_two_impulse(
            m1_scenario["r0"],
            m1_scenario["v0_minus"],
            m1_scenario["rf"],
            m1_scenario["vf_plus"],
            t_sing,
            chief_n,
        )


def test_solver_works_normally_away_from_singularity(chief_n, m1_scenario):
    # T = 1800 s is the M1/M2 representative transfer, far from any
    # nt = k*pi singularity.
    result = solve_two_impulse(
        m1_scenario["r0"],
        m1_scenario["v0_minus"],
        m1_scenario["rf"],
        m1_scenario["vf_plus"],
        1800.0,
        chief_n,
    )
    assert not result.conditioning.unsafe
    assert result.dv_total > 0.0
    assert np.isfinite(result.dv_total)
