"""M2 verification: independent numerical ODE propagation cross-check.

Covers DESIGN.md section 9 check J. scipy.integrate.solve_ivp is used
ONLY in this verification path -- never in the production propagation
code (cw.py's propagate_cw uses the closed-form STM exclusively).
"""

from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from rendezvous_cw.cw import propagate_cw

from conftest import cw_ode_rhs

_TEST_STATES = [
    np.array([100.0, -500.0, 20.0, 0.05, -0.02, 0.01]),
    np.array([0.0, -1000.0, 30.0, 0.0, 0.0, 0.0]),
    np.array([-250.0, 300.0, -40.0, -0.03, 0.04, -0.02]),
]
_TEST_TIMES = [50.0, 300.0, 1800.0, 4000.0, 5000.0]


@pytest.mark.parametrize("state0", _TEST_STATES)
@pytest.mark.parametrize("t", _TEST_TIMES)
def test_stm_matches_independent_ode_integration(chief_n, state0, t):
    n = chief_n

    sol = solve_ivp(
        cw_ode_rhs,
        t_span=(0.0, t),
        y0=state0,
        args=(n,),
        method="DOP853",
        rtol=1e-12,
        atol=1e-12,
        dense_output=False,
    )
    assert sol.success
    state_ode = sol.y[:, -1]

    state_stm = propagate_cw(state0, t, n)

    np.testing.assert_allclose(state_stm, state_ode, atol=1e-6, rtol=1e-8)
