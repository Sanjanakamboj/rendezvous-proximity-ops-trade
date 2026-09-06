"""M2 verification: chief-orbit utilities match the DESIGN.md M1 scenario."""

from __future__ import annotations

import pytest

from rendezvous_cw.orbit import (
    CHIEF_ALTITUDE,
    MU_EARTH,
    R_EARTH,
    chief_mean_motion_and_period,
    mean_motion,
    orbital_period,
    orbital_radius,
)

# DESIGN.md section 2.1 / 7.1 values.
DESIGN_MD_A_KM = 6778.137
DESIGN_MD_N_RAD_S = 1.131366654e-3
DESIGN_MD_PERIOD_S = 5553.624271


def test_constants_match_design_md():
    assert R_EARTH == pytest.approx(6378.137e3)
    assert MU_EARTH == pytest.approx(398600.4418e9)
    assert CHIEF_ALTITUDE == pytest.approx(400.0e3)


def test_orbital_radius_matches_design_md():
    a = orbital_radius(CHIEF_ALTITUDE)
    assert a == pytest.approx(DESIGN_MD_A_KM * 1e3, rel=1e-9)


def test_mean_motion_matches_design_md():
    a = orbital_radius(CHIEF_ALTITUDE)
    n = mean_motion(a)
    assert n == pytest.approx(DESIGN_MD_N_RAD_S, rel=1e-6)


def test_period_matches_design_md():
    a = orbital_radius(CHIEF_ALTITUDE)
    n = mean_motion(a)
    p = orbital_period(n)
    assert p == pytest.approx(DESIGN_MD_PERIOD_S, rel=1e-6)


def test_chief_mean_motion_and_period_convenience_wrapper():
    n, p = chief_mean_motion_and_period()
    assert n == pytest.approx(DESIGN_MD_N_RAD_S, rel=1e-6)
    assert p == pytest.approx(DESIGN_MD_PERIOD_S, rel=1e-6)
