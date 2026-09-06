"""Circular chief-orbit utilities.

Pure functions only. All quantities are SI (meters, seconds, radians)
unless a function name/docstring says otherwise. Constants match the M1
scenario documented in DESIGN.md section 2.1 (400 km circular LEO chief).
Do not change these constants without updating DESIGN.md.
"""

from __future__ import annotations

import math

#: Earth gravitational parameter, m^3/s^2 (WGS-84 value; matches DESIGN.md,
#: which quotes it in km^3/s^2 as 398600.4418).
MU_EARTH = 398600.4418e9

#: Earth mean equatorial radius, m (WGS-84 value; matches DESIGN.md's
#: 6378.137 km).
R_EARTH = 6378.137e3

#: M1 baseline chief altitude, m (400 km circular LEO -- DESIGN.md section 2.1).
CHIEF_ALTITUDE = 400.0e3


def orbital_radius(altitude: float, r_body: float = R_EARTH) -> float:
    """Circular orbit radius a = r_body + altitude, in meters."""
    return r_body + altitude


def mean_motion(a: float, mu: float = MU_EARTH) -> float:
    """Mean motion n = sqrt(mu / a^3), in rad/s, for a circular orbit of radius a."""
    return math.sqrt(mu / a**3)


def orbital_period(n: float) -> float:
    """Orbital period P = 2*pi / n, in seconds."""
    return 2.0 * math.pi / n


def chief_mean_motion_and_period(
    altitude: float = CHIEF_ALTITUDE,
    mu: float = MU_EARTH,
    r_body: float = R_EARTH,
) -> tuple[float, float]:
    """Convenience wrapper: (n, period) for a circular chief orbit at `altitude`.

    Defaults reproduce the exact M1 scenario: 400 km circular LEO,
    n = 1.131366654e-3 rad/s, period = 5553.624271 s (DESIGN.md section 2.1
    / section 7.1).
    """
    a = orbital_radius(altitude, r_body)
    n = mean_motion(a, mu)
    p = orbital_period(n)
    return n, p
