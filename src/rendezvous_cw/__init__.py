"""rendezvous_cw — Clohessy-Wiltshire / Hill-frame relative-motion analysis.

Milestone 5 status: chief-orbit utilities, the closed-form CW state-
transition matrix, Phi_rv conditioning diagnostics, the two-impulse
boundary-value solver, the transfer-time/Delta-v trade-study module
(trade.py), proximity-operations geometric screening (constraints.py),
and an independent nonlinear (inertial Cartesian two-body) validation
model with LVLH<->ECI frame kinematics (nonlinear.py) are implemented and
verified in tests/. J2/drag, finite burns, actuator limits, sensor
field-of-view, line-of-sight occultation, plume impingement, docking
dynamics, and collision-probability modeling are Milestone 6+ work and
are not implemented yet.
"""

__version__ = "0.5.0"
