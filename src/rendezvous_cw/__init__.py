"""rendezvous_cw — Clohessy-Wiltshire / Hill-frame relative-motion analysis.

Milestone 6 (final) status: chief-orbit utilities, the closed-form CW
state-transition matrix, Phi_rv conditioning diagnostics, the two-impulse
boundary-value solver, the transfer-time/Delta-v trade-study module
(trade.py), proximity-operations geometric screening (constraints.py),
and an independent nonlinear (inertial Cartesian two-body) validation
model with LVLH<->ECI frame kinematics (nonlinear.py) are implemented and
verified in tests/. See DESIGN.md section 16 for the final authoritative
result and full superseded-result hierarchy across all milestones.
J2/drag, finite burns, actuator limits, sensor field-of-view,
line-of-sight occultation, plume impingement, docking dynamics, and
collision-probability modeling are explicitly out of this project's
scope and are not implemented.
"""

__version__ = "0.6.0"
