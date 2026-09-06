"""rendezvous_cw — Clohessy-Wiltshire / Hill-frame relative-motion analysis.

Milestone 4 status: chief-orbit utilities, the closed-form CW state-
transition matrix, Phi_rv conditioning diagnostics, the two-impulse
boundary-value solver, the transfer-time/Delta-v trade-study module
(trade.py), and proximity-operations geometric screening -- a keep-out
sphere, final V-bar approach corridor, and no-chief-crossing rule
(constraints.py) -- are implemented and verified in tests/. This is
geometric screening inside the linearized CW model only: nonlinear
two-body validation, J2/drag, finite burns, actuator limits, sensor
field-of-view, line-of-sight occultation, plume impingement, docking
dynamics, and collision-probability modeling are Milestone 5+ work and
are not implemented yet.
"""

__version__ = "0.4.0"
