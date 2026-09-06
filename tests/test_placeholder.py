"""Milestone 1 placeholder test.

Verifies the package scaffolding imports cleanly. The real verification
suite (state-transition identity, ODE residual checks, STM composition,
etc. -- see DESIGN.md section 9) is added in Milestone 2 alongside the CW
solver implementation it verifies.
"""

import rendezvous_cw


def test_package_imports_and_has_version():
    assert hasattr(rendezvous_cw, "__version__")
    assert isinstance(rendezvous_cw.__version__, str)
