"""M2 diagnostic figure: the T = 1800 s representative two-impulse transfer.

This is a diagnostic plot only (per project ground rules, figures are
diagnostics until the underlying equations/tests are independently
verified -- which they are, in tests/). It is NOT the M3 portfolio-final
approach-trajectory figure or the Delta-v/time trade plot.

Run:
    python scripts/m2_plot_representative_transfer.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.orbit import chief_mean_motion_and_period
from rendezvous_cw.rendezvous import solve_two_impulse


def main() -> None:
    n, _period = chief_mean_motion_and_period()

    r0 = np.array([0.0, -1000.0, 30.0])
    v0_minus = np.array([0.0, 0.0, 0.0])
    rf = np.array([0.0, -50.0, 0.0])
    vf_plus = np.array([0.0, 0.0, 0.0])
    T = 1800.0

    result = solve_two_impulse(r0, v0_minus, rf, vf_plus, T, n)

    # Sample the post-first-burn trajectory over [0, T] using the STM
    # (the sole production propagation path -- no ODE solver here).
    t_samples = np.linspace(0.0, T, 400)
    state0 = np.concatenate([r0, result.v0_plus])
    states = np.array([propagate_cw(state0, t, n) for t in t_samples])
    x, y, z = states[:, 0], states[:, 1], states[:, 2]

    fig, (ax_xy, ax_z) = plt.subplots(1, 2, figsize=(11, 5))

    ax_xy.plot(y, x, "-", color="tab:blue", lw=1.8, label="relative trajectory")
    ax_xy.plot(r0[1], r0[0], "o", color="tab:green", ms=9, label="r0 (start)")
    ax_xy.plot(rf[1], rf[0], "s", color="tab:red", ms=9, label="rf (target)")
    ax_xy.plot(0, 0, "*", color="k", ms=12, label="chief (origin)")
    ax_xy.set_xlabel("y — along-track [m]")
    ax_xy.set_ylabel("x — radial [m]")
    ax_xy.set_title("LVLH x–y (radial vs. along-track) projection")
    ax_xy.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax_xy.axvline(0, color="0.85", lw=0.8, zorder=0)
    ax_xy.legend(loc="best", fontsize=8)
    ax_xy.set_aspect("equal", adjustable="datalim")
    ax_xy.invert_xaxis()  # y more negative (behind chief) drawn to the right, approaching origin

    ax_z.plot(t_samples, z, "-", color="tab:purple", lw=1.8)
    ax_z.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax_z.set_xlabel("time since burn 1 [s]")
    ax_z.set_ylabel("z — cross-track [m]")
    ax_z.set_title("Cross-track motion (30 m offset nulled by burn 1)")

    fig.suptitle(
        "M2 diagnostic — CW linearized relative trajectory (closed-form STM only)\n"
        f"Transfer time T = {T:.0f} s   |   "
        f"|Δv1| = {result.dv1_mag * 1000:.2f} mm/s, "
        f"|Δv2| = {result.dv2_mag * 1000:.2f} mm/s, "
        f"Δv_total = {result.dv_total * 1000:.2f} mm/s",
        fontsize=10,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))

    out_path = Path(__file__).resolve().parent.parent / "figures" / "m2_representative_transfer.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")


if __name__ == "__main__":
    main()
