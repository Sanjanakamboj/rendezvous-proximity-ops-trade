"""M4: proximity-geometry-constrained Delta-v/transfer-time trade.

Re-evaluates the M3 numerically-safe transfer-time sweep against the M4
geometric constraints (constraints.py): a 100 m keep-out sphere, a final
V-bar approach corridor, and a no-chief-crossing rule. Uses ONLY the
verified M2 solver and M3 conditioning policy -- no dynamics changed.

IMPORTANT SCOPE NOTE (see DESIGN.md M4 section and README.md): a
trajectory that passes these checks is described as "geometrically
feasible under this CW model" -- NOT "collision-free", "flight safe", or
"operationally safe". This is linear-CW geometric screening only.

Produces:
  - results/m4_constrained_trade.csv        curated representative-time table
  - results/m4_full_constrained_sweep.csv   full sweep (diagnostic/reproducibility)
  - figures/m4_constrained_dv_vs_transfer_time.png
  - figures/m4_approach_geometry.png
  - figures/m4_constraint_margin_vs_transfer_time.png

Run:
    python scripts/m4_constrained_trade.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rendezvous_cw.conditioning import phi_rv_conditioning
from rendezvous_cw.constraints import (
    CORRIDOR_X_BOUND,
    CORRIDOR_Y_MAX,
    CORRIDOR_Y_MIN,
    R_KOZ,
    classify_trajectory_constraints,
)
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.orbit import chief_mean_motion_and_period
from rendezvous_cw.rendezvous import solve_two_impulse

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"

# --- Exact M1 scenario (frozen; do not alter) -------------------------------
R0 = np.array([0.0, -1000.0, 30.0])
V0_MINUS = np.array([0.0, 0.0, 0.0])
RF = np.array([0.0, -50.0, 0.0])
VF_PLUS = np.array([0.0, 0.0, 0.0])

T_MIN_DOMAIN = 300.0
T_MAX_DOMAIN = 5100.0
SWEEP_STEP = 10.0  # seconds (M4 sweep resolution)
SAMPLE_DT = 2.0  # seconds (M4 constraint-checking dense-sample spacing)

# M3 reference results (DESIGN.md section 13 / M3 checkpoint).
M2_REFERENCE_T = 1800.0
M3_SHORT_BRANCH_MIN_T = 2539.592622
M3_LONG_BRANCH_MIN_T = 4868.130154465808


def evaluate(t: float, n: float, period: float):
    """Evaluate M3 numeric safety + M4 geometry for one transfer time.

    Returns (m3_safe, m4_geometry_safe, two_impulse_result_or_None,
    constraint_result_or_None, conditioning_diag).
    """
    diag = phi_rv_conditioning(t, n)
    if diag.unsafe:
        return False, False, None, None, diag
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, t, n)
    cc = classify_trajectory_constraints(R0, result.v0_plus, t, n, sample_dt=SAMPLE_DT)
    return True, cc.safe, result, cc, diag


def find_feasibility_boundary(n: float, lo: float, hi: float, n_iter: int = 50) -> float:
    """Bisect for the T where geometric feasibility flips between lo (safe) and hi (unsafe)."""
    for _ in range(n_iter):
        mid = (lo + hi) / 2.0
        result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, mid, n)
        cc = classify_trajectory_constraints(R0, result.v0_plus, mid, n, sample_dt=0.5)
        if cc.safe:
            lo = mid
        else:
            hi = mid
    return lo


def main() -> None:
    n, period = chief_mean_motion_and_period()
    t_sing = period / 2.0

    times = np.arange(T_MIN_DOMAIN, T_MAX_DOMAIN + 1.0, SWEEP_STEP)

    rows = []  # (t, m3_safe, m4_safe, result, cc, diag)
    for t in times:
        m3_safe, m4_safe, result, cc, diag = evaluate(t, n, period)
        rows.append((t, m3_safe, m4_safe, result, cc, diag))

    m3_safe_rows = [r for r in rows if r[1]]
    m4_safe_rows = [r for r in rows if r[1] and r[2]]

    print(f"n={n!r} rad/s, period={period!r} s")
    print(f"Sweep: {len(rows)} points, step={SWEEP_STEP} s, domain=[{T_MIN_DOMAIN},{T_MAX_DOMAIN}] s")
    print(f"M3-numerically-safe: {len(m3_safe_rows)}")
    print(f"M4-geometrically-safe (of the M3-safe set): {len(m4_safe_rows)} "
          f"({100.0 * len(m4_safe_rows) / len(m3_safe_rows):.2f}%)")

    if m4_safe_rows:
        t_feasible = [r[0] for r in m4_safe_rows]
        print(f"Feasible T window (coarse grid): [{min(t_feasible):.1f}, {max(t_feasible):.1f}] s")

    # Dominant rejection reasons among M3-safe-but-M4-rejected points.
    rejected = [r for r in rows if r[1] and not r[2]]
    reason_counts: dict[str, int] = {}
    for _, _, _, _, cc, _ in rejected:
        for reason in cc.violation_reason.split(","):
            reason_counts[reason] = reason_counts.get(reason, 0) + 1
    print(f"Rejected (M3-safe, M4-unsafe): {len(rejected)}")
    print(f"Rejection reason counts: {reason_counts}")

    # --- Find the coarse-grid best (max T, i.e. lowest dv within the fast
    # feasible window) and refine the feasibility boundary precisely. -------
    if m4_safe_rows:
        coarse_best = max(m4_safe_rows, key=lambda r: r[0])  # largest safe T -> lowest dv here
        coarse_best_t = coarse_best[0]
        idx = list(times).index(coarse_best_t)
        # Bracket: coarse_best_t (safe) to the next grid point (may be unsafe).
        hi_bracket = times[idx + 1] if idx + 1 < len(times) else coarse_best_t + SWEEP_STEP
        refined_t = find_feasibility_boundary(n, coarse_best_t, hi_bracket)
    else:
        refined_t = None

    if refined_t is not None:
        refined_result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, refined_t, n)
        refined_cc = classify_trajectory_constraints(R0, refined_result.v0_plus, refined_t, n, sample_dt=0.5)
        print(f"\nRefined M4 minimum-feasible transfer: T={refined_t!r} s "
              f"(T/P={refined_t / period:.6f})")
        print(f"  dv_total={refined_result.dv_total * 1000:.4f} mm/s, "
              f"clearance_outside_corridor={refined_cc.min_clearance_outside_corridor_m!r} m")

    # --- Check the M3 reference/minima explicitly ---------------------------
    print("\nM3 key transfer times re-evaluated under M4 geometry:")
    for label, t in [
        ("M2 reference", M2_REFERENCE_T),
        ("M3 short-branch min", M3_SHORT_BRANCH_MIN_T),
        ("M3 long-branch/global min", M3_LONG_BRANCH_MIN_T),
    ]:
        m3_safe, m4_safe, result, cc, diag = evaluate(t, n, period)
        print(f"  {label} (T={t:.3f} s): m3_safe={m3_safe} m4_geo_safe={m4_safe} "
              f"reason={cc.violation_reason if cc else 'n/a (m3-unsafe)'}")

    # ==========================================================================
    # Curated decision table
    # ==========================================================================
    representative_times = {
        "M2 reference (1800 s)": M2_REFERENCE_T,
        "M3 short-branch minimum": M3_SHORT_BRANCH_MIN_T,
        "M3 long-branch/global minimum": M3_LONG_BRANCH_MIN_T,
        "fast safe transfer (300 s)": 300.0,
        "slow safe transfer (5100 s)": 5100.0,
        "near feasibility boundary (safe side, 385 s)": 385.0,
        "near feasibility boundary (unsafe side, 390 s)": 390.0,
        "refined M4 minimum-feasible transfer": refined_t if refined_t is not None else float("nan"),
    }

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS_DIR / "m4_constrained_trade.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "label", "T_s", "T_min", "T_over_P",
                "m3_numerically_safe", "m4_geometry_safe", "cond_phi_rv",
                "dv_total_m_s", "min_distance_m", "time_of_min_distance_s",
                "min_y_m", "max_y_m", "koz_violation", "corridor_violation",
                "chief_crossing_violation", "violation_reason", "closure_error_m",
            ]
        )
        for label, t in representative_times.items():
            m3_safe, m4_safe, result, cc, diag = evaluate(t, n, period)
            writer.writerow(
                [
                    label,
                    f"{t:.6f}",
                    f"{t / 60.0:.6f}",
                    f"{t / period:.8f}",
                    m3_safe,
                    m4_safe if m3_safe else "n/a",
                    f"{diag.cond:.6e}",
                    f"{result.dv_total:.6f}" if result else "NaN",
                    f"{cc.min_distance_m:.6f}" if cc else "NaN",
                    f"{cc.time_of_min_distance_s:.6f}" if cc else "NaN",
                    f"{cc.min_y_m:.6f}" if cc else "NaN",
                    f"{cc.max_y_m:.6f}" if cc else "NaN",
                    cc.koz_violation if cc else "n/a",
                    cc.corridor_violation if cc else "n/a",
                    cc.chief_crossing_violation if cc else "n/a",
                    cc.violation_reason if cc else "m3_numerically_unsafe",
                    f"{result.terminal_position_residual_norm:.3e}" if result else "NaN",
                ]
            )
    print(f"\nWrote {csv_path}")

    full_csv_path = RESULTS_DIR / "m4_full_constrained_sweep.csv"
    with open(full_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "T_s", "T_over_P", "m3_numerically_safe", "m4_geometry_safe",
                "cond_phi_rv", "dv_total_m_s", "min_distance_m",
                "min_clearance_outside_corridor_m", "min_y_m", "max_y_m",
                "koz_violation", "corridor_violation", "chief_crossing_violation",
                "violation_reason", "closure_error_m",
            ]
        )
        for t, m3_safe, m4_safe, result, cc, diag in rows:
            writer.writerow(
                [
                    f"{t:.6f}",
                    f"{t / period:.8f}",
                    m3_safe,
                    m4_safe if m3_safe else "",
                    f"{diag.cond:.6e}",
                    f"{result.dv_total:.6f}" if result else "",
                    f"{cc.min_distance_m:.6f}" if cc else "",
                    f"{cc.min_clearance_outside_corridor_m:.6f}" if cc else "",
                    f"{cc.min_y_m:.6f}" if cc else "",
                    f"{cc.max_y_m:.6f}" if cc else "",
                    cc.koz_violation if cc else "",
                    cc.corridor_violation if cc else "",
                    cc.chief_crossing_violation if cc else "",
                    cc.violation_reason if cc else "",
                    f"{result.terminal_position_residual_norm:.3e}" if result else "",
                ]
            )
    print(f"Wrote {full_csv_path}")

    # ==========================================================================
    # Figure 1 -- headline constrained Delta-v/time trade
    # ==========================================================================
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    fig1, ax1 = plt.subplots(figsize=(9.5, 6))

    m3_only = sorted([r for r in rows if r[1] and not r[2]], key=lambda r: r[0])
    m4_ok = sorted([r for r in rows if r[1] and r[2]], key=lambda r: r[0])

    if m3_only:
        t_arr = np.array([r[0] for r in m3_only]) / 60.0
        dv_arr = np.array([r[3].dv_total for r in m3_only]) * 1000.0
        ax1.scatter(t_arr, dv_arr, s=8, color="tab:red", alpha=0.6,
                     label="M3-safe, M4-geometry-rejected")
    if m4_ok:
        t_arr = np.array([r[0] for r in m4_ok]) / 60.0
        dv_arr = np.array([r[3].dv_total for r in m4_ok]) * 1000.0
        ax1.scatter(t_arr, dv_arr, s=14, color="tab:green", alpha=0.9,
                     label="M4 geometrically feasible")

    ax1.axvline(t_sing / 60.0, color="0.3", ls="--", lw=1.2,
                 label="P/2 — conditioning-excluded region")

    m2_ref_result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M2_REFERENCE_T, n)
    m3_global_result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, M3_LONG_BRANCH_MIN_T, n)
    ax1.plot(M2_REFERENCE_T / 60.0, m2_ref_result.dv_total * 1000.0, "o", color="black", ms=9,
              label=f"M2 reference (T=1800 s, {m2_ref_result.dv_total*1000:.0f} mm/s) — geometry-rejected")
    ax1.plot(M3_LONG_BRANCH_MIN_T / 60.0, m3_global_result.dv_total * 1000.0, "X", color="tab:purple", ms=11,
              label=f"M3 unconstrained min (T={M3_LONG_BRANCH_MIN_T:.0f} s, "
                    f"{m3_global_result.dv_total*1000:.0f} mm/s) — geometry-rejected")
    if refined_t is not None:
        ax1.plot(refined_t / 60.0, refined_result.dv_total * 1000.0, "*", color="tab:green", ms=18,
                  mec="black", mew=0.7,
                  label=f"refined M4 feasible minimum (T={refined_t:.1f} s, "
                        f"{refined_result.dv_total*1000:.0f} mm/s)")

    ax1.set_xlabel("transfer time T [minutes]")
    ax1.set_ylabel("Δv_total [mm/s]")
    ax1.set_yscale("log")
    ax1.set_title(
        "M4 — Δv vs. transfer time under geometric constraints\n"
        "CW linearized dynamics; 100 m keep-out sphere + final V-bar corridor; "
        "conditioning exclusion retained",
        fontsize=10,
    )
    ax1.legend(loc="upper right", fontsize=7.5, framealpha=0.92)
    fig1.tight_layout()
    fig1_path = FIGURES_DIR / "m4_constrained_dv_vs_transfer_time.png"
    fig1.savefig(fig1_path, dpi=150)
    print(f"Wrote {fig1_path}")

    # ==========================================================================
    # Figure 2 -- approach geometry for the selected M4 transfer
    # ==========================================================================
    if refined_t is not None:
        t_samples = np.linspace(0.0, refined_t, 600)
        state0 = np.concatenate([R0, refined_result.v0_plus])
        states = np.array([propagate_cw(state0, tt, n) for tt in t_samples])
        x, y, z = states[:, 0], states[:, 1], states[:, 2]

        fig2, ax2 = plt.subplots(figsize=(8, 7))

        theta = np.linspace(0, 2 * np.pi, 200)
        ax2.plot(R_KOZ * np.cos(theta), R_KOZ * np.sin(theta), "--", color="tab:red", lw=1.3,
                  label=f"{R_KOZ:.0f} m keep-out sphere (x-y projection)")

        # Corridor box: y in [CORRIDOR_Y_MIN, CORRIDOR_Y_MAX], |x|<=CORRIDOR_X_BOUND.
        corridor_y = [CORRIDOR_Y_MIN, CORRIDOR_Y_MAX, CORRIDOR_Y_MAX, CORRIDOR_Y_MIN, CORRIDOR_Y_MIN]
        corridor_x = [-CORRIDOR_X_BOUND, -CORRIDOR_X_BOUND, CORRIDOR_X_BOUND, CORRIDOR_X_BOUND, -CORRIDOR_X_BOUND]
        ax2.plot(corridor_y, corridor_x, "-", color="tab:orange", lw=1.8,
                  label="final approach corridor (x–y bounds)")

        ax2.plot(y, x, "-", color="tab:blue", lw=1.8,
                  label=f"selected M4 trajectory (T={refined_t:.1f} s)")
        ax2.plot(R0[1], R0[0], "o", color="tab:green", ms=10, zorder=5, label="r0 (start)")
        ax2.plot(RF[1], RF[0], "s", color="tab:red", ms=10, zorder=5, label="rf (target)")
        i_min = int(np.argmin(np.linalg.norm(states, axis=1)))
        min_dist_coincides_with_target = abs(t_samples[i_min] - refined_t) < 1.0
        ax2.plot(y[i_min], x[i_min], "^", color="black", ms=13, mfc="none", mew=2.0, zorder=6,
                  label="min clearance point"
                        + (" (coincides with target rf)" if min_dist_coincides_with_target else
                           f" (t={t_samples[i_min]:.1f} s)"))
        ax2.plot(0, 0, "*", color="k", ms=14, zorder=5, label="chief (origin)")

        ax2.set_xlabel("y — along-track [m]")
        ax2.set_ylabel("x — radial [m]")
        ax2.set_title(
            "M4 — Selected feasible trajectory vs. keep-out sphere and approach corridor\n"
            "(LVLH x–y projection; passes modeled geometric constraints)",
            fontsize=10,
        )
        ax2.set_xlim(50, -180)
        ax2.set_ylim(-130, 130)
        ax2.set_aspect("equal", adjustable="box")
        ax2.legend(loc="upper center", fontsize=7.2, framealpha=0.92, ncol=1,
                    bbox_to_anchor=(0.5, -0.09))

        ax2_inset = fig2.add_axes([0.16, 0.14, 0.22, 0.16])
        ax2_inset.plot(t_samples, z, "-", color="tab:purple", lw=1.3)
        ax2_inset.axhline(0, color="0.85", lw=0.6)
        ax2_inset.set_xlabel("t [s]", fontsize=7)
        ax2_inset.set_ylabel("z [m]", fontsize=7)
        ax2_inset.set_title("cross-track vs. time", fontsize=7)
        ax2_inset.tick_params(labelsize=6)

        fig2_path = FIGURES_DIR / "m4_approach_geometry.png"
        fig2.savefig(fig2_path, dpi=150, bbox_inches="tight")
        print(f"Wrote {fig2_path}")

    # ==========================================================================
    # Figure 3 -- constraint margin vs. transfer time
    # ==========================================================================
    fig3, ax3 = plt.subplots(figsize=(9.5, 6))

    m3_rows_sorted = sorted(m3_safe_rows, key=lambda r: r[0])
    t_arr = np.array([r[0] for r in m3_rows_sorted]) / 60.0
    clearance_arr = np.array([r[4].min_clearance_outside_corridor_m for r in m3_rows_sorted])
    chief_cross_arr = np.array([r[4].chief_crossing_violation for r in m3_rows_sorted])

    ax3.plot(t_arr, clearance_arr, "-", color="tab:blue", lw=1.5,
              label="min clearance outside corridor = min(||r|| − 100 m)")
    ax3.axhline(0.0, color="tab:red", ls="--", lw=1.3, label="0 m requirement (KOZ boundary)")
    ax3.axvline(t_sing / 60.0, color="0.3", ls=":", lw=1.2, label="P/2 — conditioning-excluded")

    # Mark where chief-crossing is a *separate* failure mode not captured
    # by the clearance metric.
    if chief_cross_arr.any():
        ax3.scatter(t_arr[chief_cross_arr], clearance_arr[chief_cross_arr],
                     marker="v", color="tab:brown", s=22, zorder=5,
                     label="also chief-crossing violation (separate failure mode)")

    ax3.plot(M2_REFERENCE_T / 60.0,
              [r for r in rows if r[0] == M2_REFERENCE_T][0][4].min_clearance_outside_corridor_m
              if [r for r in rows if r[0] == M2_REFERENCE_T] else np.nan,
              "o", color="black", ms=8, zorder=6, label="M2 reference (1800 s)")

    if refined_t is not None:
        ax3.plot(refined_t / 60.0, refined_cc.min_clearance_outside_corridor_m, "*",
                  color="tab:green", ms=16, mec="black", mew=0.6, zorder=6,
                  label=f"refined M4 minimum (T={refined_t:.1f} s)")

    ax3.set_xlabel("transfer time T [minutes]")
    ax3.set_ylabel("clearance outside corridor [m]  (+safe / −violation)")
    ax3.set_title(
        "M4 — Keep-out clearance margin vs. transfer time\n"
        "(chief-crossing is a distinct failure mode, marked separately — see ▽ markers)",
        fontsize=10,
    )
    ax3.legend(loc="center right", fontsize=7.5, framealpha=0.92)
    fig3.tight_layout()
    fig3_path = FIGURES_DIR / "m4_constraint_margin_vs_transfer_time.png"
    fig3.savefig(fig3_path, dpi=150)
    print(f"Wrote {fig3_path}")


if __name__ == "__main__":
    main()
