"""M5: nonlinear two-body validation of the M4-selected constrained rendezvous,
and a local robustness refinement in the fast-transfer feasible window.

Uses ONLY the verified M2 solver (rendezvous.solve_two_impulse), M2 STM
(cw.propagate_cw), M4 constraint classification (constraints.py), and the
new, independent M5 nonlinear two-body model (nonlinear.py). No CW
equation is changed; nonlinear.py never calls cw.py.

Produces:
  - results/m5_representative_cases.csv     the 5 required representative cases
  - results/m5_local_robustness_scan.csv    365-390 s local CW+NL scan
  - results/m5_separation_sensitivity.csv   250/500/1000/1500 m sensitivity
  - figures/m5_cw_vs_nonlinear_trajectory.png
  - figures/m5_cw_error_vs_transfer_time.png
  - figures/m5_local_robustness_trade.png

Run:
    python scripts/m5_nonlinear_validation.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rendezvous_cw.constraints import CORRIDOR_X_BOUND, CORRIDOR_Y_MAX, CORRIDOR_Y_MIN, R_KOZ, classify_trajectory_constraints
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.nonlinear import classify_nonlinear_trajectory, propagate_relative_nonlinear
from rendezvous_cw.orbit import CHIEF_ALTITUDE, chief_mean_motion_and_period, orbital_radius
from rendezvous_cw.rendezvous import solve_two_impulse

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = REPO_ROOT / "results"
FIGURES_DIR = REPO_ROOT / "figures"

R0 = np.array([0.0, -1000.0, 30.0])
V0_MINUS = np.array([0.0, 0.0, 0.0])
RF = np.array([0.0, -50.0, 0.0])
VF_PLUS = np.array([0.0, 0.0, 0.0])

M2_REF_T = 1800.0
M3_SHORT_MIN_T = 2539.592622
M3_LONG_MIN_T = 4868.130154465808
M4_SELECTED_T = 386.56385551988984
M4_PRACTICAL_T = 380.0


def compare_case(t: float, n: float, period: float, a: float, n_samples: int = 80):
    """CW vs nonlinear comparison for one transfer time. Returns a dict of results."""
    result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, t, n)
    times = np.linspace(0.0, t, n_samples)
    state0 = np.concatenate([R0, result.v0_plus])
    cw_states = np.array([propagate_cw(state0, tt, n) for tt in times])
    nl = propagate_relative_nonlinear(R0, result.v0_plus, times, a)

    pos_err = np.linalg.norm(nl.lvlh_states[:, :3] - cw_states[:, :3], axis=1)
    vel_err = np.linalg.norm(nl.lvlh_states[:, 3:] - cw_states[:, 3:], axis=1)

    cw_cc = classify_trajectory_constraints(R0, result.v0_plus, t, n, sample_dt=1.0)
    nl_cc = classify_nonlinear_trajectory(R0, result.v0_plus, t, a, sample_dt=1.0)

    return {
        "t": t, "t_over_p": t / period, "result": result,
        "times": times, "cw_states": cw_states, "nl_states": nl.lvlh_states,
        "pos_err": pos_err, "vel_err": vel_err,
        "terminal_pos_err_m": float(pos_err[-1]), "terminal_vel_err_mm_s": float(vel_err[-1] * 1000.0),
        "max_pos_dev_m": float(pos_err.max()),
        "cw_cc": cw_cc, "nl_cc": nl_cc,
    }


def main() -> None:
    n, period = chief_mean_motion_and_period()
    a = orbital_radius(CHIEF_ALTITUDE)

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    # ==========================================================================
    # 1. Representative case comparison
    # ==========================================================================
    cases = {
        "M2 reference": M2_REF_T,
        "M3 short-branch minimum": M3_SHORT_MIN_T,
        "M3 long-branch/global minimum": M3_LONG_MIN_T,
        "M4 selected (CW boundary optimum)": M4_SELECTED_T,
        "M4 practical margin candidate": M4_PRACTICAL_T,
    }

    print("=" * 100)
    print("Representative case comparison (CW vs. nonlinear two-body)")
    print("=" * 100)
    case_results = {}
    csv_rows = []
    for label, t in cases.items():
        c = compare_case(t, n, period, a)
        case_results[label] = c
        print(f"\n--- {label} (T={t:.3f} s, T/P={c['t_over_p']:.4f}) ---")
        print(f"  CW dv_total={c['result'].dv_total*1000:.3f} mm/s   CW margin={c['cw_cc'].min_clearance_outside_corridor_m}")
        print(f"  terminal pos err={c['terminal_pos_err_m']:.4f} m   terminal vel err={c['terminal_vel_err_mm_s']:.4f} mm/s")
        print(f"  max CW-vs-NL pos deviation={c['max_pos_dev_m']:.4f} m")
        print(f"  NL margin={c['nl_cc'].min_clearance_outside_corridor_m}   NL safe={c['nl_cc'].safe}   NL reason={c['nl_cc'].violation_reason}")
        print(f"  CW max_y={c['cw_cc'].max_y_m:.3f} m   NL max_y={c['nl_cc'].max_y_m:.3f} m")
        csv_rows.append({
            "label": label, "T_s": t, "T_over_P": c["t_over_p"],
            "cw_dv_total_m_s": c["result"].dv_total,
            "cw_margin_m": c["cw_cc"].min_clearance_outside_corridor_m,
            "cw_safe": c["cw_cc"].safe,
            "nl_terminal_pos_err_m": c["terminal_pos_err_m"],
            "nl_terminal_vel_err_mm_s": c["terminal_vel_err_mm_s"],
            "max_cw_vs_nl_pos_dev_m": c["max_pos_dev_m"],
            "nl_margin_m": c["nl_cc"].min_clearance_outside_corridor_m,
            "nl_safe": c["nl_cc"].safe,
            "nl_violation_reason": c["nl_cc"].violation_reason,
            "cw_max_y_m": c["cw_cc"].max_y_m,
            "nl_max_y_m": c["nl_cc"].max_y_m,
        })

    with open(RESULTS_DIR / "m5_representative_cases.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)
    print(f"\nWrote {RESULTS_DIR / 'm5_representative_cases.csv'}")

    # ==========================================================================
    # 2. Local robustness scan (365-390 s)
    # ==========================================================================
    print("\n" + "=" * 100)
    print("Local robustness scan (365-390 s, 1 s step)")
    print("=" * 100)
    local_times = np.arange(365.0, 390.0 + 1.0, 1.0)
    local_rows = []
    for t in local_times:
        result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, t, n)
        cw_cc = classify_trajectory_constraints(R0, result.v0_plus, t, n, sample_dt=1.0)
        nl_cc = classify_nonlinear_trajectory(R0, result.v0_plus, t, a, sample_dt=1.0)
        local_rows.append({
            "T_s": t, "dv_total_m_s": result.dv_total,
            "cw_margin_m": cw_cc.min_clearance_outside_corridor_m, "cw_safe": cw_cc.safe,
            "nl_margin_m": nl_cc.min_clearance_outside_corridor_m, "nl_safe": nl_cc.safe,
        })
        print(f"  T={t:6.1f}  dv={result.dv_total*1000:9.2f} mm/s  "
              f"CWmargin={cw_cc.min_clearance_outside_corridor_m:9.5f}  "
              f"NLmargin={nl_cc.min_clearance_outside_corridor_m:9.5f}  NLsafe={nl_cc.safe}")

    with open(RESULTS_DIR / "m5_local_robustness_scan.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(local_rows[0].keys()))
        writer.writeheader()
        writer.writerows(local_rows)
    print(f"\nWrote {RESULTS_DIR / 'm5_local_robustness_scan.csv'}")

    # M4 boundary optimum's actual nonlinear margin.
    selected_case = case_results["M4 selected (CW boundary optimum)"]
    selected_nl_margin = selected_case["nl_cc"].min_clearance_outside_corridor_m
    print(f"\nM4 selected transfer (T={M4_SELECTED_T:.3f} s) nonlinear margin: "
          f"{selected_nl_margin:.6f} m  (PASS: {selected_case['nl_cc'].safe})")

    # --- Final M5 recommendation: T=380 s (chosen for its much larger,
    # still-modest-Delta-v-cost nonlinear margin -- see DESIGN.md section 15).
    final_t = M4_PRACTICAL_T
    final_case = case_results["M4 practical margin candidate"]
    final_result = final_case["result"]
    print(f"\nFinal M5 recommended transfer: T={final_t:.3f} s "
          f"(T/P={final_t/period:.6f})")
    print(f"  CW dv_total={final_result.dv_total*1000:.4f} mm/s")
    print(f"  CW margin={final_case['cw_cc'].min_clearance_outside_corridor_m:.6f} m")
    print(f"  NL margin={final_case['nl_cc'].min_clearance_outside_corridor_m:.6f} m")
    print(f"  NL safe={final_case['nl_cc'].safe}")
    print(f"  terminal position miss (NL vs. rf)={final_case['terminal_pos_err_m']:.6f} m")

    # ==========================================================================
    # 3. Separation sensitivity study (fixed T=1800 s)
    # ==========================================================================
    print("\n" + "=" * 100)
    print("Separation sensitivity study (T=1800 s, direction of r0 preserved)")
    print("=" * 100)
    r0_unit_mag = np.linalg.norm(R0)
    sep_rows = []
    T_sens = 1800.0
    for mag in (250.0, 500.0, 1000.0, 1500.0):
        r0_scaled = R0 * (mag / r0_unit_mag)
        result = solve_two_impulse(r0_scaled, V0_MINUS, RF, VF_PLUS, T_sens, n)
        times = np.linspace(0.0, T_sens, 60)
        state0 = np.concatenate([r0_scaled, result.v0_plus])
        cw_states = np.array([propagate_cw(state0, tt, n) for tt in times])
        nl = propagate_relative_nonlinear(r0_scaled, result.v0_plus, times, a)
        pos_err = np.linalg.norm(nl.lvlh_states[:, :3] - cw_states[:, :3], axis=1)
        sep_rows.append({
            "separation_m": mag, "r0_over_a": mag / a,
            "terminal_pos_err_m": float(pos_err[-1]), "max_pos_dev_m": float(pos_err.max()),
        })
        print(f"  |r0|={mag:6.0f} m   term_err={pos_err[-1]:9.5f} m   max_dev={pos_err.max():9.5f} m   ||r0||/a={mag/a:.3e}")

    with open(RESULTS_DIR / "m5_separation_sensitivity.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(sep_rows[0].keys()))
        writer.writeheader()
        writer.writerows(sep_rows)
    print(f"\nWrote {RESULTS_DIR / 'm5_separation_sensitivity.csv'}")

    # Check quadratic-ish scaling explicitly.
    errs = [r["terminal_pos_err_m"] for r in sep_rows]
    mags = [r["separation_m"] for r in sep_rows]
    for i in range(1, len(errs)):
        ratio_err = errs[i] / errs[i - 1]
        ratio_mag = mags[i] / mags[i - 1]
        print(f"  {mags[i-1]:.0f}->{mags[i]:.0f} m: separation ratio={ratio_mag:.3f}, "
              f"error ratio={ratio_err:.3f} (quadratic-scaling expectation: {ratio_mag**2:.3f})")

    # ==========================================================================
    # 4. Nonlinear-corrected second burn for the final recommended transfer
    # ==========================================================================
    print("\n" + "=" * 100)
    print(f"Nonlinear-corrected second-burn analysis (T={final_t:.1f} s)")
    print("=" * 100)
    v_final_nl = final_case["nl_states"][-1, 3:]
    r_final_nl = final_case["nl_states"][-1, :3]
    dv2_nl_corrected = VF_PLUS - v_final_nl
    dv2_nl_corrected_mag = float(np.linalg.norm(dv2_nl_corrected))
    total_dv_nl_corrected = final_result.dv1_mag + dv2_nl_corrected_mag
    print(f"  CW dv2 magnitude:             {final_result.dv2_mag*1000:.4f} mm/s")
    print(f"  NL-corrected dv2 magnitude:   {dv2_nl_corrected_mag*1000:.4f} mm/s")
    print(f"  difference:                   {(dv2_nl_corrected_mag-final_result.dv2_mag)*1000:.6f} mm/s")
    print(f"  total dv with NL-corrected 2nd burn: {total_dv_nl_corrected*1000:.4f} mm/s "
          f"(vs. CW-only total {final_result.dv_total*1000:.4f} mm/s)")
    print(f"  terminal position miss NOT corrected by this velocity-only burn: "
          f"{np.linalg.norm(r_final_nl - RF):.6f} m")

    # ==========================================================================
    # Figure 1 -- CW vs nonlinear selected trajectory
    # ==========================================================================
    c = final_case
    fig1, ax1 = plt.subplots(figsize=(8, 7))

    theta = np.linspace(0, 2 * np.pi, 200)
    ax1.plot(R_KOZ * np.cos(theta), R_KOZ * np.sin(theta), "--", color="tab:red", lw=1.2,
              label=f"{R_KOZ:.0f} m keep-out sphere")
    corridor_y = [CORRIDOR_Y_MIN, CORRIDOR_Y_MAX, CORRIDOR_Y_MAX, CORRIDOR_Y_MIN, CORRIDOR_Y_MIN]
    corridor_x = [-CORRIDOR_X_BOUND, -CORRIDOR_X_BOUND, CORRIDOR_X_BOUND, CORRIDOR_X_BOUND, -CORRIDOR_X_BOUND]
    ax1.plot(corridor_y, corridor_x, "-", color="tab:orange", lw=1.6, label="approach corridor")

    ax1.plot(c["cw_states"][:, 1], c["cw_states"][:, 0], "-", color="tab:blue", lw=2.0, label="CW trajectory")
    ax1.plot(c["nl_states"][:, 1], c["nl_states"][:, 0], "--", color="tab:green", lw=1.6,
              label="nonlinear two-body trajectory")
    ax1.plot(R0[1], R0[0], "o", color="tab:green", ms=9, zorder=5, label="r0 (start)")
    ax1.plot(RF[1], RF[0], "s", color="tab:red", ms=9, zorder=5, label="rf (target)")

    ax1.set_xlim(50, -180)
    ax1.set_ylim(-130, 130)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlabel("y — along-track [m]")
    ax1.set_ylabel("x — radial [m]")
    ax1.set_title(
        f"M5 — CW vs. nonlinear two-body validation (T={final_t:.0f} s, same impulsive burn 1)\n"
        "geometric constraints overlaid",
        fontsize=10,
    )
    ax1.legend(loc="upper center", fontsize=7.5, framealpha=0.92, bbox_to_anchor=(0.5, -0.09))

    ax1_inset = fig1.add_axes([0.15, 0.14, 0.24, 0.16])
    ax1_inset.plot(c["times"], c["pos_err"], "-", color="tab:purple", lw=1.3)
    ax1_inset.set_xlabel("t [s]", fontsize=7)
    ax1_inset.set_ylabel("|Δr| [m]", fontsize=7)
    ax1_inset.set_title("CW-vs-NL position error", fontsize=7)
    ax1_inset.tick_params(labelsize=6)

    fig1_path = FIGURES_DIR / "m5_cw_vs_nonlinear_trajectory.png"
    fig1.savefig(fig1_path, dpi=150, bbox_inches="tight")
    print(f"\nWrote {fig1_path}")

    # ==========================================================================
    # Figure 2 -- CW error vs transfer time (representative cases)
    # ==========================================================================
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    labels_ordered = list(cases.keys())
    ts = [cases[l] for l in labels_ordered]
    term_errs = [case_results[l]["terminal_pos_err_m"] for l in labels_ordered]
    max_devs = [case_results[l]["max_pos_dev_m"] for l in labels_ordered]
    colors = ["black", "tab:purple", "tab:red", "tab:green", "tab:blue"]

    for t, te, md, lab, col in zip(ts, term_errs, max_devs, labels_ordered, colors):
        ax2.scatter([t], [te], marker="o", s=70, color=col, zorder=5,
                     label=f"{lab} (T={t:.0f} s)")
        ax2.scatter([t], [md], marker="^", s=50, facecolors="none", edgecolors=col, zorder=5)

    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlabel("transfer time T [s]  (log scale)")
    ax2.set_ylabel("position error [m]  (log scale)")
    ax2.set_title("M5 — CW linearization error vs. transfer time", fontsize=11)
    ax2.legend(loc="upper left", fontsize=7.5, framealpha=0.92)
    ax2.grid(True, which="both", alpha=0.25)

    fig2.tight_layout(rect=(0, 0.06, 1, 1))
    fig2.text(
        0.5, 0.015,
        "● terminal position error   △ max CW-vs-nonlinear trajectory deviation "
        "(coincides exactly with the terminal error in every case shown here, so no\n"
        "separate triangle marker is visible) — M4 fast-transfer cases at left, "
        "unconstrained M3 cases at right; supporting/diagnostic figure, not a final-decision plot",
        ha="center", va="bottom", fontsize=7.5, color="0.35",
    )

    fig2_path = FIGURES_DIR / "m5_cw_error_vs_transfer_time.png"
    fig2.savefig(fig2_path, dpi=150)
    print(f"Wrote {fig2_path}")

    # ==========================================================================
    # Figure 3 -- Robustness trade near the M4 boundary
    # ==========================================================================
    fig3, (ax3a, ax3b) = plt.subplots(2, 1, figsize=(9, 7.5), sharex=True)

    t_arr = np.array([r["T_s"] for r in local_rows])
    dv_arr = np.array([r["dv_total_m_s"] for r in local_rows]) * 1000.0
    nl_margin_arr = np.array([r["nl_margin_m"] for r in local_rows])

    ax3a.plot(t_arr, dv_arr, "-", color="tab:red", lw=1.8)
    ax3a.set_ylabel("Δv_total [mm/s]")
    ax3a.set_title(
        "M5 — Local robustness trade near the M4 CW boundary optimum\n"
        "(slightly more Δv buys nonlinear-model-fidelity margin)",
        fontsize=10,
    )
    ax3a.grid(True, alpha=0.3)

    ax3b.plot(t_arr, nl_margin_arr, "-", color="tab:blue", lw=1.8, label="nonlinear geometry margin")
    ax3b.axhline(0.0, color="k", ls="--", lw=1.2, label="0 m requirement")
    ax3b.set_xlabel("transfer time T [s]")
    ax3b.set_ylabel("nonlinear margin [m]")
    ax3b.grid(True, alpha=0.3)

    for ax in (ax3a, ax3b):
        ax.axvline(M4_SELECTED_T, color="tab:purple", ls=":", lw=1.5)
        ax.axvline(final_t, color="tab:green", ls=":", lw=1.5)

    ax3a.plot(M4_SELECTED_T, final_result.dv_total * 1000 if False else
              case_results["M4 selected (CW boundary optimum)"]["result"].dv_total * 1000,
              "*", color="tab:purple", ms=16, mec="black", mew=0.6,
              label=f"M4 CW boundary optimum (T={M4_SELECTED_T:.1f} s)")
    ax3a.plot(final_t, final_result.dv_total * 1000, "*", color="tab:green", ms=16,
              mec="black", mew=0.6, label=f"M5 recommended (T={final_t:.0f} s)")
    ax3a.legend(loc="upper right", fontsize=7.5, framealpha=0.92)

    ax3b.plot(M4_SELECTED_T, selected_nl_margin, "*", color="tab:purple", ms=16, mec="black", mew=0.6)
    ax3b.plot(final_t, final_case["nl_cc"].min_clearance_outside_corridor_m, "*", color="tab:green",
               ms=16, mec="black", mew=0.6)
    ax3b.legend(loc="upper left", fontsize=7.5, framealpha=0.92)

    fig3.tight_layout()
    fig3_path = FIGURES_DIR / "m5_local_robustness_trade.png"
    fig3.savefig(fig3_path, dpi=150)
    print(f"Wrote {fig3_path}")


if __name__ == "__main__":
    main()
