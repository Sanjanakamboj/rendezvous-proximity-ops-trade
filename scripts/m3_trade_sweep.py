"""M3: Delta-v vs. transfer-time trade study for the fixed M1 scenario.

Uses ONLY the verified M2 solver (rendezvous.solve_two_impulse) and M2
conditioning policy (conditioning.py), wrapped by the M3 trade-study
module (trade.py). No CW/rendezvous equations are duplicated here.

Produces:
  - results/m3_transfer_trade.csv        curated representative-time table
  - results/m3_full_sweep.csv             full uniform sweep (diagnostic/reproducibility)
  - figures/m3_dv_vs_transfer_time.png    headline Delta-v vs. T figure
  - figures/m3_approach_trajectory_comparison.png
  - figures/m3_conditioning_vs_transfer_time.png

Run:
    python scripts/m3_trade_sweep.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from rendezvous_cw.conditioning import COND_THRESHOLD
from rendezvous_cw.cw import propagate_cw
from rendezvous_cw.orbit import chief_mean_motion_and_period
from rendezvous_cw.trade import TradeRecord, evaluate_transfer_time, find_best_safe_transfer, sweep_transfer_times

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
SWEEP_STEP = 5.0  # seconds


def build_main_sweep_times() -> np.ndarray:
    """Uniform 5 s sweep over [300, 5100] s (base resolution)."""
    return np.arange(T_MIN_DOMAIN, T_MAX_DOMAIN + 1.0, SWEEP_STEP)


def build_singularity_probe_times(t_sing: float) -> np.ndarray:
    """Dense probe of transfer times bracketing the P/2 singularity.

    Resolution (0.5 ms) is fine enough to resolve the actual unsafe
    guard band (~19 ms wide, per conditioning.py), used only to
    characterize/plot the singularity itself (Figure 3) and to
    demonstrate real, solver-derived unsafe records -- not to inflate
    the headline Delta-v sweep.
    """
    return t_sing + np.arange(-0.05, 0.05 + 1e-9, 0.0005)


def main() -> None:
    n, period = chief_mean_motion_and_period()
    t_sing = period / 2.0

    main_times = build_main_sweep_times()
    probe_times = build_singularity_probe_times(t_sing)

    # Full sweep = main grid + singularity probe, merged/sorted/deduped.
    all_times = np.unique(np.concatenate([main_times, probe_times]))
    all_records = sweep_transfer_times(R0, V0_MINUS, RF, VF_PLUS, all_times, n, period)

    safe_records = [r for r in all_records if r.safe]
    unsafe_records = [r for r in all_records if not r.safe]

    print(f"n = {n!r} rad/s, period = {period!r} s, P/2 = {t_sing!r} s")
    print(f"Total sweep points: {len(all_records)}  safe: {len(safe_records)}  unsafe: {len(unsafe_records)}")
    if unsafe_records:
        unsafe_t = sorted(r.t_s for r in unsafe_records)
        print(f"Unsafe interval observed: [{unsafe_t[0]:.6f}, {unsafe_t[-1]:.6f}] s "
              f"(width {unsafe_t[-1] - unsafe_t[0]:.6f} s), centered near P/2 = {t_sing:.6f} s")

    # --- Branch split around the singularity --------------------------------
    branch1 = [r for r in safe_records if r.t_s < t_sing]  # short transfers
    branch2 = [r for r in safe_records if r.t_s > t_sing]  # long transfers

    branch1_coarse_min = min(branch1, key=lambda r: r.dv_total)
    branch2_coarse_min = min(branch2, key=lambda r: r.dv_total)

    branch1_refined = find_best_safe_transfer(
        branch1, R0, V0_MINUS, RF, VF_PLUS, n, period, refine=True
    )
    branch2_refined = find_best_safe_transfer(
        branch2, R0, V0_MINUS, RF, VF_PLUS, n, period, refine=True
    )

    overall_refined = min(
        [branch1_refined, branch2_refined], key=lambda r: r.dv_total
    )

    print("\nBranch 1 (T < P/2, short transfers):")
    print(f"  coarse min:  T={branch1_coarse_min.t_s:.3f} s  dv_total={branch1_coarse_min.dv_total*1000:.4f} mm/s")
    print(f"  refined min: T={branch1_refined.t_s:.6f} s  dv_total={branch1_refined.dv_total*1000:.4f} mm/s")

    print("Branch 2 (T > P/2, long transfers):")
    print(f"  coarse min:  T={branch2_coarse_min.t_s:.3f} s  dv_total={branch2_coarse_min.dv_total*1000:.4f} mm/s")
    print(f"  refined min: T={branch2_refined.t_s:.6f} s  dv_total={branch2_refined.dv_total*1000:.4f} mm/s")

    print(f"\nOverall refined minimum-safe transfer over [{T_MIN_DOMAIN},{T_MAX_DOMAIN}] s: "
          f"T={overall_refined.t_s:.6f} s, dv_total={overall_refined.dv_total*1000:.4f} mm/s")

    # --- Grid-resolution convergence check (also exercised in tests) --------
    for step in (20.0, 10.0, 5.0):
        times_k = np.arange(T_MIN_DOMAIN, T_MAX_DOMAIN + 1.0, step)
        recs_k = sweep_transfer_times(R0, V0_MINUS, RF, VF_PLUS, times_k, n, period)
        safe_k = [r for r in recs_k if r.safe]
        branch2_k = [r for r in safe_k if r.t_s > t_sing]
        best_k = min(branch2_k, key=lambda r: r.dv_total)
        print(f"  grid step={step:5.1f} s -> branch-2 coarse min T={best_k.t_s:.2f} s, "
              f"dv_total={best_k.dv_total*1000:.4f} mm/s")

    # --- M2 reference comparison ---------------------------------------------
    ref_1800 = evaluate_transfer_time(R0, V0_MINUS, RF, VF_PLUS, 1800.0, n, period)
    pct_change = 100.0 * (overall_refined.dv_total - ref_1800.dv_total) / ref_1800.dv_total
    print(f"\nM2 T=1800 s reference: dv_total={ref_1800.dv_total*1000:.4f} mm/s")
    print(f"Refined minimum vs. M2 1800 s reference: {pct_change:+.2f}%")

    worst_closure = max(r.closure_residual_m for r in safe_records)
    max_safe_cond = max(r.cond_phi_rv for r in safe_records)
    print(f"Worst terminal closure residual over all safe sweep points: {worst_closure:.3e} m")
    print(f"Maximum condition number actually admitted as safe: {max_safe_cond:.3e}")

    # ==========================================================================
    # Decision table (curated representative times)
    # ==========================================================================
    representative_times = {
        "300 s": 300.0,
        "600 s": 600.0,
        "900 s": 900.0,
        "1200 s": 1200.0,
        "1800 s (M2 reference)": 1800.0,
        "2400 s": 2400.0,
        "just before singular exclusion zone": t_sing - 0.02,
        "just after singular exclusion zone": t_sing + 0.02,
        "3600 s": 3600.0,
        "4200 s": 4200.0,
        "4800 s": 4800.0,
        "5100 s": 5100.0,
        "refined minimum-dv (branch 2)": overall_refined.t_s,
    }

    table_rows = []
    for label, t in representative_times.items():
        rec = evaluate_transfer_time(R0, V0_MINUS, RF, VF_PLUS, t, n, period)
        table_rows.append((label, rec))

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = RESULTS_DIR / "m3_transfer_trade.csv"
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "label",
                "T_s",
                "T_min",
                "T_over_P",
                "safe",
                "cond_phi_rv",
                "dv1_m_s",
                "dv2_m_s",
                "dv_total_m_s",
                "closure_error_m",
            ]
        )
        for label, rec in table_rows:
            writer.writerow(
                [
                    label,
                    f"{rec.t_s:.6f}",
                    f"{rec.t_min:.6f}",
                    f"{rec.t_over_p:.8f}",
                    rec.safe,
                    f"{rec.cond_phi_rv:.6e}",
                    f"{rec.dv1_mag:.6f}" if rec.safe else "NaN",
                    f"{rec.dv2_mag:.6f}" if rec.safe else "NaN",
                    f"{rec.dv_total:.6f}" if rec.safe else "NaN",
                    f"{rec.closure_residual_m:.3e}" if rec.safe else "NaN",
                ]
            )
    print(f"\nWrote {csv_path}")

    # Full sweep CSV (diagnostic/reproducibility, not the curated table).
    full_csv_path = RESULTS_DIR / "m3_full_sweep.csv"
    with open(full_csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["T_s", "T_over_P", "safe", "cond_phi_rv", "dv1_m_s", "dv2_m_s", "dv_total_m_s", "closure_error_m"]
        )
        for rec in all_records:
            writer.writerow(
                [
                    f"{rec.t_s:.6f}",
                    f"{rec.t_over_p:.8f}",
                    rec.safe,
                    f"{rec.cond_phi_rv:.6e}",
                    f"{rec.dv1_mag:.6f}" if rec.safe else "",
                    f"{rec.dv2_mag:.6f}" if rec.safe else "",
                    f"{rec.dv_total:.6f}" if rec.safe else "",
                    f"{rec.closure_residual_m:.3e}" if rec.safe else "",
                ]
            )
    print(f"Wrote {full_csv_path}")

    # ==========================================================================
    # Figure 1 -- headline Delta-v vs. transfer time
    # ==========================================================================
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    def branch_arrays(records: list[TradeRecord]):
        records = sorted(records, key=lambda r: r.t_s)
        t = np.array([r.t_s for r in records])
        dv1 = np.array([r.dv1_mag for r in records]) * 1000.0
        dv2 = np.array([r.dv2_mag for r in records]) * 1000.0
        dvt = np.array([r.dv_total for r in records]) * 1000.0
        return t, dv1, dv2, dvt

    # Use only the base (main_times) grid for the headline curve -- the
    # dense singularity probe is for Figure 3, not for this plot's line
    # density.
    main_records = [r for r in all_records if r.t_s in set(main_times.tolist())]
    b1 = sorted([r for r in main_records if r.safe and r.t_s < t_sing], key=lambda r: r.t_s)
    b2 = sorted([r for r in main_records if r.safe and r.t_s > t_sing], key=lambda r: r.t_s)

    fig1, ax1 = plt.subplots(figsize=(9.5, 6))

    for branch, label_suffix in ((b1, ""), (b2, "")):
        if not branch:
            continue
        t_arr, dv1_arr, dv2_arr, dvt_arr = branch_arrays(branch)
        t_min_arr = t_arr / 60.0
        ax1.plot(t_min_arr, dv1_arr, color="tab:blue", lw=1.0, alpha=0.45,
                  label="|Δv1|" if branch is b1 else None)
        ax1.plot(t_min_arr, dv2_arr, color="tab:orange", lw=1.0, alpha=0.45,
                  label="|Δv2|" if branch is b1 else None)
        ax1.plot(t_min_arr, dvt_arr, color="tab:red", lw=2.2,
                  label="Δv_total" if branch is b1 else None)

    ax1.axvline(t_sing / 60.0, color="0.3", ls="--", lw=1.2,
                 label="P/2 — excluded/unsafe region (Φrv ill-conditioned)")

    ax1.plot(1800.0 / 60.0, ref_1800.dv_total * 1000.0, "o", color="black", ms=9,
              label=f"M2 reference (T=1800 s, {ref_1800.dv_total*1000:.1f} mm/s)")
    ax1.plot(overall_refined.t_s / 60.0, overall_refined.dv_total * 1000.0, "*",
              color="tab:green", ms=16, mec="black", mew=0.6,
              label=f"refined minimum (T={overall_refined.t_s:.0f} s, {overall_refined.dv_total*1000:.1f} mm/s)")

    ax1.set_xlabel("transfer time T [minutes]")
    ax1.set_ylabel("Δv [mm/s]")
    ax1.set_title(
        "M3 — Δv vs. transfer time, CW linearized dynamics, fixed M1 boundary conditions\n"
        "(unsafe ill-conditioned region near T = P/2 excluded from the search)",
        fontsize=10,
    )
    ax1.legend(loc="upper right", fontsize=8, framealpha=0.9)
    ax1.set_ylim(0, 2500)
    ax1.set_xlim(main_times.min() / 60.0, main_times.max() / 60.0)

    # The excluded interval is only ~19 ms wide (see conditioning.py / Figure 3)
    # -- far too narrow to see at this minutes-wide axis scale, even though it
    # is real and solver-derived (see figures/m3_conditioning_vs_transfer_time.png).
    # Annotate this explicitly rather than exaggerating the shaded band.
    ax1.annotate(
        "excluded region is real but only ~19 ms\nwide at this scale — see Figure 3",
        xy=(t_sing / 60.0, 2500), xytext=(t_sing / 60.0 - 24, 1650),
        fontsize=8, ha="left",
        arrowprops=dict(arrowstyle="->", color="0.3", lw=1.0),
    )

    fig1.tight_layout()
    fig1_path = FIGURES_DIR / "m3_dv_vs_transfer_time.png"
    fig1.savefig(fig1_path, dpi=150)
    print(f"Wrote {fig1_path}")

    # ==========================================================================
    # Figure 2 -- approach trajectory comparison
    # ==========================================================================
    comparison_times = {
        "fast (T=600 s)": 600.0,
        "M2 reference (T=1800 s)": 1800.0,
        f"refined minimum (T={overall_refined.t_s:.0f} s)": overall_refined.t_s,
        "slow safe (T=5100 s)": 5100.0,
    }

    fig2, ax2 = plt.subplots(figsize=(8, 7))
    colors = ["tab:blue", "black", "tab:green", "tab:purple"]
    for (label, t), color in zip(comparison_times.items(), colors):
        rec = evaluate_transfer_time(R0, V0_MINUS, RF, VF_PLUS, t, n, period)
        # Recompute via the solver directly (evaluate_transfer_time only
        # returns the summary TradeRecord; the full result object with
        # v0_plus is needed here to propagate the trajectory for plotting).
        from rendezvous_cw.rendezvous import solve_two_impulse

        full_result = solve_two_impulse(R0, V0_MINUS, RF, VF_PLUS, t, n)
        t_samples = np.linspace(0.0, t, 300)
        post_burn1 = np.concatenate([R0, full_result.v0_plus])
        states = np.array([propagate_cw(post_burn1, tt, n) for tt in t_samples])
        x, y = states[:, 0], states[:, 1]
        ax2.plot(y, x, "-", color=color, lw=1.8,
                  label=f"{label}: Δv={rec.dv_total*1000:.0f} mm/s")

    ax2.plot(R0[1], R0[0], "o", color="tab:green", ms=10, zorder=5, label="r0 (start)")
    ax2.plot(RF[1], RF[0], "s", color="tab:red", ms=10, zorder=5, label="rf (target)")
    ax2.plot(0, 0, "*", color="k", ms=14, zorder=5, label="chief (origin)")
    ax2.axhline(0, color="0.85", lw=0.8, zorder=0)
    ax2.axvline(0, color="0.85", lw=0.8, zorder=0)
    ax2.set_xlabel("y — along-track [m]")
    ax2.set_ylabel("x — radial [m]")
    ax2.set_title(
        "M3 — Approach trajectory comparison (LVLH x–y projection)\n"
        "CW linearized dynamics; lower Δv trades against transfer time / path shape",
        fontsize=10,
    )
    ax2.invert_xaxis()
    ax2.set_aspect("equal", adjustable="datalim")
    ax2.legend(loc="upper center", fontsize=7.5, framealpha=0.9, ncol=1,
                bbox_to_anchor=(0.5, -0.08))
    fig2.tight_layout()
    fig2_path = FIGURES_DIR / "m3_approach_trajectory_comparison.png"
    fig2.savefig(fig2_path, dpi=150, bbox_inches="tight")
    print(f"Wrote {fig2_path}")

    # ==========================================================================
    # Figure 3 -- conditioning near the half-period singularity
    # ==========================================================================
    probe_records = [r for r in all_records if r.t_s in set(probe_times.tolist())]
    probe_records = sorted(probe_records, key=lambda r: r.t_s)
    t_probe = np.array([r.t_s for r in probe_records])
    cond_probe = np.array([r.cond_phi_rv for r in probe_records])
    safe_mask = np.array([r.safe for r in probe_records])

    fig3, ax3 = plt.subplots(figsize=(9, 5.5))
    dt_ms = (t_probe - t_sing) * 1000.0  # ms offset from P/2
    ax3.plot(dt_ms, cond_probe, "-", color="tab:blue", lw=1.5, marker=".", ms=3)
    ax3.axhline(COND_THRESHOLD, color="tab:red", ls="--", lw=1.3,
                 label=f"cond threshold = {COND_THRESHOLD:.0e}")
    ax3.axvline(0.0, color="0.3", ls=":", lw=1.2, label="P/2 (exact singularity)")

    unsafe_dt = dt_ms[~safe_mask]
    if unsafe_dt.size:
        ax3.axvspan(unsafe_dt.min(), unsafe_dt.max(), color="0.6", alpha=0.3,
                     label="excluded/unsafe region")

    ax3.set_yscale("log")
    ax3.set_xlabel("time offset from P/2 [ms]")
    ax3.set_ylabel("cond(Φrv)  (2-norm condition number)")
    ax3.set_title(
        "M3 — Φrv conditioning near the T = P/2 half-period singularity\n"
        f"(P/2 = {t_sing:.3f} s; verification/supporting figure)",
        fontsize=10,
    )
    ax3.legend(loc="upper right", fontsize=8)
    fig3.tight_layout()
    fig3_path = FIGURES_DIR / "m3_conditioning_vs_transfer_time.png"
    fig3.savefig(fig3_path, dpi=150)
    print(f"Wrote {fig3_path}")


if __name__ == "__main__":
    main()
