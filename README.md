# Rendezvous & Proximity Operations Trade

A rigorously verified Clohessy–Wiltshire (CW) / Hill-frame relative-motion analysis
for a satellite-servicing rendezvous, culminating in:

1. an approach-trajectory plot, and
2. a Δv-vs-transfer-time trade table.

**Status: Milestone 5 complete.** Building on the unmodified, verified M2–M4 solver
and constraint machinery, M5 validates the M4-selected trajectory against an
**independent nonlinear two-body model** (inertial Cartesian propagation + rigorous
LVLH↔ECI frame kinematics — no CW equations reused) and performs a local robustness
refinement. 151 tests pass (`pytest -W error`).

**Key result: the M4 CW-boundary optimum (T≈386.56 s) survives nonlinear validation,
but only by ~27.5 mm — a razor-thin, operationally fragile margin.** Because M4's
constrained optimum sat at essentially zero (~15 µm) *CW* margin by construction,
this milestone asks whether that razor's edge holds up against true nonlinear
dynamics. It does, but barely. **M5-recommended transfer: T = 380 s** — spending only
**+1.65% more Δv** (5136.0 vs. 5052.7 mm/s) buys a **~33× larger real safety margin**
(0.90 m vs. 0.028 m, nonlinear). The mathematically-minimal M4 solution is therefore
labeled **"CW-only constrained optimum — superseded by M5 robustness validation"**;
the M4 engineering history itself is unchanged. Full detail in
[`DESIGN.md`](DESIGN.md) §15.

**Headline figures:**

![CW vs. nonlinear two-body validation](figures/m5_cw_vs_nonlinear_trajectory.png)
![Local robustness trade near the M4 boundary](figures/m5_local_robustness_trade.png)

**Supporting figure:** [`figures/m5_cw_error_vs_transfer_time.png`](figures/m5_cw_error_vs_transfer_time.png)
(CW linearization error vs. transfer time across all representative cases — grows
from ~1 cm at fast M4-regime transfers to ~2.7 m at the M3 long-branch minimum,
confirming the expected small-separation/first-order-linearization scaling).

**Result tables:** [`results/m5_representative_cases.csv`](results/m5_representative_cases.csv),
[`results/m5_local_robustness_scan.csv`](results/m5_local_robustness_scan.csv),
[`results/m5_separation_sensitivity.csv`](results/m5_separation_sensitivity.csv).

**Scope — read before citing any M5 result:** the nonlinear model here is still
**unperturbed two-body** (no J2, no drag, no finite burns, no actuator limits, no
navigation/sensor error, no collision probability). A transfer that "passes modeled
geometric constraints under nonlinear two-body validation" is **not** a claim of
"collision-free", "flight-qualified", or "operationally safe" status — this milestone
demonstrates **CW-linearization robustness**, nothing more.

**M4/M3 status (superseded selection, machinery unchanged):** the M4 keep-out/corridor
screening and M3 unconstrained trade remain fully documented and correct as
originally reported — see [`DESIGN.md`](DESIGN.md) §13–§14 and
[`results/m4_constrained_trade.csv`](results/m4_constrained_trade.csv) /
[`results/m3_transfer_trade.csv`](results/m3_transfer_trade.csv) — only the
*recommended* transfer time is updated by M5's robustness finding above.

## Model scope

Linearized relative motion about a **circular** chief orbit (Clohessy–Wiltshire / Hill
equations), in the **LVLH/Hill frame**, with **impulsive burns**. No J2, drag, finite
burns, nonlinear two-body propagation, keep-out-zone logic, docking dynamics, or
collision avoidance — see [`DESIGN.md` §1 and §10](DESIGN.md) for the full, explicit
scope statement and limitations.

## Scenario (summary — full detail in `DESIGN.md` §2)

- Chief: circular LEO, 400 km altitude (`n ≈ 1.1314e-3 rad/s`, period ≈ 92.56 min).
- Deputy starts holding 1000 m behind the chief on the V-bar, with a 30 m out-of-plane
  offset (3D scenario).
- Target: 50 m behind the chief on the V-bar, in-plane, zero relative velocity.
- Transfer time traded over `T ∈ [300 s, 5100 s]`.

## Repository layout

```
├── DESIGN.md                    Frame/sign conventions, CW equations, full STM,
│                                 two-impulse method, hand calculations, M2 solver
│                                 verification, M3 trade study, M4 geometric
│                                 constraints + feasible-minimum, verification plan,
│                                 limitations
├── README.md                     This file
├── pyproject.toml                Packaging + pytest configuration
├── src/rendezvous_cw/
│   ├── orbit.py                   Circular chief-orbit utilities (a, n, period)
│   ├── cw.py                       Closed-form CW STM blocks + propagate_cw()
│   ├── conditioning.py             Phi_rv determinant/condition-number diagnostics
│   ├── rendezvous.py                Two-impulse boundary-value solver
│   ├── trade.py                     M3 transfer-time/Δv trade-study utilities
│   │                                 (wraps the M2 solver; no equations duplicated)
│   ├── constraints.py               M4 proximity-geometry screening: keep-out
│   │                                 sphere, approach corridor, no-chief-crossing
│   └── nonlinear.py                 M5 independent nonlinear two-body model +
│                                     LVLH<->ECI frame kinematics (no CW reuse)
├── tests/                        151 tests covering the M2-M5 verification plans
├── scripts/
│   ├── m2_plot_representative_transfer.py   Generates the M2 diagnostic figure
│   ├── m3_trade_sweep.py                     Generates the M3 CSV tables + figures
│   ├── m4_constrained_trade.py               Generates the M4 CSV tables + figures
│   └── m5_nonlinear_validation.py            Generates the M5 CSV tables + figures
├── results/
│   ├── m3_transfer_trade.csv                 M3 curated decision table
│   ├── m3_full_sweep.csv                      M3 full sweep (reproducibility)
│   ├── m4_constrained_trade.csv               M4 curated decision table
│   ├── m4_full_constrained_sweep.csv          M4 full sweep (reproducibility)
│   ├── m5_representative_cases.csv            M5 CW-vs-nonlinear case comparison
│   ├── m5_local_robustness_scan.csv           M5 365-390 s local scan
│   └── m5_separation_sensitivity.csv          M5 separation sensitivity study
└── figures/
    ├── m2_representative_transfer.png        M2 diagnostic figure
    ├── m3_dv_vs_transfer_time.png             M3 headline figure
    ├── m3_approach_trajectory_comparison.png  M3 supporting figure
    ├── m3_conditioning_vs_transfer_time.png   M3 verification/supporting figure
    ├── m4_constrained_dv_vs_transfer_time.png M4 headline figure
    ├── m4_approach_geometry.png                M4 headline figure
    ├── m4_constraint_margin_vs_transfer_time.png  M4 verification/supporting figure
    ├── m5_cw_vs_nonlinear_trajectory.png       M5 headline figure
    ├── m5_cw_error_vs_transfer_time.png        M5 supporting figure
    └── m5_local_robustness_trade.png           M5 headline figure
```

## Development

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -W error
```

## Roadmap

- **M1:** scenario definition, governing equations, closed-form STM, two-impulse
  method, hand calculations, verification plan. No solver code.
- **M2:** CW STM, Φrv conditioning policy, and two-impulse solver implemented and
  verified (92 tests) against the M2 verification plan and the M1 hand calculation.
  One diagnostic figure (`figures/m2_representative_transfer.png`).
- **M3:** Δv-vs-transfer-time trade study over `T ∈ [300, 5100] s` using the
  unmodified M2 solver; refined minimum-Δv transfer identified within the searched
  domain; decision table (CSV) and three figures. 15 new tests (107 total).
- **M4:** proximity-geometry screening (100 m keep-out sphere, final V-bar approach
  corridor, no-chief-crossing) layered on the unmodified M2/M3 solver; none of the M3
  minima survive; refined minimum-Δv *geometrically feasible* transfer identified
  (T≈386.6 s, ≈31× the unconstrained M3 minimum's Δv); decision table (CSV) and three
  figures. 26 new tests (133 total). Geometric screening only.
- **M5 (this milestone):** independent nonlinear two-body validation (no CW reuse) +
  LVLH↔ECI frame kinematics; the M4 optimum passes but with only ~27.5 mm of real
  margin; local robustness refinement recommends T=380 s instead (+1.65% Δv for ~33×
  more margin); CW-linearization-error scaling confirmed vs. both transfer time and
  separation magnitude (quadratic, as expected). 18 new tests (151 total). Still no
  J2/drag/finite burns/collision probability.
- **M6 (next, pending approval):** J2/drag perturbations, finite-burn modeling, and/or
  covariance/collision-probability analysis — not yet scoped in detail.
