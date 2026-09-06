# Rendezvous & Proximity Operations Trade

A rigorously verified Clohessy–Wiltshire (CW) / Hill-frame relative-motion analysis
for a satellite-servicing rendezvous, culminating in:

1. an approach-trajectory plot, and
2. a Δv-vs-transfer-time trade table.

**Status: Milestone 2 complete.** The Clohessy–Wiltshire state-transition matrix,
Φrv conditioning/singularity diagnostics, and the two-impulse boundary-value solver
are implemented and verified (92 tests, `pytest -W error`) against the plan in
[`DESIGN.md`](DESIGN.md) §9. The representative T = 1800 s transfer matches the M1
hand-calculation target to 0.02 mm/s. **The full Δv-vs-transfer-time trade (Milestone
3) is not implemented yet.**

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
│                                 two-impulse method, hand calculations, trade plan,
│                                 verification plan, M2 implementation/verification
│                                 results, limitations
├── README.md                     This file
├── pyproject.toml                Packaging + pytest configuration
├── src/rendezvous_cw/
│   ├── orbit.py                   Circular chief-orbit utilities (a, n, period)
│   ├── cw.py                       Closed-form CW STM blocks + propagate_cw()
│   ├── conditioning.py             Phi_rv determinant/condition-number diagnostics
│   └── rendezvous.py                Two-impulse boundary-value solver
├── tests/                        92 tests covering the M2 verification plan
├── scripts/
│   └── m2_plot_representative_transfer.py   Generates the M2 diagnostic figure
└── figures/
    └── m2_representative_transfer.png       M2 diagnostic figure (not portfolio-final)
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
- **M2 (this milestone):** CW STM, Φrv conditioning policy, and two-impulse solver
  implemented and verified (92 tests) against the M2 verification plan and the M1 hand
  calculation. One diagnostic figure (`figures/m2_representative_transfer.png`).
- **M3 (next, pending approval):** Δv-vs-transfer-time trade table/script and the
  portfolio-final approach-trajectory plot.
