# Rendezvous & Proximity Operations Trade

A rigorously verified Clohessy–Wiltshire (CW) / Hill-frame relative-motion analysis
for a satellite-servicing rendezvous, culminating in:

1. an approach-trajectory plot, and
2. a Δv-vs-transfer-time trade table.

**Status: Milestone 1 complete.** Scenario, governing equations, closed-form STM, the
two-impulse method, hand-verified representative numbers, and the verification plan
are documented in [`DESIGN.md`](DESIGN.md). **No solver code exists yet** — this
milestone is documentation and scaffolding only.

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
├── DESIGN.md              Frame/sign conventions, CW equations, full STM, two-impulse
│                           method, hand calculations, trade plan, verification plan,
│                           limitations
├── README.md               This file
├── pyproject.toml          Packaging + pytest configuration
├── src/rendezvous_cw/       Package (empty placeholder in M1 — CW solver is M2)
├── tests/                   Test suite (placeholder in M1 — verification tests are M2)
├── scripts/                 Analysis/plotting scripts (M2/M3)
└── figures/                 Generated figures (M2/M3)
```

## Development

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
pytest -W error
```

## Roadmap

- **M1 (this milestone):** scenario definition, governing equations, closed-form STM,
  two-impulse method, hand calculations, verification plan. No solver code.
- **M2 (next, pending approval):** implement and unit-test the CW STM and two-impulse
  solver against the verification plan in `DESIGN.md` §9.
- **M3:** approach-trajectory plot and Δv-vs-transfer-time trade table/script.
