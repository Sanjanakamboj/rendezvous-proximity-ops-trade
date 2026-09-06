# DESIGN.md — Rendezvous & Proximity Operations Trade

Status: **Milestone 1 — scenario, equations, hand calculations, verification plan.**
No CW solver code exists yet. Every number in this document is reproducible from the
arithmetic shown inline (Python used only as a calculator; no project package code is
invoked).

---

## 1. Model scope (read this first)

This project implements **linearized relative motion about a circular chief orbit**
(Clohessy–Wiltshire / Hill equations), expressed in the **LVLH (Hill) frame**, with
**impulsive burns**. Unless a later milestone explicitly says otherwise, the following
are **out of scope** for the entire project until introduced:

- J2 or any other gravity perturbation
- Atmospheric drag
- Finite-burn / continuous-thrust dynamics, actuator limits
- Nonlinear (full two-body) relative propagation
- Keep-out-zone / approach-corridor / line-of-sight constraints
- Docking mechanics / contact dynamics
- Collision avoidance logic
- Eccentric chief orbit (CW requires circular reference orbit)
- Any regime where deputy-to-chief separation is not small relative to orbital radius

These restrictions apply to Milestone 1 (this document) and will continue to apply to
M2 (solver + tests) and M3 (trade study + plots) unless a future milestone explicitly
lifts one.

---

## 2. Servicing-mission scenario

### 2.1 Chief orbit

| Quantity | Symbol | Value |
|---|---|---|
| Chief orbit type | — | Circular LEO |
| Altitude | `alt` | 400.0 km |
| Earth GM | `mu` | 398600.4418 km³/s² (WGS-84 value) |
| Earth mean equatorial radius | `Re` | 6378.137 km |
| Chief orbit radius | `a = Re + alt` | 6778.137 km |
| Chief mean motion | `n = sqrt(mu / a^3)` | 1.131366654 × 10⁻³ rad/s |
| Chief orbital period | `T_period = 2π/n` | 5553.624 s = 92.560 min |

This is an ISS-class LEO servicing altitude — realistic and a common reference case in
the rendezvous literature (Fehse, *Automated Rendezvous and Docking of Spacecraft*, uses
a similar regime).

### 2.2 Initial deputy relative state (LVLH, at transfer start, t = 0)

The deputy (servicer) begins holding on the **V-bar** (along-track axis), 1 km behind
the chief, with a small out-of-plane offset representing realistic injection/insertion
dispersion, and zero relative velocity at the moment the transfer burn is planned:

```
r0 = (x0, y0, z0) = (   0.0 m, -1000.0 m,  30.0 m )
v0 = (vx0, vy0, vz0) = (  0.0 m/s,   0.0 m/s,   0.0 m/s )
```

Sign convention (Section 3): negative `y` is **behind** the chief along the velocity
direction, i.e. this is a classic "V-bar hold point" 1 km trailing the chief, displaced
30 m out of the orbit plane.

Note: because `z0 ≠ 0`, this initial condition is **not** a full dynamical
equilibrium of the CW equations — the cross-track component evolves as a free harmonic
oscillator `z(t) = z0 cos(nt)` even absent any burn. This is intentional: it reflects a
believable, non-trivial 3D dispersion at the start of the transfer rather than an
idealized perfectly in-plane hold, and it is why this scenario is defined as **3D**, not
planar (Section 2.4).

### 2.3 Final rendezvous/approach target state (t = T)

The final target is an in-plane approach-corridor entry point 50 m behind the chief on
the V-bar, with the deputy's cross-track offset nulled and relative velocity matched to
zero (a stationkeeping hold, staged short of contact/docking, which remains out of
scope):

```
rf = (xf, yf, zf) = ( 0.0 m, -50.0 m,  0.0 m )
vf = (vxf, vyf, vzf) = ( 0.0 m/s,  0.0 m/s,  0.0 m/s )
```

Physically: the deputy transfers 950 m closer along the V-bar, removes the 30 m
out-of-plane offset, and arrives with zero relative velocity, ready for a subsequent
(out-of-scope) final approach/docking phase.

### 2.4 Planar vs. 3D

**Baseline is 3D.** The scenario carries a nonzero initial cross-track offset (`z0 = 30
m`) that must be nulled by the transfer, in addition to the in-plane (x–y) transfer.
This is deliberate: an in-plane-only (planar, `z ≡ 0`) scenario is a valid and useful
*special case* for verification (Section 9, "planar z=0 invariance" check) but is a
degenerate trade on its own — the in-plane CW problem is decoupled from cross-track
motion, so a planar scenario would not exercise the full 3×3 STM blocks. The chosen 3D
scenario exercises all twelve nonzero STM terms and produces a physically meaningful
Δv contribution from cross-track correction.

### 2.5 Transfer-time range to trade

Transfer time `T` (impulsive burn at t=0, second impulsive burn at t=T) is traded over:

```
T ∈ [300 s, 5100 s]
```

i.e. roughly 5 minutes to ~92% of one chief orbital period (5553.6 s). This range is
chosen deliberately to:

- start well short of the orbital period, where Δv is large (short, expensive transfers),
- extend close to (but carefully short of / excluding a guard band around) one full
  orbital period, where Δv is small (slow, cheap transfers),
- straddle the **half-period singularity** of Φrv at `t = T_period/2 ≈ 2776.8 s`
  (Section 5), which must be detected and excluded rather than blindly solved,
- avoid the coincident singularity at `t → T_period ≈ 5553.6 s` by stopping the sampled
  interval at 5100 s, with the exact exclusion/rejection rule for any sampled time that
  falls in a singularity guard band specified in Section 8.

This produces a genuine, non-trivial time-vs-Δv trade (Section 7.3): Δv decreases
steeply then more gradually as T grows, blows up near the singular times, and the
"knee" of the trade curve is visible well inside the traded interval — i.e. this is not
a monotonic or trivial trade.

---

## 3. Hill / LVLH frame definition

Origin at the chief's center of mass. Right-handed orthonormal triad, co-rotating with
the chief at the chief's orbital angular rate `n` (circular orbit ⇒ constant `n`):

| Axis | Name | Direction |
|---|---|---|
| `x` | **Radial** | From Earth's center through the chief, outward (R-bar points *inward*; `+x` is anti-R-bar / zenith) |
| `y` | **Along-track** | Direction of the chief's orbital velocity (V-bar), completing the right-handed triad |
| `z` | **Cross-track** | Along the chief's orbital angular momentum vector `h = r × v` (orbit-normal, "H-bar") |

Right-handedness: `x̂ × ŷ = ẑ`.

**Sign conventions used throughout this project:**
- `+x`: deputy is farther from Earth center than the chief (radially outward / above).
- `+y`: deputy is ahead of the chief along the direction of travel.
- `+z`: deputy is on the angular-momentum side of the orbit plane (out of plane in the
  `+h` direction).
- Negative values are the opposite sense in each case (e.g. `y = -1000 m` means the
  deputy is 1000 m *behind* the chief along-track — a trailing/V-bar-behind position).
- Relative position vector: `r = (x, y, z)`. Relative velocity: `v = (ẋ, ẏ, ż)`,
  differentiated in the rotating Hill frame (not inertial-frame velocity).
- State vector convention: `s = [x, y, z, ẋ, ẏ, ż]ᵀ` (position stacked above velocity).

This is the standard convention used by Clohessy & Wiltshire (1960), Vallado, and
Curtis; it is stated explicitly here because radial/along-track sign flips are one of
the most common sources of error in independent CW implementations.

---

## 4. Governing equations (Clohessy–Wiltshire / Hill equations)

For a deputy in the linearized vicinity of a chief on a **circular** reference orbit of
mean motion `n`, the relative motion in the Hill frame is governed by:

```
ẍ − 2n ẏ − 3n² x = 0        (radial)
ÿ + 2n ẋ           = 0        (along-track)
z̈ + n² z           = 0        (cross-track)
```

Key structural properties used later as verification targets (Section 9):

- The cross-track equation is a **decoupled, undamped harmonic oscillator** of
  angular frequency exactly `n` (the chief's mean motion) — independent of the
  in-plane (x, y) motion. This gives an exact closed-form check
  (`z(t) = z0 cos(nt) + (ż0/n) sin(nt)`).
- The in-plane (x, y) equations are coupled through the Coriolis-like terms
  `-2n ẏ` and `+2n ẋ`, and the radial equation carries the tidal/gravity-gradient
  term `-3n² x`.
- These equations are the small-relative-separation, first-order Taylor expansion of
  the two-body relative equations of motion about a circular reference orbit; they are
  **not** valid for arbitrary separations or eccentric chief orbits.

---

## 5. Clohessy–Wiltshire state-transition matrix (STM)

The CW equations are linear and time-invariant in the rotating frame (for constant
`n`), so the state at time `t` is related to the state at `t = 0` by a closed-form
state-transition matrix:

```
[ r(t) ]   [ Φrr(t)  Φrv(t) ] [ r(0) ]
[ v(t) ] = [ Φvr(t)  Φvv(t) ] [ v(0) ]
```

with `r = (x, y, z)`, `v = (ẋ, ẏ, ż)`, and each `Φ··` a 3×3 block. Using shorthand
`s = sin(nt)`, `c = cos(nt)`:

### Φrr(t)  (position due to initial position)

```
Φrr = [ 4 − 3c        0     0 ]
      [ 6(s − nt)     1     0 ]
      [ 0             0     c ]
```

### Φrv(t)  (position due to initial velocity)

```
Φrv = [  s/n              2(1 − c)/n        0   ]
      [ −2(1 − c)/n       (4s − 3nt)/n      0   ]
      [  0                0                 s/n ]
```

### Φvr(t)  (velocity due to initial position)

```
Φvr = [  3n·s          0    0  ]
      [ −6n(1 − c)     0    0  ]
      [  0             0   −n·s ]
```

### Φvv(t)  (velocity due to initial velocity)

```
Φvv = [  c        2s        0 ]
      [ −2s       4c − 3    0 ]
      [  0        0         c ]
```

These are the standard closed-form CW/Hill STM blocks (Clohessy & Wiltshire 1960;
reproduced in e.g. Curtis, *Orbital Mechanics for Engineering Students*, and Vallado,
*Fundamentals of Astrodynamics and Applications*). Each entry above will be
independently unit-tested term-by-term in M2 against:
1. the `t → 0` limit (must reduce to `Φrr → I`, `Φvv → I`, `Φrv → 0`, `Φvr → 0`),
2. direct numerical integration of the CW ODEs (Section 4) for several `(r0, v0, t)`
   combinations,
3. the decoupled exact z-solution `z(t) = z0 cos(nt) + (ż0/n) sin(nt)`.

**Composition property (to be tested in M2):** `Φ(t2) · Φ(t1) = Φ(t1 + t2)` for the
full 6×6 STM (semigroup property of a linear time-invariant flow), which follows from
`n` being constant. This is an important cross-check because the closed-form blocks
above were derived directly (not by matrix exponentiation), so this identity does not
verify itself.

---

## 6. Two-impulse boundary-value method (planned for M2/M3)

Given an initial position `r0`, a desired final position `rf`, and a transfer time `T`,
the two-impulse rendezvous solves for the initial post-burn velocity that hits `rf`
exactly at `t = T`, then reconstructs both impulses:

```
v0+        = Φrv(T)⁻¹ [ rf − Φrr(T) r0 ]
Δv1        = v0+ − v0−

vT−        = Φvr(T) r0 + Φvv(T) v0+
Δv2        = vf+ − vT−

Δv_total   = ‖Δv1‖ + ‖Δv2‖
```

where `v0−` is the deputy's actual (pre-burn) relative velocity at `t = 0`, and `vf+`
is the desired relative velocity immediately after the second burn (in this project's
scenario, `vf+ = 0` — see Section 2.3; the method supports nonzero `vf+` in general,
e.g. for a matched-velocity fly-by rather than a hold).

### 6.1 Singularity / ill-conditioning of Φrv

`Φrv(T)` must be inverted, so its conditioning is central to the method's validity.

- **Exact singularity:** `det(Φrv(t)) = 0` whenever `sin(nt) = 0`, i.e.
  `nt = kπ` for integer `k ≥ 0` — equivalently `t = k · (T_period / 2)`
  (every half chief-orbit period, including `t = 0`). This was confirmed numerically
  in Section 7.2 below (`det ≈ 1.35 × 10⁻⁶` at `t = T_period/2`, versus `~10⁸` at
  points 10–50 s away — a sharp, narrow zero, not a broad ill-conditioned region).
  Physically: at `nt = kπ` every initial velocity produces a family of reachable final
  positions collapsing onto a lower-dimensional set (in particular the cross-track
  channel `Φrv,zz = sin(nt)/n` vanishes identically, decoupling `z̈` from any control
  authority via `v0` at exactly that instant), so no unique `v0+` exists that reaches
  an arbitrary `rf`.
- **Near-singular / ill-conditioned transfer times** are a real numerical hazard even
  where `det(Φrv) ≠ 0` exactly, because `Φrv` varies smoothly through zero at each
  `t = kπ/n`: transfer times within a guard band of a singular time will produce a
  well-defined but numerically huge `v0+` (and hence `Δv`) that is not a modeling
  error, but is operationally useless and numerically fragile. These must be
  **detected and rejected**, not silently solved — solving them "successfully" would
  return an answer, but one dominated by floating-point noise in the null direction of
  a near-singular matrix, which is not a trustworthy Δv.
- **Detection strategy (for M2/M3):** before inverting `Φrv(T)`, check
  `|det(Φrv(T))|` (or, more robustly, the reciprocal condition number
  `1/cond(Φrv(T))`) against a tolerance; reject/flag any `T` where the matrix is
  singular or ill-conditioned rather than returning a numerically produced but
  physically meaningless `Δv`. Exact thresholds and the rejection rule are finalized
  in M2/M3 alongside the implementation, but the *policy* — never blindly invert — is
  fixed now.
- **Desired final velocity:** in this project's baseline scenario `vf+ = 0` (a
  rendezvous with zero final relative velocity — a stationkeeping/hold arrival, not a
  fly-by). The method as written supports an arbitrary `vf+`, and that generality is
  preserved in the formulas above so a future milestone could trade a nonzero-`vf`
  case without changing the method.

---

## 7. Hand calculation (representative transfer, T = 1800 s)

All arithmetic below was performed with double-precision decimal arithmetic (shown to
the precision needed to reproduce it independently on a calculator); intermediate
values are given so each step can be checked term-by-term. This transfer will become an
M2 regression/verification target.

### 7.1 Chief mean motion and period

```
mu = 398600.4418 km^3/s^2
Re = 6378.137 km
alt = 400.0 km
a = Re + alt = 6778.137 km

n = sqrt(mu / a^3)
  = sqrt(398600.4418 / 6778.137^3)
  = sqrt(398600.4418 / 3.114256...e11)
  = sqrt(1.27987...e-6)
  = 1.131366654e-3 rad/s

T_period = 2*pi / n = 5553.624271 s = 92.560405 min
```

### 7.2 Chosen transfer time and STM evaluation

Representative transfer time: **T = 1800 s** (well inside one orbital period; nowhere
near a Φrv singularity).

```
n*T = 1.131366654e-3 * 1800 = 2.036459976 rad   (~116.68 deg)
sin(n*T) =  0.893523773
cos(n*T) = -0.449015887
```

Evaluating the STM blocks at T = 1800 s (m, m/s units; r0, rf in meters):

```
Φrr(T) = [  5.347047661       0        0          ]
         [ -6.857617219       1        0          ]
         [  0                 0       -0.449015887 ]

Φrv(T) = [   789.773828     2561.531900        0         ]
         [ -2561.531900    -2240.904687        0         ]
         [    0                0            789.773828   ]      (units: m per (m/s), i.e. seconds)

det(Φrv(T)) = 3.784310e9   (well-conditioned; far from any singularity)
```

(Full 6-decimal values, and the `Φvr`, `Φvv` blocks, are reproduced exactly by
evaluating the closed-form expressions in Section 5 at `n*T = 2.036459976176...`; they
are not hand-typed here in full to avoid transcription error — Section 5's formulas are
the source of truth, and M2's unit tests evaluate them directly.)

### 7.3 Burn reconstruction

Using `r0 = (0, -1000, 30) m`, `v0 = (0,0,0) m/s`, `rf = (0, -50, 0) m`, `vf = (0,0,0) m/s`:

```
rf − Φrr(T) r0 = (0, -50, 0) − Φrr(T)·(0,-1000,30)

Φrr(T)·r0 = ( 5.347047661*0 + 0*(-1000) + 0*30 ,
             -6.857617219*0 + 1*(-1000) + 0*30 ,
              0*0 + 0*(-1000) + (-0.449015887)*30 )
          = ( 0, -1000, -13.470477 )

rf − Φrr(T) r0 = (0 - 0, -50 - (-1000), 0 - (-13.470477))
               = (0, 950, 13.470477)         [m]

v0+ = Φrv(T)^-1 * (0, 950, 13.470477)
    ≈ (-0.5079, 0.1566, 0.01706) m/s

Δv1 = v0+ - v0- = v0+ - (0,0,0) = (-0.5079, 0.1566, 0.01706) m/s
|Δv1| = sqrt(0.5079^2 + 0.1566^2 + 0.01706^2) = 0.5317 m/s = 531.7 mm/s

vT- = Φvr(T) r0 + Φvv(T) v0+ ≈ (0.5079, 0.1566, -0.03799) m/s

Δv2 = vf+ - vT- = (0,0,0) - (0.5079, 0.1566, -0.03799)
    = (-0.5079, -0.1566, 0.03799) m/s
|Δv2| = sqrt(0.5079^2 + 0.1566^2 + 0.03799^2) = 0.5328 m/s = 532.8 mm/s

Δv_total = |Δv1| + |Δv2| ≈ 0.5317 + 0.5328 = 1.0645 m/s ≈ 1064.5 mm/s
```

**Terminal closure check:** `Φrr(T) r0 + Φrv(T) v0+` must equal `rf = (0,-50,0)` exactly
(to numerical precision). Verified: residual `< 1e-13 m` component-wise — this is the
single most important hand-check because it validates that `v0+` was solved correctly
from the linear system, independent of any Δv arithmetic.

These values (`|Δv1| ≈ 0.5317 m/s`, `|Δv2| ≈ 0.5328 m/s`, `Δv_total ≈ 1.0645 m/s`, at
`T = 1800 s`) are the **M2 regression target** for the two-impulse solver.

---

## 8. Planned Δv / transfer-time trade (M3 — not implemented yet)

- **Transfer-time interval:** `T ∈ [300 s, 5100 s]` (Section 2.5).
- **Sampling strategy:** dense uniform sampling (e.g. 1 s or a few-second step) across
  the interval, sufficient to resolve the shape of the trade curve and to visually
  locate — without solving through — the singular time at `T_period/2 ≈ 2776.8 s`.
  Exact step size is an M3 implementation choice, not fixed here.
- **Outputs to report per sampled T:** `|Δv1|`, `|Δv2|`, `Δv_total`, and
  `det(Φrv(T))` (or `1/cond(Φrv(T))`) so singular/near-singular points are visible
  alongside the Δv curve rather than hidden.
- **Singular/invalid times to reject:** any sampled `T` for which `Φrv(T)` is singular
  or numerically ill-conditioned per the detection rule in Section 6.1 — concretely,
  `T` within a guard band of `k · T_period/2` for integer `k` (this scenario's traded
  interval only contains `k = 1`, at `T ≈ 2776.8 s`, since the interval starts above
  `T = 0` and ends before `T_period ≈ 5553.6 s`). Rejected points are reported as
  invalid (e.g. `NaN`/flagged), not silently dropped or silently solved.
- Representative points computed for this document (not the final M3 sampling, just a
  sanity preview using the same formulas as Section 7): Δv_total falls from several
  m/s at `T = 300 s` to a broad minimum around `T ≈ 4000–4500 s`, then rises again
  approaching `T_period`, with a sharp spike at the `T ≈ 2776.8 s` singularity — i.e. a
  genuine, non-monotonic, non-trivial trade curve, not a flat or trivially monotonic
  one.
- No optimization (e.g. minimum-Δv search) is performed in M1 or planned before M3;
  M3 only tabulates/plots the trade.

---

## 9. Verification plan (for M2 onward)

At minimum, the following checks will be implemented as automated tests
(`pytest -W error`) before any trade results are trusted:

1. **t = 0 state-transition identity** — `Φrr(0) = I`, `Φvv(0) = I`,
   `Φrv(0) = 0`, `Φvr(0) = 0` (evaluated at a small `t → 0` limit and/or via the
   analytic limits of each closed-form term).
2. **CW differential-equation residual check** — substitute the closed-form STM
   solution `r(t) = Φrr(t) r0 + Φrv(t) v0` back into
   `ẍ − 2nẏ − 3n²x = 0`, `ÿ + 2nẋ = 0`, `z̈ + n²z = 0` (via numerical or symbolic
   differentiation) and confirm the residual is ~0 for multiple `(r0, v0)`.
3. **STM composition / propagation consistency** — `Φ(t1+t2) = Φ(t2)·Φ(t1)` for
   several `(t1, t2)` pairs, both in and out of the traded interval.
4. **Planar z=0 invariance** — with `z0 = 0`, `ż0 = 0`, confirm `z(t) ≡ 0` for all `t`
   (the in-plane and cross-track channels do not leak into each other).
5. **Exact harmonic z-motion check** — for nonzero `z0`/`ż0` (as in this project's own
   scenario), confirm `z(t) = z0 cos(nt) + (ż0/n) sin(nt)` matches the STM-based
   propagation to numerical precision.
6. **Terminal position closure** — after solving `v0+` from the two-impulse method,
   confirm `Φrr(T) r0 + Φrv(T) v0+ = rf` to numerical precision (as done by hand in
   Section 7.3).
7. **Burn reconstruction check** — confirm `Δv1 = v0+ − v0−` and
   `Δv2 = vf+ − (Φvr(T) r0 + Φvv(T) v0+)` are self-consistent with independently
   propagating the post-burn state forward and checking it lands exactly on
   `(rf, vf+)`.
8. **Zero-separation / trivial limiting case** — `r0 = rf = 0`, `v0 = vf = 0` must
   give `Δv_total = 0` exactly (no phantom burns for a no-op transfer).
9. **Near-singular Φrv detection** — confirm the solver flags/rejects transfer times
   at and near `T = k·T_period/2`, and confirm `det(Φrv(T))` computed by the code
   matches the closed-form zero-crossing behavior verified numerically in Section 6.1.
10. **Independent numerical ODE propagation cross-check** — integrate the raw CW ODEs
    (Section 4) with a numerical integrator (e.g. `scipy.integrate.solve_ivp`) for the
    same `(r0, v0, T)` and confirm agreement with the closed-form STM propagation to
    integrator tolerance — this is independent of the hand-derived STM formulas
    themselves and guards against a shared algebra error.
11. **Symmetry / time-reversal check** — propagating forward by `T` then backward by
    `T` (i.e. applying `Φ(T)` then `Φ(-T)`, or equivalently `Φ(T)^-1`) returns the
    original state, confirming `Φ(-t) = Φ(t)^-1`.

Each check will be an explicit, independently reproducible `pytest` test in M2, not a
figure or a visual inspection. Figures (M3) are treated as diagnostics only, secondary
to these numerical checks — per project ground rules, a plot is never itself
verification.

---

## 10. Limitations (restated explicitly)

- Valid only for **linearized relative dynamics** — first-order in separation over
  chief orbit radius; not valid once separation grows to be a non-negligible fraction
  of `a`.
- Valid only for a **circular chief orbit** — CW equations do not apply as-is to an
  eccentric chief (would require the time-varying Tschauner–Hempel formulation, out of
  scope).
- Assumes **small separation relative to chief orbital radius** (`|r| << a`); this
  scenario's largest separation is 1000 m against `a ≈ 6778 km`, i.e. `|r|/a ≈ 1.5e-4`,
  comfortably small.
- **Impulsive burns only** — no finite burn duration, no thrust curve, no attitude/
  pointing constraints on the burn.
- **No J2, no drag, no other perturbations** — pure two-body linearized dynamics.
- **No nonlinear two-body verification yet** — item 10 of the verification plan adds a
  numerical ODE cross-check of the *linear* CW equations, not a nonlinear two-body
  propagator comparison; that would be a separate, later validation not planned here.
- **No keep-out sphere / approach corridor / line-of-sight constraints** — the two
  burns are solved purely as a boundary-value problem; nothing in this project checks
  whether the straight-line/relative trajectory between `r0` and `rf` avoids the chief
  or respects an approach corridor.
- **No finite-thrust or actuator limits** — `Δv` values are idealized instantaneous
  velocity changes with no minimum impulse bit, thruster saturation, or duty-cycle
  constraint.

---

## 11. File/package scaffolding created in M1

No CW solver code is implemented yet. This milestone creates only:

```
rendezvous-proximity-ops-trade/
├── DESIGN.md                    (this file)
├── README.md
├── .gitignore
├── pyproject.toml               (packaging + pytest config only)
├── src/
│   └── rendezvous_cw/
│       └── __init__.py          (package docstring/version only — no CW equations)
├── tests/
│   └── test_placeholder.py      (imports the empty package; passes trivially)
├── scripts/
│   └── .gitkeep
└── figures/
    └── .gitkeep
```

The CW STM implementation, the two-impulse solver, the approach-trajectory plot, and
the Δv-vs-time trade table/script are all **Milestone 2/3 work** and are not present in
this commit.
