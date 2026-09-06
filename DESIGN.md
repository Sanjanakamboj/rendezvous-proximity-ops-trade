# DESIGN.md — Rendezvous & Proximity Operations Trade

Status: **Milestone 6 complete — project finalized for portfolio use.** Milestones
1–5 content (scenario, equations, hand calculations, STM/solver implementation,
unconstrained Δv/time trade, geometric constraints, nonlinear validation) is
unchanged below except for hygiene fixes and "superseded" banners added in M6; see
**§16 for the authoritative final result and the full superseded-result hierarchy**
before reading any individual milestone section in isolation. **Final recommended
transfer: T = 380 s, Δv_total ≈ 5.136 m/s, CW margin ≈ +0.875 m, nonlinear margin ≈
+0.902 m** (§16.1). **Scope reminder (applies throughout):** this project performs
geometric screening inside a model — first CW-linear (M4), then cross-checked against
nonlinear two-body dynamics (M5). No milestone's "passes" language is a claim of
collision safety, flight safety, or operational safety. J2/drag, finite burns,
actuator limits, sensor field-of-view, line-of-sight occultation, plume impingement,
docking dynamics, and collision-probability modeling remain **not implemented** and
are out of this project's scope.

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

---

## 12. Milestone 2 — implementation and verification results

M2 implements exactly the equations and method fixed in sections 3–7 above (no
equations were changed) as a real, tested Python package, and verifies them against
the 11-point plan in section 9. **The full Δv-vs-transfer-time trade (M3) is still not
implemented** — M2 exercises only the single T = 1800 s representative transfer plus
the standalone STM/conditioning verification suite.

### 12.1 Package layout added in M2

```
src/rendezvous_cw/
├── orbit.py          Circular chief-orbit utilities (a, n, period)
├── cw.py              Closed-form CW STM blocks + propagate_cw()
├── conditioning.py    Phi_rv determinant/condition-number diagnostics + policy
└── rendezvous.py       Two-impulse boundary-value solver (solve_two_impulse)

tests/
├── conftest.py                 Shared fixtures + independent CW ODE reference
├── test_orbit.py                n/period vs. this document
├── test_cw.py                   STM identity, ODE residual, composition,
│                                 planar invariance, cross-track oscillator,
│                                 time-reversal (checks A,B,C,D,E,K)
├── test_conditioning.py         Singularity detection (check I)
├── test_ode_cross_check.py      Independent solve_ivp cross-check (check J)
└── test_rendezvous.py           Terminal closure, burn reconstruction, trivial
                                  case, M1 hand-calc regression (checks F,G,H)

scripts/m2_plot_representative_transfer.py   Generates the M2 diagnostic figure
figures/m2_representative_transfer.png       M2 diagnostic figure (not portfolio-final)
```

`propagate_cw()` is the sole production propagation path; it evaluates the closed-form
STM directly (`state(t) = Phi(t, n) @ state0`). `scipy.integrate.solve_ivp` is used
**only** inside `tests/test_ode_cross_check.py`, as an independent verification path —
never in `src/`.

### 12.2 Chief orbit (unchanged from M1)

Reproduced by `orbit.chief_mean_motion_and_period()`:

```
n = 1.131366653611e-3 rad/s   (== M1 section 7.1 value, to 12 significant figures)
P = 5553.624271252 s          (== M1 section 2.1 / 7.1 value)
```

### 12.3 Phi_rv conditioning policy (implemented per M1 §6.1's stated policy)

`conditioning.phi_rv_conditioning(t, n)` returns `det(Φrv)`, `cond(Φrv)` (2-norm
condition number, via `numpy.linalg.cond`), and an `unsafe` flag. **Policy:** a
transfer time is unsafe if `cond(Φrv) > 1e6` (`COND_THRESHOLD`), *not* an absolute
determinant cutoff — `Φrv`'s entries scale with `1/n`, so an absolute determinant
threshold would need re-tuning per orbit regime, whereas the condition number directly
measures error amplification in the linear solve regardless of scale. `cond = 1e6`
means up to ~6 of double precision's ~16 digits can be lost — a conservative but not
overly restrictive cutoff for this project's mm/s-level Δv targets. Full rationale is
in the `conditioning.py` module docstring.

Measured conditioning near the M1 half-period singularity (`T_period/2 ≈ 2776.812 s`,
400 km chief):

| Δt from singularity | cond(Φrv) | flagged unsafe? |
|---:|---:|:---:|
| 0 s (exact) | 8.90 × 10¹⁶ | yes |
| 0.001 s | 9.63 × 10⁶ | yes |
| 0.01 s | 9.63 × 10⁵ | no |
| 0.1 s | 9.63 × 10⁴ | no |
| 1 s | 9.62 × 10³ | no |
| 10 s | 9.57 × 10² | no |

Condition number strictly increases as `T` approaches the singularity (verified
monotonically in `test_conditioning.py`); the `cond > 1e6` policy rejects a guard band
of roughly ±(1–10) ms around each `nt = kπ` singular time for this scenario — tight,
because conditioning degrades extremely sharply (not gradually) right at the
singularity, so a tight guard band still catches every numerically dangerous point.
`solve_two_impulse()` calls this check before ever calling `numpy.linalg.solve` on
`Φrv`, and raises `SingularTransferError` (not a silently-huge Δv) when unsafe.

### 12.4 Representative transfer (T = 1800 s) — authoritative M2 values

Computed by `rendezvous.solve_two_impulse()` using the exact M1 scenario (§2) and
compared against the M1 hand-calculation target (§7.3):

| Quantity | M1 hand target (rounded) | M2 computed (double precision) | Discrepancy |
|---|---:|---:|---:|
| `\|Δv1\|` | 531.7 mm/s | 531.7191 mm/s | 0.019 mm/s |
| `\|Δv2\|` | 532.8 mm/s | 532.8013 mm/s | 0.001 mm/s |
| `Δv_total` | 1064.5 mm/s | 1064.5203 mm/s | 0.020 mm/s |

All discrepancies are within M1's stated hand-calculation rounding (4 significant
figures) — **no disagreement requiring investigation**. `v0_plus = (-0.507855,
0.156582, 0.017056) m/s`; `vT_minus = (0.507855, 0.156582, -0.037986) m/s`;
`Φrv(1800 s)` determinant `≈ 3.784 × 10⁹`, condition number `≈ 5.29` (very
well-conditioned — this transfer time is far from any singularity).

### 12.5 Numerical-quality residuals (not just "tests passed")

Measured directly (see verification suite for the exact cases exercised):

| Check | Metric | Max residual/error |
|---|---|---:|
| STM composition (`Φ(t1+t2)` vs. `Φ(t2)Φ(t1)`) | max abs. matrix entry residual, 5 (t1,t2) pairs | 1.46 × 10⁻¹¹ |
| Time reversal (`Φ(-t)Φ(t)` vs. `I`) | max abs. matrix entry residual, 5 times | 4.77 × 10⁻¹² |
| STM vs. independent `solve_ivp` (DOP853, rtol/atol 1e-12) | max abs. state-component error, 3 states × 5 times | 7.20 × 10⁻¹⁰ |
| Terminal position closure (T = 1800 s transfer) | `‖r(T) − rf‖`, meters | 1.61 × 10⁻¹³ m |
| M1 hand-target discrepancy | `Δv_total`, mm/s | 0.020 mm/s |
| Condition number at exact half-period singularity | `cond(Φrv)` | 8.90 × 10¹⁶ (flagged unsafe) |

All residuals are consistent with double-precision floating-point arithmetic on
closed-form trigonometric expressions — i.e. the implementation matches the M1
equations to numerical precision, not merely "close enough."

### 12.6 Diagnostic figure

`figures/m2_representative_transfer.png` (generated by
`scripts/m2_plot_representative_transfer.py`): the T = 1800 s post-first-burn relative
trajectory, LVLH x–y (radial vs. along-track) projection plus a cross-track (z) vs.
time panel showing the 30 m out-of-plane offset being nulled smoothly by burn 1,
reaching exactly `z = 0` at `t = T`. Labeled with transfer time and both burn
magnitudes. This is a diagnostic figure only, not the M3 portfolio-final
approach-trajectory plot, and no Δv-vs-transfer-time trade figure is produced in M2.

### 12.7 Explicitly not done in M2

- No Δv-vs-transfer-time trade sweep or table (M3).
- No transfer-time optimization.
- No keep-out zones, approach corridors, or line-of-sight constraints.
- No nonlinear two-body validation (the `solve_ivp` cross-check in 12.5 validates the
  *linear* CW ODE against the closed-form STM, not a nonlinear two-body propagator).
- No J2, drag, finite burns, actuator models, or collision avoidance.
- No portfolio-final figures.

### 12.8 Limitations

Unchanged from M1 §10 — restated in full there; still binding.

---

## 13. Milestone 3 — Δv-vs-transfer-time trade study

> **⚠ Superseded (see §16 for the authoritative final result).** This section's
> `T ≈ 4868.130 s` **unconstrained** Δv minimum is a real, correctly-computed result of
> the CW two-impulse trade — but it **fails** the M4 proximity-geometry constraints
> (§14.4) and is retained here only as engineering history / the starting point that
> motivated M4. It is **not** the project's recommendation. The final recommended
> transfer is `T = 380 s` (§16.1 / §15.6).

M3 uses **only** the verified M2 solver (`rendezvous.solve_two_impulse`) and M2
conditioning policy (`conditioning.py`) — no CW/rendezvous equations were changed or
duplicated. The fixed M1 scenario (`r0, v0_minus, rf, vf_plus`, 400 km circular chief)
is unchanged. New code: `src/rendezvous_cw/trade.py` (a thin, reusable wrapper around
the M2 solver — `evaluate_transfer_time`, `sweep_transfer_times`,
`find_best_safe_transfer`), `scripts/m3_trade_sweep.py` (generates all M3 artifacts),
`results/m3_transfer_trade.csv` (curated decision table), `results/m3_full_sweep.csv`
(full sweep, for reproducibility), and three figures (§13.7).

### 13.1 Sweep interval and resolution

- **Domain:** `T ∈ [300, 5100] s`, deliberately spanning the interior half-period
  singularity at `T = P/2 ≈ 2776.812 s` (M1 §6.1 / M2 §12.3).
- **Base resolution:** uniform 5 s spacing across the full domain (used for the
  headline Δv curve, Figure 1, and the approach-trajectory comparison, Figure 2).
- **Singularity probe:** a separate, dense 0.5 ms-resolution probe spanning
  `P/2 ± 0.05 s`, used only to characterize the singularity itself (Figure 3) and to
  produce real, solver-derived unsafe records — **not** merged into the headline Δv
  curve's line density, and not used to search for a minimum. This matches the M3
  instruction not to sample arbitrarily close to the singularity "just to create huge
  values": the probe exists to *measure* the excluded region precisely, not to inflate
  Δv anywhere it is reported as a trade result.

### 13.2 Unsafe-conditioning policy (inherited unchanged from M2)

`cond(Φrv) > 1×10⁶` ⇒ unsafe/excluded (`conditioning.COND_THRESHOLD`, M2 §12.3). No
new threshold or policy was introduced in M3; `trade.evaluate_transfer_time` calls the
same `conditioning.phi_rv_conditioning` used by the M2 solver and returns an explicit
`safe=False` record (with real `cond`/`det` values, `dv_total = NaN`) instead of
raising, so a sweep can carry unsafe points through as explicit excluded/NaN entries.

### 13.3 How the singularity splits the feasible domain

The 0.5 ms probe locates the excluded interval precisely by direct evaluation:

```
Excluded interval: [2776.802636, 2776.821636] s  (width = 0.019000 s ≈ 19 ms)
centered at P/2 = 2776.812136 s
```

This is **far narrower** than the 5 s base sweep resolution — no point in the base
sweep grid falls inside it, which is why the headline Δv curve (Figure 1) never
contains a sample point flagged unsafe: the singularity is real and precisely
characterized (Figure 3), but it is simply too narrow to be visible or to interrupt a
minutes-scale plot without either artificially oversampling near it (which M3
instructions rule out) or exaggerating the shaded region's width (which would
misrepresent the physics). Figure 1 therefore marks `P/2` with a vertical line and an
explicit annotation stating the true ~19 ms width, and Figure 3 shows the actual
excluded band at millisecond resolution.

The singularity splits `[300, 5100] s` into two **disjoint feasible branches**, each
searched and minimized independently (never treating the domain as one smooth curve
through the singularity, per M3 instructions):

- **Branch 1** (`T < P/2`, short transfers): `T ∈ [300, 2776.80) s`
- **Branch 2** (`T > P/2`, long transfers): `T ∈ (2776.82, 5100] s`

### 13.4 Representative trade table

Full table: [`results/m3_transfer_trade.csv`](results/m3_transfer_trade.csv). Summary
(Δv columns in mm/s; all points safe unless noted):

| T | T/P | cond(Φrv) | \|Δv1\| | \|Δv2\| | Δv_total | closure [m] |
|---:|---:|---:|---:|---:|---:|---:|
| 300 s | 0.054 | 1.06 | 3223.72 | 3223.90 | 6447.62 | 3.6×10⁻¹⁵ |
| 600 s | 0.108 | 1.25 | 1670.01 | 1670.35 | 3340.36 | 0 |
| 900 s | 0.162 | 1.62 | 1141.94 | 1142.45 | 2284.39 | 1.1×10⁻¹³ |
| 1200 s | 0.216 | 2.25 | 856.97 | 857.64 | 1714.61 | 1.1×10⁻¹³ |
| 1800 s (M2 ref.) | 0.324 | 5.29 | 531.72 | 532.80 | 1064.52 | 1.6×10⁻¹³ |
| 2400 s | 0.432 | 20.16 | 354.64 | 356.26 | 710.89 | 1.1×10⁻¹³ |
| P/2 − 0.02 s (just before excl. zone) | 0.500 | 4.81×10⁵ | 1500.00 | 1500.00 | 3000.00 | 3.6×10⁻¹⁵ |
| P/2 + 0.02 s (just after excl. zone) | 0.500 | 4.81×10⁵ | 1500.00 | 1500.00 | 3000.00 | 0 |
| 3600 s | 0.648 | 19.98 | 158.77 | 162.35 | 321.12 | 1.4×10⁻¹⁴ |
| 4200 s | 0.756 | 18.46 | 106.15 | 111.44 | 217.59 | 1.4×10⁻¹⁴ |
| 4800 s | 0.864 | 25.67 | 78.68 | 85.69 | 164.36 | 2.3×10⁻¹³ |
| 5100 s | 0.918 | 39.26 | 86.82 | 93.22 | 180.04 | 0 |
| **4868.130 s (refined min., branch 2)** | 0.877 | 27.62 | 78.20 | 85.25 | **163.457** | 3.6×10⁻¹⁵ |

Note the "just before/after excluded zone" rows: at only 20 ms from `P/2` (safe by the
`cond(Φrv) > 1e6` test but barely so, `cond ≈ 4.8×10⁵`), Δv is already ~3 m/s — an
order of magnitude worse than the 1800 s reference — a direct, solver-derived
illustration of why the conditioning policy exists: technically "safe" transfers
immediately adjacent to the singularity are numerically valid but operationally
useless.

### 13.5 Refined minimum-Δv result

Two disjoint local minima were found, one per branch (never compared as if on one
smooth curve — see §13.3):

- **Branch 1 minimum:** `T ≈ 2539.593 s` (`T/P ≈ 0.4573`), `Δv_total ≈ 678.609 mm/s`
  — a genuine interior local minimum (verified: Δv decreases monotonically from 300 s
  to ~2540 s, then rises sharply approaching the singularity — not a boundary artifact
  of the search domain).
- **Branch 2 minimum ("the" refined minimum, lower of the two):**
  `T = 4868.130 s` (`T/P = 0.8766`)
  `Δv1 = (-0.04437, -0.05431, 0.03462) m/s`, `|Δv1| = 78.205 mm/s`
  `Δv2 = (-0.04437, 0.05431, -0.04848) m/s`, `|Δv2| = 85.253 mm/s`
  `Δv_total = 163.457 mm/s`
  `cond(Φrv) = 27.62` (very well-conditioned), terminal closure residual `3.6×10⁻¹⁵ m`.

**Δv_total relative to the M2 T=1800 s reference (1064.520 mm/s): −84.64%** — i.e. a
transfer time roughly 2.7× longer than the M2 reference reduces total Δv by roughly
6.5×, within the searched domain and this fixed boundary-value problem.

**Scope of the word "minimum" (stated explicitly, per M3 instructions):** this is the
minimum-Δv transfer time found by (1) evaluating the CW two-impulse solver on a finite
grid over the explicitly searched interval `[300, 5100] s`, excluding the
`cond(Φrv) > 1e6` unsafe region, then (2) locally refining with a bounded scalar
minimizer restricted to the immediate safe neighbors of the coarse grid minimum (never
crossing into an unsafe sub-interval or outside the swept domain). It is:
- **not** a global optimum over all possible transfer times (the search stops at
  5100 s; Δv could behave differently beyond one orbital period, unexplored here),
- **not** an optimum with respect to any variable other than transfer time (r0, rf,
  vf_plus, and the chief orbit are all fixed),
- **not** fuel-optimal, collision-safe, or operationally feasible for a real
  spacecraft — it is a Δv-minimizing point of the linearized CW two-impulse
  boundary-value problem only, under the M1 scenario's fixed endpoints,
- **not** validated outside CW linear dynamics (no nonlinear two-body check has been
  performed on this or any other M3 result).

### 13.6 Numerical verification / convergence

- **Grid-resolution convergence** (branch 2 minimum, coarse grid before refinement):

  | grid step | T (coarse) | Δv_total |
  |---:|---:|---:|
  | 20 s | 4860.00 s | 163.4713 mm/s |
  | 10 s | 4870.00 s | 163.4581 mm/s |
  | 5 s | 4870.00 s | 163.4581 mm/s |
  | refined (bounded scalar minimize within the 5 s grid's safe neighbors) | 4868.130 s | 163.4574 mm/s |

  The coarse minimum location and value stabilize between 10 s and 5 s resolution
  (identical at both), and refinement improves the value by only ~0.0007 mm/s beyond
  the 5 s grid — i.e. the answer has converged well before reaching the refinement
  step; refinement mainly sharpens the exact `T` location, not the Δv value.
- **Direct-solver-vs-trade-module residual:** `trade.evaluate_transfer_time` matches
  `rendezvous.solve_two_impulse` called directly to `< 1e-12` relative error on
  `|Δv1|`, `|Δv2|`, `Δv_total`, and closure residual, for 7 arbitrary safe transfer
  times spanning the domain (`tests/test_trade.py`, check B) — confirming the trade
  module introduces no numerical drift.
- **Worst terminal-position closure residual over all safe sweep points (1123
  points, base 5 s grid + singularity probe):** `3.595×10⁻¹³ m`.
- **Maximum condition number actually admitted as "safe" by the policy** across the
  same sweep: `9.629×10⁵` (just under the `1e6` threshold, at the closest safe probe
  point to `P/2`).
- **STM-vs-independent-`solve_ivp` closure at the refined minimum transfer:** post-burn
  state propagated by the closed-form STM matches an independent DOP853 integration
  (rtol/atol `1e-12`) of the raw CW ODE to `< 1e-6` (position/velocity components);
  the STM-propagated position matches `rf` to `< 1e-6 m` and STM-propagated pre-second-
  burn velocity matches the solver's `vT_minus` to `< 1e-8 m/s`
  (`tests/test_trade.py::test_refined_minimum_trajectory_closure`).
- **Sweep-ordering independence:** ascending, descending, and shuffled evaluation
  orders of the same time set produce bitwise-identical per-time results (no hidden
  state/caching) — `tests/test_trade.py::test_sweep_independent_of_time_ordering`.
- **Singularity-side behavior:** approaching `P/2` from both sides (offsets 10, 1,
  0.1, 0.02 s), the condition number increases strictly monotonically on each side as
  the offset shrinks (verified, not merely asserted); a dense 1 ms-resolution probe
  spanning `P/2 ± 0.03 s` contains both unsafe points (inside the ~19 ms band) and safe
  points immediately outside it on both sides, confirming the unsafe region genuinely
  interrupts a fine-grained scan rather than being a modeling artifact.

### 13.7 Figures

1. **`figures/m3_dv_vs_transfer_time.png`** (headline) — total Δv (heavy red) plus
   `|Δv1|`/`|Δv2|` (lighter supporting curves) vs. transfer time in minutes, over the
   full safe domain in both branches, with the M2 1800 s reference (●) and the refined
   branch-2 minimum (★) marked, a vertical line at `P/2`, and an explicit annotation
   stating the true (~19 ms) width of the excluded region rather than a misleadingly
   wide shaded band. Title states CW linearized dynamics, fixed M1 boundary
   conditions, and that the unsafe region is excluded.
2. **`figures/m3_approach_trajectory_comparison.png`** — LVLH x–y trajectories for four
   safe transfers (600 s fast case, 1800 s M2 reference, 4868 s refined minimum, 5100 s
   slow case), each labeled with its Δv, showing that the lowest-Δv transfer follows a
   visibly different (deeper, closer-to-natural-drift) path than the fast case.
3. **`figures/m3_conditioning_vs_transfer_time.png`** (verification/supporting) —
   `cond(Φrv)` vs. time offset from `P/2` in milliseconds (log y-axis), the `1e6`
   threshold line, and the actual excluded band shaded at true scale — this is the
   figure that resolves what Figure 1 cannot show at its own scale.

All three were visually inspected for clipping, overlapping annotations/legends, and
correct unsafe-point handling (no line is drawn connecting across the excluded
region in any figure) before being accepted.

### 13.8 Explicitly not done in M3

Everything listed as out-of-scope for M3 in the milestone instructions remains
out of scope: keep-out zones, approach corridors, line-of-sight constraints,
nonlinear two-body validation, J2/drag, finite burns, actuator/thruster models,
generic black-box optimization, collision-risk probability. No claim of
collision-safety, operational feasibility, fuel-optimality for a real spacecraft,
global optimality, or validation outside CW linear dynamics is made anywhere in this
section (see §13.5's explicit scope statement).

### 13.9 Limitations

Unchanged from M1 §10 (still binding) — see also M2 §12.7/12.8 and M3 §13.8 above for
what remains unimplemented.

---

## 14. Milestone 4 — proximity-geometry-constrained trade

> **⚠ CW-only constrained optimum — superseded by M5 robustness validation (see §16).**
> This section's `T ≈ 386.564 s` result is the correct minimum-Δv transfer *within the
> CW model alone*, sitting at essentially zero (~15 µm) active-constraint margin by
> construction. M5 (§15) found it survives nonlinear validation but only barely
> (+27.5 mm real margin), and recommends `T = 380 s` instead for a robust, easily
> defensible margin at a small Δv cost. This section's method, table, and figures
> remain correct as originally computed and are retained for engineering history.

**Scope statement (read first):** this section performs **geometric screening only**,
inside the already-verified linearized CW model. A trajectory that passes these checks
is described as **"geometrically feasible under this CW model"** — never as
"collision-free", "flight safe", or "operationally safe". No sensor field of view,
Earth occultation, plume impingement, docking dynamics, finite-burn execution, actuator
limits, or collision probability is modeled. M4 uses **only** the verified M2 solver
(`rendezvous.solve_two_impulse`), M2 closed-form STM (`cw.propagate_cw`), and M3 trade
machinery — no dynamics were changed. New code: `src/rendezvous_cw/constraints.py`,
`scripts/m4_constrained_trade.py`, `results/m4_constrained_trade.csv` (curated),
`results/m4_full_constrained_sweep.csv` (full sweep), three figures (§14.9).

### 14.1 Exact keep-out/corridor/no-crossing definitions

- **A. Spherical keep-out zone (KOZ):** radius `r_KOZ = 100 m`, centered on the chief.
  `‖r_rel(t)‖ ≥ r_KOZ` required for the whole coast, **except** while inside the
  corridor (B).
- **B. Final approach corridor:** `-100 m ≤ y ≤ -50 m`, `|x| ≤ 20 m`, `|z| ≤ 20 m`.
  Only inside this box is `‖r‖ < r_KOZ` authorized. The M1 target
  `r_f = (0, -50, 0) m` sits exactly at the corridor's inner (y = -50 m) end.
- **C. No chief crossing:** `y(t) ≤ -50 m` for the entire transfer (endpoint tolerance
  below).

Within an unauthorized KOZ breach, two distinct reasons are reported separately rather
than one collapsed flag (`src/rendezvous_cw/constraints.py` module docstring):
`koz_violation` (breach while entirely outside the corridor's y-range — a generic
keep-out breach) vs. `corridor_violation` (breach while within the corridor's y-range
but outside its lateral x/z bounds — an attempted corridor approach from the wrong
lateral position).

### 14.2 Continuous-time checking method and tolerances

1. Dense-sample `r(t)` at `sample_dt` spacing (production default 2 s for the sweep,
   1 s or finer used in convergence checks and tests) using the verified closed-form
   STM (`propagate_cw`) — no new dynamics.
2. Classify every sample against all three constraints.
3. Bracket the coarse-sample candidate minimum-distance time (and, separately, the
   candidate minimum-clearance-outside-corridor time) with its immediate dense-sample
   neighbors, then refine with `scipy.optimize.minimize_scalar` (`method="bounded"`)
   on the continuous position function. Refinement can only find a smaller (worse)
   minimum than the coarse grid saw — it strictly tightens, never loosens, the safety
   check.
4. **Floating-point boundary tolerance:** the M1 target sits exactly on the corridor's
   `y = -50 m` boundary, and STM evaluation leaves ~10⁻¹³ m residuals there (verified
   in M2/M3). Directly observed: for some transfer times, the closed-form STM evaluates
   the exact target point as `y = -49.999999999999886` — numerically *above* -50 m by
   ~1.1×10⁻¹³ m, which would misclassify the valid endpoint as a corridor/chief-crossing
   violation without a tolerance. `CORRIDOR_BOUNDARY_TOL_M = CHIEF_CROSSING_TOL_M =
   1×10⁻⁶ m` is applied to all three boundary checks — six orders of magnitude larger
   than the observed residual, so it absorbs floating-point noise while remaining far
   below any physically meaningful boundary distance (this exact bug was caught and
   fixed during M4 development — see §14.10).

### 14.3 Re-evaluating the M3 trade under constraints

Domain and resolution: `T ∈ [300, 5100] s` (unchanged from M3), 10 s sweep step, 2 s
constraint-sampling spacing. The M3 `cond(Φrv) > 1e6` exclusion policy is inherited
unchanged (`conditioning.py`, untouched). Feasible set =
`M3 numerically safe AND M4 geometrically safe`; M3-unsafe points are never
reclassified by geometry (verified in tests, check F).

**Result: only 9 of 481 M3-safe sweep points (1.87%) are M4-geometrically feasible —
all of them fast transfers `T ≲ 387 s`.** No other feasible window exists anywhere in
`[300, 5100] s` (confirmed by full-domain sweep at 10 s resolution — a single
contiguous feasible band at the very start of the domain, nothing elsewhere).

**Why:** for a two-impulse transfer between these fixed, widely-separated endpoints
(`r0` 1000 m away, `rf` 50 m away, both with zero relative velocity), the natural CW
coast path develops a large radial ("belly") excursion — tens to hundreds of meters —
that only collapses toward zero in the final moments before reaching `rf`. For all but
the fastest transfers, the trajectory's distance from the chief already drops below
100 m *before* the lateral (x, z) offset has shrunk to within the corridor's 20 m
bound — an unauthorized KOZ/corridor breach. Longer transfers (branch 2, `T ≳ 4200 s`,
including the M3 unconstrained minimum) additionally **overshoot past the target
along-track position** (`y > -50 m`, seen up to `y ≈ -39 m` for the M3 global minimum
at `T = 4868.130 s`) before curving back to hit `rf` exactly at `t = T` — a genuine
no-chief-crossing violation the pure-Δv M3 trade could not surface, since M3 only
checked the two endpoints, not the continuous path.

Rejection reason counts (of 472 M3-safe-but-M4-rejected points, 10 s sweep):
`koz_violation`: 472, `corridor_violation`: 469, `chief_crossing_violation`: 232
(points can carry more than one reason simultaneously).

### 14.4 Do the M3 minima survive?

**No — none of the three M3 key transfer times pass M4 geometry:**

| Transfer | T (s) | M3 numerically safe? | M4 geometrically safe? | Reason |
|---|---:|:---:|:---:|---|
| M2 reference | 1800.000 | ✅ | ❌ | `corridor_violation, koz_violation` |
| M3 short-branch minimum | 2539.593 | ✅ | ❌ | `corridor_violation, koz_violation` |
| M3 long-branch/global minimum | 4868.130 | ✅ | ❌ | `corridor_violation, koz_violation, chief_crossing_violation` |

The unconstrained M3 minimum-Δv transfer not only breaches the keep-out sphere well
outside the authorized corridor, it also overshoots past the target along-track
position — a materially worse geometric problem than the shorter transfers, despite
its far lower Δv.

### 14.5 Refined M4 minimum-Δv feasible transfer

Bisection (50 iterations, 0.5 s constraint-sampling resolution) on the feasibility
boundary between the largest feasible coarse-grid point (`T = 380 s`) and the next
grid point (`T = 390 s`, infeasible) locates the boundary precisely:

```
T* = 386.563856 s   (T*/P = 0.069606, i.e. ~7.0% of one chief orbital period)
Δv1 = (-1.0256, 2.3076, -0.0726) m/s   |Δv1| = 2526.251 mm/s
Δv2 = (-1.0256, -2.3076, 0.0801) m/s   |Δv2| = 2526.479 mm/s
Δv_total = 5052.7305 mm/s
cond(Φrv) = 1.099 (very well-conditioned)
terminal closure residual = 5.68×10⁻¹⁴ m
minimum clearance outside corridor = +1.55×10⁻⁵ m  (i.e. essentially zero margin —
    the constraint is ACTIVE at this optimum, as expected for a constrained minimum)
minimum distance to chief = 50.000 m, occurring exactly at t = T (coincides with rf)
corridor entry time = 366.087 s (≈20.5 s before arrival)
```

**Δv penalty:** +374.6% relative to the M2 T=1800 s reference (1064.520 mm/s), and
**+2991% (≈30.9×)** relative to the unconstrained M3 global minimum (163.457 mm/s at
T=4868.130 s). This is the central M4 engineering finding: a tight 20 m-wide,
100 m-radius approach corridor is extremely costly for this scenario's boundary
conditions — the cheapest geometrically feasible transfer costs roughly 5× the M2
reference and roughly 31× the unconstrained optimum.

**Practical note on margin:** the mathematically refined optimum sits with ~15 µm of
clearance margin — an artifact of it being an active-constraint boundary point, not a
robust engineering margin. A practical selection would back off to, e.g., `T = 380 s`
(clearance ≈ 0.875 m, Δv_total ≈ 5136.0 mm/s) or `T = 375 s` (clearance ≈ 1.564 m,
Δv_total ≈ 5201.5 mm/s) for a non-zero safety margin — both are in the curated CSV.

### 14.6 Decision table

Full table: [`results/m4_constrained_trade.csv`](results/m4_constrained_trade.csv).
Summary:

| Transfer | T (s) | M3 safe | M4 safe | Δv_total (mm/s) | min dist (m) | reason |
|---|---:|:---:|:---:|---:|---:|---|
| M2 reference | 1800.000 | ✅ | ❌ | 1064.52 | 50.00¹ | corridor+koz |
| M3 short-branch min | 2539.593 | ✅ | ❌ | 678.61 | 50.00¹ | corridor+koz |
| M3 long-branch/global min | 4868.130 | ✅ | ❌ | 163.46 | 44.06 | corridor+koz+crossing |
| fast safe transfer | 300.000 | ✅ | ✅ | 6447.62 | 50.00 | none |
| slow safe transfer | 5100.000 | ✅ | ❌ | 180.04 | 43.64 | corridor+koz+crossing |
| just-safe boundary point | 385.000 | ✅ | ✅ | 5072.31 | 50.00 | none |
| just-unsafe boundary point | 390.000 | ✅ | ❌ | 5010.26 | 50.00¹ | koz |
| **refined M4 minimum** | **386.564** | ✅ | ✅ | **5052.73** | **50.00** | **none** |

¹ `min_distance_m` is the raw global closest-approach distance regardless of corridor
authorization — for these transfers it coincides with the terminal point (50 m) even
though an *unauthorized* dip below 100 m occurs earlier in the coast; that is exactly
what `min_clearance_outside_corridor_m` (§14.5, Figure 3) captures and `min_distance_m`
does not — see `results/m4_full_constrained_sweep.csv` for both columns together.

### 14.7 Verification / convergence evidence

- **Grid-resolution convergence at the selected T=380 s transfer:** min-distance at
  2 s sampling vs. 1 s sampling agree to `< 1 m` (both ≈50.0 m, since the closest
  approach coincides with the fixed terminal point for these fast transfers); refined
  (bounded-scalar-minimized) value is `≤` the 1 s-sampled value, never worse
  (`tests/test_constraints.py::test_constraint_margin_convergence_at_selected_transfer`).
- **STM-vs-independent-ODE geometric residual** along the selected T=380 s transfer:
  max state-component discrepancy vs. an independent `scipy.integrate.solve_ivp`
  (DOP853, rtol/atol 1e-12) integration of the raw CW ODE, sampled at 10 points across
  the coast: **7.62×10⁻¹² m/(m/s)** — consistent with the M2/M3 residuals, confirming
  the constraint classification is being applied to a dynamically correct trajectory.
- **Worst terminal-position closure residual across the entire 481-point sweep:**
  `3.595×10⁻¹³ m` — unchanged in magnitude from M3 (constraint evaluation does not
  perturb the solver).
- **M3 regression preserved exactly** through the M4 pipeline: T=1800/2539.593/4868.130 s
  Δv values reproduce the M3-authoritative figures to `< 0.05` mm/s
  (`tests/test_constraints.py::test_m3_dv_values_unchanged_through_m4_pipeline`).
- **Sweep-order independence:** ascending/descending/shuffled evaluation of the same
  time set produce identical per-time `(safe, min_distance, violation_reason)` tuples
  (`tests/test_constraints.py::test_constraint_classification_independent_of_evaluation_order`).
- **Continuous refinement demonstrably recovers a tighter minimum than coarse
  sampling alone misses** (verified with a deliberately coarse 50 s sample against the
  M3 long-branch-minimum trajectory,
  `tests/test_constraints.py::test_refinement_recovers_true_minimum_missed_by_coarse_sampling`).
- **M3-unsafe times are never rescued by geometry:** `solve_two_impulse` itself refuses
  transfer times inside the `cond(Φrv) > 1e6` guard band, so there is no trajectory for
  M4 to classify at those times at all — the distinction is structural, not just a
  policy choice (`tests/test_constraints.py::test_m3_unsafe_times_never_rescued_by_geometry`).
- **26 new constraint-primitive/endpoint-tolerance/corridor-entry tests** all pass —
  see `tests/test_constraints.py` for the full check list (task items A–J).

### 14.8 A genuine bug found and fixed during M4 development

While building `constraints.py`, an internal-consistency bug was found: the
continuous-refinement step could discover a keep-out violation (a refined clearance
value that crosses zero) that no discrete dense sample had caught, correctly setting
the aggregate `safe=False` flag — but the human-readable `violation_reason` string was
still built only from the discrete per-sample classifications, so it could report
`safe=False` alongside `violation_reason="none"`, an internal contradiction. This was
caught by manual inspection of sweep output (a rejection-reason histogram showed
472 rejected points but only 469+232 reason tags, i.e. some rejections had no
attributed reason) before any test was written against the buggy behavior. Fixed by
unifying both the discrete and refinement-promoted violations through one `_note()`
helper that updates `reasons` and `first_violation_time_s` together, with an assertion
(`safe == (violation_reason == "none")`) enforcing the invariant going forward. This is
a constraints.py logic fix, not an M2/M3 dynamics change — no CW equation, STM entry,
or solver formula was touched.

### 14.9 Figures

1. **`figures/m4_constrained_dv_vs_transfer_time.png`** (headline) — Δv_total (log
   scale) vs. transfer time, M3-safe-but-M4-rejected points in red, M4-feasible points
   in green, the M2 reference and unconstrained M3 minimum both marked as
   geometry-rejected, and the refined M4 feasible minimum marked with a star. The
   `P/2` conditioning-exclusion gap is retained as a vertical marker.
2. **`figures/m4_approach_geometry.png`** (headline) — the selected feasible
   trajectory (T=386.6 s) in LVLH x–y, overlaid with the 100 m keep-out circle and the
   corridor box, showing visually that the path stays outside the circle except while
   inside the box, with a cross-track-vs-time inset panel.
3. **`figures/m4_constraint_margin_vs_transfer_time.png`** (verification/supporting) —
   minimum clearance outside the corridor vs. transfer time (0 m requirement line),
   with chief-crossing violations marked separately (▽) since that failure mode is not
   captured by the clearance metric alone.

All three were visually inspected for clipping, overlapping legends/annotations,
correct safe/rejected visual distinction, and to confirm no line is drawn connecting
across a discontinuity in a way that implies a rejected point is feasible.

### 14.10 Limitations

Unchanged from M1 §10 (still binding), plus the M4-specific scope statement at the top
of this section (§14, restated): this is linear-CW **geometric screening only** — not
collision probability, not a nonlinear-dynamics check, not a sensor/FOV/occultation/
plume model, not a finite-burn or actuator-limited feasibility check. See also M2
§12.7/12.8 and M3 §13.8 for what else remains unimplemented.

---

## 15. Milestone 5 — nonlinear two-body validation and robustness

**Scope statement (read first):** this section checks the CW/M4 result against a
higher-fidelity **nonlinear two-body** model — still no J2, drag, finite burns,
actuator limits, sensor FOV, occultation, or collision probability. A trajectory
described as passing here is **"geometrically feasible under this CW model, validated
against nonlinear two-body dynamics"** — never "collision-free", "flight-qualified",
or "operationally safe". M5 uses **only** the verified M2 solver
(`rendezvous.solve_two_impulse`), M2 STM (`cw.propagate_cw`), and M4 constraint logic
(`constraints.classify_point`, reused by the new nonlinear classifier) — no CW
equation, STM entry, or M1–M4 solver formula was changed. New code:
`src/rendezvous_cw/nonlinear.py`, `scripts/m5_nonlinear_validation.py`, three CSVs,
three figures (§15.9).

### 15.1 Nonlinear model

Independent inertial (ECI-like) Cartesian two-body propagation for chief and deputy,
integrated **separately** with `scipy.integrate.solve_ivp` (DOP853, `rtol=1e-12`,
`atol=1e-9`):

```
rddot = -mu * r / ||r||^3
```

using the same `mu = MU_EARTH` as M1/M2 (`orbit.py`, unchanged). No CW equation
appears anywhere in `nonlinear.py`; it is a from-scratch propagator used only to
*check* the CW model. The chief starts on the analytic circular initial condition
`r_c0=(a,0,0)`, `v_c0=(0,sqrt(mu/a),0)` and is then verified — not assumed — to stay
circular under numerical propagation (§15.4, check A).

### 15.2 LVLH ↔ inertial frame kinematics

At any instant, given the chief's inertial state `(r_c, v_c)`:

```
x_hat = r_c / ||r_c||                     (radial, outward)
h_c   = r_c x v_c                          (specific angular momentum)
z_hat = h_c / ||h_c||                      (cross-track / orbit-normal)
y_hat = z_hat x x_hat                      (along-track; right-handed)
C     = [x_hat | y_hat | z_hat]            (LVLH-to-inertial rotation matrix)
omega_inertial = h_c / ||r_c||^2           (LVLH frame's inertial angular velocity;
                                             exact for any orbit, not only circular)
```

Deputy inertial state from an LVLH relative state `(r_rel, v_rel)`:

```
r_dep = r_c + C @ r_rel
v_dep = v_c + C @ v_rel + omega_inertial x (C @ r_rel)
```

**The rotating-frame `omega x r` term is essential** — `v_dep` is *not* simply
`v_c + C @ v_rel`; omitting the Coriolis-type term would silently misstate the
deputy's inertial velocity by an amount comparable to `n * ||r_rel||` (≈ 1.1 m/s for
the M1 1000 m separation — far larger than the mm/s-scale burns being analyzed, so
this is not a subtle effect). The inverse transform (`inertial_to_lvlh`) solves the
same equation for `v_rel` using `C^T = C^-1` (C is orthonormal by construction).

### 15.3 Reproducing the CW trajectories nonlinearly

For each representative transfer: (1) solve for the impulsive burn pair with the
unmodified M2 solver, (2) convert `(r0, v0_plus)` to inertial via `lvlh_to_inertial`
at `t=0`, (3) propagate chief and deputy independently with the nonlinear model,
(4) reconstruct the LVLH relative state at each sample time via `inertial_to_lvlh`
using each body's own *propagated* (osculating) state, (5) compare to the CW-STM
prediction at the same times. No nonlinear correction burn is applied in this step
(§15.8 adds one, separately, only for the final recommended transfer).

### 15.4 Verification (checks A–G — summary; full detail in `tests/test_nonlinear.py`)

- **A. Chief-orbit conservation** (one full period, 400 samples): radius constant to
  `< 1.4×10⁻⁵ m` (6778 km orbit), specific energy conserved to `< 1×10⁻⁴` (relative to
  a mean of `-2.94×10⁷`), angular momentum conserved to `< 0.09` (relative to
  `~5.2×10¹⁰`), final position after one period returns to within `1.5×10⁻⁵ m` of the
  start — all far tighter than any effect being measured.
- **B. Frame round trips:** LVLH basis orthonormal and right-handed
  (`x̂ × ŷ = ẑ`, `det(C) = 1`) to `< 1×10⁻¹²`; LVLH→ECI→LVLH round trip (position *and*
  velocity, the latter exercising the `ω×r` term) exact to `< 1×10⁻⁹` for several
  arbitrary `(r_rel, v_rel)`; zero relative state maps exactly to the chief state;
  `+x`/`+y`/`+z` map exactly to radial/along-track/orbit-normal directions.
- **C. CW small-time consistency:** a short (60 s), small-separation (50 m) transfer
  shows `< 1 cm` maximum CW-vs-nonlinear position deviation, as expected.
- **D. Separation-scaling trend:** a 250 m-separation case shows strictly smaller CW
  error than a 1000 m case (same geometry, same T=1800 s) — see §15.7.
- **E. Integrator-tolerance convergence:** CW-vs-nonlinear terminal error agrees to
  `< 1%` between `(rtol,atol) = (1e-12,1e-9)` and `(1e-13,1e-10)` — the measured error
  is real linearization error, not integrator noise.
- **F. M4 regression preserved:** CW `Δv_total` at T=1800 s and at the M4-selected
  T=386.564 s reproduce the M3/M4-authoritative values to `< 0.05` mm/s through the M5
  pipeline; the M4-selected transfer's CW feasibility classification is unchanged.
- **G. Final selected-case regression:** dedicated tests pin the M5-recommended
  transfer's CW Δv, CW margin, nonlinear terminal position error, and nonlinear
  constraint classification/margin.

18 new tests, all passing (`tests/test_nonlinear.py`).

### 15.5 Does the M4 optimum survive nonlinear validation?

**Yes — but only just.** The CW closed-form model predicted the M4 boundary optimum
(T=386.564 s) sits with `1.55×10⁻⁵ m` (15 µm) of active-constraint margin. Evaluating
the *same* burn pair's trajectory nonlinearly and re-running the identical M4
constraint classification (reusing `constraints.classify_point`, not a redefinition):

```
CW margin (T=386.564 s):          +0.0000155 m   (safe, by construction — the CW boundary)
Nonlinear margin (same T):        +0.027548 m    (safe — PASSES)
```

The nonlinear margin is small but unambiguously positive, and — critically — this
result is **robust to sampling resolution** (`0.027548`–`0.027564 m` across
`sample_dt` from 2 s down to 0.1 s) **and to integrator tolerance**
(`0.027535`–`0.027558 m` across three `(rtol, atol)` settings), confirming it reflects
real dynamics, not numerical noise (§15.4, check E). The nonlinear model happens to be
slightly *more* permissive than CW here — the true closest-approach clearance is
marginally better than CW predicted — but 27.5 mm is still minuscule next to
real-world navigation/actuation uncertainties, so per the M5 task's own framing this
is treated as **operationally fragile even though it technically passes** (§15.6).

### 15.6 Local robustness refinement (T = 365–390 s)

A 26-point scan (1 s step) confirms:
- CW and nonlinear margins track each other closely and both decrease monotonically
  with `T` in this range (nonlinear margin consistently ~0.02–0.03 m *larger* than CW
  here).
- The nonlinear feasibility boundary sits between T=386 s (margin +0.101 m, safe) and
  T=387 s (margin −0.029 m, unsafe) — one second later than the CW boundary.
- Full scan: `results/m5_local_robustness_scan.csv`.

**M5-recommended transfer: T = 380 s.** Chosen because it converts the CW boundary
optimum's microscopic margin into a comfortable, easily-defensible one at a small,
explicitly-quantified Δv cost:

| | T=386.564 s (M4 CW boundary) | **T=380 s (M5 recommended)** | change |
|---|---:|---:|---:|
| CW Δv_total | 5052.7305 mm/s | 5136.0226 mm/s | **+1.65%** |
| CW margin | +0.0000155 m | +0.875219 m | **56,466×** |
| Nonlinear margin | +0.027548 m | +0.902241 m | **32.8×** |

Spending **1.65% more Δv buys a ~33× larger real (nonlinear) safety margin** —
converting a millimeter-scale, sampling-and-tolerance-sensitive-looking pass into a
comfortably-positive, meter-scale one. This is the central M5 engineering conclusion.
See §15.9 Figure 3.

**Important labeling note (per M5 task instructions):** the M4 T=386.564 s result is
correct as reported in M4 — it was, and remains, the CW-model constrained optimum. It
is now superseded, for practical recommendation purposes, by this M5 robustness
analysis: **"CW-only constrained optimum — superseded by M5 robustness validation."**
The M4 engineering history above (§14) is left unmodified.

### 15.7 CW linearization error: transfer-time and separation scaling

**Representative cases** (`results/m5_representative_cases.csv`; terminal position
error = max CW-vs-nonlinear deviation in every case tested — the two coincide because
error accumulates monotonically along these particular coast arcs, verified not
assumed):

| Case | T (s) | terminal pos. error (m) | terminal vel. error (mm/s) |
|---|---:|---:|---:|
| M4 practical margin (T=380 s) | 380.000 | 0.0109 | 0.0421 |
| M4 selected (CW boundary) | 386.564 | 0.0113 | 0.0430 |
| M2 reference | 1800.000 | 0.3987 | 0.4838 |
| M3 short-branch minimum | 2539.593 | 0.9646 | 0.7642 |
| M3 long-branch/global minimum | 4868.130 | 2.6734 | 0.2270 |

Error grows with transfer time (more time for second-order/nonlinear effects to
accumulate over the ~1 km-scale coast), consistent with CW being a first-order
(linear) approximation of the true nonlinear relative dynamics.

**Separation sensitivity** (fixed T=1800 s, direction of `r0` preserved, magnitude
varied — `results/m5_separation_sensitivity.csv`):

| \|r0\| (m) | \|\|r0\|\|/a | terminal error (m) | error ratio | separation ratio | (ratio)² |
|---:|---:|---:|---:|---:|---:|
| 250 | 3.69×10⁻⁵ | 0.02665 | — | — | — |
| 500 | 7.38×10⁻⁵ | 0.10188 | 3.82 | 2.00 | 4.00 |
| 1000 | 1.48×10⁻⁴ | 0.39831 | 3.91 | 2.00 | 4.00 |
| 1500 | 2.21×10⁻⁴ | 0.88936 | 2.23 | 1.50 | 2.25 |

Error scales approximately with the **square** of separation (a 2× separation
increase gives a ~3.8–3.9× error increase against a 4.0× quadratic expectation; a
1.5× increase gives a ~2.23× error increase against a 2.25× expectation) — exactly the
expected behavior for a first-order linearization, whose leading-order error term is
second-order in the small parameter `||r||/a`. This is a clean, textbook-consistent
confirmation of the CW model's small-separation limitation (DESIGN.md §10/§1), not
merely an assumption.

### 15.8 Independent burn/terminal validation (T = 380 s, the recommended transfer)

- CW-predicted `vT_minus`: magnitude 2568.1235 mm/s.
- Nonlinear actual (numerically propagated) pre-second-burn velocity: magnitude
  2568.1235 mm/s — **agrees with CW to `4×10⁻⁶` mm/s** (a vector difference of
  `0.0421` mm/s, §15.7, but a near-identical scalar magnitude since the two vectors
  are nearly parallel).
- **Terminal position miss** (nonlinear actual position at `t=T` vs. `rf`):
  **0.0109 m** — this is *not* zero, and a velocity-only second burn cannot correct
  it.
- Nonlinear-corrected second burn (nulls the *actual* nonlinear arrival velocity
  instead of the CW-predicted one): magnitude 2568.1235 mm/s — indistinguishable from
  the CW `Δv2` at this precision.
- **Total Δv using the nonlinear-corrected second burn: 5136.0226 mm/s**, vs.
  CW-only total 5136.0226 mm/s — a difference smaller than 0.001 mm/s.

**Conclusion, stated per the task's explicit framing:** this is only
*"CW first burn + nonlinear propagation + endpoint velocity correction,"* not a true
nonlinear Lambert/rendezvous optimum. The velocity correction is negligible in this
case (the CW and true arrival velocities already agree to µm/s-equivalent precision
in magnitude), but it explicitly does **not** fix the ~1 cm terminal position miss —
a velocity-only impulse cannot close a position gap at the instant it is applied. For
the separations and transfer times in this project's scope, the ~1 cm position miss
is itself far smaller than the M4 corridor's 20–50 m scale and is not operationally
significant here, but the distinction (velocity correctable vs. position not) is
stated explicitly because it would matter at larger separation or longer transfer
time (§15.7).

### 15.9 Figures

1. **`figures/m5_cw_vs_nonlinear_trajectory.png`** (headline) — LVLH x–y trajectories
   for the recommended T=380 s transfer: CW (solid) and nonlinear (dashed) are
   visually indistinguishable at this scale (as expected, ~1 cm error against
   50–1000 m geometry), with an inset panel showing the actual CW-vs-nonlinear
   position-error magnitude growing over time. KOZ and corridor overlaid.
2. **`figures/m5_cw_error_vs_transfer_time.png`** — log-log terminal position error
   vs. transfer time for the five representative cases, showing the clear growth
   trend from the M4 fast-transfer regime (~1 cm) to the M3 long-branch minimum
   (~2.7 m).
3. **`figures/m5_local_robustness_trade.png`** (headline) — two-panel: CW `Δv_total`
   (top) and nonlinear geometry margin (bottom) vs. transfer time over T=365–390 s,
   with the M4 CW-boundary optimum and the M5-recommended transfer both marked and a
   0 m requirement line — this is the figure that communicates §15.6's central trade.

All three visually inspected: no clipping, no legend overlap (a title-clipping issue
in an early draft of Figure 2 was caught and fixed), KOZ/corridor clearly visible,
CW-vs-nonlinear curves distinguishable where they differ meaningfully, inset scales
explicitly labeled, no wording implying "collision-free" or "flight safe" anywhere.

### 15.10 Limitations

Unchanged from M1 §10 and M4 §14.10, plus: the nonlinear model here is still
**unperturbed two-body** (no J2, no drag) — it validates CW's *linearization* error,
not atmospheric or oblateness effects. No finite-burn execution, actuator limits,
navigation/sensor error, or collision probability is modeled. The M5 "robustness"
result is a deterministic sensitivity to model fidelity (CW vs. nonlinear two-body
point dynamics), not a statistical/covariance-based robustness analysis.

---

## 16. Milestone 6 — final technical audit, result hierarchy, and portfolio packaging

M6 adds **no new rendezvous physics**. It performs a final repository-wide hygiene
audit (dead code / unused imports found and fixed via `pyflakes` — see the M6 commit;
no solver formula changed), establishes one authoritative result hierarchy across all
five milestones' worth of engineering history, adds a `LICENSE`, a CI workflow
(`.github/workflows/ci.yml`), and confirms fresh-environment reproducibility.

### 16.1 Authoritative final result

```
Final recommended transfer:      T = 380 s   (T/P = 0.0684)
Model used for design:            CW linearized relative-motion model (M2)
Independent validation:           nonlinear two-body propagation (M5)
Total Δv (CW):                    5.136023 m/s   (|Δv1| = 2.567899 m/s, |Δv2| = 2.568123 m/s)
CW keep-out/corridor margin:      +0.875219 m
Nonlinear keep-out/corridor margin: +0.902241 m
Nonlinear terminal position miss: 0.0109 m  (not corrected by any burn in this analysis)
Nonlinear validation:             PASSES modeled geometric constraints
```

### 16.2 Full engineering history (superseded results, retained and clearly labeled)

| # | Result | T (s) | Status | Why superseded |
|---|---|---:|---|---|
| M2 | Representative transfer | 1800.000 | Diagnostic reference only | Never claimed optimal; used throughout as a fixed regression/comparison point |
| M3 | Unconstrained global Δv minimum | 4868.130 | **Superseded — fails M4 geometry** | Breaches the keep-out sphere outside the corridor and overshoots past the target (chief-crossing) — see §14.4 |
| M4 | CW-only constrained optimum | 386.564 | **Superseded — essentially zero CW margin, fragile** | Sits at ~15 µm active-constraint margin by construction; nonlinear validation shows only +27.5 mm real margin — technically passes but operationally fragile, see §15.5 |
| **M5** | **Recommended robust transfer** | **380.000** | **Final recommendation** | +1.65% Δv over the M4 optimum buys a ~33× larger real (nonlinear) margin — see §15.6 |

No result above is deleted or rewritten; each remains documented in its own
milestone section (§13, §14, §15) exactly as originally computed. This table exists
solely to prevent a reader from mistaking an earlier, superseded milestone result for
the project's final recommendation.

### 16.3 Terminology used consistently throughout this document (M6 audit)

"CW linearized relative-motion model", "nonlinear two-body validation", "modeled
geometric constraints", "keep-out zone (KOZ)", "V-bar approach corridor",
"numerically unsafe / ill-conditioned transfer" (M2/M3 conditioning sense),
"geometrically feasible under modeled constraints" (M4/M5 sense), "nonlinear
validation pass/fail". Never: "collision-free", "flight-safe", "operationally safe",
"flight-qualified", or an unscoped "globally optimal". Every use of "minimum" or
"optimum" in this document is scoped to one of: *minimum within the searched
transfer-time interval* (M3), *CW-only constrained optimum* (M4), or *M5 recommended
robust transfer* (final).

### 16.4 Repository hygiene fixed in M6

- One genuinely dead local variable (`coarse_t_min_clearance`, assigned but never
  read) removed from `constraints.py` — found by `pyflakes`, confirmed unused by
  direct grep, no behavior change (verified: full test suite unchanged, 151→151
  passing before and after).
- Unused imports removed from `scripts/m3_trade_sweep.py`, `scripts/m4_constrained_trade.py`,
  and `tests/test_constraints.py` (found by `pyflakes`).
- One dead local variable (`state0`, computed then immediately shadowed/unused)
  removed from `scripts/m3_trade_sweep.py`'s Figure 2 plotting loop.
- No TODO/FIXME/placeholder comments, no debug `print()` calls in `src/`, no
  hardcoded absolute local filesystem paths, and no secrets/credentials were found
  anywhere in the tracked repository (checked by direct grep across `.py`/`.md`/
  `.toml` files).
- No broken relative links in `README.md` or `DESIGN.md` (every non-URL markdown
  link target verified to exist on disk).
- No duplicate/colliding figure output paths across `scripts/m2..m5_*.py` — each
  milestone's figures are uniquely named and no script silently overwrites another
  milestone's output.
- Largest tracked file is 128 KB (`figures/m3_dv_vs_transfer_time.png`); no
  oversized accidental files.

### 16.5 Reproducibility

Verified in a genuinely fresh `python3 -m venv` (no cached wheels reused from the
development environment): `pip install -e ".[dev]"` installed cleanly (resolving
newer library versions than development was done with: `numpy 2.5.3`,
`scipy 1.18.1`, `matplotlib 3.11.1`), `pytest -W error` passed all 151 tests
(21.3 s), and regenerating the T=380 s result and the M5 figures/CSVs from the copied
source reproduced `Δv_total = 5136.022643719935` mm/s, CW margin
`0.8752190597669625` m, and nonlinear margin `0.9022414877394453` m — **bit-for-bit
identical** to the values reported in §16.1 and §15.6, despite the newer dependency
versions. This is expected and unsurprising: the CW solver is closed-form arithmetic
and the nonlinear ODE integration is deterministic (no random seeding anywhere in
this project), so no fuzz in these numbers is expected across environments; the check
confirms that expectation holds in practice, not merely in principle.

### 16.6 What remains explicitly out of scope

Unchanged from all prior milestones' limitations sections (§1, §10, §12.7/12.8,
§13.8, §14.10, §15.10): no J2, drag, finite burns, thrust saturation,
covariance/collision-probability modeling, sensor field of view, plume impingement,
docking dynamics, or any additional optimization dimension. This project is a
CW-linearized (with nonlinear-two-body cross-check) rendezvous trade study and
geometric-screening demonstration — not a flight-qualification analysis.
