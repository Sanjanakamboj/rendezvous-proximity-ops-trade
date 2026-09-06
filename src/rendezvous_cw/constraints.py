"""Milestone 4: proximity-operations geometric screening (keep-out sphere,
final V-bar approach corridor, no-chief-crossing) for a CW coast trajectory.

Scope (restated explicitly, per DESIGN.md M4 section): this module performs
GEOMETRIC screening only, inside the linearized CW model already verified
in M2/M3. It does not model collision probability, sensor fields of view,
plume impingement, Earth occultation, finite burns, or any nonlinear
two-body effect. A trajectory that "passes" these checks is described as
"geometrically feasible under this CW model" -- never as "collision-free",
"flight safe", or "operationally safe".

Three constraints, applied to the post-first-burn coast trajectory
r(t) = (x(t), y(t), z(t)) in the LVLH/Hill frame (DESIGN.md section 3):

  A. Spherical keep-out zone (KOZ), radius r_KOZ = 100 m, centered on the
     chief, EXCEPT that entry is authorized while inside the final
     approach corridor (B).
  B. Final approach corridor: -100 m <= y <= -50 m, |x| <= 20 m,
     |z| <= 20 m. Only inside this box is `||r|| < r_KOZ` authorized.
  C. No chief crossing: y(t) <= -50 m for the entire transfer (with a
     small documented floating-point tolerance at the endpoint, since the
     M1 target itself sits exactly at y = -50 m).

Within an unauthorized keep-out breach (constraint A/B), two distinct
failure reasons are distinguished rather than collapsed into one flag:
  - `koz_violation`: the trajectory enters `||r|| < r_KOZ` while outside
    the corridor's along-track (y) range entirely -- a generic keep-out
    breach unrelated to any approach-corridor attempt.
  - `corridor_violation`: the trajectory enters `||r|| < r_KOZ` while
    within the corridor's y-range but outside its lateral (x, z) bounds --
    i.e. attempting to approach at the right along-track distance but
    through the wrong lateral position.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.optimize import minimize_scalar

from .cw import propagate_cw

#: Keep-out sphere radius, meters (DESIGN.md M4 section, constraint A).
R_KOZ = 100.0

#: Final approach corridor bounds (DESIGN.md M4 section, constraint B).
CORRIDOR_Y_MIN = -100.0
CORRIDOR_Y_MAX = -50.0
CORRIDOR_X_BOUND = 20.0
CORRIDOR_Z_BOUND = 20.0

#: No-chief-crossing limit (constraint C) and its floating-point tolerance.
#: The M1 target itself sits exactly at y = -50 m, and M2/M3 closure
#: residuals are ~1e-13 m, so a 1e-6 m tolerance comfortably absorbs
#: floating-point noise at the endpoint while still catching any physically
#: meaningful crossing.
CHIEF_CROSSING_Y_LIMIT = -50.0
CHIEF_CROSSING_TOL_M = 1.0e-6

#: Boundary tolerance applied to the corridor's own y/x/z limits, for the
#: same reason as CHIEF_CROSSING_TOL_M: the M1 target r_f = (0, -50, 0)
#: sits EXACTLY on the corridor's y = -50 m boundary, and closed-form STM
#: evaluation leaves ~1e-13 m floating-point residuals there (verified in
#: M2/M3). Without this tolerance, the valid terminal point itself is
#: occasionally (depending on transfer time, hence STM argument) computed
#: as e.g. y = -49.999999999999886 -- a ~1.1e-13 m residual that is
#: numerically ABOVE -50 m and would be misclassified as outside the
#: corridor. 1e-6 m comfortably absorbs this while remaining far below any
#: physically meaningful corridor-boundary distance.
CORRIDOR_BOUNDARY_TOL_M = 1.0e-6

#: Default dense-sampling spacing for continuous-time constraint checking
#: (DESIGN.md M4 section; "recommended production spacing: <= 1 s").
DEFAULT_SAMPLE_DT_S = 1.0


def inside_lateral_corridor(
    x: float, z: float,
    x_bound: float = CORRIDOR_X_BOUND, z_bound: float = CORRIDOR_Z_BOUND,
    tol: float = CORRIDOR_BOUNDARY_TOL_M,
) -> bool:
    """True if (x, z) is within the corridor's lateral (radial/cross-track) box.

    `tol` (default CORRIDOR_BOUNDARY_TOL_M) absorbs floating-point residuals
    at the box boundary -- see that constant's docstring.
    """
    return abs(x) <= x_bound + tol and abs(z) <= z_bound + tol


def inside_y_corridor(
    y: float, y_min: float = CORRIDOR_Y_MIN, y_max: float = CORRIDOR_Y_MAX,
    tol: float = CORRIDOR_BOUNDARY_TOL_M,
) -> bool:
    """True if y is within the corridor's along-track range (+/- `tol`)."""
    return (y_min - tol) <= y <= (y_max + tol)


def inside_approach_corridor(
    r_vec: np.ndarray,
    y_min: float = CORRIDOR_Y_MIN, y_max: float = CORRIDOR_Y_MAX,
    x_bound: float = CORRIDOR_X_BOUND, z_bound: float = CORRIDOR_Z_BOUND,
) -> bool:
    """True if r_vec = (x, y, z) is inside the full final approach corridor box."""
    x, y, z = r_vec
    return inside_y_corridor(y, y_min, y_max) and inside_lateral_corridor(x, z, x_bound, z_bound)


def keep_out_violation(r_vec: np.ndarray, r_koz: float = R_KOZ) -> bool:
    """Raw keep-out check (ignores the corridor exception): True if ||r|| < r_koz."""
    return float(np.linalg.norm(r_vec)) < r_koz


@dataclass(frozen=True)
class PointClassification:
    """Classification of a single relative-position sample."""

    distance_m: float
    inside_corridor: bool
    koz_violation: bool
    corridor_violation: bool
    chief_crossing_violation: bool

    @property
    def safe(self) -> bool:
        return not (self.koz_violation or self.corridor_violation or self.chief_crossing_violation)


def classify_point(
    r_vec: np.ndarray,
    r_koz: float = R_KOZ,
    y_min: float = CORRIDOR_Y_MIN, y_max: float = CORRIDOR_Y_MAX,
    x_bound: float = CORRIDOR_X_BOUND, z_bound: float = CORRIDOR_Z_BOUND,
    chief_crossing_limit: float = CHIEF_CROSSING_Y_LIMIT,
    chief_crossing_tol: float = CHIEF_CROSSING_TOL_M,
) -> PointClassification:
    """Classify one relative-position sample against all three M4 constraints."""
    r_vec = np.asarray(r_vec, dtype=float)
    x, y, z = r_vec
    distance = float(np.linalg.norm(r_vec))
    in_corridor = inside_approach_corridor(r_vec, y_min, y_max, x_bound, z_bound)

    koz = False
    corridor_viol = False
    if distance < r_koz and not in_corridor:
        if inside_y_corridor(y, y_min, y_max):
            # Right along-track range, wrong lateral placement.
            corridor_viol = True
        else:
            # Keep-out breach with no relation to the approach corridor.
            koz = True

    chief_crossing = y > (chief_crossing_limit + chief_crossing_tol)

    return PointClassification(
        distance_m=distance,
        inside_corridor=in_corridor,
        koz_violation=koz,
        corridor_violation=corridor_viol,
        chief_crossing_violation=chief_crossing,
    )


@dataclass(frozen=True)
class TrajectoryConstraintResult:
    """Result of continuous-time M4 constraint screening for one coast trajectory."""

    t: float  # transfer time, s
    sample_dt_s: float

    safe: bool

    min_distance_m: float
    time_of_min_distance_s: float

    min_y_m: float
    max_y_m: float

    # Minimum of (||r(t)|| - r_KOZ) over samples NOT inside the corridor;
    # NaN if the trajectory is inside the corridor at every sample.
    # Positive = safe margin, negative = KOZ violation, per DESIGN.md.
    min_clearance_outside_corridor_m: float

    koz_violation: bool
    corridor_violation: bool
    chief_crossing_violation: bool

    first_violation_time_s: Optional[float]
    violation_reason: str

    corridor_entry_time_s: Optional[float]  # first t where inside_approach_corridor


def _dense_sample_times(t: float, sample_dt: float) -> np.ndarray:
    n_samples = max(2, int(np.ceil(t / sample_dt)) + 1)
    return np.linspace(0.0, t, n_samples)


def trajectory_constraint_metrics(
    r0: np.ndarray, v0_plus: np.ndarray, t: float, n: float,
    sample_dt: float = DEFAULT_SAMPLE_DT_S,
) -> tuple[np.ndarray, np.ndarray]:
    """Dense-sample the post-burn-1 coast trajectory; return (times, positions).

    positions has shape (N, 3). Uses only the verified M2 propagate_cw
    (closed-form STM) -- no new dynamics.
    """
    times = _dense_sample_times(t, sample_dt)
    state0 = np.concatenate([np.asarray(r0, dtype=float), np.asarray(v0_plus, dtype=float)])
    states = np.array([propagate_cw(state0, tt, n) for tt in times])
    return times, states[:, :3]


def classify_trajectory_constraints(
    r0: np.ndarray, v0_plus: np.ndarray, t: float, n: float,
    sample_dt: float = DEFAULT_SAMPLE_DT_S,
    refine: bool = True,
) -> TrajectoryConstraintResult:
    """Continuous-time M4 constraint screening for one post-burn-1 coast trajectory.

    Method (DESIGN.md M4 section):
      1. Dense-sample r(t) at `sample_dt` spacing (default 1 s) using the
         verified closed-form CW STM (propagate_cw).
      2. Classify every sample against the three constraints.
      3. Identify the candidate minimum-distance sample and, separately,
         the candidate minimum-clearance-outside-corridor sample.
      4. If `refine`, bracket each candidate with its immediate dense-sample
         neighbors and refine with a bounded scalar minimizer
         (`scipy.optimize.minimize_scalar`, method="bounded") on the
         continuous position function -- this can only find a smaller
         (worse) minimum than the coarse grid saw, never a larger one, so
         it strictly improves (never weakens) the safety check.
    """
    times, positions = trajectory_constraint_metrics(r0, v0_plus, t, n, sample_dt)

    classifications = [classify_point(p) for p in positions]
    distances = np.array([c.distance_m for c in classifications])
    ys = positions[:, 1]

    # Coarse candidates.
    i_min_dist = int(np.argmin(distances))
    coarse_t_min_dist = float(times[i_min_dist])
    coarse_min_dist = float(distances[i_min_dist])

    outside_corridor_mask = np.array([not c.inside_corridor for c in classifications])
    if outside_corridor_mask.any():
        clearance_outside = distances[outside_corridor_mask] - R_KOZ
        i_outside = np.where(outside_corridor_mask)[0][np.argmin(clearance_outside)]
        coarse_min_clearance = float(distances[i_outside] - R_KOZ)
    else:
        coarse_min_clearance = float("nan")

    def state_at(tt: float) -> np.ndarray:
        state0 = np.concatenate([np.asarray(r0, dtype=float), np.asarray(v0_plus, dtype=float)])
        return propagate_cw(state0, tt, n)[:3]

    def distance_at(tt: float) -> float:
        return float(np.linalg.norm(state_at(tt)))

    # Large finite sentinel (not inf -- scipy's bounded Brent solver can
    # produce spurious NaN internal comparisons when the objective returns
    # inf) used to exclude in-corridor samples from the outside-corridor
    # clearance minimization.
    _NOT_A_CANDIDATE = 1.0e9

    def clearance_outside_at(tt: float) -> float:
        r = state_at(tt)
        if inside_approach_corridor(r):
            return _NOT_A_CANDIDATE
        return float(np.linalg.norm(r)) - R_KOZ

    min_distance_m = coarse_min_dist
    time_of_min_distance_s = coarse_t_min_dist
    min_clearance_outside_corridor_m = coarse_min_clearance
    refined_clearance_violation_time_s: Optional[float] = None

    if refine:
        lo = float(times[max(0, i_min_dist - 1)])
        hi = float(times[min(len(times) - 1, i_min_dist + 1)])
        if hi > lo:
            res = minimize_scalar(distance_at, bounds=(lo, hi), method="bounded")
            if res.fun < min_distance_m:
                min_distance_m = float(res.fun)
                time_of_min_distance_s = float(res.x)

        if outside_corridor_mask.any() and np.isfinite(coarse_min_clearance):
            lo_c = float(times[max(0, i_outside - 1)])
            hi_c = float(times[min(len(times) - 1, i_outside + 1)])
            if hi_c > lo_c:
                res_c = minimize_scalar(clearance_outside_at, bounds=(lo_c, hi_c), method="bounded")
                if res_c.fun < _NOT_A_CANDIDATE and res_c.fun < min_clearance_outside_corridor_m:
                    min_clearance_outside_corridor_m = float(res_c.fun)
                    refined_clearance_violation_time_s = float(res_c.x)

    min_y_m = float(np.min(ys))
    max_y_m = float(np.max(ys))

    # Aggregate violation flags. Start from the dense per-sample
    # classifications, then -- separately, and just as importantly --
    # promote a violation discovered only by the continuous refinement
    # step (a smaller-than-sampled clearance that crosses zero between
    # two samples). Both sources feed the SAME `reasons`/timing
    # bookkeeping below, so the reported `violation_reason` string can
    # never disagree with the aggregate `safe` flag (unlike an earlier
    # version of this function, which could report `safe=False` with
    # `violation_reason="none"` when only the refinement step found the
    # violation -- fixed here).
    koz_violation = any(c.koz_violation for c in classifications)
    corridor_violation = any(c.corridor_violation for c in classifications)
    chief_crossing_violation = bool(max_y_m > (CHIEF_CROSSING_Y_LIMIT + CHIEF_CROSSING_TOL_M))

    refinement_found_koz_violation = (
        np.isfinite(min_clearance_outside_corridor_m) and min_clearance_outside_corridor_m < 0.0
    )
    if refinement_found_koz_violation:
        koz_violation = True

    first_violation_time_s = None
    reasons: list[str] = []

    def _note(reason: str, tt: float) -> None:
        nonlocal first_violation_time_s
        if reason not in reasons:
            reasons.append(reason)
        if first_violation_time_s is None or tt < first_violation_time_s:
            first_violation_time_s = tt

    for c, tt in zip(classifications, times):
        if c.koz_violation:
            _note("koz_violation", float(tt))
        if c.corridor_violation:
            _note("corridor_violation", float(tt))
        if c.chief_crossing_violation:
            _note("chief_crossing_violation", float(tt))

    if refinement_found_koz_violation and refined_clearance_violation_time_s is not None:
        _note("koz_violation", refined_clearance_violation_time_s)

    if chief_crossing_violation and "chief_crossing_violation" not in reasons:
        idx = int(np.argmax(ys > (CHIEF_CROSSING_Y_LIMIT + CHIEF_CROSSING_TOL_M)))
        _note("chief_crossing_violation", float(times[idx]))

    violation_reason = "none" if not reasons else ",".join(reasons)
    safe = not (koz_violation or corridor_violation or chief_crossing_violation)
    assert safe == (violation_reason == "none"), (
        f"internal inconsistency: safe={safe} but violation_reason={violation_reason!r}"
    )

    corridor_entry_time_s = None
    for c, tt in zip(classifications, times):
        if c.inside_corridor:
            corridor_entry_time_s = float(tt)
            break

    return TrajectoryConstraintResult(
        t=t,
        sample_dt_s=sample_dt,
        safe=safe,
        min_distance_m=min_distance_m,
        time_of_min_distance_s=time_of_min_distance_s,
        min_y_m=min_y_m,
        max_y_m=max_y_m,
        min_clearance_outside_corridor_m=min_clearance_outside_corridor_m,
        koz_violation=koz_violation,
        corridor_violation=corridor_violation,
        chief_crossing_violation=chief_crossing_violation,
        first_violation_time_s=first_violation_time_s,
        violation_reason=violation_reason,
        corridor_entry_time_s=corridor_entry_time_s,
    )
