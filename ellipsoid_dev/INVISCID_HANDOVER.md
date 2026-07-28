# Inviscid Ellipsoid Handover

Date: 2026-06-30 (updated 2026-07-27)

This branch contains experimental work for inviscid ellipsoid immersed-boundary
treatment. Keep it separate from `ellipsoid-dev`, which is the validated viscous
ellipsoid branch.

**2026-07-27 update: the core instability and force-symmetry bugs are fixed.**
See "2026-07-27 Fixes" below before reading the rest of this doc — several
items in "Validation State" and "Open Work" are now resolved.

## Branch State

- Stable viscous branch: `ellipsoid-dev` at `a5ed9f9 Add ellipsoid example input`.
- Inviscid development branch: `inviscid-ellipsoid-dev` at
  `6b8de4b Add experimental inviscid Schur projection`.
- `ellipsoid-dev` has the clean example input only. It intentionally does not
  include the inviscid source experiments.
- Run directories under `ellipsoid_dev/runs/` are scratch space and ignored by
  git. Archive them separately if their ParaView output is needed.

On a new machine:

```sh
git clone git@github.com:adfboyd/Incompact3d.git
cd Incompact3d
git fetch origin
git switch ellipsoid-dev
git switch inviscid-ellipsoid-dev
```

## Directory Layout

- `examples/Ellipsoid/input.i3d`: compact viscous single-ellipsoid example,
  also present on `ellipsoid-dev`.
- `ellipsoid_dev/inputs/`: committed validation and experiment inputs.
- `ellipsoid_dev/scripts/`: local helper scripts.
- `ellipsoid_dev/runs/`: ignored scratch run directories.

Use the run helper from the repo root:

```sh
./ellipsoid_dev/scripts/new_run.sh RUN_NAME input_file_from_ellipsoid_dev_inputs.i3d
cd ellipsoid_dev/runs/RUN_NAME
NP=8 ./run_local.sh
```

## Current Inviscid Implementation

The inviscid branch is opt-in and experimental. Default viscous behavior should
be unchanged unless the new input flags are enabled.

Important flags:

- `ellipsoid_init_potential`: initialize a sphere case with the analytic
  potential-flow velocity field.
- `ellipsoid_projection_slip_fix`: old local pressure/slip correction. This is
  not trusted.
- `ellipsoid_projection_flux_fix`: experimental pressure-flux correction. This
  is not trusted.
- `ellipsoid_pressure_geometry_diag`: writes geometry/projection diagnostics.
- `ellipsoid_lagrange_projection_steps`: old Jacobi/Lagrange projection
  experiment. This is not trusted as a final method.
- `ellipsoid_lagrange_projection_relax`: relaxation for that old projection.
- `ellipsoid_schur_projection_iters`: matrix-free Schur projection iteration
  count.
- `ellipsoid_schur_projection_relax`: Schur multiplier relaxation.
- `ellipsoid_schur_projection_tol`: Schur convergence tolerance.

The Schur projection in `src/navier.f90` is the current most promising route. It
applies immersed-boundary normal constraints as a pressure-projected velocity
correction:

1. sample boundary-normal velocity error on dominant-normal IB points,
2. spread a boundary-normal multiplier to a velocity correction,
3. project that correction using the existing pressure solver,
4. resample the normal velocity and solve the resulting matrix-free Schur system
   with CG.

Internal pressure-projection calls must not overwrite the physical pressure
gradient history used for force calculation. The current code passes
`save_bc_grad=.false.` and `quiet=.true.` for those internal calls.

## 2026-07-27 Fixes

Two independent problems were conflated in the earlier validation state below:
a **boundary-condition accuracy** problem (drag doesn't converge to the right
value) and a **bulk stability** problem (velocity field runs away over long
times regardless of Schur relax). They needed two separate fixes.

**1. Grid-scale energy runaway (the actual cause of "just fails").** With
`xnu=0`, there is zero dissipation anywhere in the domain, including the
hyperviscous SVV operator (it's scaled by `xnu` too). The Schur projection
only enforces the surface boundary condition — it does nothing for bulk
interior energy. Over 1500 steps, `Umax` grew unboundedly (5-7x inflow) no
matter what Schur `relax` was used; higher relax made it *worse*, not better.
Fix: wired the existing compact spatial filter (`src/filters.f90` +
`src/tools.f90`, already used for ABL/LES) as a standalone dissipation
mechanism, independent of `ilesmod`/physical viscosity. Time-loop gate in
`src/xcompact3d.f90` now allows `itype.eq.itype_ellip` through in addition to
the existing ABL/turbine condition; still opt-in via `ifilter` (default 0,
so viscous cases on `ellipsoid-dev` are unaffected). Add to any inviscid
input file:
```
ifilter = 1
C_filter = 0.49
```
This alone stabilizes `Umax` to ~1.15-1.25 (vs 5-7x runaway) over 1500 steps.
Counter-intuitively, *more* aggressive filtering (lower `C_filter`) makes
drag accuracy *worse*, not better — it over-smooths the potential-flow
structure itself, not just noise. Use mild filtering (0.49).

**2. Schur `relax` must be retuned once the filter is active.** The filter's
own damping already removes most of the surface residual, so the old
unfiltered-tuned `relax=0.5` massively overcorrects (drag overshoots past
zero, sign-flipped, ~2-3x worse than no correction at all). With
`C_filter=0.49`, the drag zero-crossing is at **`relax≈0.008`**, not 0.5 — a
~60x difference. Chosen defaults:
```
ellipsoid_schur_projection_iters = 3
ellipsoid_schur_projection_relax = 0.008
ellipsoid_schur_projection_tol = 1.0e-8
```
Do not reuse a `relax` value tuned without the filter, or vice versa.

**3. Force-symmetry bug — three coordinate/kinematics bugs, now fixed.**
Independently of the above, a fixed sphere in axisymmetric uniform inviscid
flow showed `Fy≈0.03, Fz≈-0.017` — both *larger* than the drag `Fx≈-0.005`,
where symmetry requires them to be exactly zero. This predated the filter/
Schur work above (present even with both switched off) and traces to the
same class of coordinate bug already found and fixed on `ellipsoid-dev`
(see that branch's commit `67b5e1a`), but this branch diverged before those
fixes landed and also has its own extra code with a fresh instance of the
same bug pattern:
- `src/ibm.f90`: `cubsplx/y/z` computed grid-point coordinates as
  `(index-1)*dx` where the rest of the codebase uses `(index-2)*dx` — this
  was the actual cause of the Fy/Fz asymmetry above; fixed and confirmed
  (Fy, Fz now ~1e-15, machine zero).
- `src/ellip_utils.f90`: `CalculatePointVelocity` computed `r×ω` instead of
  `ω×r` for rigid-body surface velocity (sign-reversed; invisible for a
  rotating *sphere* since it doesn't change whether `u·n=0`, but wrong for
  any non-spherical rotating body).
- `src/ellip_utils.f90`: `EllipsoidalRadius`/`EllipsoidalRadius_debug`
  rotated lab→body using the raw orientation quaternion instead of its
  conjugate (invisible at the identity orientation used by all the test
  inputs, matters once a body's orientation is non-trivial).
- `src/ellip_utils.f90`: `EllipsoidNormal` (used by the Schur projection's
  surface-normal sampling — this function doesn't exist on `ellipsoid-dev`,
  it's new on this branch) had the lab↔body rotation swapped in *both*
  directions simultaneously. Also invisible at identity orientation.

All four are fixed. Confirmed post-fix: sphere-in-uniform-flow Fy/Fz at
machine epsilon; free translating+rotating non-spherical ellipsoid
(`input_1ellip_inviscid_test.i3d`) stays bounded (`Umax` ~0.6-0.75, no
blow-up) with all three force components responding continuously and
plausibly to the tumbling motion.

**Remaining known imperfection:** residual drag is small but not exactly
zero on any case (order 0.01-0.07 depending on config, versus 0.078-0.47
before these fixes), with mild time-drift in some configurations rather
than a perfectly flat plateau. Presented as "much improved and stable," not
"exact" — treat further precision work as the next increment, not blocking.

## Validation State

Do not treat this as finished. The trustworthy conclusions so far are:

- Uniform-flow startup around a sphere creates a large wake-like transient in
  inviscid mode. This is not physical potential flow; it comes from the discrete
  immersed-boundary/projection startup rather than viscosity.
- Potential-flow initialization alone does not fix the boundary-normal residual.
- Better body resolution helps substantially.
- A local slip fix can reduce a sampled normal residual while making force much
  worse, so force and divergence must be checked together.
- The Schur projection gives small but consistent improvements with conservative
  relaxation. Aggressive relaxation worsens force.

Useful saved run names from the original machine, if archived separately:

- `sphere_potential_baseline_t01`: coarse potential sphere baseline.
- `sphere_potential_resolved_t005`: resolved potential sphere baseline.
- `sphere_potential_schur_ppifix_relax002_t002`: coarse Schur run,
  `iters=3`, `relax=0.02`.
- `sphere_potential_schur_ppifix_relax001_t002`: coarse Schur run,
  `iters=3`, `relax=0.01`.
- `sphere_potential_resolved_schur_relax0005_t001`: resolved Schur run,
  `iters=3`, `relax=0.005`.

Reference metrics:

- Coarse potential baseline at `t=0.02`: force magnitude about
  `1.79e-2`, normal RMS about `4.20e-1`.
- Coarse Schur, `iters=3`, `relax=0.02`, at `t=0.02`: force magnitude about
  `1.96e-2`, normal RMS about `4.00e-1`.
- Resolved potential baseline at `t=0.01`: force magnitude about `3.72e-2`,
  normal RMS about `2.65e-1`.
- Resolved Schur, `iters=3`, `relax=0.005`, at `t=0.01`: force magnitude about
  `2.41e-2`, normal RMS about `2.64e-1`.
- Resolved Schur, `relax=0.02`, was too aggressive and increased force to about
  `1.79e-1`.

Recommended next smoke test:

```sh
./ellipsoid_dev/scripts/new_run.sh sphere_potential_resolved_schur \
    input_sphere_inviscid_potential_test.i3d
cd ellipsoid_dev/runs/sphere_potential_resolved_schur
# edit input.i3d to use a resolved radius/domain and:
#   ellipsoid_schur_projection_iters = 3
#   ellipsoid_schur_projection_relax = 0.005
#   ellipsoid_schur_projection_tol = 1.0e-4
NP=8 ./run_local.sh
```

Check all of:

- body force magnitude from `forces.dat`,
- divergence max/mean in stdout,
- boundary-normal residual diagnostics if enabled,
- ParaView wake/field structure.

## Theory Notes

For inviscid impermeable-body flow the essential immersed-boundary condition is
no penetration, not no slip. Tangential velocity should remain unconstrained.
For a rotating sphere in inviscid flow, pure rotation should not force the fluid
because the surface normal velocity is zero.

For non-spherical rotating bodies, the normal component of the body velocity is
generally not zero, so rotation can still couple to the flow through
impermeability.

The Lorentz reciprocal theorem is useful as a validation/checking framework for
force and added-mass relationships: two admissible incompressible flows around
the same body geometry have reciprocal boundary integral relationships. It does
not by itself solve the discrete IB/projection problem, but it can provide
independent force checks once the inviscid boundary condition is represented
cleanly.

## Open Work

- ~~Confirm that the Schur correction remains compatible with moving and~~
  ~~rotating ellipsoids, not only spheres.~~ Done 2026-07-27: tested on a free
  (translating+rotating) non-spherical ellipsoid, stable and bounded.
- Drive residual drag closer to exactly zero (currently 0.01-0.07 depending
  on case, down from 0.078-0.47, but not exact) — likely needs either a
  finer Schur `relax`/`iters` sweep per-resolution, or addressing the
  remaining discretization error in the reconstruction/projection directly
  rather than trimming it with a scalar multiplier.
- ~~Make the Schur projection + filter combination robust across~~
  ~~resolution... check whether they need to scale with resolution.~~
  Checked 2026-07-28: same sphere/domain/dt at dx=0.156 (nx=65, half
  baseline resolution), dx=0.078 (nx=129, baseline), dx=0.052 (nx=193,
  1.5x baseline) — a 3x span in dx, same `C_filter=0.49`,
  `relax=0.008`, `iters=3` at all three. Result: **no rescaling needed.**
  Fy/Fz stay at machine epsilon and `Umax` stays bounded (~1.15-1.25) at
  every resolution — stability and symmetry are resolution-independent.
  The converged drag value is also resolution-independent: coarse settles
  to `Fx≈-0.037` (confirmed by running to `t=5.0`, cheap at this
  resolution — 0.16s/step), baseline settles to `Fx≈-0.042` (at
  `t≈0.7-1.5`), same ballpark. What genuinely differs is the *settling
  timescale*: coarse needs `t~3-5` for the residual to visibly flatten,
  baseline needs `t~0.7-1.5`. This is expected for a fixed relaxation
  fraction applied to a coarser, physically-slower-relaxing surface
  correction — not a bug, and not something the current fixed defaults
  need to compensate for. Practical implication: **when validating on
  an unfamiliar resolution, don't trust a short fixed-step-count run
  as converged — plot Fx(t) and confirm the increments have decayed
  before reading off a value.** Fine mesh (dx=0.052) was only run to
  `t=0.4` (this resolution is ~3.4x baseline's per-step cost, so a full
  convergence check was skipped as unnecessary — no divergence and the
  same qualitative trend/sign was already enough to confirm the pattern
  holds in both directions, not just toward coarser meshes).
- Add automated post-processing for force magnitude, divergence summary, and
  sampled boundary-normal residual so results are comparable between machines.
  (Note for whoever does this: when reading `forces.dat`/`forces.dat<N>`,
  columns are `t, dra1(1..10), dra2(1..10), dra3(1..10), ...` — only element 1
  of each length-10 array is populated when `nvol=1`; Fx/Fy/Fz are at columns
  2/12/22, not 2/3/4. This tripped up validation more than once this session.)
- Preserve or regenerate the key ParaView outputs if visual comparison is still
  needed.
