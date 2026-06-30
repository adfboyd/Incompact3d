# Inviscid Ellipsoid Handover

Date: 2026-06-30

This branch contains experimental work for inviscid ellipsoid immersed-boundary
treatment. Keep it separate from `ellipsoid-dev`, which is the validated viscous
ellipsoid branch.

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

- Make the Schur projection robust across resolution, time step, and body shape.
- Confirm that the Schur correction remains compatible with moving and rotating
  ellipsoids, not only spheres.
- Add automated post-processing for force magnitude, divergence summary, and
  sampled boundary-normal residual so results are comparable between machines.
- Preserve or regenerate the key ParaView outputs if visual comparison is still
  needed.
