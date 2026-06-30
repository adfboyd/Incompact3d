# Running Xcompact3d

This note gives the practical steps for building and running this repository.
It also records the extra setup needed for the ellipsoid immersed-boundary case.

The examples below assume commands are run from the repository root unless a
different working directory is shown.

## Build

Requirements:

- CMake 3.20 or newer.
- A modern Fortran compiler.
- MPI and a matching MPI Fortran wrapper, for example `mpifort` or `mpif90`.

Configure and build:

```sh
export FC=mpifort
cmake -S . -B build -DCMAKE_BUILD_TYPE=RELEASE
cmake --build build -j 8
```

The executable is written to:

```text
build/bin/xcompact3d
```

If CMake finds the wrong MPI launcher, configure with:

```sh
cmake -S . -B build -DMPIEXEC_EXECUTABLE=/path/to/mpirun
```

To run the default tests:

```sh
ctest --test-dir build --output-on-failure
```

## Run A Case

Use a separate run directory for each case. The solver writes output into the
current working directory, so running directly in `examples/` or the repository
root will mix generated files with source files.

Example:

```sh
mkdir -p runs/my_case
cp examples/TGV-Taylor-Green-vortex/input.i3d runs/my_case/input.i3d
cd runs/my_case
mpirun --oversubscribe -np 8 ../../build/bin/xcompact3d input.i3d
```

On an HPC system, replace the local `mpirun` command with the site launcher
requested by the job scheduler.

The input file is a Fortran namelist file. Important general controls are:

- `nx`, `ny`, `nz`: grid size.
- `xlx`, `yly`, `zlz`: domain size.
- `nclx1`, `nclxn`, `ncly1`, `nclyn`, `nclz1`, `nclzn`: velocity boundary
  conditions. `0` is periodic, `1` is free-slip, and `2` is Dirichlet.
- `p_row`, `p_col`: MPI pencil decomposition. If both are `0`, the code chooses
  a decomposition automatically. Otherwise `p_row*p_col` must equal the number
  of MPI ranks.
- `dt`, `ifirst`, `ilast`: time-step size and iteration range.
- `ioutput`: snapshot interval.
- `nvisu`: visualization downsampling stride.
- `icheckpoint`: checkpoint interval.
- `irestart`: set to `1` to restart from an existing checkpoint.

For non-periodic directions, remember that the grid convention differs from
periodic directions. In practice this often means adding one point, for example
using `nx=129` where the equivalent periodic case would use `nx=128`.

## Output Files

The code creates output relative to the run directory:

- `data/`: visualization snapshots. With the default MPI-IO/XDMF output, open
  files such as `data/snapshot-0.xdmf` in ParaView.
- `out/`: case-specific text output.
- `probes/`: probe output if probes are enabled.
- `restart*` and `restart.info`: checkpoint/restart state.
- `forces.dat*`: control-volume force output when `iforces=1`.

Use `ioutput` and `nvisu` to control ParaView output volume. For a quick
visualization run, use a small `ioutput` and `nvisu=1`. For long production
runs, increase `ioutput` or set it beyond `ilast` if snapshots are not needed.

## Restarting

To restart a run:

1. Keep the previous run directory, including `restart*`, `restart.info`, and
   any case-specific restart files.
2. Set `irestart=1`.
3. Set `ifirst` and `ilast` for the continuation interval.
4. Run from the same directory, or copy all restart files into a new run
   directory with the updated `input.i3d`.

Ellipsoid runs also write `body_state.dat`. That file is required on restart
and the `nbody` value in `input.i3d` must match the value stored in
`body_state.dat`.

## Ellipsoid Case

The ellipsoid case is selected with:

```fortran
&BasicParam
itype = 16
iibm = 3
...
/End
```

The committed example is:

```text
examples/Ellipsoid/input.i3d
```

Run it from a separate directory:

```sh
mkdir -p runs/ellipsoid_example
cp examples/Ellipsoid/input.i3d runs/ellipsoid_example/input.i3d
cd runs/ellipsoid_example
mpirun --oversubscribe -np 8 ../../build/bin/xcompact3d input.i3d
```

This example is a viscous moving-ellipsoid case. It uses the standard IBM
no-slip treatment.

### Ellipsoid Geometry And Motion

Ellipsoid setup is in the `&ibmstuff` namelist:

```fortran
&ibmstuff
imove=1
ce=1.5,1.75,2.0
sh=0.36,0.30,0.246
ori=1.0,0.0,0.0,0.0
lv=1.0,0.0,0.0
av=1.77,1.77,0.0
rho_s=4.0
ra=0.15
nraf=10
nobjmax=1
iforces=1
nvol=1
bodies_fixed=0
torques_flag=1
orientations_free=1
torq_flip=1
cvl_scalar=1.2
nbody=1
/End
```

Key fields:

- `ce`: body centre.
- `sh`: unnormalised ellipsoid axis ratios in the body frame. The code
  normalizes these by their geometric mean to get the stored shape.
- `ra`: body size scale applied to the normalized shape. For a sphere, use
  equal `sh` components and set the radius with `ra`.
- `ori`: orientation quaternion, ordered as `w,x,y,z`.
- `lv`: linear velocity.
- `av`: angular velocity.
- `rho_s`: body density.
- `nbody`: number of bodies.
- `nobjmax`: maximum number of immersed objects used by the IBM geometry
  search. For ellipsoid cases, set it at least as large as `nbody`.
- `nraf`: refined geometry sampling factor for the immersed-boundary mask.
- `imove`: enables moving-body updates.
- `bodies_fixed`: keeps body state fixed when set.
- `orientations_free`: enables orientation updates.
- `torques_flag`: enables torque calculation and angular-velocity updates.
- `torq_flip`: sign convention switch for torque response used by the current
  ellipsoid implementation.
- `cvl_scalar`: scales the automatic control volume around each ellipsoid.

For multiple ellipsoids, supply body-major entries and set `nbody`
accordingly. For example, use `ce=x1,y1,z1,x2,y2,z2` for two bodies, and
`ori=w1,x1,y1,z1,w2,x2,y2,z2` for two orientation quaternions. Use one scalar
entry per body for `rho_s` and `ra`.

### Forces And Torques

Set:

```fortran
iforces=1
nvol=1
torques_flag=1
```

The `&ForceCVs` namelist must be present when `iforces=1`:

```fortran
&ForceCVs
xld(1) = 1.2
xrd(1) = 1.8
yld(1) = 1.45
yud(1) = 2.05
zld(1) = 1.7
zrd(1) = 2.3
/End
```

For `itype=16`, the code initializes and updates force-control-volume bounds
from the ellipsoid position and `cvl_scalar`. Keep `&ForceCVs` in the input so
the namelist read succeeds, but expect the runtime bounds to follow the body.

Output files:

- `forces.dat1`, `forces.dat2`, ...: force history for each control volume.
- `torques.dat1`, `torques.dat2`, ...: torque history when torque output is
  enabled.
- `body.dat1`, `body.dat2`, ...: body position/orientation history.
- `body_state.dat`: latest body state for restart.

### Boundary Conditions For Uniform Flow

For a periodic box, use periodic boundary flags:

```fortran
nclx1 = 0
nclxn = 0
ncly1 = 0
nclyn = 0
nclz1 = 0
nclzn = 0
```

For inflow/outflow in `x`, use non-periodic x boundary conditions:

```fortran
nclx1 = 2
nclxn = 2
```

Then increase the non-periodic grid count by one compared with the equivalent
periodic FFT-friendly size. For example, use `nx=129` rather than `nx=128`.

The ellipsoid case applies its inlet condition through `u1`, `u2`, and
`inflow_noise`, and its outflow update in `Case-Ellipsoid.f90`.

### Practical Checks

For a new ellipsoid run, check:

- The stdout banner says `Simulating Ellipsoid`.
- The requested MPI rank count is compatible with `p_row` and `p_col`.
- `data/`, `out/`, and force/body history files are written in the run
  directory, not the repository root.
- `body_state.dat` exists before attempting a restart.
- The body remains inside the domain and the control volume does not intersect
  another body or the domain boundary.
- For moving bodies, `dt` is small enough that the body does not move too far in
  one step relative to the grid and IBM interpolation width.

## Development Branches

The stable ellipsoid branch is intended for the viscous no-slip ellipsoid case.
Experimental inviscid work should remain on `inviscid-ellipsoid-dev`, whose
handover notes are in `ellipsoid_dev/INVISCID_HANDOVER.md` on that branch.
