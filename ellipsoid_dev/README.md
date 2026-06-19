# Ellipsoid Development Checks

This directory holds local validation cases and helper scripts for the ellipsoid
IBM development work.

- `inputs/`: input files for reproducible ellipsoid, inviscid-slip, and viscous
  comparison checks.
- `scripts/`: local post-processing or batch helper scripts.
- `runs/`: scratch run directories. Contents are ignored by git.

Recommended run pattern:

```sh
./ellipsoid_dev/scripts/new_run.sh rotating_sphere input_sphere_inviscid_rotating_test.i3d
cd ellipsoid_dev/runs/rotating_sphere
./run_local.sh
```

Keeping runs under `ellipsoid_dev/runs/` prevents generated files such as
`data/`, `out/`, `forces.dat*`, and `torques.dat*` from cluttering the project
root.
