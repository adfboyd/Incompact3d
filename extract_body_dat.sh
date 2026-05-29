#!/bin/bash
# Reconstruct body.dat{N} files from xcompact3d stdout log files.
# Supply log files in chronological order; duplicate timesteps at restart
# boundaries are dropped automatically (first occurrence wins).
#
# Usage:
#   ./extract_body_dat.sh [--nbody N] run1.log [run2.log ...]
#
# Output: body.dat1 ... body.dat{N} written to the current directory.
# Default nbody=2; override with --nbody if you have more/fewer bodies.

set -euo pipefail

nbody=2
logs=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --nbody) nbody="$2"; shift 2 ;;
        *) logs+=("$1"); shift ;;
    esac
done

if [[ ${#logs[@]} -eq 0 ]]; then
    echo "Usage: $0 [--nbody N] log1.log [log2.log ...]" >&2
    exit 1
fi

# Pre-clear output files so a re-run doesn't append to stale data
for ((i = 1; i <= nbody; i++)); do
    : > "body.dat${i}"
done

# Parse log files.
#
# body.dat column order (from Case-Ellipsoid.f90 write statement):
#   t  x y z  q1 q2 q3 q4  vx vy vz  wx wy wz  Fx Fy Fz  Tx Ty Tz  eek
#
# Note: the write statement uses angularVelocity(i,2:4), skipping index 1.
# The log prints all 4 angular velocity components; $5 $6 $7 give indices 2-4.
# eek (fluid kinetic energy, masked by ep1) is printed once per timestep
# before the per-body blocks; the same value is written to all body.dat files.
#
# Field positions in each log line after leading-whitespace stripping by awk:
#   "Time step = N/ M, Time unit = T"  →  $NF = T
#   "Kinetic Energy = eek"              →  $NF = eek
#   "Body  N"                           →  $1="Body"  $2=N  NF==2
#   "Position = x y z"                 →  $3 $4 $5
#   "Orientation = q1 q2 q3 q4"        →  $3 $4 $5 $6
#   "Linear velocity = vx vy vz"       →  $4 $5 $6
#   "Angular velocity = w1 w2 w3 w4"   →  $5 $6 $7   (skip w1 = index 1)
#   "Linear Force = Fx Fy Fz"          →  $4 $5 $6
#   "Torque = Tx Ty Tz"                →  $3 $4 $5

awk -v nbody="$nbody" '
/Time step =.*Time unit =/ {
    current_t = $NF + 0
    in_body = 0
}

$1 == "Kinetic" && $2 == "Energy" {
    eek = $NF + 0
}

$1 == "Body" && NF == 2 {
    in_body = $2 + 0
}

in_body > 0 && $1 == "Position" {
    x = $3; y = $4; z = $5
}

in_body > 0 && $1 == "Orientation" {
    q1 = $3; q2 = $4; q3 = $5; q4 = $6
}

in_body > 0 && $1 == "Linear" && $2 == "velocity" {
    vx = $4; vy = $5; vz = $6
}

in_body > 0 && $1 == "Angular" && $2 == "velocity" {
    wx = $5; wy = $6; wz = $7
}

in_body > 0 && $1 == "Linear" && $2 == "Force" {
    fx = $4; fy = $5; fz = $6
}

in_body > 0 && $1 == "Torque" {
    tx = $3; ty = $4; tz = $5
    if (in_body <= nbody && current_t > last_t[in_body] + 1e-14) {
        print current_t, x, y, z, q1, q2, q3, q4, vx, vy, vz, wx, wy, wz, fx, fy, fz, tx, ty, tz, eek \
            > ("body.dat" in_body)
        last_t[in_body] = current_t
    }
    in_body = 0
}
' "${logs[@]}"

echo "Done. Written body.dat1 ... body.dat${nbody}" >&2
