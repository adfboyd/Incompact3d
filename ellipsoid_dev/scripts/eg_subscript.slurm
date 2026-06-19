#!/bin/bash
#SBATCH --job-name=x3d-ellipsoids
#SBATCH --nodes=32
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --account=e643
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -e
set -o pipefail

module load PrgEnv-gnu
module load cray-mpich

cd "$SLURM_SUBMIT_DIR"

# ---- Lustre striping ----
# data/ is a directory of snapshot files — stripe the dir; new files inherit
mkdir -p data
lfs setstripe -c -1 -S 128M data 2>/dev/null || true

# checkpoint / restart-forces are FILES in MPIIO mode — pre-stripe if absent
[ -e checkpoint    ] || lfs setstripe -c -1 -S 128M ./checkpoint
[ -e restart-forces ] || lfs setstripe -c -1 -S 128M ./restart-forces

# Warn if an existing file has wrong stripe (doesn't auto-fix — user's call)
for f in checkpoint restart-forces; do
  if [ -e "$f" ]; then
    s=$(lfs getstripe -S "$f" 2>/dev/null)
    if [ "$s" != "134217728" ]; then
      echo "WARNING: $f stripe size is $s, expected 134217728 (128M). Re-stripe manually if intended."
    fi
  fi
done

# ---- MPI-IO hints ----
export MPICH_MPIIO_HINTS="\
checkpoint:romio_cb_write=enable:cb_nodes=$SLURM_NNODES:cray_cb_nodes_multiplier=1:cray_cb_write_lock_mode=1:cb_buffer_size=134217728,\
restart-forces:romio_cb_write=disable:romio_ds_write=disable,\
data/*:romio_cb_write=disable:romio_ds_write=disable"
export MPICH_MPIIO_HINTS_DISPLAY=1
export MPICH_MPIIO_STATS=1
echo "=== MPICH_MPIIO_HINTS = $MPICH_MPIIO_HINTS"

EXEC=/work/e643/e643/shared/Incompact3d_git/Incompact3d/build/bin/xcompact3d
RUN_DIR=$SLURM_SUBMIT_DIR


cd $RUN_DIR

srun $EXEC input.i3d