#!/bin/bash
#SBATCH --job-name=pvrender
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=4
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --account=e643           # <-- your ARCHER2 budget code
#SBATCH --hint=nomultithread
#SBATCH --distribution=block:block
#SBATCH --array=0-15%16
#SBATCH --output=pvrender-%j.out
#SBATCH --error=pvrender-%j.err

# Render snapshots with paraview_render.py
#
# Submit with:   sbatch render.slurm
# Dry-run check: sbatch --test-only render.slurm
#
# For the short QoS instead, replace the two #SBATCH lines above with:
#   #SBATCH --partition=standard
#   #SBATCH --qos=short
# (short: <=20 min, <=32 nodes; serial: longer walltime, no MPI)

set -euo pipefail

module load paraview

# Run inside the directory containing the snapshots.
cd "${SLURM_SUBMIT_DIR}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
mkdir -p frames

# Round-robin shard the snapshot list across array tasks.
ALL_SNAPS=( data/snapshot-*.xdmf )
NSNAPS=${#ALL_SNAPS[@]}
NTASKS=${SLURM_ARRAY_TASK_COUNT}
TASK_ID=${SLURM_ARRAY_TASK_ID}

MY_SNAPS=()
for i in "${!ALL_SNAPS[@]}"; do
    if (( i % NTASKS == TASK_ID )); then
        MY_SNAPS+=( "${ALL_SNAPS[$i]}" )
    fi
done

echo "Array task ${TASK_ID}/${NTASKS}: ${#MY_SNAPS[@]} of ${NSNAPS} snapshots"
if (( ${#MY_SNAPS[@]} == 0 )); then
    echo "Nothing to do for this task."
    exit 0
fi

srun --unbuffered \
    pvbatch paraview_render.py \
        "${MY_SNAPS[@]}" \
        --outdir frames \
        --rotate-offset ${TASK_ID} \
        --rotate-stride ${NTASKS} \
        --total-frames ${NSNAPS}
       

echo "Task ${TASK_ID} done."