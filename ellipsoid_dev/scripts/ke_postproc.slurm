#!/bin/bash
#SBATCH --job-name=ke-postproc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
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
module load cray-python   # or: module load python; adjust to what's available

cd "$SLURM_SUBMIT_DIR"

python3 -u compute_ke.py \
    . \
    --input input.i3d \
    --chunk 100000000 \
    --out ke.csv
