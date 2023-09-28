#!/bin/bash
#SBATCH --job-name jl_lbl
#SBATCH -o ./labelF-%A.out
#SBATCH -p work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=400GB
#SBATCH --time=24:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/label_events.jl