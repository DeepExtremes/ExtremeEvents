#!/bin/bash
#SBATCH --job-name jl_rchk
#SBATCH -o ./rechunk-%A.out
#SBATCH -p work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ --heap-size-hint=7GB ../scripts/Rechunk_PET.jl