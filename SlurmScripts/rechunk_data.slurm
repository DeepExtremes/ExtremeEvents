#!/bin/bash
#SBATCH --job-name rchk_mn
#SBATCH -o ./rechunk-mn-%A.out
#SBATCH -p work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/Rechunk_Data.jl