#!/bin/bash
#SBATCH --job-name rank_jl
#SBATCH -o ./rank-%A.out
#SBATCH -p work
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/smooth_events.jl