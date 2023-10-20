#!/bin/bash
#SBATCH --job-name jl_st
#SBATCH -o ./testjob-%A.out
#SBATCH -p work
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/stats_extremes.jl