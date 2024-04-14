#!/bin/bash
#SBATCH --job-name jl_18_bck
#SBATCH -o ./ind18bk-%A.out
#SBATCH -p work
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=24:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/plot_ind_annual_backext.jl