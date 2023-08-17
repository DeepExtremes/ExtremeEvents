#!/bin/bash
#SBATCH --job-name julia_ind
#SBATCH -o ./ind-%A.out
#SBATCH -p work
#SBATCH --ntasks=5
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=24:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ --procs=${SLURM_NTASKS} --heap-size-hint=8GB ../scripts/plot_ind_annual.jl