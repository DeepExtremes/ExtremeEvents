#!/bin/bash
#SBATCH --job-name julia
#SBATCH -o ./testjob-%A.out
#SBATCH -p work
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40GB

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../Rechunk_Fluxcom.jl