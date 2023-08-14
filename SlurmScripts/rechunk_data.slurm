#!/bin/bash
#SBATCH --job-name jl_rchk
#SBATCH -o ./rechunk-%A.out
#SBATCH -p work
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ --procs=${SLURM_NTASKS} --heap-size-hint=25GB ../scripts/Rechunk_Data.jl