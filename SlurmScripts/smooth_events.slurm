#!/bin/bash
#SBATCH --job-name julia_smooth
#SBATCH -o ./testjob-%A.out
#SBATCH -p work
#SBATCH --ntasks=20
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=12:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ --procs=${SLURM_NTASKS} ../scripts/smooth_events.jl