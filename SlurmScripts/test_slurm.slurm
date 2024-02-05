#!/bin/bash
#SBATCH --job-name jl_tst
#SBATCH -o ./test-%A.out
#SBATCH -p work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --time=120:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/test_slurm.jl