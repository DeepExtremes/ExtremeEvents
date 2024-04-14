#!/bin/bash
#SBATCH --job-name pei_jl
#SBATCH -o ./pei-%A.out
#SBATCH -p work
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/compute_pei.jl