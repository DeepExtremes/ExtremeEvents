#!/bin/bash
#SBATCH --job-name julia_qdoy
#SBATCH -o ./detect-%A.out
#SBATCH -p work
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20GB
#SBATCH --time=24:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/qdoy.jl