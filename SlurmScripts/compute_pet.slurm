#!/bin/bash
#SBATCH --job-name pet2_jl18
#SBATCH -o ./pet22-%A.out
#SBATCH -p work
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=48:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../  ../scripts/compute_pet.jl