#!/bin/bash
#SBATCH --job-name pet_jl10
#SBATCH -o ./pet-%A.out
#SBATCH -p work
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=2-00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../  ../scripts/compute_pet.jl 2023