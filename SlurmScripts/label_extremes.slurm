#!/bin/bash
#SBATCH --job-name jl_lbl
#SBATCH -o ./label-%A.out
#SBATCH -p big
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=60GB
#SBATCH --time=120:00:00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/label_events.jl