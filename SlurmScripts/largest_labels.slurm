#!/bin/bash
#SBATCH --job-name jl_llbls
#SBATCH -o ./llbls-%A.out
#SBATCH --mail-type FAIL,END
#SBATCH -p work
#SBATCH --ntasks=8
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=7GB
#SBATCH --time=5-00

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/largest_labels.jl