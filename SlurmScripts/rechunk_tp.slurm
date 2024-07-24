#!/bin/bash
#SBATCH --job-name rchk_tp
#SBATCH -o ./rechunk-htp-%A.out
#SBATCH --partition work
#SBATCH --mail-type FAIL,END
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=30GB
#SBATCH --time=72:00:00
#SBATCH --nodes=1

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/Rechunk_htp.jl