#!/bin/bash
#SBATCH --job-name pei_jl
#SBATCH -o ./pei-%A.out
#SBATCH -p work
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10GB
#SBATCH --time=48:00:00
#SBATCH --mail-type FAIL,END

export JULIA_NUM_THREADS=${SLURM_CPUS_PER_TASK}

julia --project=../ ../scripts/compute_pei_add.jl 2023 2023

OLDCUBE="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr"
NEWCUBE="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube_2023.zarr"
srun /Net/Groups/BGI/scratch/mweynants/miniconda3/bin/python -c "import add2cube as ac;ac.add2cube(\"${OLDCUBE}\", \"${NEWCUBE}\")"

echo "Success!"