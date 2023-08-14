# Workflow for DeepExtremes large scale event detection v3

## Preprocessing

1. Rechunk ERA5 data, agregating from hourly to daily: t2m (min, mean, max), tp (sum), ssrd (sum).

```
sbatch rechunk_data.slurm
```

2. Compute PET on hourly data (t2m, tp , u10, v10, sp, snr (calc), vpd_cf (calc)) and aggregate to daily.

```
sbatch compute_pet.slurm
```

3. Rechunk PET and add to ER5cube
```
sbatch rechunk_pet.slurm
```

4. Consolidate metadata of ERA5Cube

```
source consolidate_metadata.sh
```

5. Compute SPEI.
