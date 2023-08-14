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

The daily PET is rechunked to match the chunk size of the ERA5cube and added to it.
```
sbatch rechunk_pet.slurm
```

4. Consolidate metadata of ERA5Cube

Metadata from the different variables in the ERA5cube are consolidated, so as to reduce the number of read operations on the backend store.
```
source consolidate_metadata.sh
```

5. Compute PEI

The precipitation evaporation index (PEI) is a moving average of the water balance between daily potential evapotranspiration and precipitation. The moving window is 30, 90 or 180 days.

```
sbatch compute_pei.slurm
```

## Temporal ranking and spatial smoothing

The time series of the four indiactors: t2mmax, PEI_30, PEI_90 and PEI_180 are rescaled between 0 and 1. A convolutional spatial filter is run on the results to smoothe the edges of the extreme events.

```
sbatch smooth_events.slurm
```