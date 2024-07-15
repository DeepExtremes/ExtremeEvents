# Workflow for DeepExtremes large scale event detection v3

## Preprocessing

1. Rechunk ERA5 data, agregating from hourly to daily: t2m (min, mean, max), tp (sum), ssrd (sum).

```
sbatch rechunk_data.slurm
```

*Output*: ERA5Cube.zarr

2. Compute PET on hourly data (t2m, tp , u10, v10, sp, snr (calc), vpd_cf (calc)) and aggregate to daily.

```
sbatch compute_pet.slurm
```

*Output*: PET yearly cubes

3. Rechunk PET and add to ER5cube

The daily PET is rechunked to match the chunk size of the ERA5cube and added to it.
```
sbatch rechunk_pet.slurm
```

*Output*: updated ERA5Cube.zarr

4. Consolidate metadata of ERA5Cube

Metadata from the different variables in the ERA5cube are consolidated, so as to reduce the number of read operations on the backend store.
```
source consolidate_metadata.sh
```

*Output*: consolidated ERA5Cube.zarr

5. Compute PEI

The precipitation evaporation index (PEI) is a moving average of the water balance between daily potential evapotranspiration and precipitation. The moving window is 30, 90 or 180 days.

```
sbatch compute_pei.slurm
```

*Output*: PEICube.zarr


## Processing
### Temporal rescaling and spatial smoothing

The time series of the four indicators: t2mmax, PEI_30, PEI_90 and PEI_180 are rescaled between 0 and 1. No convolutional spatial filter is run on the results to smoothe the extent of the extreme events.

```
sbatch smooth_events.slurm
```

*Output*: tmax_ranked.zarr and pei_ranks.zarr

### Compute extremes
A pass over threshold is applied to the rescaled indicators and they are combined into a Byte integer (Int8), with one bit for each indicator and an extra bit encoding for non extremes. The smallest bit (little endian) encodes the maximum temperature extremes.

```
sbatch compute_events.slurm
```

*Output*: EventCube.zarr

### Label extreme events
Unique labels are assigned to blobs of co-occurrent hot and dry extremes, i.e. where values are uneven (t2mmax extremes) and larger than one (PEI extremes), connected in space and time. A filter is applied before runnning the connected component analysis: temperature extremes must last at least three consecutive days.

*Note*: this will reduce the total number of tmax extremes in the cube...

Because the connected component analysis requires to load the full cube into memory and the algorithm is greedy, the analysis was split into seven tasks, covering each 13 years, with three years overlap between two successive periods.

```
sbatch label_events.slurm
```

*Output*: 
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_1950_1962.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_1960_1972.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_1970_1982.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_1980_1992.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_1990_2002.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_2000_2012.zarr
- labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022.zarr

The labels are then merged into a single mergedlabels cube.


### Compute statistics
Statistics for all labelled events are computed and gathered in a unique table.

```
sbatch stats_extremes.slurm
```

*Output*: 
- MergedEventStats_landonly.csv

```julia
plot_stats.jl
```
Where are the labelled events in the period 2010-2022?

<img src="../v3/fig/nolabel_ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022.png" alt="plot_event" width="800"/>

![alt text](../v3/fig/nolabel_ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022.png)

## Postprocessing

### Sanity check


### Annual continental summary

