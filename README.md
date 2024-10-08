# Workflow for building and analysing Dheed v3
Dheed v3 is an ERA5 based global database dry and hot extreme events from 1950 to 2022, developed in the context of ESA funded project [DeepExtremes](https://eo4society.esa.int/projects/deep-extremes/).
The workflow runs in Julia 1.10.0, except for the consolidation of the data cubes, which is run in python. 
Most steps of the workflow were run on the [MPI BGC-jena](https://bgc-jena.mpg.de) cluster. The input data are hourly ERA5 data retrieved from the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/) and stored on a local server as netcdf files. Some derived variables were calculated prior to the processing presented here.

## Preprocessing

1. Rechunk ERA5 data, agregating from hourly to daily: t2m (min, mean, max), tp (sum), ssrd (sum).

```
cd SlurmScripts
sbatch rechunk_data.slurm
```

*Output*: ERA5Cube.zarr

2. Compute PET on hourly data (t2m, tp , u10, v10, sp, snr (calc), vpd_cf (calc)) and aggregate to daily.

```
cd SlurmScripts
sbatch compute_pet.slurm
```

*Output*: PET yearly cubes

3. Rechunk PET and add to ER5cube

The daily PET is rechunked to match the chunk size of the ERA5cube and added to it.

```
cd SlurmScripts
sbatch rechunk_pet.slurm
```

*Output*: updated ERA5Cube.zarr

4. Consolidate metadata of ERA5Cube

Metadata from the different variables in the ERA5cube are consolidated, so as to reduce the number of read operations on the backend store.

```
python -c import zarr; g = zarr.open_group(path2cube); zarr.consolidate_metadata(g.store)
```

*Output*: consolidated ERA5Cube.zarr

5. Compute PEI

The precipitation evaporation index (PEI) is a moving average of the water balance between daily potential evapotranspiration and precipitation. The moving window is 30, 90 or 180 days.

```
cd SlurmScripts
sbatch compute_pei.slurm
```

*Output*: PEICube.zarr


## Processing
### Temporal analysis

The time series of the four indicators: t2mmax, PEI_30, PEI_90 and PEI_180 are ranked transformed between 0 and 1. No convolutional spatial filter is run on the results to smoothe the extent of the extreme events.

```
cd SlurmScripts
sbatch smooth_events.slurm
```

*Output*: tmax_ranked.zarr and pei_ranks.zarr

### Compute extremes
A pass over threshold is applied to the rescaled indicators and they are combined into a Byte integer (Int8), with one bit for each indicator and an extra bit encoding for non extremes. The first bit (little end) encodes the maximum temperature extremes.

```
cd SlurmScripts
sbatch compute_events.slurm
```

*Output*: EventCube.zarr

### Label extreme events
Unique labels are assigned to blobs of co-occurrent hot and dry extremes, i.e. where values are uneven (t2mmax extremes) and larger than one (PEI extremes), connected in space and time. A filter is applied before runnning the connected component analysis: temperature extremes must last at least three consecutive days.

*Note*: this will reduce the total number of tmax extremes in the cube...

Because the connected component analysis requires to load the full cube into memory and the algorithm is greedy, the analysis was split into seven tasks, covering each 13 years, with three years overlap between two successive periods.

```
cd SlurmScripts
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

The labels are then merged into a single mergedlabels cube with `scripts/merge_labels.jl`.

*Output*: 
- mergedlabels.zarr

### Compute statistics
Statistics for all labelled events are computed and gathered in a unique table.

```
cd SlurmScripts
sbatch stats_extremes_merged.slurm
```

*Output*: 
- MergedEventStats_landonly.csv

## Figures and Postprocessing

### fig01_workflow.png

Flowchart designed in ppt. 

### fig02_workflow-plot-33.png

Example of dry and hot extreme event detection workflow over the 2003 summer heatwave in Europe.

See `scripts/fig4dheed.jl`

### fig03: timeseries

Timeseries of indicators for 2 contrasted locations: Jena (timeseries_11.59_50.92_2012_2023.png) and Niamey (timeseries_2.1254_13.5116_2012_2023.png).

See `scripts/plot_t2mmax.jl`

### fig04 - fig06: Trends -- Annual global/continental summary

Trends in annual global/continental indicators.

Compute annual statistics of indicators, globally and by continent.

See `scripts/plot_ind_annual.jl`

Compute annual statistics of EventCube at global and continental scale, by summing all extremes by type.

See `scripts/hist_EventType.jl`

Plot the results of the EventCube analysis.

```
julia --project="ExtremeEvents.toml" plot_EventType.jl
```

### fig07: events_stats_ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022_landonly_1970.png

Extract largest and longest events from `MergedEventStats_landonly` with `largest_labels.jl` and plot statistics with `plot_stats.jl`.

### fig08: largest_ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022_landonly_1970.png

A map of the spatial footprint of the largest events is generated with `plot_stats.jl`.

### fig09: Validation

Compare MergedEventStats_landonly with table of reported events compiled *a priori* with `SanityCheck.jl`. 

### fig10: Lytton.png

Extract timeseries at single locations with `plot_city.jl`

