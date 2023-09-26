using SlurmClusterManager, Distributed
#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere begin
    using EarthDataLab, YAXArrays, Zarr, DiskArrays
end

using YAXArrayBase

outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr"

# need to find a way to 
# - rename axis Time to time
# - rename cubes from :layer to :pet
map(1950:2022) do yr
    mv("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/layer", "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/pet")
    # mv("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/Time", "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/time")
    # mv("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/time", "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr/Time")
end
ds = YAXArrays.Datasets.open_mfdataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/*.zarr")
# # ERROR: Dimension ranges overlap
# # check cubes individually and identify which ones are problematic
# for yr in 1950:2022
#     println(yr)
#     @show ds = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr")
# end
## problem comes from the definition of the time axis. it wants Time and doesn't recognise time. So revert rename Time axis to time.

## Renaming ds time axis
# renameaxis!(ds, :Time => :time)
# ERROR: MethodError: no method matching renameaxis!(::Dataset, ::Pair{Symbol, Symbol})
# Closest candidates are:
# renameaxis!(::YAXArray, ::Pair)
## Do it at a still later stage
ds = setchunks(ds,(lon=60,lat=60,time=5844))

savedataset(ds, path=outpath,append=true,max_cache=1e9)
println("Done!")

# rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/Time
# rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/layer



# consolidate metadata