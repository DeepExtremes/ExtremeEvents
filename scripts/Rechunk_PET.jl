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

ds = YAXArrays.Datasets.open_mfdataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/*.zarr")
ds = setchunks(ds,(lon=60,lat=60,time=5844))
# need to find a way to rename cubes from :layer to :pet

savedataset(ds, path=outpath,append=true,max_cache=1e9)