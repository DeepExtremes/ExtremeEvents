@show (startyear, endyear) = (parse(Int, ARGS[1]), parse(Int, ARGS[2]))

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
    using YAXArrays, Zarr, DiskArrays
end

using YAXArrayBase

outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube_$startyear.zarr"

# input

filelist = [joinpath("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/", "$yr.zarr") for yr in startyear:endyear]

ds = YAXArrays.Datasets.open_mfdataset(filelist)
ds_chunked = setchunks(ds,(lon=60,lat=60,time=5844))
# add missing properties
era = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr")
ds_chunked.pet.properties["units"] = era.pet.properties["units"]
ds_chunked.pet.properties["long_name"] = era.pet.properties["long_name"]
ds_chunked.pet.properties["aggfun"] = era.pet.properties["aggfun"]
# rename time axis
renameaxis!(ds_chunked.pet, :Ti => :time)

savedataset(ds_chunked, path=outpath,append=true,max_cache=1e9)
println("Done!")


# consolidate metadata



