using SlurmClusterManager, Distributed

using EarthDataLab, YAXArrays, Zarr, RollingFunctions, DiskArrays

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
ds = open_dataset(zg)
# seems that missing are considered missing even if fill_as_missing = false

# tp is in m/day, while pet is in mm/day
# downward fluxes are >0
diffcube = map((i,j)-> i*1e3+j,ds.tp,ds.pet)

windowsizes = [30,90,180]
windowax = CategoricalAxis("Variable",map(ws->string("pei_",ws),windowsizes))
indims = InDims("Time")
outdims = OutDims("Time",windowax,
    chunksize = :input, 
    path="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr",
    overwrite=true,
)

@everywhere function compute_spei(xout, xin, windows)
    for i in eachindex(windows)
        xout[:,i] = runmean(xin, windows[i])
    end
end

spei = mapCube(compute_spei,diffcube,windowsizes; indims, outdims, max_cache=1e9, showprog = true)