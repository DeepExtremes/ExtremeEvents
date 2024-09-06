# The precipitation evaporation index (PEI) is a moving average of the water balance between daily 
# potential evapotranspiration and precipitation.
# The moving window is 30, 90 or 180 days.

@show (startyear, endyear) = (parse(Int, ARGS[1]), parse(Int, ARGS[2]))

using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
        addprocs(1)
        sleep(0.001)
    end
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using YAXArrays, Zarr, RollingFunctions, DiskArrays

zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
ds = open_dataset(zg)[time = Date(startyear) .. Date(endyear,12,31,)]
# seems that missing are considered missing even if fill_as_missing = false

# tp is in m/day, while pet is in mm/day
# downward fluxes are >0
diffcube = map((i,j)-> i*1e3+j,ds.tp,ds.pet)

windowsizes = [30,90,180]
windowax = Dim{:Variable}(map(ws->string("pei_",ws),windowsizes))
indims = InDims("time")
outdims = OutDims("time",windowax,
    chunksize = :input, 
    path="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube_$startyear.zarr",
    overwrite=true,
)

@everywhere function compute_pei(xout, xin, windows)
    for i in eachindex(windows)
        xout[:,i] = runmean(xin, windows[i])
    end
end

pei = mapCube(compute_pei,diffcube,windowsizes; indims, outdims, max_cache=1e9, showprog = true)

# python xarray.to_zarr
# outdims are now :Ti instead of :time => trouble to add data with xr.to_zarr
# => manually edit .zattrs for Ti, pei_xx + rename Ti to time.
