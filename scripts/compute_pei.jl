# The precipitation evaporation index (PEI) is a moving average of the water balance between daily 
# potential evapotranspiration and precipitation.
# The moving wiindow is 30, 90 or 180 days.
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
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
ds = open_dataset(zg)
pet = ds.pet
# seems that missing are considered missing even if fill_as_missing = false

zg1 = zopen("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/ERA5Cube.zarr")
ds1 = open_dataset(zg1)
tp = ds1.tp

# tp is in m/day, while pet is in mm/day
# downward fluxes are >0
diffcube = map((i,j)-> i*1e3+j,tp,pet)

outpath="/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/PEICube.zarr"


windowsizes = [30,90,180]
windowax = Dim{:Variable}(map(ws->string("pei_",ws),windowsizes))
indims = InDims("time")
outdims = OutDims("time",windowax,
    chunksize = :input, 
    path = outpath,
    overwrite=true,
)

@everywhere function compute_pei(xout, xin, windows)
    for i in eachindex(windows)
        xout[:,i] = runmean(xin, windows[i])
    end
end

pei = mapCube(compute_pei,diffcube,windowsizes; indims, outdims, max_cache=1e9, showprog = true)

# consolidate_metadata
using CondaPkg; CondaPkg.add("xarray"); CondaPkg.add("zarr")
using  PythonCall
zr = pyimport("zarr")
g = zr.open_group(outpath)
g.consolidate_metadata()
