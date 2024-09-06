# annual number of events by type and by grid cell (lat, lon)
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

@everywhere using Zarr, YAXArrays, DimensionalData, Dates

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

# Event Cube
trial = "ranked_pot0.01_ne0.1"
eventcube = Cube(path * "EventCube_$trial.zarr")#[time = Date(2001) .. Date(2010), latitude = 45.0 .. 46.0, longitude = 15.0 .. 16.0]
outpath = "$(path)AnnualEventCount_$trial.zarr"

@everywhere begin
# number of events by type
function net(x)
    nhe = sum(skipmissing(x .& UInt8(1)))
    nde30 = sum(skipmissing((x .& UInt8(2)) .== UInt8(2)))
    nde90 = sum(skipmissing((x .& UInt8(4)) .== UInt8(4)))
    nde180 = sum(skipmissing( (x .& UInt8(8)) .== UInt8(8)))
    ndhe = sum(skipmissing((x .> 1) .& .!(iseven.(x))))
    return [nhe, nde30, nde90, nde180, ndhe]
end

# number of events by type and by year
function net!(xout, xin, indices)
    ne = map(indices) do idxi
        net(xin[idxi])
    end
    # xout has size (numberofyears, numberofeventtypes)
    xout[:] = stack(ne, dims=1)[:]
    return xout
end
end # begin

ti = lookup(eventcube, :Ti, )

# create indices
indices = map(Dates.year(ti[1]) : Dates.year(ti[end])) do y
    findall(Dates.year.(ti) .== y)
end

# Input dimension : Time
indims = InDims(:Ti)
# Output chunking by map
outdims = OutDims( 
    Dim{:year}(Dates.year(ti[1]) : Dates.year(ti[end])), 
    Dim{:net}(["h", "d30", "d90", "d180", "dh"]);
        outtype = Int16, # maximum value is 366 if all days in the year are extremes
        chunksize = Dict("longitude" => 1440, "latitude" => 421, "year" => 1, "net" => 5), # one chunk is 6kB
        path = outpath,
        overwrite = true,
        backend = :zarr,
)

@time aec = mapCube(net!, eventcube, indices; indims = indims, outdims = outdims)
# 1151.221801 seconds (21.95 M allocations: 1.460 GiB, 0.08% gc time, 1.33% compilation time)

@show aec 
# aec = 73×5×1440×721 DiskArrayTools.CFDiskArray{Union{Missing, Int16}, 4, Int16, ZArray{Int16, 4, Zarr.BloscCompressor, DirectoryStore}, Int16}


