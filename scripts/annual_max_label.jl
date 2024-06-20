# annual number of events by type and by grid cell (lat, lon)
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:haskey(ENV, "SLURM_NTASKS")
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
    patho = "https://s3.bgc-jena.mpg.de:9000/deepextremes/v3/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
    patho = "/Net/Groups/BGI/work_1/scratch/s3/deepextremes/v3/"
end

# Label Cube
labelcube = Cube("$(patho)mergedlabels.zarr/")#[time = Date(2001) .. Date(2010), latitude = 45.0 .. 46.0, longitude = 15.0 .. 16.0]
outpath = "$(path)AnnualMaxLabel_cmp_S1_T3.zarr"

@everywhere begin
# max label by year
function mxlb!(xout, xin, indices)
    ml = map(indices) do idxi
        maximum(skipmissing(xin[idxi]))
    end
    # xout has size (numberofyears, numberofeventtypes)
    xout[:] = stack(ml, dims=1)[:]
    return xout
end
end # begin

ti = lookup(labelcube, :Ti, )

# create indices
indices = map(Dates.year(ti[1]) : Dates.year(ti[end])) do y
    findall(Dates.year.(ti) .== y)
end

# Input dimension : Time
indims = InDims(:Ti)
# Output chunking by map
outdims = OutDims( 
    Dim{:year}(Dates.year(ti[1]) : Dates.year(ti[end])), 
        # outtype = Int64, # same as input
        # chunksize = Dict("longitude" => 1440, "latitude" => 421, "year" => 1), 
        path = outpath,
        overwrite = true,
        backend = :zarr,
)

@time aml = mapCube(mxlb!, labelcube, indices; indims = indims, outdims = outdims)
# 2103.752381 seconds (21.82 M allocations: 1.450 GiB, 0.05% gc time, 0.79% compilation time)


@show aml 
# aml = 73×1440×721 DiskArrayTools.CFDiskArray{Union{Missing, Int64}, 3, Int64, ZArray{Int64, 3, Zarr.BloscCompressor, DirectoryStore}, Int64}



