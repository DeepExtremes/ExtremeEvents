# Compute summary statistics on labelled extrene events
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

@everywhere using Zarr, YAXArrays, OnlineStats, WeightedOnlineStats, NetCDF, DimensionalData
@everywhere using DataFrames
import Tables
import CSV

@everywhere begin

# stats run on land only
landonly = "_landonly"

# include functions to compute stats and create DataFrame
include("../src/stats.jl")

end # begin

if occursin("/Users", pwd())
    path = "https://s3.bgc-jena.mpg.de:9000/deepextremes/v3/"
    patho = "./"
else
    path = "/Net/Groups/BGI/work_1/scratch/s3/deepextremes/v3"
    patho = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

# Open all Data Cubes
pei = open_dataset("$(path)PEICube.zarr")

zg = zopen("$(path)ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)

eventspath = "$(path)EventCube_ranked_pot0.01_ne0.1.zarr"
eventcube = open_dataset(eventspath)

# LandSeaMask
lsmask = open_dataset("$(path)lsm.1440.721.static.nc")
lsmask_notime = lsmask[time = At(DateTime("2019-01-01T13:00:00"))]

# # Fluxcom carbon fluxes (gC m^(-2) year^(-1))
# # gross primary productivity
# gpp = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_GPP_final.zarr"))
# renameaxis!(gpp, "Time" => "time")
# # net ecosystem exchange
# nee = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_NEE_final.zarr"))
# renameaxis!(nee, "Time" => "time")
# # terrestrial ecosystem respiration
# ter = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_TER_final.zarr"))
# renameaxis!(ter, "Time" => "time")

# ranks
pei_ranks = open_dataset("$(path)pei_ranks.zarr")
tmax_ranks = open_dataset("$(path)tmax_ranked.zarr")

labelpath = "$(path)mergedlabels.zarr/"
labels = open_dataset(labelpath)

ti = lookup( labels.labels, :Ti)
period = ti[1] .. ti[end]

# create iterable table with data cube label
# each chunk can be read in memory
tab = CubeTable(
    label    = labels.labels[time=period], # [region=region,time=period] # [time=period]
    pei_30  = pei.pei_30[time=period], 
    pei_90  = pei.pei_90[time=period], #
    pei_180 = pei.pei_180[time=period], #
    t2mmax   = era.t2mmax[time=period], 
    t2m      = era.t2m[time=period],
    t2mmin   = era.t2mmin[time=period],
    tp       = era.tp[time=period],
    pet      = era.pet[time=period],
    event    = eventcube.layer[time=period],
    # gpp      = gpp[ time=period],
    # nee      = nee[ time=period],
    # ter      = ter[ time=period],
    landmask = lsmask_notime.lsm,
    # rank
    rt = tmax_ranks.layer[time=period, ],
    rd30 = pei_ranks.pei_30[time=period,],
    rd90 = pei_ranks.pei_90[time=period,],
    rd180 = pei_ranks.pei_180[time=period,],
    )

#tab[1]

# compute all stats on CubeTable
@time res = fitalldata1(tab);
# 123757.268823 seconds (111.79 G allocations: 55.172 TiB, 5.85% gc time, 0.01% compilation time)

# delete empty lines 
# sort results by decreasing Volume of events
sort!(res, by=i->i[end].v, rev=true);

# convert res to data frame
# res is a vector of unnamed tuples
# first convert tuple to named tuple
df = toDF(res)
# write DataFrame out to CSV file
outname = "$(patho)MergedEventStats$landonly.csv"
CSV.write(outname, df)
println(outname)
println("done!")
# df = CSV.read(outname, DataFrame)
# convert back to NamedTuple
# Tables.rowtable(df)

# close workers
rmprocs(workers())
