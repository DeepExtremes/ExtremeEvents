
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere include("../src/MyWeightedVar.jl")

@everywhere mergefun!(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fitc(df)
    dfg = groupby(df,[:cont])####### !!!
    continents = range(0,8)
    allstats = Dict("$(k)" => WeightedSum() for  k = continents)
    for k = keys(dfg)
        if !ismissing(k[1])
            dfs = dfg[k] ### !!!
            # @show dfs
            fit!(allstats["$(k[1])"], ones(length(dfs.cont)) , cosd.(dfs.latitude) .* (dfs.lsm .> 0.5) )
        end
    end
    return allstats
end

cont = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/era5-continents/output/continents.zarr")
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))

function ct(;lon, lat)
    t = CubeTable(
        cont = cont.cont[longitude = lon, latitude = lat],
        lsm = lsm_notime.lsm[longitude = lon, latitude = lat]
    )
    return t
end

latbound = (-90.0,90.0)
lonbound = (0.0,360.0)

# t1 = ct(lon = (26.0,27.0), lat = (41.0,42.0));
# @time r1 = fitc(DataFrame(t1[1]))

t = ct(lon = lonbound, lat = latbound, )
@time stats = pmapreduce(mergefun!,t) do tab
    fitc(DataFrame(tab))
end

using JLD
save("$path2v/land_wstats_continents.jld",stats)
rmprocs(workers())

allres = vec([value(stats["$k"]) for k = range(0,8)])

import CSV
df = DataFrame(
    x = range(0,8), 
    cont = ["Null", "Africa", "Asia", "Australia", "North America", "Oceania", "South America", "Antarctica", "Europe"],
     value = allres) 
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/land_wstats_continents.csv"
CSV.write(outname, df)
