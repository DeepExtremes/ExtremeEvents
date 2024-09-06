# Compute area by year, by continent and by event type
println(Base.active_project())
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
else
    addprocs(8)
end

@everywhere begin
    using Pkg
    if !isnothing(findfirst("Slurm", @__DIR__)) || !isnothing(findfirst("scripts", @__DIR__))
        Pkg.activate("$(@__DIR__)/..")
    else
        Pkg.activate("$(@__DIR__)")
    end
end

@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    dfg = groupby(df, [:y, :cont])
    continents = range(0,8)
    allhists = Dict("$(i).$(k)" => WeightedHist(-0.5:1.0:16.5) for i in 1950:2022, k = continents)
    for k in keys(dfg)
        if !ismissing(k[2])
            dfs = dfg[k]
            fit!(allhists["$(k[1]).$(k[2])"], dfs.event, cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
        end
    end
    allhists
end

pot = 0.01
ne = 0.1
path2v = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3"
ds = open_dataset("$path2v/EventCube_ranked_pot$(pot)_ne$(ne).zarr/")
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))
cont = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/era5-continents/output/continents.zarr")

# sds = subsetcube(ds, time = (1998:1998), longitude = (325.0,326.0), latitude = (45.0,46.0))
# slsm = subsetcube(lsm_notime, longitude = (325.0,326.0), latitude = (45.0,46.0))
function ct(;lon = (0.0,360.0), lat = (-90.0,90.0), tim = (1950:2022))
    t = CubeTable(
        event = ds.layer[time = tim, longitude = lon, latitude = lat],
        lsm = lsm_notime.lsm[longitude = lon, latitude = lat],
        cont = cont.cont[longitude = lon, latitude = lat]
    )
    return t
end

# t1 = ct(lon = (26.0,27.0), lat = (41.0,42.0), tim = (1950:1951))
# r1 = fit1(DataFrame(t1[1]))
# sum(value(r1["1950.8"]).y)
# sum(value(r1["1951.8"]).y)
# r1["1951.8"]
# r1["1950.2"]
# r1["1951.2"]

t = ct()

annualstats = pmapreduce(mergefun,t) do tab
    fit1(DataFrame(tab))
end
# sum(value(annualstats["1950.8"]).y)
# sum(value(annualstats["1951.8"]).y)
# sum(value(annualstats["1967.8"]).y)

using JLD
save("$path2v/YearlyEventType_ranked_pot" * string(pot) * "_ne" * string(ne) * "_land.jld",annualstats)
rmprocs(workers())