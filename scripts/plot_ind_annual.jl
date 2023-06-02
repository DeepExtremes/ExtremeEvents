# plot annual and monthly avergae of 
# t2mmax
# PEI_30, _90, _180

using YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF, Distributed

addprocs(10)
@everywhere begin
    using Pkg
    Pkg.activate("/Net/Groups/BGI/scratch/mweynants/ExtremeEvents/")
    # Pkg.instantiate()
end
@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    dfg = groupby(df,:y)
    variable_names = ["t2mmax", "pei_30", "pei_90", "pei_180"]
    allstats = Dict(i*"."*s*"."*j => ifelse(s=="m", WeightedMean(), WeightedVariance()) for i = string.(1950:2021), j = variable_names, s=["m","v"])
    for k = keys(dfg), j = variable_names, s=["m","v"]
        dfs = dfg[k]
        fit!(allstats[string(k[1])*"."*s*"."*j], dfs[!,j], cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
    end
    return allstats
end
# i'm not sure that a good idea: will each chunk be touched only once? Should be 

pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))

# try on small subset
# sera = subsetcube(era, time = (1998:1998), longitude = (25.0,26.0), latitude = (45.0,46.0))
# spei = subsetcube(pei, time = (1998:1998), longitude = (25.0,26.0), latitude = (45.0,46.0))
# slsm = subsetcube(lsm_notime, longitude = (25.0,26.0), latitude = (45.0,46.0))

# t = CubeTable(t2mmax = sera.t2mmax,
#             pei_30 = spei.pei_30,
#             pei_90 = spei.pei_90,
#             pei_180 = spei.pei_180,
#             lsm = slsm.lsm)

t = CubeTable(t2mmax = era.t2mmax,
            pei_30 = pei.pei_30,
            pei_90 = pei.pei_90,
            pei_180 = pei.pei_180,
            lsm = lsm_notime.lsm)
# @time r1 = fit1(DataFrame(t[1]))
# # 84.553397 seconds (24.68 M allocations: 14.434 GiB, 2.81% gc time, 10.86% compilation time)
# printnl(length(t))
# # 85 * 3605 => about 85 hours, divided by nb of workers * merge time => >10 hours

@time annualstats = pmapreduce(mergefun,t) do tab
    fit1(DataFrame(tab))
end
# 16500.125021 seconds (8.60 M allocations: 454.161 MiB, 0.00% gc time, 0.02% compilation time)
rmprocs(workers())

allres = vec([value(annualstats[i*"."*s*"."*j]) for i = string.(1950:2021), j = ["t2mmax", "pei_30", "pei_90", "pei_180"], s=["m","v"]])
# NaNs introduced. Why? Probably there are missing values, but there shouldn't!
# would be OK if I used skipmissing, but how to pass it on the weights?
allx = vec([i*"."*s*"."*j for i = string.(1950:2021), j = ["t2mmax", "pei_30", "pei_90", "pei_180"], s=["m","v"]])

import CSV
df = DataFrame(x = allx, value = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:yr, :stat, :variable]))
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/indicators_annual_wstats_land.csv"
CSV.write(outname, df)

print("done.")

# do plots
df = CSV.read(outname, DataFrame)
