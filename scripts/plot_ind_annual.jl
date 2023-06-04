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
@everywhere include("/Net/Groups/BGI/scratch/mweynants/ExtremeEvents/src/MyWeightedVar.jl")
@everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    dfg = groupby(df,:y)
    variable_names = ["t2mmax", "pei_30", "pei_90", "pei_180"]
    allstats = Dict(i*"."*j => MyWeightedVariance() for i = string.(1950:2021), j = variable_names)
    for k = keys(dfg), j = variable_names
        dfs = dfg[k]
        fit!(allstats[string(k[1])*"."*j], dfs[!,j], cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
    end
    return allstats
end
# i'm not sure that a good idea: will each chunk be touched only once? Should be 

pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))

# # try on small subset
# sera = subsetcube(era, time = (1950:1950), longitude = (12.0,13.0), latitude = (45.0,46.0))
# spei = subsetcube(pei, time = (1950:1950), longitude = (12.0,13.0), latitude = (45.0,46.0))
# slsm = subsetcube(lsm_notime, longitude = (12.0,13.0), latitude = (45.0,46.0))

# ts = CubeTable(t2mmax = sera.t2mmax,
#             pei_30 = spei.pei_30,
#             pei_90 = spei.pei_90,
#             pei_180 = spei.pei_180,
#             lsm = slsm.lsm)
# @time myrs = fit1(DataFrame(ts[1]))

# myrs["1950.t2mmax"].µ
# myrs["1950.t2mmax"].σ2
# myrs["1950.t2mmax"].n

# # with WeightedVariance
# @time rs = fit1(DataFrame(ts[1]))

# rs["1950.t2mmax"].µ
# rs["1950.t2mmax"].σ2
# rs["1950.t2mmax"].n # n is larger but µ and σ2 are the same


t = CubeTable(t2mmax = era.t2mmax,
            pei_30 = pei.pei_30,
            pei_90 = pei.pei_90,
            pei_180 = pei.pei_180,
            lsm = lsm_notime.lsm)
# @time r1 = fit1(DataFrame(t[1]))
# # 84.553397 seconds (24.68 M allocations: 14.434 GiB, 2.81% gc time, 10.86% compilation time)
# printnl(length(t))
# # 85 * 3605 => about 85 hours, divided by nb of workers * merge time => >10 hours

# # t[1] goes from 1950 to 1965. All of them get values of NaN
# r1["1950.t2mmax"] # value=NaN
# r1["1950.t2mmax"].µ
# r1["1950.t2mmax"].σ2
# r1["1950.t2mmax"].n
# # NaNs introduced. Why? Probably there are missing values, but there shouldn't!
# # I don't understand how and where NaNs are introduced. When I do
# sum(t[1].t2mmax .* cosd.(t[1].latitude).* (t[1].lsm .> 0.5))./sum(cosd.(t[1].latitude).* (t[1].lsm .> 0.5))
# sum(cosd.(t[1].latitude).* (t[1].lsm .> 0.5))
# # there is also a problem.
# # but if I skip the wi = 0 with MyWeightedVariance, it should work.


@time annualstats = pmapreduce(mergefun,t) do tab
    fit1(DataFrame(tab))
end
# 16500.125021 seconds (8.60 M allocations: 454.161 MiB, 0.00% gc time, 0.02% compilation time)
rmprocs(workers())

allres = vec([fn(annualstats[i*"."*j]) for fn = (mean, var), i = string.(1950:2021), j = ["t2mmax", "pei_30", "pei_90", "pei_180"]])
allx = vec([i*"."*s*"."*j for  s=["mean","var"], i = string.(1950:2021), j = ["t2mmax", "pei_30", "pei_90", "pei_180"]])

import CSV
df = DataFrame(x = allx, value = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:yr, :stat, :variable]))
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/indicators_annual_wstats_land.csv"
CSV.write(outname, df)

print("done.")

# do plots
# df = CSV.read(outname, DataFrame)
