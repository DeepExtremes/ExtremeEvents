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
    variable_names = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"]
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


t = CubeTable(t2mmin = era.t2mmin,
            t2m = era.t2m,
            t2mmax = era.t2mmax,
            pei_30 = pei.pei_30,
            pei_90 = pei.pei_90,
            pei_180 = pei.pei_180,
            ssrd = era.ssrd,
            tp = era.tp,
            pet = era.pet,
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

allres = vec([fn(annualstats[i*"."*j]) for fn = (mean, var), i = string.(1950:2021), j = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"]])
allx = vec([i*"."*s*"."*j for  s=["mean","var"], i = string.(1950:2021), j = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"]])

import CSV
df = DataFrame(x = allx, value = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:yr, :stat, :variable]))
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/indicators_annual_wstats_land.csv"
CSV.write(outname, df)

print("done.")

# do plots
using DataFrames
df = CSV.read(outname, DataFrame)
# tp is in m/day while pet and pei are in mm/day => *1e3
df1 = df |>
    (df -> DataFrames.filter(:variable => ==("tp"), df)) |>
    (df -> DataFrames.select(df, [:value, :stat] => ByRow((x,y) -> ifelse(y == "var", x * 1e6, x * 1e3)) => :ve3))

df[df.variable .== "tp", :value] = df1[:, :ve3]

using StatsPlots
tdf = df |> 
    (df -> filter(:stat => ==("mean"), df)) |>
    (df -> filter(:variable => x -> contains(x,"t2m"), df))
p = @df tdf scatter(
        :yr, :value, group = :variable, smooth = true,
        legend = :outerright, lw = 1,
        xlabel = "Year", 
        ylabel = "Yearly global average of \n temperature over land",
        size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
        colour = [colorant"#ffab56", colorant"#ff9223", colorant"#ffd1a2"]',
        # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
        )
savefig(p, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/t2mmax_mean_annual_land.png")

pdf = df |> 
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => x -> contains(x,"pei"), df))
p1 = @df pdf scatter(
            :yr, :value, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            xlabel = "Year", 
            ylabel = "Yearly global average of \nPrecipitation-Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#002D5A", colorant"#A6C5E8", colorant"#4C7FB8", ]',
            # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
            )
savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/pei_mean_annual_land.png")

edf = df |> 
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => (x -> x=="tp" || x=="pet"), df)) 
p_e = edf[edf.variable .== "tp" , :value]  .+ edf[edf.variable .== "pet", :value]   

p2 = @df edf scatter(
            :yr, :value, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            xlabel = "Year", 
            ylabel = "Yearly global average of \nPrecipitation and Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
            # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
            )
savefig(p2, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/tp_pet_mean_annual_land.png")
       
scatter( 1950:2021, p_e, smooth = true, label = "tp+pet", color = colorant"#4C7FB8",
    xlabel = "Year", ylabel = "Difference of gloabl annual averages of \n total precipitation and potential evapotranspiration",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/tp_pet_mean_annual_land.png")
scatter(1950:2021, df[df.stat .== "mean" .&& df.variable .== "tp", :value], smooth = true, color = colorant"#002D5A",
    label = "tp", xlabel = "Year", ylabel = "Gloabl annual average of \n total precipitation",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/tp_mean_annual_land.png")
#hline!([sum(df[df.stat .== "mean" .&& df.variable .== "tp", :value])./(2021-1950+1)], label = "time series mean")
scatter(1950:2021, df[df.stat .== "mean" .&& df.variable .== "pet", :value], smooth = true, color = colorant"#A6C5E8",
    label = "pet", xlabel = "Year", ylabel = "Gloabl annual average of \n potential evapotranspiration",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/pet_mean_annual_land.png")



tvdf = df |> 
    (df -> filter(:stat => ==("var"), df)) |>
    (df -> filter(:variable => x -> contains(x,"t2m"), df)) |>
    (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p3 = @df tvdf scatter(
        :yr, :std, group = :variable, smooth = true,
        legend = :outerright, lw = 1,
        xlabel = "Year", 
        ylabel = "Yearly global standard deviation of \n temperature over land",
        size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
        colour = [colorant"#ffab56", colorant"#ff9223", colorant"#ffd1a2"]',
        # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
        )
savefig(p3, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/t2mmax_std_annual_land.png")

pvdf = df |> 
        (df -> filter(:stat => ==("var"), df)) |>
        (df -> filter(:variable => x -> contains(x,"pei"), df)) |>
        (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p4 = @df pvdf scatter(
            :yr, :std, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            xlabel = "Year", 
            ylabel = "Yearly global standard deviation of \nPrecipitation-Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#002D5A",colorant"#A6C5E8", colorant"#4C7FB8", ]',
            # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
            )
savefig(p4, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/pei_std_annual_land.png")

pvdf = df |> 
        (df -> filter(:stat => ==("var"), df)) |>
        (df -> filter(:variable => (x -> x=="tp" || x=="pet"), df)) |>
        (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p5 = @df pvdf scatter(
            :yr, :std, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            xlabel = "Year", 
            ylabel = "Yearly global standard deviation of \nPrecipitation and Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
            # xrotation = 45.0, xtickfontsize = 6,xlims = (0,(2021-1950+1)),xticks=(.5:5:(2021-1950+1),string.(1950:5:2021))
            )
savefig(p5, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/tp_pet_std_annual_land.png")
        
scatter(1950:2021, sqrt.(df[df.stat .== "var" .&& df.variable .== "tp", :value]), smooth = true, color = colorant"#002D5A",
    label = "tp", xlabel = "Year", ylabel = "Gloabl annual standard deviation of \n total precipitation",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/tp_std_annual_land.png")
#hline!([sum(df[df.stat .== "mean" .&& df.variable .== "tp", :value])./(2021-1950+1)], label = "time series mean")
scatter(1950:2021, sqrt.(df[df.stat .== "var" .&& df.variable .== "pet", :value]), smooth = true, color = colorant"#A6C5E8",
    label = "pet", xlabel = "Year", ylabel = "Gloabl annual standard deviation of \n potential evapotranspiration",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/pet_std_annual_land.png")
