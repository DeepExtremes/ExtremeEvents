# plot annual averages of 
# t2mmax
# PEI_30, _90, _180

using Pkg
Pkg.instantiate()
using Distributed

# addprocs(10) # workers called from slurm script.
@everywhere begin
    using Pkg
    Pkg.activate("/Net/Groups/BGI/scratch/mweynants/ExtremeEvents/")
    Pkg.instantiate()
end
@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere include("/Net/Groups/BGI/scratch/mweynants/ExtremeEvents/src/MyWeightedVar.jl")
@everywhere mergefun!(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    dfg = groupby(df,[:y, :cont])####### !!!
    variable_names = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"]
    continents = range(0,8)
    allstats = Dict("$(i).$(j).$(k)" => MyWeightedVariance() for i = string.(1950:1978), j = variable_names, k = continents)
    for k = keys(dfg), j = variable_names
        if !ismissing(k[2])
            dfs = dfg[k] ### !!!
            fit!(allstats["$(k[1]).$(j).$(k[2])"], dfs[!,j], cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
        end
    end
    return allstats
end
# i'm not sure that a good idea: will each chunk be touched only once? Should be 

pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)

cont = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/era5-continents/output/continents.zarr")

lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))

# # try on small subset

function ct(;lon, lat, tim)
    t = CubeTable(t2mmin = era.t2mmin[time = tim, longitude = lon, latitude = lat],
        t2m = era.t2m[time = tim, longitude = lon, latitude = lat],
        t2mmax = era.t2mmax[time = tim, longitude = lon, latitude = lat],
        pei_30 = pei.pei_30[time = tim, longitude = lon, latitude = lat],
        pei_90 = pei.pei_90[time = tim, longitude = lon, latitude = lat],
        pei_180 = pei.pei_180[time = tim, longitude = lon, latitude = lat],
        ssrd = era.ssrd[time = tim, longitude = lon, latitude = lat],
        tp = era.tp[time = tim, longitude = lon, latitude = lat],
        cont = cont.cont[longitude = lon, latitude = lat],
        pet = era.pet[time = tim, longitude = lon, latitude = lat],
        lsm = lsm_notime.lsm[longitude = lon, latitude = lat])
    return t
end

# t1 = ct(lon = (12.0,13.0), lat = (45.0,46.0), tim = (1950:1950));
# t2 = ct(lon = (10.0,10.5), lat = (45.0,46.0), tim = (1950:1950));

# @time r1 = fit1(DataFrame(t1[1]))
# r1["1950.t2m"]
# @time r2 = fit1(DataFrame(t2[1]))
# r2["1950.t2m"]
# @time r = mergefun!(r1,r2)
# r["1950.t2m"]
# # # 84.553397 seconds (24.68 M allocations: 14.434 GiB, 2.81% gc time, 10.86% compilation time)
# # printnl(length(t))
# # # 85 * 3605 => about 85 hours, divided by nb of workers * merge time => >10 hours

# # NaNs introduced. if I skip the wi = 0 with MyWeightedVariance, it works (no division by 0 anymore)
# latbound = (45.5,46.0)
# lonbound = [(359.5,360.0), (0.0,0.5)]
latbound = (-90.0,90.0)
lonbound = (0.0,360.0)
timebound = (1950:1978)

t = ct(lon = lonbound, lat = latbound, tim = timebound)
@time annualstats = pmapreduce(mergefun!,t) do tab
    fit1(DataFrame(tab))
end

# 16500.125021 seconds (8.60 M allocations: 454.161 MiB, 0.00% gc time, 0.02% compilation time)
rmprocs(workers())

### !! Need to aggregate all continents before saving to csv

allres = vec([fn(annualstats["$i.$j.$k"]) for fn = (mean, var), i = string.(1950:1978), j = ["t2mmin", "t2m", "t2mmax", "tp","pet", "pei_30", "pei_90", "pei_180" ], k = range(0,8)]) # 
allx = vec(["$i.$s.$j.$k" for  s=["mean","var"], i = string.(1950:1978), j = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"], k = range(0,8)]) # 

import CSV
df = DataFrame(x = allx, value = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:yr, :stat, :variable, :continent]))
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/indicators_annual_wstats_continents_backext.csv"
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

continents = Dict(
        1 => "Africa", 
        2 => "Asia", 
        3 => "Australia", 
        4 => "North America",
        5 => "Oceania", 
        6 => "South America", 
        7 => "Antarctica",  
        8 => "Europe",
        )
for (rf, rn) in continents
    println("$rn has number $rf")
    region = rn
tdf = df |> 
    (df -> filter(:continent => ==(rf), df)) |>
    (df -> filter(:stat => ==("mean"), df)) |>
    (df -> filter(:variable => x -> contains(x,"t2m"), df))
p = @df tdf scatter(
        :yr, :value, group = :variable, smooth = true,
        legend = :outerright, lw = 1,
        title = region,
        xlabel = "Year", 
        ylabel = "Yearly continental average of \n temperature over land",
        size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
        colour = [colorant"#ffab56", colorant"#ff9223", colorant"#ffd1a2"]',
        )
# display(p)
savefig(p, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/t2mmax_mean_annual_land_backext_$(region).png")

# PEI
pdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => x -> contains(x,"pei"), df))
p1 = @df pdf scatter(
            :yr, :value, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \nPrecipitation-Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#002D5A", colorant"#A6C5E8", colorant"#4C7FB8", ]',
            )
savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/pei_mean_annual_land_backext$region.png")

# TP
df1 = df |> 
(df -> filter(:continent => ==(rf), df)) |>
(df -> filter(:stat => ==("mean"), df)) |>
(df -> filter(:variable => x -> contains(x,"tp"), df))

p2 = @df df1 scatter(
            :yr, :value, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \nPrecipitation over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
            )
# display(p2)


# TP and PET
edf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => x -> contains(x,"tp") || contains(x, "pet"), df)) 
p_e = edf[edf.variable .== "tp" , :value]  .+ edf[edf.variable .== "pet", :value]   

p2 = @df edf scatter(
            :yr, :value, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \nPrecipitation and Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
            )
savefig(p2, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_pet_mean_annual_land_backext$region.png")
       
scatter( 1950:1978, p_e, smooth = true, label = "tp+pet", color = colorant"#4C7FB8",
    xlabel = "Year", ylabel = "Difference of continental annual averages of \n total precipitation and potential evapotranspiration",
    title = rn,
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_pet_mean_annual_land_backext$region.png")

scatter(1950:1978, df[df.stat .== "mean" .&& df.variable .== "tp" .&& df.continent .== rf, :value], smooth = true, color = colorant"#002D5A",
    title = region,
    label = "tp", xlabel = "Year", ylabel = "Continental annual average of \n total precipitation",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_mean_annual_land_backext$region.png")
#hline!([sum(df[df.stat .== "mean" .&& df.variable .== "tp", :value])./(2021-1950+1)], label = "time series mean")
scatter(1950:1978, df[df.stat .== "mean" .&& df.variable .== "pet" .&& df.continent .== rf, :value], smooth = true, color = colorant"#A6C5E8",
    title = region,
    label = "pet", xlabel = "Year", ylabel = "Continental annual average of \n potential evapotranspiration",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/pet_mean_annual_land_backext$region.png")

tvdf = df |> 
    (df -> filter(:continent => ==(rf), df)) |>
    (df -> filter(:stat => ==("var"), df)) |>
    (df -> filter(:variable => x -> contains(x,"t2m"), df)) |>
    (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p3 = @df tvdf scatter(
        :yr, :std, group = :variable, smooth = true,
        legend = :outerright, lw = 1,
        title = region,
        xlabel = "Year", 
        ylabel = "Yearly continental standard deviation of \n temperature over land",
        size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
        colour = [colorant"#ffab56", colorant"#ff9223", colorant"#ffd1a2"]',
        )
savefig(p3, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/t2mmax_std_annual_land_backext$region.png")

pvdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("var"), df)) |>
        (df -> filter(:variable => x -> contains(x,"pei"), df)) |>
        (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p4 = @df pvdf scatter(
            :yr, :std, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental standard deviation of \nPrecipitation-Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#002D5A",colorant"#A6C5E8", colorant"#4C7FB8", ]',
            )
savefig(p4, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/pei_std_annual_land_backext$region.png")

pvdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("var"), df)) |>
        (df -> filter(:variable => (x -> x=="tp" || x=="pet"), df)) |>
        (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
p5 = @df pvdf scatter(
            :yr, :std, group = :variable, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental standard deviation of \nPrecipitation and Evapotranspiration over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
            )
savefig(p5, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_pet_std_annual_land_backext$region.png")
        
scatter(1950:1978, sqrt.(df[df.stat .== "var" .&& df.variable .== "tp" .&& df.continent .== rf, :value]), smooth = true, color = colorant"#002D5A",
    title = region,
    label = "tp", xlabel = "Year", ylabel = "Continental annual standard deviation of \n total precipitation",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_std_annual_land_backext$region.png")
#hline!([sum(df[df.stat .== "mean" .&& df.variable .== "tp", :value])./(2021-1950+1)], label = "time series mean")
scatter(1950:1978, sqrt.(df[df.stat .== "var" .&& df.variable .== "pet" .&& df.continent .== rf, :value]), smooth = true, color = colorant"#A6C5E8",
    title = region,
    label = "pet", xlabel = "Year", ylabel = "Gloabl annual standard deviation of \n potential evapotranspiration",
    size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/pet_std_annual_land_backext$region.png")

end