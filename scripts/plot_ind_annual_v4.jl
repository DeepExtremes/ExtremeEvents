# plot annual averages of 
# t2mmin t2m t2mmax
# tp
# pet
# PEI_30, _90, _180
# by continent (/Net/Groups/BGI/scratch/mweynants/DeepExtremes/era5-continents/output/continents.zarr)
println(Base.active_project())
using SlurmClusterManager, Distributed

using YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF

# region = "" # World: "" # "Europe"
# if region == ""
    latbound = (-90.0,90.0)
    lonbound = (0.0,360.0)
# elseif region == "Europe" # 35°N to 72°N and the longitudes 25°W to 65°E.
#     latbound = (35.0,72.0)
#     lonbound = [(360.0 - 25.0, 360.0),(0.0,65.0)]
# else 
#     # abort
#     error("region not set")
# end

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
        addprocs(1)
        sleep(0.001)
    end
end

@everywhere begin
    startyear = 1950
    endyear = 2023
end
timebound = (startyear:endyear)

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere include("../src/MyWeightedVar.jl")
@everywhere mergefun!(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.Ti)
    dfg = groupby(df,[:y, :continent])####### !!!
    variable_names = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"]
    continents = range(0,8)
    allstats = Dict("$(i).$(j).$(k)" => MyWeightedVariance() for i = string.(startyear:endyear), j = variable_names, k = continents)
    for k = keys(dfg), j = variable_names
        if !ismissing(k[2])
            dfs = dfg[k] ### !!!
            fit!(allstats["$(k[1]).$(j).$(k[2])"], dfs[!,j], cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
        end
    end
    return allstats
end
# i'm not sure that a good idea: will each chunk be touched only once? Should be 

pei = open_dataset("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/PEICube.zarr")

eraold = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr")
eranew = open_dataset("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/ERA5Cube.zarr")
cont = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/era5-continents/output/continents.zarr")
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = lsm[time = At(DateTime("2019-01-01T13:00:00"))]

# # try on small subset

function ct(;lon, lat, tim)
    t = CubeTable(
        t2mmin = eranew.t2mmin[time = tim, longitude = lon, latitude = lat],
        t2m = eranew.t2m[time = tim, longitude = lon, latitude = lat],
        t2mmax = eranew.t2mmax[time = tim, longitude = lon, latitude = lat],
        pei_30 = pei.pei_30[time = tim, longitude = lon, latitude = lat],
        pei_90 = pei.pei_90[time = tim, longitude = lon, latitude = lat],
        pei_180 = pei.pei_180[time = tim, longitude = lon, latitude = lat],
        # ssrd = era.ssrd[time = tim, longitude = lon, latitude = lat],
        tp = eranew.tp[time = tim, longitude = lon, latitude = lat],
        continent = cont.cont[longitude = lon, latitude = lat],
        pet = eraold.pet[time = tim, longitude = lon, latitude = lat],
        lsm = lsm_notime.lsm[longitude = lon, latitude = lat])
    return t
end

# t1 = ct(lon = (26.0,27.0), lat = (41.0,42.0), tim = (1950:1951));
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

# if typeof(lonbound) <: Vector
#     # split into 2 CubeTables, then merge
#     @time annualstats = mapreduce(mergefun!, lonbound) do lonbox
#         t = ct(lon = lonbox, lat = latbound, tim = timebound)
#         @time annualstats = pmapreduce(mergefun!, t) do tab
#             fit1(DataFrame(tab))
#         end
#     end
# else
    t = ct(lon = lonbound[1] .. lonbound[end], lat = latbound[1] .. latbound[end], tim = Date(startyear) .. Date(endyear, 12, 31))
    @time annualstats = pmapreduce(mergefun!,t) do tab
        fit1(DataFrame(tab))
    end
# end
# 16500.125021 seconds (8.60 M allocations: 454.161 MiB, 0.00% gc time, 0.02% compilation time)
# 48781.671483 seconds (7.85 M allocations: 404.246 MiB, 0.00% gc time, 0.01% compilation time)
using JLD
save("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/indicators_annual_wstats_continents.jld", annualstats)

# do not need workers any more
rmprocs(workers())
annualstats = load("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/indicators_annual_wstats_continents.jld")
allres = vec([fn(annualstats["$i.$j.$k"]) for fn = (mean, var), i = string.(startyear:endyear), j = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"], k = range(0,8)])
allx = vec(["$i.$s.$j.$k" for  s=["mean","var"], i = string.(startyear:endyear), j = ["t2mmin", "t2m", "t2mmax", "tp", "pet", "pei_30", "pei_90", "pei_180"], k = range(0,8)])

import CSV
df = DataFrame(x = allx, value = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:yr, :stat, :variable, :continent]))
outname = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/indicators_annual_wstats_continents.csv"
CSV.write(outname, df)

println("done.")

# do plots
using DataFrames
df = CSV.read(outname, DataFrame)
# tp is in m/day while pet and pei are in mm/day => *1e3
df1 = df |>
    (df -> DataFrames.filter(:variable => ==("tp"), df)) |>
    (df -> DataFrames.select(df, [:value, :stat] => ByRow((x,y) -> ifelse(y == "var", x * 1e6, x * 1e3)) => :ve3))

df[df.variable .== "tp", :value] = df1[:, :ve3]

using CairoMakie

# Theil-Sen
include("../src/stats.jl")
# m, b = theilsen(dfpp[!, :Year], dfpp[!,:Area_pc_sum])

### compare with Mann Kendall implementation by @mixstam1821
# reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(dfpp.Year, dfpp.Area_pc_sum)


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
rf=1;rn="Africa";
for (rf, rn) in continents
    println("$rn has number $rf")
    region = rn
tdf = df |> 
    (df -> DataFrames.filter(:continent => ==(rf), df)) |>
    (df -> DataFrames.filter(:stat => ==("mean"), df)) |>
    (df -> DataFrames.filter(:variable => x -> contains(x,"t2m"), df)) |>
    (df -> DataFrames.transform!(df,:value => ByRow(x -> x - 273.15) => :value, ))
f = Figure(size=(900,600));
ax = Axis(f[1,1],
        title = region,
        xlabel = "Year", 
        ylabel = "Yearly continental average of \n temperature over land"
    );
t2mcols = [colorant"#ffd1a2", colorant"#ffab56", colorant"#ff9223", ]'
for (i,variable) in zip(1:3,unique(tdf.variable)) 
    tdfv = DataFrames.filter(:variable => ==(variable), tdf)
    reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(tdfv.yr, tdfv.value)
    scatter!(ax, tdfv.yr, tdfv.value, color = t2mcols[i], label = variable)
    lines!(ax, tdfv.yr,  (slope .* tdfv.yr .+ intercept), color = t2mcols[i], label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
end
Legend(f[2,1], ax, orientation = :horizontal, nbanks = 2, framevisible = false)
f
# display(p)
save("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/t2mmax_mean_annual_land_$(region).png", f)

# PEI
pdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => x -> contains(x,"pei"), df))
f = Figure(size=(900,600));
ax = Axis(f[1,1],
        title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \nPrecipitation-Evapotranspiration over land",
        );
peicols = [colorant"#002D5A", colorant"#A6C5E8", colorant"#4C7FB8", ]'
for (i,variable) in zip(1:3,unique(pdf.variable)) 
    pdfv = DataFrames.filter(:variable => ==(variable), pdf)
    reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(pdfv.yr, pdfv.value)
    scatter!(ax, pdfv.yr, pdfv.value, color = peicols[i], label = variable)
    lines!(ax, pdfv.yr,  (slope .* pdfv.yr .+ intercept), color = peicols[i], label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
end
Legend(f[2,1], ax, orientation = :horizontal, nbanks = 2, framevisible = false)
f
save("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/pei_mean_annual_land_$region.png", f)

# TP
pdf = df |> 
    (df -> filter(:continent => ==(rf), df)) |>
    (df -> filter(:stat => ==("mean"), df)) |>
    (df -> filter(:variable => x -> contains(x,"tp"), df))
f = Figure(size=(900,600));
ax = Axis(f[1,1],
        title = region,
        xlabel = "Year", 
        ylabel = "Yearly continental average of \nPrecipitation over land",
    );
peicols = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]'
for (i,variable) in zip(1:3,unique(pdf.variable)) 
    pdfv = DataFrames.filter(:variable => ==(variable), pdf)
    reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(pdfv.yr, pdfv.value)
    scatter!(ax, pdfv.yr, pdfv.value, color = peicols[i], label = variable)
    lines!(ax, pdfv.yr,  (slope .* pdfv.yr .+ intercept), color = peicols[i], label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
end
Legend(f[2,1], ax, orientation = :horizontal, nbanks = 2, framevisible = false)
f


# TP and PET
edf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => x -> contains(x,"tp") || contains(x, "pet"), df)) 
p_e = edf[edf.variable .== "tp" , :value]  .+ edf[edf.variable .== "pet", :value]   
f = Figure(size=(900,600));
ax = Axis(f[1,1],
        title = region,
        xlabel = "Year", 
        ylabel = "Yearly continental average of \nPrecipitation and Evapotranspiration over land",
    );
pecols = [colorant"#4C7FB8", colorant"#C40000"]'
for (i,variable) in zip(1:3,unique(edf.variable)) 
    edfv = DataFrames.filter(:variable => ==(variable), edf)
    reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(edfv.yr, edfv.value)
    barplot!(ax, edfv.yr, edfv.value, color = pecols[i], label = variable)
    lines!(ax, edfv.yr,  (slope .* edfv.yr .+ intercept), color = pecols[i], label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
end
Legend(f[2,1], ax, orientation = :horizontal, nbanks = 2, framevisible = false)
f
save("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/tp_pet_mean_annual_land_$region.png", f)
       
# tvdf = df |> 
#     (df -> filter(:continent => ==(rf), df)) |>
#     (df -> filter(:stat => ==("var"), df)) |>
#     (df -> filter(:variable => x -> contains(x,"t2m"), df)) |>
#     (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
# p3 = @df tvdf scatter(
#         :yr, :std, group = :variable, smooth = true,
#         legend = :outerright, lw = 1,
#         title = region,
#         xlabel = "Year", 
#         ylabel = "Yearly continental standard deviation of \n temperature over land",
#         size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#         colour = [colorant"#ffab56", colorant"#ff9223", colorant"#ffd1a2"]',
#         )
# savefig(p3, "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/t2mmax_std_annual_land_$region.png")

# pvdf = df |> 
#         (df -> filter(:continent => ==(rf), df)) |>
#         (df -> filter(:stat => ==("var"), df)) |>
#         (df -> filter(:variable => x -> contains(x,"pei"), df)) |>
#         (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
# p4 = @df pvdf scatter(
#             :yr, :std, group = :variable, smooth = true,
#             legend = :outerright, lw = 1,
#             title = region,
#             xlabel = "Year", 
#             ylabel = "Yearly continental standard deviation of \nPrecipitation-Evapotranspiration over land",
#             size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#             colour = [colorant"#002D5A",colorant"#A6C5E8", colorant"#4C7FB8", ]',
#             )
# savefig(p4, "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/pei_std_annual_land_$region.png")

# pvdf = df |> 
#         (df -> filter(:continent => ==(rf), df)) |>
#         (df -> filter(:stat => ==("var"), df)) |>
#         (df -> filter(:variable => (x -> x=="tp" || x=="pet"), df)) |>
#         (df -> DataFrames.transform(df, :value => (x -> sqrt.(x)) => :std))
# p5 = @df pvdf scatter(
#             :yr, :std, group = :variable, smooth = true,
#             legend = :outerright, lw = 1,
#             title = region,
#             xlabel = "Year", 
#             ylabel = "Yearly continental standard deviation of \nPrecipitation and Evapotranspiration over land",
#             size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#             colour = [colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]',
#             )
# savefig(p5, "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig/tp_pet_std_annual_land_$region.png")
        
# scatter(1950:1978, sqrt.(df[df.stat .== "var" .&& df.variable .== "tp" .&& df.continent .== rf, :value]), smooth = true, color = colorant"#002D5A",
#     title = region,
#     label = "tp", xlabel = "Year", ylabel = "Continental annual standard deviation of \n total precipitation",
#     size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
# savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/tp_std_annual_land_$region.png")
# #hline!([sum(df[df.stat .== "mean" .&& df.variable .== "tp", :value])./(2021-1950+1)], label = "time series mean")
# scatter(1950:1978, sqrt.(df[df.stat .== "var" .&& df.variable .== "pet" .&& df.continent .== rf, :value]), smooth = true, color = colorant"#A6C5E8",
#     title = region,
#     label = "pet", xlabel = "Year", ylabel = "Gloabl annual standard deviation of \n potential evapotranspiration",
#     size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),)
# savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/pet_std_annual_land_$region.png")

end