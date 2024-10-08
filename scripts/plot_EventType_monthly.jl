# Plot area by year by event type
using Revise
using YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF, Distributed

addprocs(8)
@everywhere begin
    using Pkg
    Pkg.activate("/Net/Groups/BGI/scratch/mweynants/ExtremeEvents/")
end
@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates, NetCDF
@everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    df.m = month.(df.time)
    dfg = groupby(df, [:y, :m])
    allhists = Dict(i=>WeightedHist(-0.5:1.0:32.5) for i in reduce(hcat,broadcast(.+, (1950:2021).*100, [1:12]')))
    for k in keys(dfg)
        dfs = dfg[k]
        fit!(allhists[k[1]*100+k[2]], dfs.event, cosd.(dfs.latitude) .* (dfs.lsm .> 0.5))
    end
    allhists
end

ds = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr/")
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))
t = CubeTable(event = ds.layer,
            lsm = lsm_notime.lsm)

# sds = subsetcube(ds, time = (1998:1998), longitude = (25.0,26.0), latitude = (45.0,46.0))
# slsm = subsetcube(lsm_notime, longitude = (25.0,26.0), latitude = (45.0,46.0))
# t = CubeTable(event = sds.layer,
#             lsm = slsm.lsm)


monthlystats = pmapreduce(mergefun,t) do tab
    fit1(DataFrame(tab))
end
rmprocs(workers())

moyr = reduce(hcat,broadcast(.+, (1950:2021).*100, [1:12]'))'
allx = [value(monthlystats[yr]).x[i] for i in 1:17, yr in moyr]
allres = [value(monthlystats[yr]).y[i] for i in 1:17, yr in moyr]
using Plots
allres_norm = allres./sum(allres,dims=1)
labels = map([UInt8(i) for _ in 1:1, i in 1:16]) do i
    n = ""
    ((i & 0x01) > 0) && (n = join((n,"HW"),"_"))
    ((i & 0x02) > 0) && (n = join((n,"D30"),"_"))
    ((i & 0x04) > 0) && (n = join((n,"D90"),"_"))
    ((i & 0x08) > 0) && (n = join((n,"D180"),"_"))
    lstrip(n,'_')
end
p = plot(allres_norm[1:end-1,:,1]',
    labels=labels, legend = :outerright, lw=1, 
    xlabel = "Month-Year",
    xticks = (0.0:120:840,string.(1950:10:2021)),
    ylabel = "Monthly fraction of land area and days affected",
    size=(800,450), dpi=300)

savefig(p,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/n_extremes_month_land.png")

# save to csv
import CSV
using DataFrames
# not sure about this one
df = DataFrame(dropdims(allres,dims=3), :auto)
pot = 0.01
ne = 0.1
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/MonthlyEventType_ranked_pot" * string(pot) * "_ne" * string(ne) * "_land.csv"
CSV.write(outname, df)
df = CSV.read(outname, DataFrame)
rename!(df, vec(string.(moyr)))
df.ev = 0:16

using StatsPlots
# land Area x Days by Int8
v = [
        RGBA(1,1,1,0), # 0x00 # 0
        RGBA(1,0,0,1), # 0x01 # 1 !!
        RGBA(0,0,.8,1), # 0x02 # 2 
        RGBA(1,0,.8,1), # 0x03 # 3 !!
        RGBA(0,0,.6,1), # 0x04 # 4 
        RGBA(1,0,.6,1), # 0x05 # 5 !!
        RGBA(0,0,.6,1), # 0x06 # 6 
        RGBA(1,0,.6,1), # 0x07 # 7 !!
        RGBA(0,0,.4,1), # 0x08 # 8 
        RGBA(1,0,.4,1), # 0x09 # 9 !!
        RGBA(0,0,.4,1), # 0x0a # 10
        RGBA(1,0,.4,1), # 0x0b # 11!!
        RGBA(0,0,.4,1), # 0x0c # 12
        RGBA(1,0,.4,1), # 0x0d # 13!!
        RGBA(0,0,.4,1), # 0x0e # 14
        RGBA(1,0,.4,1), # 0x0f # 15!!
        RGBA(.7,.7,.7,1), # 0x10 # 16
    ]
# ChatGPT suggestion of colourblind friendly colourscale
    # Dark Blue: #002D5A
    # Medium Blue: #4C7FB8
    # Light Blue: #A6C5E8
    # Light Orange: #FFB86F
    # Medium Purple (Dark): #65498C
    # Medium Purple (Medium): #8464A5
    # Medium Purple (Light): #A386BB
cols = [:white,  
    colorant"#FFB86F", # 1 Light Orange
    colorant"#A6C5E8", # 2 Light Blue
    colorant"#A386BB", # 3 Medium Purple (Light)
    colorant"#4C7FB8", # 4 Medium Blue
    colorant"#8464A5", # 5 Medium Purple (Medium)
    colorant"#4C7FB8", # 6 Medium Blue
    colorant"#8464A5", # 7 Medium Purple (Medium)
    colorant"#002D5A", # 8 Dark Blue
    colorant"#65498C", # 9 Medium Purple (Dark)
    colorant"#002D5A", # 10 Dark Blue
    colorant"#65498C", # 11 Medium Purple (Dark)
    colorant"#002D5A", # 12 Dark Blue
    colorant"#65498C", # 13 Medium Purple (Dark)
    colorant"#002D5A", # 14 Dark Blue
    colorant"#65498C", # 15 Medium Purple (Dark)
    colorant"#65498C", # 16
    ] 
@df df |> 
    (df -> stack(df, Not(:ev), variable_name = :YearMonth, value_name = :Area)) |>
    (df -> DataFrames.transform(df, :YearMonth => ByRow(x -> collect((x[1:4],x[5:6]))) => [:Year, :Month])) |>
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100 .* (2021-1950+1)) => :Area_pc) |>
    (df -> subset(df, :ev => x-> x.>0 .&& x.<16)) 
    ) plot(:YearMonth, :Area_pc,
     group = :ev, legend = :outerright, lw = 1,
    xlabel = "Year", 
    ylabel = "Percentage of monthly days and land area affected",
    size=(800,500), dpi=300, left_margin = (5, :mm),
    colour = cols[2:end]',)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/monthly_landArea_by_Int8.png")


function macrotype(x)
    if x == 1
        y = "Only hot"
    elseif iseven(x)
        y = "Only dry (any)"
    else
        y = "Hot and dry"
    end
    return y
end

# prepare data for plotting
dfp = df |> 
    # pivot_longer
    (df -> stack(df, Not(:ev), variable_name = :YearMonth, value_name = :Area)) |>
    # :YearMonth => :Year, :Month
    (df -> DataFrames.transform(df, :YearMonth => ByRow(x -> collect((x[1:4],x[5:6]))) => [:Year, :Month])) |>
    # relative land area
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100 ) => :Area_pc)) |>
    # filter
    (df -> subset(df, :ev => x-> x.>0 .&& x.<16)) |>
    # :ev to UInt8
    (df -> DataFrames.transform(df, :ev => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # macro type
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> map(macrotype,x)) => :MacroType)) |>
    # group by year, month and MacroType
    (df -> groupby(df, [:YearMonth, :MacroType])) |>
    # sum area
    (gdf -> combine(gdf, :Area_pc => sum))

    # plot
@df dfp plot(:YearMonth, :Area_pc_sum, group = :MacroType, 
    legend = :top, lw = 1,
    xlabel = "Year", 
    ylabel = "Percentage of monthly days and land area affected",
    size=(800,450), dpi=300, left_margin = (5, :mm),
    colour = [colorant"#65498C" colorant"#002D5A" colorant"#FFB86F" ])
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/monthly_landArea_by_macroType.png")

for type in levels(dfp.MacroType)
    println(sort(subset(dfp,:MacroType => x -> x .== type), :Area_pc_sum, rev = true)[1:20,:])
end

# group by indicator
dfp1 = df |>
    # pivot_longer
    (df -> stack(df, Not(:ev), variable_name = :YearMonth, value_name = :Area)) |>
    # relative land area
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100 ) => :Area_pc)) |>
    # filter
    (df -> subset(df, :ev => x-> x.>0 .&& x.<16)) |>
    # :ev to UInt8
    (df -> DataFrames.transform(df, :ev => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # add variables
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x01) > 0)) => :heat)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x02) > 0)) => :d30)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x04) > 0)) => :d90)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x08) > 0)) => :d180))

function single_event_plot(df::DataFrame, var::Symbol, g::Symbol)
    df = df |> 
        (df -> subset(df, var)) |>
        (df -> groupby(df, g)) |> 
        (gdf -> combine(gdf, :Area_pc => sum))
    @df df  plot(g, :Area_pc_sum, label = String(var))
    print("Sum area over time for " * String(var) * ": ")
    print(sum(df.Area_sum))
end

function single_event_plot!(p::Plots.Plot, df::DataFrame, var::Symbol, g::Symbol, color::Color)
    df = df |> 
        (df -> subset(df, var)) |>
        (df -> groupby(df, g)) |> 
        (gdf -> combine(gdf, :Area_pc => sum))
    @df df plot!(p, :YearMonth, :Area_pc_sum, 
        label = String(var), legend=:top,
        color = color
        )
    print("Sum area over time for " * String(var) * ": ")
    println(sum(df.Area_pc_sum))
    println(sort(df, :Area_pc_sum, rev = true)[1:10,:])
    return p
end
p = plot(xlabel = "YearMonth", ylabel = "Percentage of monthly days and land area affected",
    size=(800,450), dpi=300, left_margin = (5, :mm),
    );
color_palette = [colorant"#FFB86F", colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]
for (var, color) in zip([:heat, :d30, :d90, :d180], color_palette)
    single_event_plot!(p, dfp1, var, :YearMonth, color)
end
p
savefig(p, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/monthly_landArea_by_EventType.png")

# sum area over time is correctly always the same
# there seems to be something wrong with the ERA5 extension (1950-1978)

# map plots of extreme months. Do they match analyses from EDO?

function mapplot(data;kwargs...)
    data = data[[721:1440;1:720],end:-1:1]
    heatmap(data';kwargs...)
end

lsmmask = lsm_notime.lsm[:,:];

hwyear = ds.layer[time=1952:1952][:,:,:];
m = sum(==(0x01),hwyear,dims=1);
mapplot(m[1,:,:] .* lsmmask, title = "Hot days on land in 1952", dpi = 300)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/HotDays1952.png")


dyear = ds.layer[time=1951:1951][:,:,:];
md = sum(x-> x >=(0x01) && x < (0x10),dyear,dims=1);
mapplot(md[1,:,:] .* lsmmask, title = "Dry days in 1951", dpi = 300)
savefig("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/DryDays1951.png")




