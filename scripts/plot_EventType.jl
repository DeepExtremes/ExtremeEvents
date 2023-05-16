# Plot area by year by event type
using Revise
using YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates

addprocs(8)
@everywhere begin
    using Pkg
    Pkg.activate("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/ExtremeEvents/")
end
@everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates
@everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
@everywhere function fit1(df)
    df.y = year.(df.time)
    dfg = groupby(df,:y)
    allhists = Dict(i=>WeightedHist(-0.5:1.0:32.5) for i in 1950:2021)
    for k in keys(dfg)
        dfs = dfg[k]
        fit!(allhists[k[1]], dfs.event, cosd.(dfs.latitude) .* lsm > 0.5)
    end
    @show myid()
    allhists
end

ds = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr/")
lsm = open_dataset("Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsm_notime = subsetcube(lsm, time = DateTime("2019-01-01T13:00:00"))

sds = subsetcube(ds, time = (2021:2021), longitude = (), latitude = ())

t = CubeTable(event = ds.layer,
            lsm = lsm_notime.lsm)

annualstats = pmapreduce(mergefun,t) do tab
    fit1(DataFrame(tab))
end
rmprocs(workers())


allx = [value(annualstats[yr]).x[i] for i in 1:17, yr in 1950:2021]
allres = [value(annualstats[yr]).y[i] for i in 1:17, yr in 1950:2021]
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
p = plot(1950:2021,allres_norm[1:end-1,:]',labels=labels,lw=1,size=(800,400),dpi=300)

savefig(p,"n_extremes_land.png")

# save to csv
import CSV
using DataFrames
# not sure about this one
df = DataFrame(allres)
pot = 0.01
ne = 0.1
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/YearlyEventType_ranked_pot" * string(pot) * "_ne" * string(ne) * "_fabian_land.csv"
CSV.write(outname, df)
df = CSV.read(outname)

function macrotype(x)
    if x == 1
        y = "Heatwave"
    elseif iseven(x)
        y = "Drought"
    else
        y = "compound"
    end
    return y
end

# prepare data for plotting
dfp = df |> 
    # pivot_longer
    (df -> stack(df, Not(:ev), variable_name = :Year, value_name = :Area)) |>
    # filter
    (df -> subset(df, :ev => x-> x.>0 .&& x.<16)) |>
    # :ev to UInt8
    (df -> transform(df, :ev => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # macro type
    (df -> transform(df, :EventType => ByRow(x -> map(macrotype,x)) => :MacroType)) |>
    # group by year and MacroType
    (df -> groupby(df, [:Year, :MacroType])) |>
    # sum area
    (gdf -> combine(gdf, :Area => sum))

    # plot
using StatsPlots
@df dfp plot(:Year, :Area_sum, group = :MacroType, legend = :top)

# group by indicator
dfp1 = df |>
    # pivot_longer
    (df -> stack(df, Not(:ev), variable_name = :Year, value_name = :Area)) |>
    # filter
    (df -> subset(df, :ev => x-> x.>0 .&& x.<16)) |>
    # :ev to UInt8
    (df -> transform(df, :ev => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # add variables
    (df -> transform(df, :EventType => ByRow(x -> ((x & 0x01) > 0)) => :heat)) |>
    (df -> transform(df, :EventType => ByRow(x -> ((x & 0x02) > 0)) => :d30)) |>
    (df -> transform(df, :EventType => ByRow(x -> ((x & 0x04) > 0)) => :d90)) |>
    (df -> transform(df, :EventType => ByRow(x -> ((x & 0x08) > 0)) => :d180))

function single_event_plot(df::DataFrame, var::Symbol)
    df = dfp1 |> 
        (df -> subset(df, var)) |>
        (df -> groupby(df, :Year)) |> 
        (gdf -> combine(gdf, :Area => sum))
    @df df  plot(:Year, :Area_sum, label = String(var))
    print("Sum area over time for " * String(var) * ": ")
    print(sum(df.Area_sum))
end
single_event_plot(dfp1, :heat)
single_event_plot(dfp1, :d30)
single_event_plot(dfp1, :d90)
single_event_plot(dfp1, :d180)

# sum area over time is correctly always the same