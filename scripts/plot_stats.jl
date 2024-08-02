using YAXArrays, OnlineStats, WeightedOnlineStats, Zarr, NetCDF
# using EarthDataLab
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import Statistics #StatsBase
# import Plots
# using Measures
using CairoMakie, GeoMakie

include("../src/plots.jl")

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

trial = "ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022"
landonly = "_landonly"
# events = CSV.read(path * "EventStats_" * trial * landonly * ".csv", DataFrame);
ev = CSV.read(path*"MergedEventStats"*landonly*".csv", DataFrame)
# look for intersection between spatial and temporal range of events from the table or directly in the labelcube
# labelpath = path * "labelcube_$trial.zarr"
labelpath = "/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/MergeLabelCubes/mergedlabels.zarr"
labels = Cube(labelpath) # labels = open_dataset(labelpath )# 
# labels_all = open_dataset(path * "labelcube_ranked_pot0.01_ne0.1_cmp_2016_2021.zarr")
eventcube = Cube(path * "EventCube_ranked_pot0.01_ne0.1.zarr")
lsm = Cube("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")[
    time = At(DateTime("2019-01-01T13:00:00")),
    # region = region,
    ]

# Plots.heatmap(lsm.data'[end:-1:1,:])
# lsmask = lsm.data .> 0.5;

f, ax, h = hist(log10.(ev.volume[1:1000]), 
    color = :tomato,
    label = "land only",
    figure = (; size = (400, 400)),
    axis = (;xlabel = "log10 (Volume)", ylabel = "# events")
    )
save(path * "/fig/events_1000_volume_$trial$landonly.png", f)

function myhexbin(ev::DataFrame)
    # 2d histogram (hexbin)
    
    f = Figure(size = (800,400))
    ax1 =  Axis(f[1,1],
        xlabel = "log10 (Area)",
        ylabel = "duration (days)",)
    hb = hexbin!(ax1, log10.(ev.area), ev.d, cellsize = (0.1, 1.2),
        colorscale=log10)
    Colorbar(f[2,1], hb,
        label = "Number of labelled events",
        vertical = false,
        # height = Relative(0.5)
    )

    # volume over time
    ax2 =  Axis(f[1,2],
        xlabel = "Year of event onset",
        ylabel = "log10(Volume)",)
    hb = hexbin!(ax2, year.(ev.start_time), log10.(ev.volume), cellsize = (1, 0.15),
        colorscale=log10,
        )
    Colorbar(f[2, 2], hb,
        label = "Number of labelled events",
        vertical = false,
        # height = Relative(0.5)
    )
    return f
end
f = myhexbin(ev)
save(path * "fig/events_stats_$trial$(landonly)_1950.png", f,)

# filter from 1970
fev = ev |>
    (df -> transform(df, :duration => (x -> parse.(Int, replace.(x, r" day(s)?" => ""))) => :d)) |>
    (df -> filter([:start_time, :d] => (t, d) -> t .>= DateTime(1970,1,1) .&& d .> 2, df))
f1 = myhexbin(fev)
save(path * "fig/events_stats_$trial$(landonly)_1970.png", f1,)
# extract stats
print(Statistics.quantile(fev.d, [0.05, 0.25, 0.5, 0.75, 0.95]))
# [4.0, 4.0, 4.0, 5.0, 9.0]
print(Statistics.quantile(fev.area, [0.05, 0.25, 0.5, 0.75, 0.95]))
# [0.4146932363510132, 0.8429072499275208, 0.9996573328971863, 2.985333263874054, 24.471434984888276]
print(Statistics.quantile(fev.volume, [0.05, 0.25, 0.5, 0.75, 0.95]))
# [1.7510369420051575, 3.587491035461426, 5.430251181125641, 15.98317289352417, 164.24727816581685]

# volume quantiles by year
gdf = fev |>
    (df -> transform(df, :start_time => (x -> year.(x)) => :yr)) |>
    (df -> DataFrames.groupby(df, :yr)) |>
    (gdf -> combine(gdf) do sdf
        (q05, q25, q50, q75, q95) = Statistics.quantile(sdf.volume, [0.05, 0.25, 0.5, 0.75, 0.95])
        DataFrame(yr = sdf.yr[1], q05 = q05, q25 = q25, q50 = q50, q75 = q75, q95 = q95)
    end)

fig = Figure();
ax = Axis(fig[1, 1], yscale = log10)
for q in [:q05, :q25, :q50, :q75, :q95]
    lines!(ax, gdf[:,:yr], gdf[:, q])
end

fig

# test hm
ls = Cube("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
f,ax = hm(ls.data[:,:,:], ls.axes);
f

sublabels = labels[latitude = -90 .. 90];
# 10158,2018-07-13T00:00:00.0,2018-08-10T00:00:00.0,4.5,57.0,50.0,73.0 # european heatwave
period = Date("2018-07-01") .. Date("2018-09-30")
sublabels = labels[time = period, latitude = -90 .. 90];
### 
# plot where lsmask AND labels always == 0 (max(labels) == 0)
f,ax = hm(labels.data[:,:,:], labels.axes;
   fn = x -> maximum(x) == 0, 
   reduced = :Ti, 
   colormap = :reds,
   );
# # coastlines
# cl=lines!(ax, 
#         GeoMakie.coastlines(),
#         color = :black, linewidth=0.85)
lsmask = (lsm.data .< 0.5)[:,end:-1:1];
# substitute 0 by NaN
rtmp = convert(Array{Float64},lsmask);
replace!(rtmp, 0 => NaN);
x = lookup(lsm.axes, :longitude);
y = lookup(lsm.axes,:latitude)[end:-1:1];
heatmap!(ax, x, y, rtmp; 
    colormap = cgrad(:greys, rev=true),
    # colorbar=:none,
    )
save(path * "fig/nolabel_$trial.png", f)

# period = Date("2018-01-01") .. Date("2018-12-31")
# region = "Germany"
# reg_lon = EarthDataLab.known_regions[region][1] .. EarthDataLab.known_regions[region][3]
# reg_lat = EarthDataLab.known_regions[region][2] .. EarthDataLab.known_regions[region][4]

# DimensionalData.dim2key(subevents.layer.axes)
