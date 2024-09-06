using YAXArrays, OnlineStats, WeightedOnlineStats, Zarr, NetCDF
# using EarthDataLab
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
using StatsBase # Statistics #
using QuantileRegressions
# import Plots
# using Measures
using CairoMakie, GeoMakie
using Printf: @sprintf

include("../src/plots.jl")
include("../src/stats.jl")

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

# f, ax, h = hist(log10.(ev.volume[1:1000]), 
#     color = :tomato,
#     label = "land only",
#     figure = (; size = (400, 400)),
#     axis = (;xlabel = "log10 (Volume)", ylabel = "# events")
#     )
# save(path * "/fig/events_1000_volume_$trial$landonly.png", f)

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
    return f, ax1, ax2
end
f = myhexbin(ev)
save(path * "fig/events_stats_$trial$(landonly)_1950.png", f,)

# filter from 1970
fev = ev |>
    (df -> transform(df, :duration => (x -> parse.(Int, replace.(x, r" day(s)?" => ""))) => :d)) |>
    (df -> filter([:start_time, :d] => (t, d) -> t .>= DateTime(1970,1,1) .&& d .> 2, df))
f1,ax1,ax2 = myhexbin(fev)

# extract stats
print(quantile(fev.d, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
# [4.0, 4.0, 4.0, 4.0, 4.0, 5.0, 7.0, 9.0, 15.0]
print(quantile(fev.area, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
# [0.22920039296150208, 0.4146932363510132, 0.5628049373626709, 0.8429072499275208, 0.9996573328971863, 2.985333263874054, 10.421100580692293, 24.471434984888276, 115.2952259035902]
print(quantile(fev.volume, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]))
# [0.9846131801605225, 1.7510369420051575, 2.4072601795196533, 3.587491035461426, 5.430251181125641, 15.98317289352417, 63.19659821987154, 164.24727816581685, 1082.8372119736412]

# volume quantiles by year
gdf = fev |>
    (df -> transform(df, :start_time => (x -> year.(x)) => :yr)) |>
    (df -> DataFrames.groupby(df, :yr)) |>
    (gdf -> combine(gdf) do sdf
        (q01, q05, q10, q25, q50, q75, q90, q95, q99) = quantile(sdf.volume, [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99])
        DataFrame(yr = sdf.yr[1], q01 = q01, q05 = q05, q10 = q10, q25 = q25, q50 = q50, q75 = q75, q90 = q90, q95 = q95, q99 = q99)
    end)

# Quantile Regression
qrres = Dict()
for q in [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]
    qrres["q$(@sprintf("%0.2f",q))"] = 
        qreg(@formula(log10(volume)~yr), fev |> 
            (df -> transform(df, :start_time => (x -> year.(x)) => :yr)),
            q, IP())
end
DataPlot = reduce(vcat, [[coeftable(qrres["q$(@sprintf("%0.2f",q))"]).cols[2][2] coeftable(qrres["q$(@sprintf("%0.2f",q))"]).cols[2][1]] for q in [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99]])

# Bands
# band!(ax2,  gdf.yr, DataPlot[1,1] .* gdf.yr .+ DataPlot[1,2], DataPlot[7,1] .* gdf.yr .+ DataPlot[7,2], color = (:grey85, 0.3), label = "[p1, p99]")
# band!(ax2,  gdf.yr, DataPlot[2,1] .* gdf.yr .+ DataPlot[2,2], DataPlot[6,1] .* gdf.yr .+ DataPlot[6,2], color = (:grey90, 0.3), label = "[p5, p95]")
# band!(ax2,  gdf.yr, DataPlot[3,1] .* gdf.yr .+ DataPlot[3,2], DataPlot[5,1] .* gdf.yr .+ DataPlot[5,2], color = (:grey95, 0.3), label = "[p25, p75]")
# lines!(ax2, gdf.yr, DataPlot[4,1] .* gdf.yr .+ DataPlot[4,2], color = (:white, 0.5), label = "p50")
# Legend(f1[3,2],ax2, orientation = :horizontal, backgroundcolor = colorant"#7ad151ff", framecolor = colorant"#7ad151ff")
# # fig
# f1

# # lines
# lines!(ax2, gdf.yr, DataPlot[1,1] .* gdf.yr .+ DataPlot[1,2], color = (:grey80, 0.8), linestyle = (:dash, :loose))
# lines!(ax2, gdf.yr, DataPlot[2,1] .* gdf.yr .+ DataPlot[2,2], color = (:grey80, 0.8), linestyle = (:dash, :normal))
# lines!(ax2, gdf.yr, DataPlot[3,1] .* gdf.yr .+ DataPlot[3,2], color = (:grey80, 0.8), linestyle = (:dash, :dense))
# lines!(ax2, gdf.yr, DataPlot[4,1] .* gdf.yr .+ DataPlot[4,2], color = (:grey80, 0.8), linestyle = :solid)
# lines!(ax2, gdf.yr, DataPlot[5,1] .* gdf.yr .+ DataPlot[5,2], color = (:grey80, 0.8), linestyle = (:dash, :dense))
# lines!(ax2, gdf.yr, DataPlot[6,1] .* gdf.yr .+ DataPlot[6,2], color = (:grey80, 0.8), linestyle = (:dash, :normal))
# lines!(ax2, gdf.yr, DataPlot[7,1] .* gdf.yr .+ DataPlot[7,2], color = (:grey80, 0.8), linestyle = (:dash, :loose))
# f1

# Largest events
lev = fev |>
    (df -> sort(df, :volume, rev = true)) |>
    (df -> select(df, [:label, :start_time, :end_time, :d, :area, :volume])) |>
    (df -> first(df,5))
text!(ax1, log10.(lev.area), lev.d, text = " <-- " .* string.(lev.label), 
    align = (:left, :center),
    fontsize = 8, 
    color = colorant"#22a884ff",
    rotation = π/6,
    label = "largest events",
    )

levd = fev |>
    (df -> sort(df, :d, rev = true)) |>
    (df -> select(df, [:label, :start_time, :end_time, :d, :area, :volume])) |>
    (df -> first(df,5))
text!(ax1, log10.(levd.area), levd.d, text = " <-- " .* string.(levd.label), 
    align = (:left, :center), 
    fontsize = 8, 
    color = colorant"#414487ff",
    rotation = π/6,
    label = "longest events")
# Legend(f1[3,1], ax2, orientation = :horizontal)
# plot point to increase limits
scatter!(ax1, 3.6,63, color = :white)
f1
save(path * "fig/events_stats_$trial$(landonly)_1970.png", f1,)

# export table
# largest volume
df = fev |>
           (df -> sort(df, :volume, rev = true)) |>
           (df -> first(df,10));
df.start_date = Date.(df.start_time);
df.end_date = Date.(df.end_time);
show(stdout, MIME("text/latex"), select(df, [:label, :start_date, :end_date, :longitude_min,  :longitude_max,  :latitude_min,  :latitude_max, :duration, :area, :volume]))
largest = df.label'
# largest = [42561  51252  24092  55983  55755  25632  18958  44770  36790  53015]

# longest duration
df = fev |>
           (df -> sort(df, :d, rev = true)) |>
           (df -> first(df,10));
df.start_date = Date.(df.start_time);
df.end_date = Date.(df.end_time);
show(stdout, MIME("text/latex"), select(df, [:label, :start_date, :end_date, :longitude_min,  :longitude_max,  :latitude_min,  :latitude_max, :duration, :area, :volume]))
longest = df.label'
# longest = [31767  31866  49981  44843  18998  24109  42561  28223  50340  54071]

# plot largest events
# llabels = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/largest_longest_idx_labels.zarr")
lla = convert(Array{Float64},llabels.largest);
replace!(lla, 0 => NaN);
lon = lookup(llabels.largest, :longitude);
lat = lookup(llabels.longest, :latitude);
fig = Figure();
gax = GeoAxis(fig[1,1], 
        # limits=(extrema(lon), extrema(lat)), 
        source="+proj=latlong +datum=WGS84", # src CRS
        dest="+proj=eqearth", # destination CRS, in which you want to plot
        # coastlines = true # plot coastlines from Natural Earth, as a reference.
    )
h = heatmap!(gax, lon, lat, lla; 
    colormap = cgrad(:tab10, categorical = true),
    colorrange = (1,10),
    # colormap = Makie.Categorical(:viridis)
    )
# coastlines
cl=lines!(gax, 
    GeoMakie.coastlines(),
    color = :black, linewidth=0.85)
translate!(cl, 0, 0, 1000)

cbar = Colorbar(fig[2,1], h, 
    label = "Largest dry and hot events (1970 - 2022)",
    vertical = false,)
nlb = length(largest)
cbar.limits = (1,nlb)
cbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), [string(largest[1,i]) for i in 1:10])
            

# remove gridlines
gax.xgridcolor[] = colorant"transparent";
gax.ygridcolor[] = colorant"transparent";
gax.xticklabelsvisible = false;
gax.yticklabelsvisible = false;
fig
save(path * "fig/largest_$trial$(landonly)_1970.png", fig, size = (800, 494))

# plot longest events
llo = convert(Array{Float64},llabels.longest);
replace!(llo, 0 => NaN);
fig = Figure();
gax = GeoAxis(fig[1,1], 
        # limits=(extrema(lon), extrema(lat)), 
        source="+proj=latlong +datum=WGS84", # src CRS
        dest="+proj=eqearth", # destination CRS, in which you want to plot
        # coastlines = true # plot coastlines from Natural Earth, as a reference.
    )
h = heatmap!(gax, lon, lat, llo; 
    colormap = cgrad(:tab10, categorical = true),
    colorrange = (1,10),
    # colormap = Makie.Categorical(:viridis)
    )
# coastlines
cl=lines!(gax, 
    GeoMakie.coastlines(),
    color = :black, linewidth=0.85)
translate!(cl, 0, 0, 1000)

cbar = Colorbar(fig[2,1], h, 
    label = "Longest dry and hot events (1970 - 2022)",
    vertical = false,)
nlb = length(longest)
cbar.limits = (1,nlb)
cbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), [string(longest[1,i]) for i in 1:10])

# remove gridlines
gax.xgridcolor[] = colorant"transparent";
gax.ygridcolor[] = colorant"transparent";
gax.xticklabelsvisible = false;
gax.yticklabelsvisible = false;
fig
save(path * "fig/longest_$trial$(landonly)_1970.png", fig, size = (800, 494))


# # test hm
# f,ax = hm(lsm.data[:,:], lsm.axes);
# f

# # 10158,2018-07-13T00:00:00.0,2018-08-10T00:00:00.0,4.5,57.0,50.0,73.0 # european heatwave
# period = Date("2018-07-01") .. Date("2018-09-30")
# sublabels = labels[time = period, latitude = -90 .. 90];
# ### 
# # plot where lsmask AND labels always == 0 (max(labels) == 0)
# f,ax = hm(labels.data[:,:,:], labels.axes;
#    fn = x -> maximum(x) == 0, 
#    reduced = :Ti, 
#    colormap = :reds,
#    );
# # # coastlines
# # cl=lines!(ax, 
# #         GeoMakie.coastlines(),
# #         color = :black, linewidth=0.85)
# lsmask = (lsm.data .< 0.5)[:,end:-1:1];
# # substitute 0 by NaN
# rtmp = convert(Array{Float64},lsmask);
# replace!(rtmp, 0 => NaN);
# x = lookup(lsm.axes, :longitude);
# y = lookup(lsm.axes,:latitude)[end:-1:1];
# heatmap!(ax, x, y, rtmp; 
#     colormap = cgrad(:greys, rev=true),
#     # colorbar=:none,
#     )
# save(path * "fig/nolabel_$trial.png", f)

# # period = Date("2018-01-01") .. Date("2018-12-31")
# # region = "Germany"
# # reg_lon = EarthDataLab.known_regions[region][1] .. EarthDataLab.known_regions[region][3]
# # reg_lat = EarthDataLab.known_regions[region][2] .. EarthDataLab.known_regions[region][4]

# # DimensionalData.dim2key(subevents.layer.axes)
