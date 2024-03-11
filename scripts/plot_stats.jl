using YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import StatsBase
# import Plots
using Measures
using CairoMakie, GeoMakie

include("../src/plots.jl")

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

trial = "ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022"
landonly = "_landonly"
events = CSV.read(path * "EventStats_" * trial * landonly * ".csv", DataFrame);
# look for intersection between spatial and temporal range of events from the table or directly in the labelcube
labelpath = path * "labelcube_$trial.zarr"
# labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_$trial.zarr"
labels = Cube(labelpath) # labels = open_dataset(labelpath )# 
# labels_all = open_dataset(path * "labelcube_ranked_pot0.01_ne0.1_cmp_2016_2021.zarr")
eventcube = open_dataset(path * "EventCube_ranked_pot0.01_ne0.1.zarr")
lsm = Cube("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")[
    time = At(DateTime("2019-01-01T13:00:00")),
    # region = region,
    ]

# Plots.heatmap(lsm.data'[end:-1:1,:])
# lsmask = lsm.data .> 0.5;

f, ax, h = hist(log10.(events.volume[1:1000]), 
    color = :tomato,
    label = "land only",
    figure = (; size = (400, 400)),
    axis = (;xlabel = "log10 (Volume)", ylabel = "# events")
    )
save(path * "/fig/events_1000_volume_$trial$landonly.png", f)


# days = parse.(Int, replace.(events.duration, r" day(s)?" => "")); days[1:10]
# p1 = Plots.scatter(log10.(events.area[1001:end]), days[1001:end],
#     label = "smaller events",
#     xlabel = "log10 (Area)",
#     ylabel = "duration",
#     color = :grey60,
#     markershape = :circle,
#     markersize = 2,
#     markerstrokecolor = :grey60,
#     )
# Plots.scatter!(log10.(events.area[1:1000]), days[1:1000],
#     label="1000 largest events in Volume", 
#     color = :grey30,
#     markershape = :circle,
#     markersize = 2,
#     markerstrokecolor = :grey30,)
# savefig(p1,path * "fig/events_area_vs_duration_$trial$landonly.png")

# Plots.scatter(log10.(events.area[1:1000]), days[1:1000],color=:red,label="1000 largest events in Volume")

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

period = Date("2018-01-01") .. Date("2018-12-31")
region = "Germany"
reg_lon = EarthDataLab.known_regions[region][1] .. EarthDataLab.known_regions[region][3]
reg_lat = EarthDataLab.known_regions[region][2] .. EarthDataLab.known_regions[region][4]

subevents = eventcube[Ti = period]#, longitude = reg_lon, latitude = reg_lat]
DimensionalData.dim2key(subevents.layer.axes)

# plot where eventcube is always ==16 or at least where it is never 1 to 15 (always 0 or 16) + mask out ocean
function noev!(output,x) 
    output .= all(skipmissing((x .== 0) .| (x .== 16)))
end

# plot yearly number of event days (any event type)
function nev!(output,x) 
    output .= sum(skipmissing((x .!= 0) .& (x .!= 16)))
end

Indims = InDims(:Ti)
Outdims = OutDims()

# @time no_events = mapCube(noev!, subevents.layer, indims = Indims, outdims = Outdims)
# (f,ax) = heatmap(lookup(no_events, :longitude), lookup(no_events, :latitude)[end:-1:1], no_events.data[:,end:-1:1]);
# f

@time n_events = mapCube(nev!, subevents.layer, indims = Indims, outdims = Outdims)
lon = lookup(n_events, :longitude);
lat = lookup(n_events, :latitude)[end:-1:1];

fig, ax, plt = heatmap(lon, lat, n_events.data[:,end:-1:1]; colormap = cmap,
    axis = (; aspect=DataAspect()),
    figure = (; size = (1200,600), fontsize=24))
cbar = Colorbar(fig[1,2], plt, label = "Number of event days")
fig

# projected figure
axs = modaxs(n_events.axes)
# modify n_events: shift lon
shifts = getshifts(axs)
n_events1 = circshift(n_events.data, shifts);
# δlon = (lon[2]-lon[1])/2
# nlon = lon .- 180 .+ δlon;
lon_dim =  dimnum(axs, :longitude);
nlon = collect(axs[lon_dim]);

# cmap = cgrad(:viridis, 100, categorical = true);
cmap = vcat(colorant"grey90", resample_cmap(cgrad(:viridis, 100, categorical = true), 100));
f = Figure(;size=(1200,600));
ax = GeoAxis(f[1,1],
        title = string(period),
        );
h = surface!(ax, nlon, lat, n_events1[:,end:-1:1],
    colormap = cmap,
    colorrange=(0,100),
    highclip=cmap[end],
    shading=false,
    );
cbar = Colorbar(f[1,2], h, label = "Number of event days")
cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :white, linewidth=0.85,
)
translate!(cl, 0, 0, 1000);
# remove gridlines
ax.xgridcolor[] = colorant"transparent";
ax.ygridcolor[] = colorant"transparent";
# ax.xticklabelsvisible = false;
# ax.yticklabelsvisible = false;
f
save(path * "fig/event_count_$(trial)_year$(string(period)[1:4]).png", f)

# Labelled events in given year
sublabels = labels[time = period]
# max label
function mlev!(output,x) 
    output .= maximum(skipmissing(x))
end
@time ml_events = mapCube(mlev!, sublabels, indims = Indims, outdims = Outdims)
data = replace!(ml_events.data[:,:], 0 => missing);
f1,ax1,h = heatmap(lon, lat, data[:,:]; 
    colormap = :Paired_12,
    axis = (; aspect=DataAspect()),
    figure = (; size = (1200,600), fontsize=24));
Colorbar(f1[1,2], h, label = "max event label");
f1
DimensionalData.dim2key(sublabels.axes)

# projected
# shifts are the same
data1 = circshift(data, shifts);
f1 = Figure(;size=(1200,600));
ax = GeoAxis(f1[1,1],
        title = string(period),
        );
s = surface!(ax, nlon, lat, data1[:,end:-1:1],
    colormap = :Paired_12,
    shading=false,
    );
cbar = Colorbar(f1[1,2], s, label = "Labelled events")
cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :grey, linewidth=0.85,
)
translate!(cl, 0, 0, 1000);
# remove gridlines
ax.xgridcolor[] = colorant"transparent";
ax.ygridcolor[] = colorant"transparent";
# ax.xticklabelsvisible = false;
# ax.yticklabelsvisible = false;
f1
save(path * "fig/max_label_$(trial)_year$(string(period)[1:4]).png", f1)
