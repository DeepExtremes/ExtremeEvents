using YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using DataFrames, Dates
import CSV
import StatsBase
using Plots, Measures

include("../src/plots.jl")

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/"
end

trial = "ranked_pot0.01_ne0.1_cmp_2016_2021"
landonly = "_landonly"
events = CSV.read(path * "EventStats_" * trial * landonly * ".csv", DataFrame)
# look for intersection between spatial and temporal range of events from the table or directly in the labelcube
labelpath = path * "labelcube_$trial.zarr"
# labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_$trial.zarr"
labels = Cube(labelpath) # labels = open_dataset(labelpath )# 
# labels_all = open_dataset(path * "labelcube_ranked_pot0.01_ne0.1_cmp_2016_2021.zarr")
eventcube = open_dataset(path * "EventCube_ranked_pot0.01_ne0.1.zarr")
lsm = subsetcube(
    Cube("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc"),
    time=DateTime("2019-01-01T13:00:00"),
    # region = region,
    )
Plots.heatmap(lsm.data'[end:-1:1,:])
lsmask = lsm.data .> 0.5;

p = Plots.histogram(log10.(events.volume[1:1000]), 
    # yscale = :log10,
    xlabel = "log10 (Volume)",
    ylabel = "# events",
    label = "land only"
    )
savefig(path * "events_1000_volume_$trial$landonly.png")


days = parse.(Int, replace.(events.duration, r" day(s)?" => "")); days[1:10]
p1 = Plots.scatter(log10.(events.area[1001:end]), days[1001:end],
    label = "smaller events",
    xlabel = "log10 (Area)",
    ylabel = "duration",
    color = :grey60,
    markershape = :circle,
    markersize = 2,
    markerstrokecolor = :grey60,
    )
Plots.scatter!(log10.(events.area[1:1000]), days[1:1000],
    label="1000 largest events in Volume", 
    color = :grey30,
    markershape = :circle,
    markersize = 2,
    markerstrokecolor = :grey30,)
savefig(p1,path * "events_area_vs_duration_$trial$landonly.png")

Plots.scatter(log10.(events.area[1:1000]), days[1:1000],color=:red,label="1000 largest events in Volume")

sublabels = subsetcube(labels, latitude = (90, -90))

# plot where lsmask AND labels always == 0 (max(labels) == 0)
p = hm(sublabels.data[:,:,:];
   axs = sublabels.axes, 
   fn = x -> maximum(x) == 0, 
   reduced = "tim",
   c = palette(:reds),
   );
lsmask = (lsm.data .< 0.5)[:,:]'[end:-1:1,:];
# substitute 0 by NaN
rtmp = convert(Array{Float64},lsmask);
replace!(rtmp, 0 => NaN)
x = sublabels.axes[2][:];
y = sublabels.axes[3][end:-1:1];
Plots.heatmap!(x, y, rtmp, 
    c = cgrad(:inferno, rev=true),
    colorbar=:none,)
savefig(path * "nolabel_$trial.png")

period = (Date("2019-07-20"), Date("2019-07-30"))
region = "Luxembourg"

# plot where eventcube is always ==16 or at least where it is never 1 to 15 (always 0 or 16) + mask out ocean
subevents = subsetcube(eventcube, variable="layer", time = period, region = region)
function noev(output,x) 
    output = all((x .== 0) .| (x .== 16))
end
Indims = InDims("time")
Outdims = OutDims()
# lst_monthly_high = mapCube(median_by_index, lst, indims = Indims, outdims = Outdims; index_list = index_in_cube, showprog = true)
no_events = mapCube(noev, subevents, indims = Indims, outdims = Outdims)
no_events = mapslices(noev, subevents, dims = "time")

# plot labelcube flat layer of timestep count and number of events
# for each lon, lat, reduce over time: 
# list of unique labels -> length
# count -> sum
# # maybe easier to do it through a tableCube and have results as a nested df
# ct = CubeTable(label = labels[time = period, region = region],
#     landmask = lsm[region = region]
#     )
# ct[1]

sublabels = labels[time = period, region = region]

function freqlabels(itr, latitude)
    vals = unique(itr)
    if length(vals) == 1
        return vals[1]
    end
    freq = map(vals) do v
        count(itr .== v) .* cosd(latitude)
    end
end

Indims = InDims("time")
Outdims = OutDims()
# lst_monthly_high = mapCube(median_by_index, lst, indims = Indims, outdims = Outdims; index_list = index_in_cube, showprog = true)
no_events = mapCube(noev, subevents)

# testing
zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
t2m = subsetcube(era, var = "t2m", region="Luxembourg", time=period)
tmp = mapslices(mean ∘ skipmissing, t2m, dims="time")

function mymean(output, pixel)
    # @show pixel
    @show output = (mean ∘ skipmissing)(pixel)
end
indims = InDims("time")
outdims = OutDims()
tmp1 = mapCube(mymean, t2m, indims=indims, outdims=outdims)
tmp1.data[:,:]
# ????
# 4×3 Matrix{Union{Missing, Float64}}:
#  missing  missing  missing
#  missing  missing  missing
#  missing  missing  missing
#  missing  missing  missing

print(events[1:10,1:10])
print(events[1:10,11:20])
print(events[1:10,21:end])