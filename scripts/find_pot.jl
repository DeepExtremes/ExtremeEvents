# Determine min pot to detect historical events
# for each in a list of historic events, 
# subset smoothed cube in time and space
# extract histogram

using YAXArrays, EarthDataLab
using DataFrames
import CSV
using Plots
using Statistics

inpath_t = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr"
inpath_pei = "/Net/Groups/BGI/scratch/fgans/DeepExtremes/smoothed_pei_ranks.zarr"
r_t = open_dataset(inpath_t)
r_pei = open_dataset(inpath_pei)
inputs = (r_t.layer,
        r_pei.pei_30, 
        r_pei.pei_90,
        r_pei.pei_180)


# load list of events
obs1 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart1_csv_Ltime.csv", DataFrame)
# remove spaces in columns names. DataConvenience.cleannames! from DataConvenience.jl
# import DataConvenience
# DataConvenience.cleannames!(obs)
sort!(obs1, :when_from, rev=true)
# select drought and heatwave (drop floods)
filter!(:Event_type=>!=("flood"), obs1)


# for each event documented in the literature, subset smoothed_cube according to time and spatial extent, 
# compute weighted "volume" of each type of event, (min, mean, max) for t2mmax, 

# run stats.jl first
include("stats.jl")

# for tres in [0.001, 0.0025, 0.005, 0.01]
#     obs_check = sanity_check(obs1, tres);

#     # join data frames df and obs
#     obs1.rowid = rownumber.(eachrow(obs1));
#     obs_with_check = leftjoin(obs1, obs_check, on=:rowid=>:label )
#     # export to csv
#     CSV.write("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart1_SanityCheck_$tres.csv", obs_with_check)

# end

# obs_with_check[:,[:Name,:Event_type, :Year, :heat, :drought30, :drought90, :drought180, :compound]]

### EventPart2
obs = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart2.csv", DataFrame; header=3)
# subset cube
obs = dropmissing!(obs)
obs.Start .= replace.(obs.Start, r"\." => "-")
obs.End .= replace.(obs.End, r"\." => "-")
for obs_event = 1:(size(obs)[1])
    period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End]))
    #### !!! need to fix longitude to match 0-360 !!!
    lat = (obs[obs_event,:North], obs[obs_event,:South])
    lon = (obs[obs_event,:West], obs[obs_event,:East])
    # subset smoothed indicators
    subt = subsetcube(r_t, time=period, latitude=lat, longitude=lon)
    subpei = subsetcube(r_pei, time=period, latitude=lat, longitude=lon)
    # histogram
    t =  reshape( subt.layer[:,:,:], :);
    d30 = reshape(subpei.pei_30[:,:,:], :);
    h1 = Plots.histogram(t, title="Smoothed t2mmax in bbox Longitude" * string(lon) * " Latitude " * string(lat) * " Time " * string(period))
    h2 = Plots.histogram(d30)

end

    



