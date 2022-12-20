# Determine min pot to detect historical events
# for each in a list of historic events, 
# subset smoothed cube in time and space
# extract histogram

using YAXArrays, EarthDataLab
using DataFrames
import CSV
using Plots
using PlotlyJS
using Statistics
using StatsBase
using ImageMorphology
using ImageFiltering

inpath_t = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr"
inpath_pei = "/Net/Groups/BGI/scratch/fgans/DeepExtremes/smoothed_pei_ranks.zarr"
s_t = open_dataset(inpath_t)
s_t40 = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed_40_2016.zarr")
s_t80 = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed_80_2016.zarr")
s_pei = open_dataset(inpath_pei)
s_pei40 = open_dataset("/Net/Groups/BGI/scratch/fgans/DeepExtremes/smoothed_pei_ranks_40.zarr")
r_t = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr")
r_pei = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")
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
obs_event = 5# for obs_event = 1:(size(obs)[1])
    period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End]))
    #### !!! need to fix longitude to match 0-360 !!!
    lat = (obs[obs_event,:North], obs[obs_event,:South])
    lon = (obs[obs_event,:West], obs[obs_event,:East])
    # subset smoothed indicators
    subst = subsetcube(s_t, time=period, latitude=lat, longitude=lon)
    subst40 = subsetcube(s_t40, time=period, latitude=lat, longitude=lon)
    subst80 = subsetcube(s_t80, time=period, latitude=lat, longitude=lon)
    subspei = subsetcube(s_pei, time=period, latitude=lat, longitude=lon)
    subspei40 = subsetcube(s_pei40, time=period, latitude=lat, longitude=lon)
    # histogram smoothed
    st =  reshape( subst.layer[:,:,:], :);
    st40 =  reshape( subst40.layer[:,:,:], :);
    sd30 = reshape(subspei.pei_30[:,:,:], :);
    s40d30 = reshape(subspei40.pei_30[:,:,:], :);
    sd90 = reshape(subspei.pei_90[:,:,:], :);
    s40d90 = reshape(subspei40.pei_90[:,:,:], :);
    sd180 = reshape(subspei.pei_180[:,:,:], :);
    h1 = Plots.histogram(st, xlabel = "Smoothed t2mmax", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    h2 = Plots.histogram(sd30, xlabel = "Smoothed PEI_30", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))    
    h3 = Plots.histogram(sd90, xlabel = "Smoothed PEI_90", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    h4 = Plots.histogram(sd180, xlabel = "Smoothed PEI_180", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    # subset ranked indicators
    subrt = subsetcube(r_t, time=period, latitude=lat, longitude=lon)# readcubedata(subsetcube(Cube(r_t)
    subrpei = subsetcube(r_pei, time=period, latitude=lat, longitude=lon)
    # histogram smoothed
    rt =  reshape( subrt.layer[:,:,:], :);
    rd30 = reshape(subrpei.pei_30[:,:,:], :);
    rd90 = reshape(subrpei.pei_90[:,:,:], :);
    rd180 = reshape(subrpei.pei_180[:,:,:], :);
    h1 = Plots.histogram(rt, xlabel = "Ranked t2mmax", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    h2 = Plots.histogram(rd30, xlabel = "Ranked PEI_30", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))    
    h3 = Plots.histogram(rd90, xlabel = "Ranked PEI_90", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    h4 = Plots.histogram(rd180, xlabel = "Ranked PEI_180", title="bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    # compare distribution of ranked and smoothed
    # s1 = PlotlyJS.plot(PlotlyJS.histogram2d(x=rt,y=st), Layout(xaxis_title = "Ranked t2mmax", yaxis_title = "Smoothed t2mmax", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2])))
    s1 = Plots.histogram2d(rt,st)
    # s1f = PlotlyJS.plot(PlotlyJS.histogram2d(x=rt,y=st40), Layout(xaxis_title = "Ranked t2mmax", yaxis_title = "Smoothed t2mmax lbord 40", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2])))
    s1f = Plots.histogram2d(rt,st40, xlabel= "Ranked t2mmax", ylabel = "Smoothed t2mmax lbord 40") 
    # s1c = PlotlyJS.plot(PlotlyJS.histogram2d(x=st,y=st40), Layout(xaxis_title = "Smoothed t2mmax lbord 20", yaxis_title = "Smoothed t2mmax lbord 40", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2])))
    #savefig(s1,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/hist2d_rt_st.png")
    h2d20 = Plots.histogram2d(rd30,sd30, nbins = 100, xlabel = "Ranked PEI_30", ylabel = "Smoothed PEI_30", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    h2d40 = Plots.histogram2d(rd30,s40d30, nbins = 100, xlabel = "Ranked PEI_30", ylabel = "Smoothed PEI_30 (lbord = 40)", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2]))
    savefig(h2d20,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/hist2d_rd30_s20d30.png")
    savefig(h2d40,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/hist2d_rd30_s40d30.png")
# end

# Smoothing is too strong
# we loose the skewness of histogram during extremes. And too many values that were > 0.1 in ranked are <0.01 in smoothed.
# would then be better to work on ranked, but set more conditions to labelling
# e.g. element could only be labelled if > 50% neighbouring pixels are 1
# that would mean creating another spatial filter
# opening?
# tmpt =  [0 0 1 1 0 0 0;
#         1 1 0 0 1 0 0;
#         0 1 1 1 1 0 1;
#         0 0 0 1 1 1 1;
#         0 0 0 0 1 0 1;]
# tmp0 = reshape(subrt.layer[:,:,:],(:, 1, size(subrt.time)[1])) .<= 0.01;
# rtmp0 = reduce(+, tmp0, dims=1);
# Plots.bar(subrt.time[:],reshape(rtmp0, (:)))
# tmp = subrt.layer[:, :, 40] .<= 0.01;
# Plots.heatmap(tmp, title="Ranked <=0.01: no filter")
# tmps = subst.layer[:, :, 40] .<= 0.01;
# Plots.heatmap(tmps)
# # erode(tmp, strel_diamond(tmp))
# tmp1 = opening(tmp, strel_diamond(tmp)) .* tmp;
# Plots.heatmap(tmp1, title="Ranked <=0.01: opening")
# sum(tmp) # dilate(erode(tmp))
# sum(tmp1)
# # thinning(iszero.(tmp)) # 
# f = function(img)
#     #@show img
#     diamond = [
#         0 1 0;
#         1 1 1;
#         0 1 0;
#     ]
#     img1 = copy(img) .* diamond
#     # @show img1
#     v = img1[2,2]
#     # @show v
#     s = sum(img1) >= 3 
#     # @show s
#     return s && v
# end
# tmp2 = mapwindow(f, tmp, (3,3));
# sum(tmp2)
# Plots.heatmap(tmp2, title="Ranked <=0.01 : mapwindow diamond >= 3")
# import StatsBase
# tmp3 = mapwindow(StatsBase.mode, tmp, (3,3)) .* tmp;
# Plots.heatmap(tmp3, title="Ranked <=0.01 : mapwindow mode")

## with time dimension on top....
axes_rt = subrt.layer.axes
tmp = subrt.layer[:, :, :] .<= 0.01;


hmr = hm(tmp, title = "Ranked <=0.01: original")
Plots.savefig(hmr, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax")
# smoothed
axes_st = subst40.layer.axes
tmps40 = subst40.layer[:,:,:] .<= 0.01;
hms40 = hm(tmps40, axs = axes_st, title = "Smoothed (lbord = 40) <= 0.01")
Plots.savefig(hms40, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_smoothed40_tmax")
tmps80 = subst80.layer[:,:,:] .<= 0.01;
hms80 = hm(tmps80, axs = axes_st, title = "Smoothed (lbord = 80) <= 0.01")
Plots.savefig(hms80, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_smoothed80_tmax")
# I'm not sure I trust my hm function to do exactly what I want....

# opening
tmp1 = opening(tmp, strel_diamond(tmp)) .* tmp;
hm1 = hm(tmp1, title="Ranked <=0.01: diamond opening")
Plots.savefig(hm1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax_opening")
# diamond 60%

f3 = function(img)
    diamond = get_diamond((5,3,3))
    img1 = copy(img) .* diamond
    v = img1[2,2,3]
    s = sum(img1) >= sum(diamond)*.6 
    return s && v
end
tmp2 = mapwindow(f3, tmp, (5,3,3));
hm2 = hm(tmp2, title="Ranked <=0.01: mapwindow diamond (5,3,3) >= 60%")
Plots.savefig(hm2, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax_diamond60.png")
# mode
tmp3 = mapwindow(StatsBase.mode, tmp, (5,3,3)) .* tmp;
hm3 = hm(tmp3, title="Ranked <=0.01 : mapwindow mode (5,3,3)")
Plots.savefig(hm3, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax_mode335")
# diamond in space 60% + at least 3 contiguous times
f4 = function(img)
    img = permutedims(img, (2,3,1))
    # @show img
    # img has size (3,3,5)
    # central value
    v = copy(img[2,2,3])
    # @show typeof(v)
    # 2 d diamond
    diamond =  get_diamond(3);
    #@show diamond
    # apply diamond spatially on central slice
    img1 = copy(img[:,:,3]) .& diamond
    #@show img1
    s = sum(img1) >= 3
    # @show typeof(s)
    # 3rd D 
    img2 = copy(img[2,2,:]);
    d = all(img[1:3]) || all(img2[2:4]) || all(img[3:5])
    # @show d
    return v && s && d
end
tmp4 = mapwindow(f4, tmp, (5,3,3));
hm4 = hm(tmp4, title="Ranked <=0.01: mapwindow spatial diamond >= 3 \n AND 3 contiguous times")
Plots.savefig(hm4, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax_spdiam3_t3")
# diamond in space 60% + at least 5 contiguous times
f5 = function(img)
    img = permutedims(img, (2,3,1))
    # img has size (3,3,7)
    # central value
    v = copy(img[2,2,4])
    # @show typeof(v)
    # 2 d diamond
    diamond =  get_diamond(3);
    #@show diamond
    # apply diamond spatially on central slice
    img1 = copy(img[:,:,3]) .& diamond
    #@show img1
    s = sum(img1) >= 3
    # @show typeof(s)
    # 3rd D 
    img2 = copy(img[2,2,:]);
    d = all(img[1:5]) || all(img2[2:6]) || all(img[3:7])
    # @show d
    return v && s && d
end

tmp5 = mapwindow(f5, tmp, (7,3,3));
hm5 = hm(tmp5, title="Ranked <=0.01: mapwindow spatial diamond (3,3) >= 3 \n AND 5 contiguous times")
Plots.savefig(hm5, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/heatmap_EUsummer2018_ranked_tmax_spdiam3_t5.png")

# trend of tmax_smoothed
jena_tmax = subsetcube(r_t, lat = 50.9, lon = 11.5)
time_axis = jena_tmax.time[:]
p1 = PlotlyJS.plot(PlotlyJS.histogram2d(x=time_axis,y=jena_tmax.layer[:]))


# label components tests
temp = zeros(3,3,5);
temp[2,2,[1,5]].=1;
# temp[2,:,[2,4]].=1;
# temp[:,2,[2,4]].=1;
temp[:,:,3] .= 1;
temp
label_components(temp, strel_diamond((3,3,5)))
# components that do not touch but are in the the same strucyiral elements are connected.
