## connected components
using YAXArrays, EarthDataLab
using Statistics, DiskArrays
using ImageMorphology#, ImageFiltering


pot=0.01
ne=0.1
period = 1950:2022 #2016:2022 # (Date("2018-06-01"), Date("2018-09-16"))#could be added to file name "*"_"*replace("$period", ":" => "_")*"
aperiod = "_1950_2022" #"_2016_2021"
compound_events = true
cmp = compound_events ? "_cmp" : ""
filter_events = false
filter = filter_events ? "_Sdiam3_T5" : ""
filter_land = false
land = filter_land ? "_land" : ""
# region = "Italy"

# set window and t accordingly
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_sreld335.zarr"
outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * cmp * filter * aperiod * land * ".zarr"
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"

inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/EventCube_ranked_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
# inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
c = Cube(inpath)

clastyears = subsetcube(c, time=period)
# create binary array flagging all events (regardless of type of event) and load to memory
# keep out non-extremes, i.e. where event is < tresne or > 1-tresne
if compound_events
    # compound events are events where both temperature and any of the drought indicators are extremes
    # i.e. where the first bit (coding for temp extreme) is 1 and any other bit is also 1
    # function indsel(x)
    #     x > 0x01 && isodd(x)
    # end
    # @time maskarray = permutedims(
    #     map(indsel, clastyears.data)[:,:,:],
    #     (2,3,1)
        # );
    @time maskarray = permutedims(
        ((broadcast(x->isodd(x),clastyears.data)) .& (clastyears.data .> 0x01))[:,:,:],
        (2,3,1)
        );
    # 136.463942 seconds (14.02 M allocations: 28.601 GiB, 0.68% gc time, 2.64% compilation time)
    # 1950-2021: 2084.692204 seconds (16.20 M allocations: 317.017 GiB, 0.35% gc time, 0.21% compilation time)

else
    @time maskarray = permutedims(
        ((clastyears.data .> 0x00) .& (clastyears.data .< 0x10))[:,:,:],
        (2,3,1)
        );
end

# import Plots
# red = reshape(reduce(+,maskarray, dims=3), size(maskarray)[1:2]);
# size(red)
# Plots.heatmap(red)

if filter_land
    lsm = subsetcube(
    Cube(open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")),
    time=DateTime("2019-01-01T13:00:00"),
    # region = region,
    )
    lsmask = lsm.data .> 0.5;
    @time for t in 1: size(maskarray)[3]
        maskarray[:,:,t] = ifelse.(lsmask, maskarray[:,:,t], 0)
    end
end
# Plots.heatmap(lsmask)
# red1 = reshape(reduce(+,maskarray, dims=3), size(maskarray)[1:2]);
# Plots.heatmap(red1)

if filter_events
    # filter maskarray to remove small events
    include("../src/detection.jl")
    # diamond in space 60% + at least 3 contiguous times
    # myfilter function needs following global variables
    # set window (time window: Tx: x*2-1)
    window = (3,3,9)
    # central get_diamond_indices
    Nh = map(x->Int((x+1)/2), window);
    # indices of 2D diamond
    diamondindices = get_diamond_indices(window[1])
    t = 0.6 * length(diamondindices);
    # use ImageFiltering.mapwindow to apply filter
    @time filteredarray = mapwindow(myfilter, maskarray, window);
    # !!! 32027.975212 seconds (98.20 G allocations: 3.299 TiB, 1.51% gc time, 0.00% compilation time)
else
    filteredarray = maskarray;
end
# axs = clastyears.axes[[2,3,1]]
# include("../src/plots.jl");
# hmf = hm(filteredarray, axs=axs)

# use ImageMorphology.label_components to label the connected blobs of events
# use Fabian's version to wrap around the globe
include("../src/label.jl")
# set d to longitude dimension 
size(filteredarray)
# (lon, lat, time)
d = 1;
@time r = label_components(filteredarray, wrapdims=(d,));
# 8.853066 seconds (159.63 k allocations: 16.965 GiB, 1.89% gc time, 1.12% compilation time)
# 1950-2021 166.081424 seconds (162.74 k allocations: 203.452 GiB, 0.01% gc time, 0.11% compilation time)

# # look at results
# Plots.heatmap(r[:,:,1]'[end:-1:1,:])
# Plots.heatmap(r[:,:,10]'[end:-1:1,:])
# Plots.heatmap(r[:,:,20]'[end:-1:1,:])
# Plots.heatmap(r[:,:,30]'[end:-1:1,:])
# Plots.heatmap(r[:,:,40]'[end:-1:1,:])
# Plots.heatmap(r[:,:,50]'[end:-1:1,:])
# Plots.heatmap(r[:,:,60]'[end:-1:1,:])
# Plots.heatmap(r[:,:,70]'[end:-1:1,:])
# Plots.heatmap(r[:,:,80]'[end:-1:1,:])
# Plots.heatmap(r[:,:,90]'[end:-1:1,:])
# Plots.heatmap(r[:,:,100]'[end:-1:1,:])

# define YAXArray
# TO DO : find a way to give the variable a name!
@time labelcube = YAXArray(caxes(clastyears)[[2, 3, 1]], r, chunks = DiskArrays.GridChunks(r,(120,120,90)))
# @time labelcube = YAXArray(caxes(clastyears), permutedims(r,(3,1,2)), chunks = DiskArrays.GridChunks(r,(90,120,120))) # ERROR: OutOfMemoryError()
# write cube
@time labelcube = savecube(labelcube,outpath,overwrite=true)
print("done!")
# check
# using Plots
# labelcube = Cube(outpath)
# include("../src/plots.jl")
# simpleplot(labelcube,100, 2018, 4, colours = cgrad(:darkterrain, 50, categorical = true), replacement = 0=>5e3)
