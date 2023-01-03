## connected components
using YAXArrays, EarthDataLab
using Statistics, DiskArrays
using ImageMorphology, ImageFiltering


pot=0.01
ne=0.1
period = 2016:2022 # (Date("2018-06-01"), Date("2018-09-16"))#could be added to file name "*"_"*replace("$period", ":" => "_")*"
compound_events = true
cmp = compound_events ? "_cmp" : ""
filter = "_Sdiam3_T5"
# set window and t accordingly
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_sreld335.zarr"
outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * cmp * filter * ".zarr"
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"

inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_ranked_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
# inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
c = Cube(inpath)

clastyears = subsetcube(c, time=period)
# create binary array flagging all events (regardless of type of event) and load to memory
# keep out non-extremes, i.e. where event is < tresne or > 1-tresne
if compound_events
    # compound events are events where both temperature and any of the drought indicators are extremes
    # i.e. where the first bit (coding for temp extreme) is 1 and any other bit is also 1
    function indsel(x)
        x > 0x01 && isodd(x)
    end
    maskarray = permutedims(
        map(indsel, clastyears.data)[:,:,:],
        (2,3,1)
        );
else
    maskarray = permutedims(
        ((clastyears.data .> 0x00) .& (clastyears.data .< 0x10))[:,:,:],
        (2,3,1)
        );
end
# filter maskarray to remove small events

include("../src/detection.jl")

# filter function
# set window
# time window: Tx: x*2-1
window = (3,3,9)
Nh = map(x->Int((x+1)/2), window);
diamondindices = get_diamond_indices(window[1])
t = 0.6;
# diamond in space 60% + at least 3 contiguous times
# 2 d diamond
# diamond =  get_diamond(3);
# get_diamond_indices(diamond) = findall(diamond)
# diamondindices = get_diamond_indices(diamond)

# filter_sp60_t3 = function(img)
#     # img has size (3,3,5)
#     # central value
#     v = img[2,2,3]
#     # @show v
#     # apply diamond spatially on central slice
#     s = sum(diamondindices) do ind
#         view(img,:,:,3)[ind]
#     end
#     # @show s
#     # 3rd D 
#     d = all(view(img,2,2,1:3)) || all(view(img,2,2,2:4)) || all(view(img,2,2,3:5))
#     # @show d
#     return v && s >=3 && d
# end
# filter_sp60_t5 = function(img)
#     # img has size (3,3,7)
#     # central value
#     v = img[2,2,4]
#     # @show v
#     # apply diamond spatially on central slice
#     s = sum(diamondindices) do ind
#         view(img,:,:,4)[ind]
#     end
#     # @show s
#     # 3rd D 
#     d = all(view(img,2,2,1:5)) || all(view(img,2,2,2:6)) || all(view(img,2,2,3:7))
#     # @show d
#     return v && s >=3 && d
# end
# img = rand(10,10,10).>0.5;
# size(img)
# typeof(img)

# @time mapwindow(filter_sp60_t3, img, (3,3,5));
# @time mapwindow(f4, permutedims(img,(3,1,2)), (5,3,3));

# apply filter
# @time filteredarray = mapwindow(filter_sp60_t3, maskarray, (3,3,5));

@time filteredarray = mapwindow(myfilter, maskarray, window);

# axs = clastyears.axes[[2,3,1]]
# include("../src/plots.jl");
# hmf = hm(filteredarray, axs=axs)

# use ImageMorphology.label_components to label the connected blobs of events
# use Fabian's version to wrap around the globe
include("../src/label.jl")
# set d to longitude dimension 
size(filteredarray)
# (2192, 1440, 721) # (time, lon, lat)
d = 1;
@time r = label_components(filteredarray, wrapdims=(d,));
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

# import PlotlyJS
# tmp = reshape(r, (:));
# PlotlyJS.plot(PlotlyJS.histogram(x=tmp[tmp.>0]))
# define YAXArray
labelcube = YAXArray(caxes(clastyears), permutedims(r,(3,1,2)), chunks = DiskArrays.GridChunks(r,(90,120,120)))
# write cube
labelcube = savecube(labelcube,outpath,overwrite=true)

# check
# using Plots
# labelcube = Cube(outpath)
# include("../src/plots.jl")
# simpleplot(labelcube,100, 2018, 4, colours = cgrad(:darkterrain, 50, categorical = true), replacement = 0=>5e3)
