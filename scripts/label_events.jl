## connected components
using YAXArrays, EarthDataLab
using Statistics, DiskArrays
using ImageMorphology

pot=0.01
ne=0.1
period = 2016:2021 # could be added to file name "*"_"*replace("$period", ":" => "_")*"
outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_sreld335.zarr"
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"

inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_ranked_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
# inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
c = Cube(inpath)

clastyears = subsetcube(c, time=period)
#clastyears = subsetcube(c, time=(Date("2019-01-01"), Date("2019-12-31")), region="Europe")
# create binary array flagging all events (regardless of type of event) and load to memory
# keep out non-extremes, i.e. where event is < tresne or > 1-tresne
maskarray = ((clastyears.data .> 0x00) .& (clastyears.data .< 0x10))[:,:,:];
# filter maskarray to remove small events

# use ImageMorphology.label_components to label the connected blobs of events
r = label_components(maskarray, strel_diamond((3,3,5)));
# import PlotlyJS
# tmp = reshape(r, (:));
# PlotlyJS.plot(PlotlyJS.histogram(x=tmp[tmp.>0]))
# define YAXArray
labelcube = YAXArray(caxes(clastyears), r, chunks = DiskArrays.GridChunks(r,(90,120,120)))
# write cube
labelcube = savecube(labelcube,outpath,overwrite=true)

# check
using Plots
labelcube = Cube(outpath)
include("../src/plots.jl")
simpleplot(labelcube,182, 2019, 4, colours = cgrad(:darkterrain, 50, categorical = true), replacement = 0=>1e5)
