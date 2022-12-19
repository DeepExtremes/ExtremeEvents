## connected components
using YAXArrays, EarthDataLab
using Statistics, DiskArrays
using ImageMorphology

pot=0.01
ne=0.1
period = 2016:2021 # could be added to file name "*"_"*replace("$period", ":" => "_")*"
#outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_sreld335.zarr"
outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_Sdiam_T3.zarr"
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
# filter function
f4 = function(img)
    # img has size (3,3,5)
    # 2 d diamond
    diamond =  [
        false true false;
        true  true true;
        false true false;
    ];
    img1 = copy(img) .* diamond
    # @show img1
    v = img1[2,2,3]
    # @show typeof(v)
    s = sum(img1) >= 3 
    # @show typeof(s)
    # 3rd D 
    img2 = copy(img[2,2,:]);
    d = all(img2[2:4]) || all(img[1:3]) || all(img[3:5])
    # @show d
    return v && s && d
end
# 
get_diamond = function(dim1::Int)
    Nh = Int((dim1+1)/2);
    range_vec = cat(1:Nh, Nh-1:-1:1, dims = 1);
    out = (range_vec .+ transpose(range_vec)) .> Nh
end
get_diamond = function(dim::Tuple)
    Nh = map(x->Int((x+1)/2), dim);
    range_vec = map(x->cat(1:x, x-1:-1:1, dims = 1), Nh);
    ind = convert(Vector{Int},ones(length(dim)))
    fun = function(y) 
        x = range_vec[y];
        indc = copy(ind); indc[y] = dim[y]; 
        return reshape(x, tuple(indc...))
    end
    vecs = map(fun, 1:length(dim))
    out = reduce(.+, vecs) .> reduce(.+, Nh) - minimum(Nh) #### this is not right. I don't understand how to define the threshold
end


filter_array = function(img)
    dim = size(img);
    centre = map(x->Int(ceil(x/2)), dim)
    # 2D diamond
    diamond = get_diamond(dim[1:2])

end
# apply filter
filteredarray = mapwindow(filter_array, maskarray, (3,3,5));

# use ImageMorphology.label_components to label the connected blobs of events
# use Fabian's version to wrap around the globe
include("../src/label.jl")
# set d to longitude dimension 
# size(maskarray)
# (2192, 1440, 721) # (time, lon, lat)
d = 2;
label_components(filteredarray,wrapdims=(d,));
#r = label_components(filteredarray, strel_diamond((3,3,3)));

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
