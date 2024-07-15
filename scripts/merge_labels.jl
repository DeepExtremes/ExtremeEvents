using YAXArrays, Zarr, Dates, ImageMorphology
using DiskArrays: ConcatDiskArray, DiskArrays
using ImageMorphology: DisjointMinSets
using CSV


select_decade(layer,decade) = if decade == 2010
    layer[time=Date(2010,1,1)..Date(2022,12,31)]
else 
    layer[time=Date(decade,1,1)..Date(decade+9,12,31)]
end
function open_decade(decade)
    filename = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/labelcube_ranked_pot0.01_ne0.1_cmp_S1_T3_$(decade)_$(decade+12).zarr"
    arr = open_dataset(filename)
    select_decade(arr.layer,decade)
end

"""
    merge_components(components, nlabels, concat_dim)

Given a list of `arrays` which are the result of a connected component analysis using 
`label_components` merges the components of each array by investigating the edges
and returning a new `minlabel` vector for re-labelling all arrays.
"""
function merge_components(components, nlabels, ::Val{concat_dim}) where concat_dim

    label_offsets = [0;cumsum(nlabels)]
    
    dms = DisjointMinSets{Int}(last(label_offsets))

    for i in 1:length(components)-1
        println("Labelling interface $i and $(i+1)")
        a1 = components[i]
        a2 = components[i+1]
        ind1 = ntuple(i->i==concat_dim ? (size(a1,i):size(a1,i)) : Colon(),ndims(a1))
        ind2 = ntuple(i->i==concat_dim ? (1:1) : Colon(),ndims(a2))
        slice1 = add_labeloffset.(view(a1,ind1...),label_offsets[i])
        slice2 = add_labeloffset.(view(a2,ind2...),label_offsets[i+1])
        slicemerged = cat(Array(slice1), Array(slice2), dims=concat_dim)
        conn_merged = zeros(Int,size(slicemerged))
        label_components!(conn_merged, (slicemerged.>0))
        println("Filling results")
        fill_conn_merged(conn_merged,slicemerged,dms)
    end
    dms,label_offsets
end
add_labeloffset(x,offset) = x==0 ? x : x+offset

function fill_conn_merged(conn_merged,slicemerged,dms)
    conn_merged_roots = Int[]
    for j in 1:length(conn_merged)
        val = conn_merged[j]
        if val > 0
            oldcomp = slicemerged[j]
            if val > length(conn_merged_roots)
                push!(conn_merged_roots,oldcomp)
            else
                if conn_merged_roots[val] != oldcomp
                    ImageMorphology.union!(dms,conn_merged_roots[val],oldcomp)
                end
            end
        end
    end
end

function relabel_component!(labels, relabels, minlabels, offset::Int)
    for i in eachindex(labels)
        oldlabel = labels[i]
        if oldlabel > 0
            labels[i] = minlabels[ImageMorphology.find_root!(relabels, oldlabel+offset)]
        end
    end
end

function compute_nlabels(decades,allars)
    nlabels = map(decades, allars) do dec, ar
        println("Computing max label of $dec")
        maximum(ar)
    end
    CSV.write("nlabels.csv", [(;label = l) for l in nlabels])
end

function create_outdataset(outpath,allars)
    newtimes = YAXArrays.DD.Ti(Date(1950,1,1):Day(1):Date(2022,12,31))
    outar = ConcatDiskArray(reshape(allars,1,1,length(allars)))
    mergedyax = YAXArray((allyax[1].longitude, allyax[1].latitude,newtimes),outar)
    mergedyax = setchunks(mergedyax,(120,120,90))
    savedataset(Dataset(labels = mergedyax),skeleton=true,overwrite = true, path = outpath)
end


decades = 1950:10:2010
outpath = "./mergedlabels.zarr"

allyax = open_decade.(decades)
allars = map(i->i.data,allyax)

compute_nlabels(decades, allars)

create_outdataset(outpath,allars)

ds = open_dataset(zopen(outpath,"w"))
outyax = ds.labels
concat_dim = Val(3)
nlabels = [row.label for row in CSV.Rows("nlabels.csv",types=Int)]

println("Merging components")
relabels,offsets = merge_components(allars, nlabels, concat_dim)

minlabels = ImageMorphology.minlabel(relabels)

foreach(allars,offsets,decades) do arr,offset,dec
    println("Writing decade $dec")
    foreach(DiskArrays.eachchunk(arr)) do cc
        labels = arr[cc...]
        relabel_component!(labels, relabels, minlabels, offset)
        outview = select_decade(outyax,dec)
        outview[cc...] = labels
    end
end

using Makie, CairoMakie

dsready = open_dataset(outpath)

ar = dsready.labels[ti = At(Date(2010,7,15))][:,:]
heatmap(ar)