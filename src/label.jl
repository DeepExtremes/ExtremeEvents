using ImageMorphology: ImageMorphology, DisjointMinSets, _maybe_build_symmetric_strel, label_components, label_components!, 
  is_symmetric, strel_split, minlabel, find_root!, union!

function ImageMorphology.label_components!(out::AbstractArray{T}, A::AbstractArray, se; bkg=zero(eltype(A)), wrapdims = ()) where {T<:Integer}
    axes(out) == axes(A) || throw_dmm(axes(out), axes(A))
    se = _maybe_build_symmetric_strel(se) # compat patch
    is_symmetric(se) || throw(ArgumentError("Non-symmetric structuring element is not supported yet"))
    upper_se, _ = strel_split(CartesianIndex, se)
    fill!(out, zero(T))
    sets = DisjointMinSets{T}()#DisjointSets{T}()#
    sizehint!(sets.parents, floor(Int, sqrt(length(A))))
    @inbounds for i in CartesianIndices(A)
        val = A[i]
        val == bkg && continue
        label = typemax(T)    # sentinel value
        for Δi in upper_se
            ii = i + Δi
            checkbounds(Bool, A, ii) || continue
            if A[ii] == val
                newlabel = out[ii]
                label = if ((label == typemax(T)) | (label == newlabel))
                    newlabel
                else
                    # @show sets, label, newlabel
                    ImageMorphology.union!(sets, label, newlabel)
                end
            end
        end
        if label == typemax(T)
            label = push!(sets)
        end
        out[i] = label
    end
    #Now wrap the labels along the wrapdims
    foreach(wrapdims) do d
        wrap_labels(sets,selectdim(out,d,1),selectdim(out,d,size(out,d)))
    end
    # Now parse sets to find the labels
    newlabel = minlabel(sets)
    @inbounds for i in eachindex(A, out)
        if A[i] != bkg
            out[i] = newlabel[find_root!(sets, out[i])]
        end
    end
    return out
end

function wrap_labels(sets,o1,o2)
    for (x1,x2) in zip(o1,o2)
        if !iszero(x1) && !iszero(x2)
            ImageMorphology.union!(sets,x1,x2)
        end
    end
end



# From here, some tests start:
I1 = (2:8,2:3,1:10);
I2 = (1:2,5:8,12:20);
I3 = (19:20,7:9,20:23);
I4 = (19:20,10,30);

function make_x()
    x = falses(20,10,30);
    x[I1...] .= true; #Non-overlapping event
    x[I2...] .= true;#Second event that should be wrapped
    x[I3...] .= true; #Third blob that should be wrapped with second
    x[I4...] .= true; # Another event that should not be wrapped
    x
end

x = make_x();


using Test
r = label_components(x);

@test all(==(1),r[I1...])
@test all(==(2),r[I2...])
@test all(==(3),r[I3...])
@test all(==(4),r[I4...])

rw = label_components(x,wrapdims = (1,));

@test unique(rw) == 0:3
@test all(==(1),rw[I1...])
@test all(==(2),rw[I2...])
@test all(==(2),rw[I3...])
@test all(==(3),rw[I4...])


