"""
    rescale(inputcube::YAXArray, outputpath::String; <keyword arguments>))

Rescale time series between 0 and 1
### Arguments
* `inputcube`: Data cube with Time dimension
* `outputpath`: Path to write result to disk
### Keyword arguments
* `backend`: Backend for writing cube to disk. Defaults to :zarr
* `overwrite`: If false, output will only be written to outputpath if no such file exists. Defaults to true.
* `chunksize`: Size of chunks of output cube.
* `max_cache`: Maximum cache size. Defaults to 1e9.
"""
function rescale(
    inputcube::YAXArray,
    outputpath::AbstractString,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    chunksize::Dict{String, Int64} = Dict("Time"=>1000),
    max_cache::Float64 = 1e9
)  

"""
    rank_transform!(xout,xin)

Transforms input into values between 0 and 1
"""
    function rank_transform!(xout,xin)
        N=length(xout)
        for (i,si) in enumerate(sortperm(xin))
            xout[si]=i/N
        end
        xout
    end

    # apply rank_transform to time series
    indims = InDims("Time")
    outdims = OutDims(
        "Time", 
        path=outputpath, 
        backend=:zarr, 
        overwrite=true, 
        chunksize=chunksize,
        )

    mapCube(
        rank_transform!,
        inputcube;
        indims=indims,
        outdims=outdims,
        max_cache=1e9
        )
end

# should write another method for input::Tuple{YAXArray}

"""
    smooth(inputcube::YAXArray, outputpath::AbstractString; <keyword arguments>)
apply low-pass spatial filter
## Arguments
* `inputcube`: Input data cube (of type YAXArray) with Lon and Lat dimensions, on which to apply the low-pass spatial filter.
* `outputpath`: Path to write result to disk
## Keyword Arguments
* `lbord`: parameter of low-pass filter. Defaults to 20
* `width`: parameter of low-pass filter. Defaults to 2
* `backend`: Backend for writing cube to disk. Defaults to :zarr
* `overwrite`: If false, output will only be written to outputpath if no such file exists. Defaults to true.
* `max_cache`: Maximum cache size. Defaults to 5e8.
"""
function smooth(
    inputcube::YAXArray,
    outputpath::AbstractString;
    lbord::Integer = 20,
    width::Integer = 2,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 5e8
)
    function applylowpass(xout,xin;lbord = 20, width=2)
        im2 = [xin[j,i] for i in size(xin,2):-1:1,j in 1:(size(xin,1)-1)]
        imfiltered = SphericalConvolutions.lowpass(Float64.(im2),lbord = lbord,width=width)
        xout[1:end-1,end:-1:1] = permutedims(imfiltered)
        xout[end,:] = xout[end-1,:]
        xout
    end

    # apply low pass filter to spatial data
    indims = InDims("Lon","Lat")
    outdims = OutDims(
        "Lon","Lat",
        path=outputpath,
        backend=backend,
        overwrite=overwrite)

    mapCube(
        applylowpass,
        inputcube,
        indims=indims,
        outdims=outdims,
        max_cache=5e8)
end


"""
    compute_extremes(
    inputs::Tuple{YAXArray}, 
    tres::Float64,
    outputpath::AbstractString;
    tresne::Union{Float64,Missing}=missing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9
)

Compute extreme events as Peaks-over-Threshold (POT). If values have been scaled between 0 and 1, threshold is the quantile *under* which an event is considered an extreme.
The output is a single layer cube with integer values computed as bitwise or (|) of layers encoded each on one bit.
E.G. if `inputs` has 3 layers, extreme events in each layer will be encoded respectively as 1, 2 and 4. A combined extreme event of all 3 variables will have a value of 7.
Optionally, an extra bit encodes non extreme events as events comprised between tresne and 1-tresne.

"""


function getextremes!(xout, xin; tres=0.01, tresne = nothing)
    ninputs = size(xin)[1]
    ints = map(x -> UInt8(2^x), 0:ninputs)
    # first determine thresold exceedance, then allocate to bits
    bytear = broadcast(*, xin .< tres, ints[1:ninputs]) 
    # reduce along variable dimension
    jointar = reduce(|, bytear, dims=1, init = 0x00) 
    # non extremes
    if isnothing(tresne)
        xout[:] .= jointar[:]
        return xout
    else # assign value non extreme events
        nee = reduce(&, (xin .> (tresne)) .& (xin .< (1-tresne)), dims = 1, init = true) * ints[ninputs+1]
        jointar = nee .| jointar
    end
    #@show size(xin),size(xout), size(jointar[:])
    xout[:] .= jointar[:]
end


# with input as a tupple of YAXArrays
function compute_extremes(
    c::YAXArray, 
    tres::Float64,
    outputpath::String;
    tresne = nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9
)   
    indims = InDims("Variable","Time")
    @show indims
    outdims = OutDims("Time",
        outtype=UInt8,
        chunksize = :input, 
        path = outputpath,
        overwrite=overwrite,
        backend=backend,
    )

    mapCube(getextremes!,c; indims=indims, outdims=outdims, max_cache=max_cache, tres = tres, tresne = tresne)

end