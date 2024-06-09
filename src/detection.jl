using YAXArrays, EarthDataLab, DiskArrays
"""
    rescale(inputcube::YAXArray, outputpath::String; <keyword arguments>))

Rescale time series between 0 and 1
### Arguments
* `inputcube`: Data cube with Time dimension
* `outputpath`: Path to write result to disk
### Keyword arguments
* `multiplier`: multiplies input before scaling. Typically -1 if data need to be reversed for extremes detection. Defaults to `nothing` for no multiplication.
* `backend`: Backend for writing cube to disk. Defaults to :zarr
* `overwrite`: If false, output will only be written to outputpath if no such file exists. Defaults to true.
* `chunksize`: A Dict specifying the chunksizes for the output dimensions of the cube, or `:input` (default) to copy chunksizes from input cube axes or `:max` to not chunk the inner dimensions
* `max_cache`: Maximum cache size. Defaults to 1e9.
"""
function rescale(
    inputcube::YAXArray,
    outputpath::AbstractString;
    multiplier::Union{Int64,Nothing} = nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    chunksize::Any = :input, #Dict("Time"=>1000), # eachchunk(inputcube)
    max_cache::Float64 = 1e9
)  

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
        max_cache=1e9,
        multiplier=multiplier
        )
end

# should write another method for input::Tuple{YAXArray}


"""
    rank_transform!(xout,xin)

Transforms input into values between 0 and 1
"""
function rank_transform!(xout, xin; multiplier = nothing)
        # modify the xin columns for which sign need to be changed
        if !isnothing(multiplier)
            #@show xin[1:10]
            xin .*= multiplier
            #@show xin[1:10]
        end
        N=length(xout)
        for (i,si) in enumerate(sortperm(xin))
            xout[si]=i/N
        end
        xout
end


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
    function applylowpass(xout,xin;lbord=lbord,width=width)
        # @show size(xin)
        im2 = [xin[j,i] for i in size(xin,2):-1:2,j in 1:(size(xin,1)-1)]
        # @show im2
        # @show lbord
        # @show width
        imfiltered = SphericalConvolutions.lowpass(Float64.(im2),lbord = lbord,width=width)
        # @show imfiltered
        xout[1:end-1,end:-1:2] = permutedims(imfiltered)
        xout[end,:] = xout[end-1,:]
        xout[:,1] = xout[:,end]
        # @show size(xout)
        # @show sum(ismissing(xout))
        xout
    end

    # apply low pass filter to spatial data
    indims = InDims("Lon","Lat")
    outdims = OutDims(
        "Lon","Lat",
        path=outputpath,
        backend=backend,
        overwrite=overwrite)
    
    # limit number of threads to 1 per worker. (lowpass doesn't work with multiple threads)
    ntr = Dict(w=>1 for w in workers())
    @show ntr

    mapCube(
        applylowpass,
        inputcube,
        indims=indims,
        outdims=outdims,
        max_cache=max_cache,
        nthreads=ntr;
        lbord=lbord,
        width=width)
end


# with input as a tupple of YAXArrays, or rather anything that is not a YAXArray
"""
    compute_extremes(
    input::Any, 
    tres::Float64,
    outputpath::AbstractString;
    tresne::Union{Float64,Nothing} = nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9)

Compute extreme events as Peaks-over-Threshold (POT). If values have been scaled between 0 and 1, threshold is the quantile *UNDER* which an event is considered an extreme. For example, to detect high temperature extremes, the input should be (- T) or (1 - T_scaled) or argument `multiplier` must be provided as a 1x(number of Variables). The inputs will be processed as [input_1 ... input_n] .* multiplier. Values should be {-1, 1}.

`input` can be a datacube of type YAXArray or a tuple of YAXArray Datasets. In case of a tuple, each Dataset must have the same dimensions. Dimension Time is mandatory. The maximum number of input layers is 7.

The output is a single layer cube with UInt8 values computed as bitwise OR (|) of layers encoded each on one bit.
E.G. if `input` has 3 layers, extreme events in each layer will be encoded respectively as 1, 2 and 4. A combined extreme event of all 3 variables will have a value of 7.

Optionally, an extra bit encodes non extreme events as events comprised between `tresne` and `1-tresne` for all layers simultaneously (bitwise AND (&)).
"""
function compute_extremes(
    input::Any,
    tres::Float64,
    outputpath::String;
    tresne::Union{Float64,Nothing} = nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9
)   
    indims = ntuple(_->InDims("Time"),length(input))
    #@show indims
    outdims = OutDims("Time",
        outtype = UInt8,
        chunksize = :input, 
        path = outputpath,
        overwrite=overwrite,
        backend=backend,
    )
    mapCube(getextremes!, input; indims=indims, outdims=outdims, max_cache=max_cache, tres = tres, tresne = tresne)
end
# method for YAXArray
"""
    compute_extremes(
    input::YAXArray, 
    tres::Float64,
    outputpath::AbstractString;
    tresne::Union{Float64,Nothing}=nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9
)
Method for YAXArray data cube, with dimensions Time and Variable
"""
function compute_extremes(
    c::YAXArray, 
    tres::Float64,
    outputpath::String;
    tresne::Union{Float64,Nothing} = nothing,
    backend::Symbol = :zarr,
    overwrite::Bool = true,
    max_cache::Float64 = 1e9
)   
    indims = InDims("Time","Variable",)
    #@show indims
    outdims = OutDims("Time",
        outtype = UInt8,
        chunksize = :input, 
        path = outputpath,
        overwrite = overwrite,
        backend = backend,
    )

    mapCube(getextremes!,c; indims=indims, outdims=outdims, max_cache=max_cache, tres = tres, tresne = tresne)

end

function getextremes!(xout, xin...; tres=0.01, tresne = nothing)
    # reduce tuple of vectors to matrix dim(Time) x dim(Variable)
    xin = reduce(hcat, xin)
    # then run the method for matrix
    xout = getextremes!(xout, xin; tres, tresne)
end

function getextremes!(xout, xin; tres=0.01, tresne = nothing)
    ninputs = size(xin)[2]
    ints = map(x -> UInt8(2^x), 0:ninputs)
    # first determine thresold exceedance, then allocate to bits
    bytear = broadcast(*, xin .< tres, ints[1:ninputs]') 
    # reduce along variable dimension
    jointar = reduce(|, bytear, dims=2, init = 0x00) 
    # non extremes
    if isnothing(tresne)
        xout[:] .= jointar[:]
        return xout
    else # assign value non extreme events
        nee = reduce(&, (xin .> (tresne)) .& (xin .< (1-tresne)), dims = 2, init = true) * ints[ninputs+1]
        jointar = nee .| jointar
    end
    #@show size(xin),size(xout), size(jointar[:])
    xout[:] .= jointar[:]
end

# filters
# diamond matrix/array
function get_diamond(dim1::Int)
    Nh = Int((dim1+1)/2);
    range_vec = cat(1:Nh, Nh-1:-1:1, dims = 1);
    out = (range_vec .+ transpose(range_vec)) .> Nh
end


function get_diamond(dim::Tuple)
    Nh = map(x->Int((x+1)/2), dim);
    MNh = maximum(Nh)
    range_vec = map(x->cat(1:x, x-1:-1:1, dims = 1), Nh);
    # @show range_vec
    ind = convert(Vector{Int},ones(length(dim)))
    fun = function(y) 
        x = range_vec[y];
        indc = copy(ind); indc[y] = dim[y]; 
        return reshape(x, tuple(indc...))
    end
    vecs = map(fun, 1:length(dim))
    # @show vecs
    # @show reduce(.+,Nh)/length(dim)*(length(dim)-1)# maximum(Nh) * (length(dim)-1) #
    out = reduce(.+, vecs) .> reduce(.+,Nh)/length(dim)*(length(dim)-1) #    
end

function get_diamond_indices(window)
    diamond =  get_diamond(window);
    diamondindices = findall(diamond)
end

myfilter = function(img)
    # img has size = window
    # window = size(img)
    # @show window
    # # central get_diamond_indices
    # Nh = map(x->Int((x+1)/2), window)
    # # @show Nh
    # # indices of 2D diamond
    # diamondindices = get_diamond_indices(window[1])
    # t = length(diamondindices); # 0.6 * length(diamondindices);
    # # @show t
    # central value
    v = img[Nh[1],Nh[2],Nh[3]]
    # @show v
    # apply diamond spatially on central slice
    s = sum(diamondindices) do ind
        view(img,:,:,Nh[3])[ind]
    end
    # @show s
    # time (3rd) dimension
    timewindow = function()
        d = false;
        for i in 0:(Nh[3]-1)
            ind = (Nh[3]-i):((Nh[3]*2)-i-1)
            d = d || all(view(img,Nh[1],Nh[2],ind))
            d ? break : d
        end
        return d
    end
    d = timewindow()
    # @show d
    return v && d && s >= t
end

mytimefilter = function(img)
    # # central get_diamond_indices
    # Nh = map(x->Int((x+1)/2), window)
    # Nh = Int((size(img,3) + 1) /2 )
    # # time (3rd) dimension
    timewindow = function()
        d = false;
        Nh2 = (Nh*2)-1
        for i in 0:(Nh-1)
            ind = (Nh-i):(Nh2-i)
            d = d || all(view(img,1,1,ind))
            d ? break : d
        end
        return d
    end
    return timewindow()
end

####### Towards anomalies

"""
    qdoy(c::YAXArray, outputpath::String; ref::NTuple{Int,Int} = (1971,2000), w::Int = 15, 
    q::Vector = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99], 
    overwrite = true, backend = :zarr)

Function to compute quantiles q over reference period ref for each day of the year using moving window 2*w+1 centered on day of the year.
    For leap years, the missing day of the year, the window ignores the missing 29th February (only 2*w values).
    Days from the end (beginning) of the time series are added for dates <= (>) w.

"""
#
function qdoy(
    c::YAXArray,
    outputpath::String;
    ref::Tuple{Int,Int} = (1971,2000), 
    w::Int = 15, 
    q::Vector = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
    )

    # assert that c has a time axis with years ref
    s = try
        # c[time = Date(ref[1])..Date(ref[2]+1) - Day(1), latitude = At(Lytton[2], atol=0.25), longitude = At(Lytton[1], atol=0.25)]
        # c[time = Date(ref[1])..Date(ref[2]+1) - Day(1), latitude = 45.0 .. 46.0, longitude = 15.0 .. 17.0]
        c[time = Date(ref[1])..Date(ref[2]+1) - Day(1)]
    catch
        @error "Input YAXArray c should have a time axis containing reference years $ref"
    end

    # day of year for a leap year
    doy = map(x->Dates.dayofyear(x), Dim{:Time}(Date("2000-01-01"):Day(1):Date("2000-12-31")))
    
    # map quantile computation over each grid cell
    indims = InDims(:Ti;
        # artype = YAXArray
        artype = DimArray
        )
 
    outdims = OutDims(Dim{:doy}(doy), Dim{:quantiles}(q),
        outtype = Float32,
        # chunksize = :input, 
        path = outputpath,
        overwrite = overwrite,
        backend = backend,
    )
    # output YAXArray should have time dimension as day of year
    # time axis
    ti = lookup(s, :Ti)
    mapCube(getquantiles!, s, ti, q, w; indims = indims, outdims = outdims)

end

# from YAXArrays tuto
function getquantiles!(xout, xin, ti, q, w)
    # even if artype is set to YAXArray, in is passed as a SubArray (not DimArray). 
    # Hence, add ti as argument
    c = xin
    ndays = 365
    ## filterig by month-day
    monthday = map(x->Dates.format(x, "u-d"), ti)
    # # for a leap year
    # datesid = map(x->Dates.format(x, "u-d"), Dim{:Time}(Date("2000-01-01"):Day(1):Date("2000-12-31")))
    # ## number of years
    NpY = size(monthday,1)÷ndays
    idx = Array{Union{Missing, Int64}}(missing, ndays+1, NpY)
    ## leap years
    leapyears = map(x -> isleapyear(Date(x)), Dates.year(ti[1]):Dates.year(ti[end]))
    ## get the day-month indices for data subsetting
    for i in 1:(ndays+1)
        if i != 60
            idx[i,:] = Int.(findall(x-> x == datesid[i], monthday))
        else
            # Feb-29
            idx[i,leapyears] = Int.(findall(x-> x == datesid[i], monthday))
        end
    end
    # version below only for each dayofyear without window.
    quantiles = map(eachslice(idx, dims=1)) do x
        values = map(skipmissing(x)) do y
            if y <= w
                o = c[1 : y+w]
                # add end of time series
                cat(o,c[end+(y-w) : end], dims = 1)
            elseif y > length(ti)-w
                o = c[y-w : end]
                # add beginning of time series
                cat(c[1 : y+w-length(ti)], o, dims = 1)
            else
                c[y-w : y+w]
            end
        end
        # typeof(values)
        # Vector{YAXArray{Union{Missing, Float64}, 1, Vector{Union{Missing, Float64}}, D, Dict{String, Any}} where D} (alias for Array{YAXArray{Union{Missing, Float64}, 1, Array{Union{Missing, Float64}, 1}, D, Dict{String, Any}} where D, 1})
        # how could I get just a YAXArray?
        # if I cat(values, dims=:Ti), Ti axis is empty. but maybe not a problem as long as I can access the data
        v = cat(values..., dims=1)
        Statistics.quantile(skipmissing(v), q)
    end
  
    # # from YAXArrays
    # InDims(YAXArrays.MovingWindow("Time", 15, 15); window_oob_value = missing)
    xout[:] = stack(quantiles, dims=1)[:]
    return xout
end

