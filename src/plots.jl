# using Plots
# import PlotUtils
using Dates
using DimensionalData
using DimensionalData.LookupArrays

"""
    simpleplot(
        dc::YAXArray,
        d::Int, 
        year::Int, 
        nlayer::Int;
        variable::Union{nothing, AbstractString} = nothing,
        colours::Union{nothing, PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient} = nothing, 
        replacement::Union{nothing, Pair{::Any, ::Any} = nothing
        )
    Plots a heatmap of cube `dc` on day `d` of Year `year`.
"""
function simpleplot(dc::YAXArray, 
    d::Int, 
    year::Int, 
    nlayers::Int; 
    variable::Union{Nothing, AbstractString} = nothing,
    replacement::Any = Pair(nothing,nothing),
    kwargs...
    )
    if isnothing(getAxis("Variable", dc))
            sdc = subsetcube(dc, time = year:year)
    else
        if isnothing(variable)
            ErrorException("Please specify `variable` to be plotted")
        else
            sdc = subsetcube(dc, time = year:year, variable = variable)
        end
    end
    # if isnothing(colours)
    #     colours = palette(:darkterrain, nlayers^2)
    # else
    #     if typeof(colours) in [PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient]
    #         ErrorException("colours should be of type ColorGradient")
    #     end
    # end
    # transpose and count backwards to get the map correctly
    plotdata = sdc.data[:,:,d]'[end:-1:1,:];
    # replace 0 
    if !isa(replacement, Pair{Nothing, Nothing})
        replace!(plotdata, replacement);
    end
    @show 2^nlayers
    # Plots.heatmap(sdc.data[d,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [0,2^nlayers], title = Date("$year") + Day(d))
    Plots.heatmap(plotdata, zlims = [0,2^nlayers], title = Date("$year") + Day(d), ylabel="latitude", xlabel="longitude", kwargs...)
end

"""
    simpleplot(
        dc::YAXArray,
        date::Date, 
        variable::Union{nothing, AbstractString} = nothing,
        replacement::Union{nothing, Pair{::Any, ::Any} = nothing,
        kwargs...
        )
    Plots a heatmap of cube `dc` on Date `date`.
"""
function simpleplot(dc::YAXArray, 
    date::Date; 
    variable::Union{Nothing, AbstractString} = nothing,
    replacement::Any = Pair(nothing,nothing),
    kwargs...
    )
    if isnothing(getAxis("Variable", dc))
            sdc = subsetcube(dc, time = date)
    else
        if isnothing(variable)
            ErrorException("Please specify `variable` to be plotted")
        else
            sdc = subsetcube(dc, time = date, variable = variable)
        end
    end
    # if isnothing(colours)
    #     colours = palette(:darkterrain, nlayers^2)
    # else
    #     if typeof(colours) in [PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient]
    #         ErrorException("colours should be of type ColorGradient")
    #     end
    # end
    # transpose and count backwards to get the map correctly
    plotdata = sdc.data[:,:]'[end:-1:1,:];
    # replace 0 
    if !isa(replacement, Pair{Nothing, Nothing})
        replace!(plotdata, replacement);
    end
    # Plots.heatmap(sdc.data[d,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [0,2^nlayers], title = Date("$year") + Day(d))
    Plots.heatmap(plotdata; title = date, ylabel="latitude", xlabel="longitude", kwargs...)
end

# lon, lat heatmap with reduction over time
function axnms(axes)
    nms = DimensionalData.dim2key(axes)
end

function modaxs(axs;lon = :longitude)
    # axs_nms = DimensionalData.Dimensions.dim2key(axs)
    lon_dim =  dimnum(axs, lon)
    axsl = map(x -> x > 180 ? x-360 : x, axs[lon_dim]) |> sort
    axsv = [i for i in axs]
    axsv[lon_dim] = Dim{lon}(axsl)
    axst = Tuple(axsv)
    return axst
end

function getshifts(axs::Tuple{Vararg{DimensionalData.Dimensions.Dimension}};lon = :longitude)
    axs_nms = DimensionalData.Dimensions.dim2key(axs)
    lon_dim =  dimnum(axs, lon) #lon_dim =  findfirst(axs_nms .== lon[1:3])
    # lon, lat time shifts
    shifts = [sum(axs[lon_dim] .< 0), 0, 0]
    dims = ["lo", "la", "ti"];
    i = map(x -> findfirst(dims .== lowercase(string(x))[1:2]), [i for i in axs_nms])
    return shifts[i]
end

function prephm(tmp,axs,fn;reduced = :Ti)
    # axs_nms = axnms(axs)
    red_dim = dimnum(axs, reduced)
    time_dim = dimnum(axs, :Ti)
    lon_dim =  dimnum(axs, :longitude)
    lat_dim =  dimnum(axs, :latitude)
    rtmp = mapslices(fn, tmp, dims=red_dim);
    # @show size(rtmp)
    # convert to Float to discard 0 in plot
    rtmp = convert(Array{Float64},rtmp);
    replace!(rtmp, 0 => NaN)
    # replace!(rtmp, missing => NaN)
    # define x, y, z
    if reduced  == :Ti
        x = lookup(axs, lon_dim);
        y = lookup(axs, lat_dim)[end:-1:1];
        # z = permutedims(dropdims(rtmp,dims=red_dim), (time_dim < lat_dim ? lat_dim-time_dim : lat_dim, time_dim < lon_dim ? lon_dim - time_dim : lon_dim))[end:-1:1,:];   
        z = dropdims(permutedims(rtmp,(lon_dim, lat_dim, time_dim)), dims=3)[:,end:-1:1];
    elseif reduced == :longitude
        # x = convert(Vector{Date}, lookup(axs, time_dim));
        x = lookup(axs, time_dim);
        y = lookup(axs, lat_dim)[end:-1:1];
        z = dropdims(permutedims(rtmp, (time_dim, lat_dim, lon_dim)),dims=3)[:,end:-1:1];   
    elseif reduced == :latitude
        # 
        x = lookup(axs, lon_dim);
        y = lookup(axs, time_dim);
        z = dropdims(permutedims(rtmp, (time_dim, lon_dim, lat_dim)),dims=3);   
    else
        error("Reduction with $reduced is not defined")
    end

    return x,y,z
end

# function hm(tmp::BitArray{3}; title = missing, axs = axes_rt, fn = sum, reduced = :Ti, kwargs...)
#     x,y,z = prephm(tmp,axs,fn;reduced)
#     sum(map(x->!isnan(x),z))
#     heatmap(x, y, z; axis = (title = title), kwargs...)
# end

# function hm(tmp::Array{Bool, 3},args...;kwargs...)
#     tmp = convert(BitArray{3}, tmp)
#     hm(tmp,args...;kwargs...)
# end

# axs::Tuple{Vararg{DimensionalData.Dimensions.Dimension}
function hm!(ax, tmp::BitArray{3}, axs::Tuple{Vararg{DimensionalData.Dimensions.Dimension}} ; fn = sum, reduced = :Ti, kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    h = heatmap!(ax, x, y, z; kwargs...)
    return h
end

function hm!(ax,tmp::Array{Bool, 3},args...;kwargs...)
    tmp = convert(BitArray{3}, tmp)
    h = hm!(ax,tmp,args...;kwargs...)
    return h
end

function hm!(ax, tmp::Array{Int64, 3}, axs; fn = sum, reduced = :Ti, kwargs...)
    x,y,z = prephm(tmp, axs, fn; reduced)
    h = heatmap!(ax, x, y, z;kwargs...)
    return h
end

function hm(tmp::Any, axs::Tuple{Vararg{DimensionalData.Dimensions.Dimension}} ; fn = sum, reduced = :Ti, kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    f = Figure()
    ax = Axis(f[1, 1])
    h = heatmap!(ax, x, y, z; kwargs...)
    return f, ax, h
end

function myfig!(fig, lon, lat, data, q, units)
    ax = GeoAxis(fig[1,1], 
        limits=(extrema(lon), extrema(lat)), 
        source="+proj=latlong +datum=WGS84", # src CRS
        dest="+proj=eqearth", # destination CRS, in which you want to plot
        coastlines = true # plot coastlines from Natural Earth, as a reference.
    )
    s = surface!(ax, lon, lat, data; 
        colorrange=(-maximum(abs.(q)), maximum(abs.(q))),
        # highclip=:black,
        # lowclip=:grey8,
        #colorscale = sc,
        colormap, nan_color=:grey80,
        shading=NoShading,
    )
    # # coastlines
    # cl=lines!(ax, 
    #     # GeoMakie.coastlines(),
    #     x1,y1,
    #     color = :black, linewidth=0.85)
    # translate!(cl, 0, 0, 1000)
    Colorbar(fig[1,2], s, label = units)
    # remove gridlines
    ax.xgridcolor[] = colorant"transparent";
    ax.ygridcolor[] = colorant"transparent";
    ax.xticklabelsvisible = false;
    ax.yticklabelsvisible = false;
    return(fig)
end

function num2col(i, lbls, cols)
    j = findfirst(x -> x == i, lbls)
    return Base.parse(Colorant, cols[j])
    # if i == 0
    #     return colorant"white"
    # elseif i == 1
    #     return colorant"green"
    # elseif i == 2
    #     return colorant"red"
    # elseif i == 3
    #     return colorant"black"
    # end
end

function plotEvent(df, labelcube, label_row)
    periodl =( df[label_row,:start_time], df[label_row,:end_time]+Day(1))
    latl = (df[label_row,:latitude_max], df[label_row,:latitude_min])
    lonl = (df[label_row,:longitude_min], df[label_row,:longitude_max])
    sublabels = Cube(subsetcube(labelcube, time=periodl, latitude=latl, longitude=lonl))
    @show sublabels
    # load to memory and flag pixels equal to label
    label = df[label_row, :label];
    sublabels1 = DimArray((sublabels.data .== label)[:,:,:], sublabels.axes)#(sublabels.data .== label)[:,:,:];
    @show size(sublabels1)
    if any(lonl.>180)
        # modify axes
        axs = modaxs(sublabels.axes)
        # modify sublabels1: shift lon
        shifts = getshifts(axs)
        sublabels1 = circshift(sublabels1, shifts)
    else
        axs = sublabels.axes
    end
    # # drop longitudes that are 0 and modify axes accordingly
    # axs_nms = axnms(axs)
    # no_lon_dim =  findall(axs_nms .!= "lon")
    # indlon = mapslices(sublabels1,dims=no_lon_dim) do x
    #     any(x .!= 0)
    # end
    # @show (size(indlon))
    # # i am stuck here

    # view over time
    p = hm(sublabels1, 
        axs = axs, 
        reduced = "tim", 
        xlab = "Longitude", 
        ylab = "Latitude", 
        c = cgrad(:inferno, categorical = true),
        title = "Event $label \n from " * string(Date(df[label_row,:start_time])) * " to " * string(Date(df[label_row,:end_time]))
        )  
end

# function to extract cube subset
function getsubcube(cube, periodt, lato, lon)
    subcube = cube[time=periodt[1]..periodt[2], latitude=lato[1]..lato[2], longitude=lon[1]..lon[2]]
    subcubedata = subcube.data[:,:,:]
    if any(lon.>180)
        # modify axes
        axs = modaxs(subcube.axes)
        # modify subcube shift lon
        shifts = getshifts(axs)
        subcubedata = circshift(subcubedata, shifts)
    else
        axs = subcube.axes
    end
    return subcubedata, axs
end

# function to plot a heatmap (default) or other function of a cube subset reduced over time with mode (default) or other function, on existsing Makie Axis
function cubeplot!(ax, cube, periodt, lato, lon; plotfn = heatmap!, reducefn = mode, kwargs...)
    # subset cube
    if typeof(lon) <: Vector
        # load 2 parts and join them
        subp1, axsp1 = getsubcube(cube, periodt, lato, lon[1]);
        # subp1, axsp1 = getsubcube(tmax, periodt, lato, lon[1]);
        subp2, axsp2 = getsubcube(cube, periodt, lato, lon[2]);
        # subp2, axsp2 = getsubcube(tmax, periodt, lato, lon[2]);
        dimlon = dimnum(axsp1, :longitude)
        subcube = cat(subp1, subp2, dims = dimlon);
        axs = (dimlon == 1 ? Dim{:longitude}(vcat(axsp1[1].val, axsp2[1][1]:0.25:axsp2[1][end])) : axsp1[1], 
            dimlon == 2 ? Dim{:longitude}(vcat(axsp1[2].val, axsp2[2][1]:0.25:axsp2[2][end])) : axsp1[2],
            dimlon == 3 ? Dim{:longitude}(vcat(axsp1[3].val, axsp2[3][1]:0.25:axsp2[3][end])) : axsp1[3]
            )
    else
        subcube, axs = getsubcube(cube, period, lato, lon)
    end
    x,y,z = prephm(subcube, axs, reducefn)
    h = plotfn(ax, x, y, z; kwargs...)
    return ax, h
end

# plot labels on existing Makie Axis
function labelplot!(ax, labels, period, lat, lon, lblt; kwargs...)
    # subset labelcube
    if typeof(lon) <: Vector
        # load 2 parts and join them
        sublabelsp1, axsp1 = getsublabels(labels, periodt, lato, lon[1], lblt)
        sublabelsp2, axsp2 = getsublabels(labels, periodt, lato, lon[2], lblt)
        sublabels1 = cat(sublabelsp1, sublabelsp2, dims = 2)
        axs = (axsp1[1], Dim{:longitude}(vcat(axsp1[2].val, axsp2[2][1]:0.25:axsp2[2][end])), axsp1[3])
    else
        sublabels1, axs = getsublabels(labels, period, lat, lon, lblt)
    end
    # cols = Makie.Categorical(:tab20)
    # fg = Figure(); ax1 = Axis(fg)
    x,y,z = prephm(sublabels1, axs, mode)
    res = Dict()
    foreach((x, y) -> push!(res, x => Float64(y)), lblt, eachindex(lblt))
    # res
    replace!(z, Tuple(res)...);
    # heatmap(x, y, z; colormap = cols[1:length(lblt)])
    h = heatmap!(ax, x, y, z; kwargs...)
    # h = hm!(ax, sublabels1, axs; fn = mode, kwargs...)
    return ax, h
end

function getsublabels(labels, period, lat, lon, lblt)
    sublabels = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
        # load to memory and set all other values to 0 (so that they will be set to NaN by prephm)
        # small events have been discarded from the plot
        sublabels1 = map(x -> x in lblt ? x : 0, (sublabels.data)[:,:,:]);
    
        if any(lon.>180)
            # modify axes
            axs = modaxs(sublabels.axes)
            # modify sublabels1: shift lon
            shifts = getshifts(axs)
            sublabels1 = circshift(sublabels1, shifts)
        else
            axs = sublabels.axes
        end
        return sublabels1, axs
    end
