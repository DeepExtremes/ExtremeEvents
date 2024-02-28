using Plots
# import PlotUtils
using Dates

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
# function axnms(axes)
#     types = map(x -> string(typeof(x)), axes)
#     indnms = map(x -> findfirst(":", x)[end] + 1, types)
#     nms = map((x) -> types[x][indnms[x]:(indnms[x]+2)], 1:length(types))
# end

function modaxs(axs;lon = :longitude)
    axs_nms = DimensionalData.Dimensions.dim2key(axs)
    lon_dim =  dimnum(axs, lon)
    axsl = map(x -> x > 180 ? x-360 : x, axs[lon_dim]) |> sort
    axs[lon_dim] = Dim{lon}(axsl)#RangeAxis(lon, range(axsl[1] , axsl[end], length = length(axsl)))
    return axs
end

function getshifts(axs::Tuple{Vararg{DimensionalData.Dimensions.Dimension}};lon = :longitude)
    axs_nms = DimensionalData.Dimensions.dim2key(axs)
    lon_dim =  dimnum(axs, lon) #lon_dim =  findfirst(axs_nms .== lon[1:3])
    shifts = [sum(axs[lon_dim] .< 0), 1, 1]
    dims = ["lon", "lat", "ti"];
    i = map(x -> findfirst(dims .== x), axs_nms)
    return shifts[i]
end


# function prephm(tmp,axs,fn)
#     axs_nms = axnms(axs)
#     time_dim = findfirst(axs_nms .== "tim")
#     lon_dim =  findfirst(axs_nms .== "lon")
#     lat_dim =  findfirst(axs_nms .== "lat")
#     rtmp = dropdims(reduce(fn, tmp, dims=time_dim),dims=time_dim);
#     # @show size(rtmp)
#     # convert to Float to discard 0 in plot
#     rtmp = convert(Matrix{Float64},rtmp);
#     replace!(rtmp, 0 => NaN)
#     # replace!(rtmp, missing => NaN)
#     x = axs[lon_dim][:];
#     y = axs[lat_dim][end:-1:1];
#     # @show x, y
#     z = permutedims(rtmp, (time_dim < lat_dim ? lat_dim-time_dim : lat_dim, time_dim < lon_dim ? lon_dim - time_dim : lon_dim))[end:-1:1,:];   
#     return x,y,z
# end

function prephm(tmp,axs,fn;reduced="tim")
    # axs_nms = axnms(axs)
    red_dim = findfirst(axs_nms .== reduced)
    time_dim = findfirst(axs_nms .== "tim")
    lon_dim =  findfirst(axs_nms .== "lon")
    lat_dim =  findfirst(axs_nms .== "lat")
    rtmp = mapslices(fn, tmp, dims=red_dim)
    # @show size(rtmp)
    # convert to Float to discard 0 in plot
    rtmp = convert(Array{Float64},rtmp);
    replace!(rtmp, 0 => NaN)
    # replace!(rtmp, missing => NaN)
    # define x, y, z
    if reduced  == "tim"
        x = axs[lon_dim][:];
        y = axs[lat_dim][end:-1:1];
        z = permutedims(dropdims(rtmp,dims=red_dim), (time_dim < lat_dim ? lat_dim-time_dim : lat_dim, time_dim < lon_dim ? lon_dim - time_dim : lon_dim))[end:-1:1,:];   
    elseif reduced == "lon"
        # 
        x = convert(Vector{Date}, axs[time_dim][:]);
        y = axs[lat_dim][end:-1:1];
        z = dropdims(permutedims(rtmp, (lat_dim, time_dim, lon_dim)),dims=3)[end:-1:1,:];   
    else # reduced == "lat"
        # 
        x = axs[lon_dim][:];
        y = axs[time_dim];
        z = dropdims(permutedims(rtmp, (time_dim, lon_dim, lat_dim)),dims=3)[end:-1:1,:];   
    end

    return x,y,z
end

function hm(tmp::BitArray{3}; title = missing, axs = axes_rt, fn = sum, reduced = "tim", kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    @show sum(map(x->!isnan(x),z))
    Plots.heatmap(x, y, z; title = title, kwargs...)
end

function hm(tmp::Array{Bool, 3},args...;kwargs...)
    tmp = convert(BitArray{3}, tmp)
    hm(tmp,args...;kwargs...)
end

function hm!(tmp::BitArray{3}; axs = axes_rt, fn = sum, reduced = "tim", xlab = "longitude", ylab = "latitude", kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    Plots.heatmap!(x, y, z; kwargs...)
end

function hm!(tmp::Array{Bool, 3},args...;kwargs...)
    tmp = convert(BitArray{3}, tmp)
    hm!(tmp,args...;kwargs...)
end

function hm!(tmp::Array{Int64, 3}; axs = axes_rt, fn = sum, reduced = "tim", kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    Plots.heatmap!(x, y, z;kwargs...)
end

function hm(tmp::Any; axs = axes_rt, fn = sum, reduced = "tim", kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    Plots.heatmap(x, y, z;kwargs...)
end

function cf!(tmp::Array{Int64, 3}; axs = axes_rt, fn = sum, reduced = "tim", kwargs...)
    x,y,z = prephm(tmp,axs,fn;reduced)
    Plots.contourf!(x, y, z; kwargs...)
end

# import PlotUtils
function getColours(colours::Union{Nothing, PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient};n=10)
    if isnothing(colours)
        colours = palette(:darkterrain, n)
    else
        if typeof(colours) in [PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient]
            ErrorException("colours should be of type ColorGradient")
        end
    end
end
# cols = Tuple((Light1 = "#28828F",
# Dark2 = "#6E6E6E",
# Light2 = "#9E9E9E",
# Accent1 = "#C8C8C8", 
# Accent2 = "#366570",
# Accent3 = "#8C8C8C",
# Accent4 = "#57A9BA",
# Accent5 = "#FFD966",
# Accent6 = "#EAF1F3"))

# DeepAIcols = map(x -> Base.parse(Colorant, x)),
#               cols
# )

# using DataFrames, VegaLite
# function labelplot(tmp::Array{Int64, 3}; axs = axes_rt, fn = sum, reduced = "tim", xlab = "longitude", ylab = "latitude", title= "", colours = nothing)
#     x,y,z = prephm(tmp,axs,fn;reduced)
#     df = DataFrame(x=repeat(x,inner=size(z)[1]), y=repeat(y, size(z)[2]), z=z[:]);
#     df |> @vlplot(
#         :rect, 
#         x={"x:o", title=xlab}, 
#         y={"y:o", title=ylab}, 
#         color={
#             "z:n",
#             legend={title="Event label"}, 
#             # scale={domain=[], range=[]},
#             },
#         title=title,
#         )
# end

# using GeoMakie, CairoMakie
# function makielabel(tmp::Array{Int64, 3}; axs = axes_rt, fn = sum, reduced = "tim", xlab = "longitude", ylab = "latitude", title= "", colours = nothing)
#     x,y,z = prephm(tmp,axs,fn;reduced)
#     fig = Figure()
#     ax = GeoAxis(fig[1,1],
#         source = "+proj=longlat +datum=WGS84", 
#         dest = "+proj=eqearth",
#         lonlims=(x[1],x[end]),
#         latlims=(y[1],y[end]),
#         coastlines=true)
#     GeoMakie.surface!(ax, x, y, z; shading = false)
#     # GeoMakie.heatmap!(ax, x, y, z; shading = false)
#     fig
# end

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
#
# map(x -> num2col(x, ulbls, cols), ulbls)
# m = rand(ulbls, (10,10));
# colours = map(x -> num2col(x, ulbls, cols),m)
# Plots.heatmap(1:10,1:10,m,c=colours,yflip=true)
# colGRAD = cgrad(collect(cols), categorical=true)
# Plots.heatmap(1:10,1:10,m,c=colGRAD,yflip=true)
# plot(colours)


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