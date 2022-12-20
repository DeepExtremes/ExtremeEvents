using Plots
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
    colours::Union{Nothing, PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient} = nothing, 
    replacement::Any = Pair(nothing,nothing)
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
    if isnothing(colours)
        colours = palette(:darkterrain, nlayers^2)
    else
        if typeof(colours) in [PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient]
            ErrorException("colours should be of type ColorGradient")
        end
    end
    # transpose and count backwards to get the map correctly
    plotdata = sdc.data[d,:,:]'[end:-1:1,:];
    # replace 0 
    if !isa(replacement, Pair{Nothing, Nothing})
        replace!(plotdata, replacement);
    end
    @show 2^nlayers
    # Plots.heatmap(sdc.data[d,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [0,2^nlayers], title = Date("$year") + Day(d))
    Plots.heatmap(plotdata, c = colours, zlims = [0,2^nlayers], title = Date("$year") + Day(d), xlabel="latitude", ylabel="longitude")
end

# lon, lat heatmap with reduction over time
axnms = function(axes)
    types = map(x -> string(typeof(x)), axes)
    indnms = map(x -> findfirst(":", x)[end] + 1, types)
    nms = map((x) -> types[x][indnms[x]:(indnms[x]+2)], 1:length(types))
end

hm = function(tmp::BitArray{3}; title = missing, axs = axes_rt, fn = +)
    axs_nms = axnms(axs)
    time_dim = findfirst(axs_nms .== "tim")
    lon_dim =  findfirst(axs_nms .== "lon")
    lat_dim =  findfirst(axs_nms .== "lat")
    rtmp = convert(Matrix{Float64},dropdims(reduce(fn, tmp, dims=time_dim),dims=time_dim));
    # @show size(rtmp)
    replace!(rtmp, 0 => NaN)
    x = axs[lon_dim][:];
    y = axs[lat_dim][end:-1:1];
    # @show x, y
    z = permutedims(rtmp, (time_dim < lat_dim ? lat_dim-time_dim : lat_dim, time_dim < lon_dim ? lon_dim - time_dim : lon_dim))[end:-1:1,:];   
    # @show size(z)
    Plots.heatmap(x, y, z, title = title, xlabel = "longitude", ylabel = "latitude")
end

hm = function(tmp::Array{Bool, 3}; title = missing, axs = axes_rt, fn = +)
    axs_nms = axnms(axs)
    time_dim = findfirst(axs_nms .== "tim")
    lon_dim =  findfirst(axs_nms .== "lon")
    lat_dim =  findfirst(axs_nms .== "lat")
    rtmp = convert(Matrix{Float64},dropdims(reduce(fn, tmp, dims=time_dim),dims=time_dim));
    @show typeof(rtmp)
    replace!(rtmp, 0 => NaN)
    @show typeof(rtmp)
    x = axs[lon_dim][:];
    y = axs[lat_dim][end:-1:1];
    @show typeof(x), typeof(y)
    z = permutedims(rtmp, (time_dim < lat_dim ? lat_dim-time_dim : lat_dim, time_dim < lon_dim ? lon_dim - time_dim : lon_dim))[end:-1:1,:];   
    @show typeof(z)
    Plots.heatmap(x, y, z, title = title, xlabel = "longitude", ylabel = "latitude")
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