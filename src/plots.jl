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
        colours = palette(:darkterrain, 16)
    else
        if typeof(colours) in [PlotUtils.CategoricalColorGradient, PlotUtils.ContinuousColorGradient]
            ErrorException("colours should be of type ColorGradient")
        end
    end
    # transpose and count backwards to get the map correctly
    plotdata = copy(sdc.data[d,:,:]'[end:-1:1,:]);
    # replace 0 
    if !isa(replacement, Pair{Nothing, Nothing})
        replace!(plotdata, replacement);
    end
    # Plots.heatmap(sdc.data[d,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [0,2^nlayers], title = Date("$year") + Day(d))
    Plots.heatmap(plotdata, c = colours, zlims = [0,2^nlayers], title = Date("$year") + Day(d))
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