using Plots
using Dates

"""
    simpleplot(
        dc::YAXArray,
        d::Int, 
        year::Int, 
        nlayer::Int;
        variable::Union{nothing, AbstractString}
        )
    Plots a heatmap of cube `dc` on day `d` of Year `year`.
"""
function simpleplot(dc::YAXArray, d::Int, year::Int, nlayers::Int; variable = nothing)
    if isnothing(getAxis("Variable", dc))
            sdc = subsetcube(dc, time = year:year)
    else
        if isnothing(variable)
            ErrorException("Please specify `variable` to be plotted")
        else
            sdc = subsetcube(dc, time = year:year, variable = variable)
        end
    end
    # v = repeat([RGB(1,1,1)],2^nlayers + 1)
    # v[end] = RGB(.7,.7,.7)
    
    #     RGBA(1,0,0,1), # 0x01 # 1 !!

    # transpose and count backwards to get the map correctly
    Plots.heatmap(sdc.data[d,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [0,2^nlayers], title = Date("$year") + Day(d))
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