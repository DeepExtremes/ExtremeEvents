
using JLD, DataFrames
import CSV
using WeightedOnlineStats
# using Plots
# using StatsPlots
using CairoMakie

if occursin("/Users", pwd())
    path2v = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/v3"
else
    path2v = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3"
end
# path2v = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3"
pot = 0.01
ne = 0.1

annualstats = load("$path2v/YearlyEventType_ranked_pot" * string(pot) * "_ne" * string(ne) * "_land.jld")
allx = vec(["$i.$yr.$k" for i in 0:16, yr in 1950:2022, k in range(1,8)]);
allres = vec([value(annualstats["$(yr).$k"]).y[i] for i in 1:17, yr in 1950:2022, k in range(1,8)]);


# not sure about this one
df = DataFrame(x = allx, Area = allres) |>
    (df -> DataFrames.transform(df, :x => ByRow(x -> split(x, ".")) => [:Type, :Year, :cont])) |>
    (df -> subset(df, :cont => x -> x .> "0"))
df.Type = parse.(Int8, df.Type);
df.Year = parse.(Int, df.Year);
continents = Dict(
        # "0" => "Null",
        "1" => "Africa", 
        "2" => "Asia", 
        "3" => "Australia", 
        "4" => "North America",
        "5" => "Oceania", 
        "6" => "South America", 
        "7" => "Antarctica",  
        "8" => "Europe",
        )
df.Continent = map(x -> continents[x], df.cont)
outname = "$path2v/YearlyEventType_ranked_pot$(pot)_ne$(ne)_land.csv"
# CSV.write(outname, df)
df = CSV.read(outname, DataFrame)

cadf = CSV.read("$path2v/land_wstats_continents.csv", DataFrame)
sum(cadf.value)
# 1.9e5 NOT e7!!! e7 comes from the time dimension: 365


allres_norm = allres./sum(allres,dims=1)
labels = map([UInt8(i) for _ in 1:1, i in 0:16]) do i
    n = ""
    ((i & 0x01) > 0) && (n = join((n,"HW"),"_"))
    ((i & 0x02) > 0) && (n = join((n,"D30"),"_"))
    ((i & 0x04) > 0) && (n = join((n,"D90"),"_"))
    ((i & 0x08) > 0) && (n = join((n,"D180"),"_"))
    lstrip(n,'_')
end
# p = plot(1950:2022, permutedims(reshape(allres, (17,:,9))[2:end-1,:,9], (2,1)),
#     labels=labels, legend = :outerright, lw=1, 
#     xlabel = "Year",
#     ylabel = "Annual fraction of land area and days affected",
#     size=(1000,400), dpi=300)

# savefig(p,"n_extremes_land.png")

# include("mytheme.jl")
# theme(:mytheme)
# data check: every year, the total (land area) x (days in the year) should be the same, except for leap years
# Hence, all bars should have the same height!!!
barplot(df.Year, df.Area)
# this is not the case AND there are 2 outliers: 1959 and 1979?
sum(value(annualstats["1959.2"]).y)
sum(value(annualstats["1960.2"]).y)
sum(value(annualstats["1961.2"]).y)
filter([:Year, :cont] => (x,y) -> x in range(1959, 1962) &&  y in range(2,3), df)

gdf = df |> (df -> DataFrames.groupby(df, :Year)) |> (gdf -> combine(gdf, :Area => sum))
# @df df |> (df -> DataFrames.groupby(df, :Year)) |> (gdf -> combine(gdf, :Area => sum)) plot(:Year, :Area_sum)
plot(gdf.Year, gdf.Area_sum)
# small variations on top of leap years. probably due to precision errors. 
# at least no outliers. maybe bar chart has some bug?!

# if we group by event Type:
gdf = df |> (df -> DataFrames.groupby(df, [:Year, :Type])) |> (gdf -> combine(gdf, :Area => sum))
barplot(gdf.Year, gdf.Area_sum, color = gdf.Type, stack = gdf.Type)
# this one looks fine, but actually, the values are wrong. Should get to around 6.9e7

# if we group by continent
gdf = df |> (df -> DataFrames.groupby(df, [:Year, :Continent, :cont])) |> (gdf -> combine(gdf, :Area => sum))
barplot(gdf.Year, gdf.Area_sum, color = gdf.cont, stack = gdf.cont,)

# group by Continent and Type
gdf = df |> (df -> DataFrames.groupby(df, [:Year, :Type, :Continent, :cont])) |> (gdf -> combine(gdf, :Area => sum))
barplot(gdf.Year, gdf.Area_sum, color = gdf.cont *100 + gdf.Type, stack = gdf.cont *100 + gdf.Type)
# looks OK, but too many levels

# so if I split by continent, type and year, it should be fine

# value check
df |> 
    # group by continent
    (df -> DataFrames.groupby(df, :cont)) |>
    # relative area over all types and years
    (gdf -> DataFrames.transform(gdf, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
    # select only heatwaves
    (df -> subset(df, :Type => X -> map(x -> !iseven(x), X))) |>
    # combine over years
    (df -> DataFrames.groupby(df, :cont)) |>
    (df -> combine(df, :Area_pc => sum))
# I have 0.997637 % Area over all years. So that is correct!
# pei_30
df |> 
    # group by continent
    (df -> DataFrames.groupby(df, :cont)) |>
    # relative area over all types and years
    (gdf -> DataFrames.transform(gdf, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
    # select only pei_30
    (df -> subset(df, :Type => X -> map(x -> ((x & 0x02) > 0), X))) |>
    # combine over years
    (df -> DataFrames.groupby(df, :cont)) |>
    (df -> combine(df, :Area_pc => sum))
# exact same results. 
# Isn't it suspicious that all continents and all extremes have the exact same percentage to the sixth decimal?
# no: it is expected, given that for each pixel, over the full time series, there is always the exact same number of extremes.

# land Area x Days by Int8
v = [
        RGBA(1,1,1,0), # 0x00 # 0
        RGBA(1,0,0,1), # 0x01 # 1 !!
        RGBA(0,0,.8,1), # 0x02 # 2 
        RGBA(1,0,.8,1), # 0x03 # 3 !!
        RGBA(0,0,.6,1), # 0x04 # 4 
        RGBA(1,0,.6,1), # 0x05 # 5 !!
        RGBA(0,0,.6,1), # 0x06 # 6 
        RGBA(1,0,.6,1), # 0x07 # 7 !!
        RGBA(0,0,.4,1), # 0x08 # 8 
        RGBA(1,0,.4,1), # 0x09 # 9 !!
        RGBA(0,0,.4,1), # 0x0a # 10
        RGBA(1,0,.4,1), # 0x0b # 11!!
        RGBA(0,0,.4,1), # 0x0c # 12
        RGBA(1,0,.4,1), # 0x0d # 13!!
        RGBA(0,0,.4,1), # 0x0e # 14
        RGBA(1,0,.4,1), # 0x0f # 15!!
        RGBA(.7,.7,.7,1), # 0x10 # 16
    ]
# ChatGPT suggestion of colourblind friendly colourscale
    # Dark Blue: #002D5A
    # Medium Blue: #4C7FB8
    # Light Blue: #A6C5E8
    # Light Orange: #FFB86F
    # Medium Purple (Dark): #65498C
    # Medium Purple (Medium): #8464A5
    # Medium Purple (Light): #A386BB
cols = [colorant"#ffffff",  
    colorant"#FFB86F", # 1 Light Orange
    colorant"#A6C5E8", # 2 Light Blue
    colorant"#A386BB", # 3 Medium Purple (Light)
    colorant"#4C7FB8", # 4 Medium Blue
    colorant"#8464A5", # 5 Medium Purple (Medium)
    colorant"#4C7FB8", # 6 Medium Blue
    colorant"#8464A5", # 7 Medium Purple (Medium)
    colorant"#002D5A", # 8 Dark Blue
    colorant"#65498C", # 9 Medium Purple (Dark)
    colorant"#002D5A", # 10 Dark Blue
    colorant"#65498C", # 11 Medium Purple (Dark)
    colorant"#002D5A", # 12 Dark Blue
    colorant"#65498C", # 13 Medium Purple (Dark)
    colorant"#002D5A", # 14 Dark Blue
    colorant"#65498C", # 15 Medium Purple (Dark)
    colorant"#BBBBBB", # 16
    ] 
# colours from https://personal.sron.nl/~pault/
cont_cols = [colorant"#77aadd", # light Blue
             colorant"#99ddff", # light cyan]
             colorant"#44bb99", # mint
             colorant"#bbcc33", # pear
             colorant"#aaaa00", # olive
             colorant"#eedd88", # light yellow
             colorant"#ee8866", # orange
             colorant"#ffaabb", # pink
            #  colorant"#dddddd", #light gray
            ]


function mygroupedbar(dfp, grouping; startyear = 1950, endyear = 2022, kwargs...) # for whatever reason I don't understand, I can't pass kwargs...
    f = Figure(size=(800,460));
    ax = Axis(f[1,1],
        xlabel = "Year", 
        ylabel = "Percentage of annual days and land area",
    );
    barplot!(ax, 
        dfp[!,:Year], 
        dfp[!,:Area_pc_sum], 
        color = dfp[!,grouping], 
        stack = dfp[!,grouping]; 
        # colormap = colours,
        kwargs...
    )
    # Legend(f[1,2])

    # , dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
         
    # xrotation = 45.0, xtickfontsize = 6, 
    ax.xticks = (startyear:5:endyear, string.(startyear:5:endyear))
    xlims!(startyear-1, endyear+1)
    # bar_position = bar_position,
    # bar_position = :stack,
    # kwargs...
    return f
end
            
df1=  df |> 
    # (df -> stack(df, Not(:ev), variable_name = :Year, value_name = :Area)) |>
    (df -> DataFrames.groupby(df, [:Year])) |> 
    # total extremes of 1 type =1% over 73 years, divided by n years (73) ~= 0.0137 %
    # but on average should be 1 % per year
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
    (df -> subset(df, :Type => x-> x.>0 .&& x.<16))
# aggregate area by continent
df2 = df1  |>
    (df -> DataFrames.groupby(df, [:Year, :Continent, :cont])) |>
    (df -> combine(df, :Area_pc => sum)) 
p = mygroupedbar(df2, :cont, colormap = cgrad(cont_cols, categorical = true), colorrange = (1,8))
save("$path2v/fig/landArea_by_Cont.png",p)
p = mygroupedbar(df2|> (df -> subset(df, :Year => x -> x .>= 1970)), :cont; startyear = 1970, colormap = cgrad(cont_cols, categorical = true), colorrange = (1,8))
ecbar = Colorbar(p[2,1],
    # p.content[1].scene.plots[1],)
    colormap = cgrad(cont_cols, categorical = true),
    limits = (0.5,8.5),
    halign = :left,
    spinewidth = 0,
    ticksvisible = true,
    vertical = false,
        )
    ecbar.ticks = (
        1:8, 
        [
            continents["1"],
            continents["2"],
            continents["3"],
            continents["4"],
            continents["5"],
            continents["6"],
            continents["7"],
            continents["8"],
        ],
    )
save("$path2v/fig/landArea_by_Cont_1970_2022.png",p)

# # relative area by continent
# df22 = df |> 
#     # (df -> stack(df, Not(:ev), variable_name = :Year, value_name = :Area)) |>
#     (df -> DataFrames.groupby(df, [:Year, :Continent, :cont])) |> 
#     # total extremes of 1 type =1% over 73 years, divided by n years (73) ~= 0.0137 %
#     # but on average should be 1 % per year
#     (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
#     (df -> subset(df, :Type => x-> x.>0 .&& x.<16)) |>
#     (df -> DataFrames.groupby(df, [:Year, :Continent, :cont])) |>
#     (df -> combine(df, :Area_pc => sum)) 
# # p = mygroupedbar(df22 |> (df -> subset(df, :Year => x -> x .>= 1970)), :Continent, cont_cols', startyear = 1970, bar_position = :dodge)
# p = mygroupedbar(df22|> (df -> subset(df, :Year => x -> x .>= 1970)), :cont; startyear = 1970, dodge = df.cont, colormap = cgrad(cont_cols, categorical = true), colorrange = (1,8))
# ecbar = Colorbar(p[2,1],
#     # p.content[1].scene.plots[1],)
#     colormap = cgrad(cont_cols, categorical = true),
#     limits = (0.5,8.5),
#     halign = :left,
#     spinewidth = 0,
#     ticksvisible = true,
#     vertical = false,
#         )
#     ecbar.ticks = (
#         1:8, 
#         [
#             continents["1"],
#             continents["2"],
#             continents["3"],
#             continents["4"],
#             continents["5"],
#             continents["6"],
#             continents["7"],
#             continents["8"],
#         ],
#     )
# save("$path2v/fig/landArea_by_ContDodge_1970_2022.png",p)



# aggregate area by Event type
df3 = df1  |>
    (df -> DataFrames.groupby(df, [:Year, :Type])) |>
    (df -> combine(df, :Area_pc => sum)) 
# p = mygroupedbar(df3, :Type, cols[2:end]')
p = mygroupedbar(df3|> (df -> subset(df, :Year => x -> x .>= 1970)), :Type; startyear = 1970, colormap = cgrad(cols[2:end-1], categorical = true), colorrange = (1,15))
ecbar = Colorbar(p[2,1], 
            colormap = cgrad(cols[[2,9,5,3,10,6,4]], categorical=true),
            limits = (-0.5,6.5),
            halign = :left,
            ticksvisible = false,
            spinewidth = 0,
            vertical = false,
        )
    ecbar.ticks = (
        [0,1,4],#0:4, 
        [
            "only hot",
            "only dry",
            "dry and hot",
        ],
    )
save("$path2v/fig/landArea_by_Int8_1970_2022.png",p)

function macrotype(x; ot::Type = Int)
    if x == 1
        # y = "Only hot"
        y = ot == Int ? 1 : "Only hot"
    elseif iseven(x)
        # y = "Only dry (any)"
        y = ot == Int ? 2 : "Only dry (any)"
    else
        # y = "Hot and dry"
        y = ot == Int ? 3 : "Hot and dry"
    end
    return y
end

# prepare data for plotting
dfp = df1 |> 
    # :ev to UInt8
    # (df -> DataFrames.transform(df, :ev => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # macro type
    (df -> DataFrames.transform(df, :Type => ByRow(x -> map(macrotype,x)) => :MacroType)) |>
    # group by year and MacroType
    (df -> DataFrames.groupby(df, [:Year, :MacroType])) |>
    # sum area
    (gdf -> combine(gdf, :Area_pc => sum))

    # plot
colours = [colorant"#65498C", colorant"#4C7FB8", colorant"#FFB86F" ]
p = mygroupedbar(dfp, :MacroType, colormap = cgrad(colours, categorical = true, rev = true))
save("$path2v/fig/landArea_by_macroType.png", p)

p = mygroupedbar(dfp|> (df -> subset(df, :Year => x -> x .>= 1970)), :MacroType; startyear = 1970, colormap = cgrad(colours, categorical = true, rev = true))
ecbar = Colorbar(p[2,1], 
            colormap = cgrad(colours, categorical=true, rev = true),
            limits = (0.5,3.5),
            halign = :left,
            ticksvisible = false,
            spinewidth = 0,
            vertical = false,
        )
    ecbar.ticks = (
        1:3, 
        [
            "only hot",
            "only dry (any)",
            "dry and hot",
        ],
    )
save("$path2v/fig/landArea_by_macroType_1970_2022.png", p)


# p = @df dfp groupedbar(:Year, :Area_pc_sum, group = :MacroType, 
#     legend = :top, lw = 0,
#     xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     size=(800,460), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#     colour = colours,
#     bar_position = :dodge,
#     xrotation = 45.0, xtickfontsize = 6,
#     xlims = (1949,2023),xticks=(1950.5:5:2022+1,string.(1950:5:2022))
#     )
# for (type,i) in zip(levels(dfp.MacroType), 1:3)
#     tops = sort(subset(dfp,:MacroType => x -> x .== type), :Area_pc_sum, rev = true)[1:20,:]
#     println(tops)
#     hline!(p,[tops[20,:Area_pc_sum]], color = colours[i], label = "top 20 years")
# end
# p
# savefig(p, "$path2v/fig/landArea_by_macroType_top20.png")

# # compound
# p = @df dfp|> 
#     (df -> filter(:MacroType => ==("Hot and dry"), df)) bar(:Year, :Area_pc_sum, 
#     legend = :top, lw = 0,
#     xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     size=(800,460), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#     colour = colours,
#     bar_position = :stack,
#     label = "Hot and dry",
#     xrotation = 45.0, xtickfontsize = 6,
#     xlims = (1950-1,2022+1),xticks=(1950.5:5:(2022+1),string.(1950:5:2022)));

# tops = sort(subset(dfp,:MacroType => x -> x .==("Hot and dry")), :Area_pc_sum, rev = true)[1:20,:]
# hline!(p,[tops[20,:Area_pc_sum]], color = colours, label = "top 20 years")
# savefig(p, "$path2v/fig/landArea_by_hotndry_top20.png")

# p = @df dfp|> 
#     (df -> filter(:MacroType => ==("Hot and dry"), df)
#     ) scatter(:Year, :Area_pc_sum, 
#     legend = :top,
#     xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     size=(800,460), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#     colour = colours,
#     # smooth = true,
#     label = "Hot and dry",
#     xrotation = 45.0, xtickfontsize = 6,);
# hline!(p,[tops[20,:Area_pc_sum]], color = colours, label = "top 20 years")
# savefig(p, "$path2v/fig/landArea_by_hotndry_scatter.png")

# since 1970
dfpp = dfp|> 
    (df -> filter(:MacroType => ==(3), df)) |>
    (df -> subset(df, :Year => x -> x .>=1970)
    ) ;
# data check
dfpp |>
    (df -> combine(df, :Area_pc_sum => mean))
dfpp |>
    (df -> filter(:Year => x -> x .> 2000, df)) |>
    (df -> combine(df, :Area_pc_sum => mean))

# Theil-Sen
include("../src/stats.jl")
m, b = theilsen(dfpp[!, :Year], dfpp[!,:Area_pc_sum])

### compare with Mann Kendall implementation by @mixstam1821
reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(dfpp.Year, dfpp.Area_pc_sum)
# slope is the same. intercept differs from 3rd decimal.
# # other implementation https://github.com/mmhs013/MannKendall.jl/
# using MannKendall
# mk_original_test(dfpp.Area_pc_sum)

# p = @df dfpp scatter(:Year, :Area_pc_sum, 
#     legend = :top, lw = 1,
#     xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     size=(800,460), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#     colour = colours,
#     # smooth = true,
#     label = "Hot and dry",
#     xrotation = 45.0, xtickfontsize = 6,);
# plot!(p, dfpp[!, :Year],  m .* dfpp[!, :Year] .+ b, colour = :grey, label = "Theil-Sen estimator: $(round(m; sigdigits = 2)) * Year + $(round(b; sigdigits = 2))")
# png(p, "$path2v/fig/landArea_by_hotndry_scatter_1970.png")

# p = @df dfpp bar(:Year, :Area_pc_sum, 
#     legend = :top, lw = 0,
#     xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     size=(800,460), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
#     colour = colours,
#     bar_position = :stack,
#     label = "Hot and dry",
#     xrotation = 45.0, xtickfontsize = 6,
#     xlims = (1970-1,2022+1),xticks=(1970.5:5:(2022+1),string.(1970:5:2022)));
f = Figure(size = (800,460));
ax = Axis(f[1,1])
p = barplot!(ax,dfpp.Year, dfpp.Area_pc_sum, color = colours[1], label = "Hot and dry")
ax.xlabel = "Year"
ax.ylabel = "Percentage of annual days and land area"
xlims!(1970-1,2022+1)
ax.xticks=(1970:5:(2022),string.(1970:5:2022))
l = lines!(ax, dfpp.Year,  (slope .* dfpp.Year .+ intercept), color = :grey, label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
Legend(f[2,1], ax, orientation = :horizontal, nbanks = 1, framevisible = false)
f
save("$path2v/fig/landArea_by_hotndry_1970.png",f)


# group by indicator
dfp1 = df |>
    # relative land area
    (df -> DataFrames.groupby(df, :Year)) |>
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
    # filter
    (df -> subset(df, :Type => x-> x.>0 .&& x.<16)) |>
    # :Type to UInt8
    (df -> DataFrames.transform(df, :Type => ByRow(x -> map(UInt8, x)) => :EventType)) |>
    # add variables
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x01) > 0)) => :heat)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x02) > 0)) => :d30)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x04) > 0)) => :d90)) |>
    (df -> DataFrames.transform(df, :EventType => ByRow(x -> ((x & 0x08) > 0)) => :d180))

dfp2 = dfp1 |> 
    (df -> subset(df, :Year => x -> x .>= 1970))

# function single_event_plot!(p::Plots.Plot, df::DataFrame, var::Symbol, color::Color)
#     df = df |> 
#         (df -> subset(df, var)) |>
#         (df -> groupby(df, :Year)) |> 
#         (gdf -> combine(gdf, :Area_pc => sum))
#     @df df plot!(p, :Year, :Area_pc_sum, 
#         label = String(var), legend=:top,
#         color = color
#         )
#     print("Sum area over time for " * String(var) * ": ")
#     println(sum(df.Area_pc_sum))
#     return p
# end
# p = plot(xlabel = "Year", ylabel = "Percentage of annual days and land area affected",
#     size=(800,450), dpi=300, left_margin = (5, :mm),
#     );
color_palette = [colorant"#FFB86F", colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]
# for (var, color) in zip([:heat, :d30, :d90, :d180], color_palette)
#     single_event_plot!(p, dfp1, var, color)
# end
# p
# savefig(p, "$path2v/fig/landArea_by_EventType.png")
# # sum area over time is correctly always the same
# # is there something wrong with the ERA5 extension (1950-1978)? Too many droughts?

# # since 1970
# p = plot(xlabel = "Year", ylabel = "Percentage of annual days and land area affected",
#     size=(800,450), dpi=300, left_margin = (5, :mm),
#     );
# for (var, color) in zip([:heat, :d30, :d90, :d180], color_palette)
#     single_event_plot!(p, dfp2, var, color)
# end
# p


# try with layout with 3 subplots 
# function single_event_bar(p::Tuple, df::DataFrame, var::Symbol, color::Color; slope=false)
#     df = df |> 
#         (df -> subset(df, var)) |>
#         (df -> groupby(df, :Year)) |> 
#         (gdf -> combine(gdf, :Area_pc => sum))
#     b = @df df bar(:Year, :Area_pc_sum, 
#         label = String(var), 
#         legend=:outerright,
#         color = color,
#         ylim = (0.0,3.1),
#         lw=0,
#         size=(1000,400), dpi=300, 
#         # left_margin = (5, :mm), bottom_margin = (5, :mm),
#         )
#     if slope
#         m, b1 = theilsen(df[!, :Year], df[!,:Area_pc_sum])
#         plot!(b, df[!, :Year],  m .* df[!, :Year] .+ b1, colour = :grey, label = "Theil-Sen estimator: $(round(m; sigdigits = 2)) * Year + $(round(b1; sigdigits = 2))")
#     end
#     return p = (p..., b)
# end
# l = @layout [a;b;c;d]
# p = ()
# for (var, color) in zip([:heat, :d30, :d90, :d180], color_palette)
#     p = single_event_bar(p, dfp1, var, color)
# end
# plot(p... , layout = l)

# Makie
function single_event_bar!(ax::Makie.Axis, df::DataFrame, var::Symbol, color; slope=false)
    df = df |> 
        (df -> subset(df, var)) |>
        (df -> groupby(df, :Year)) |> 
        (gdf -> combine(gdf, :Area_pc => sum))
    barplot!(ax, df.Year, df.Area_pc_sum,
        color = color, label = String(var))
        if slope
            reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(df.Year, df.Area_pc_sum)
            lines!(ax, df.Year,  (slope .* df.Year .+ intercept), color = color, label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
        end
    return ax
end

figtypsub = function(df)
    f = Figure(size = (800,660));
    for (var, color, i) in zip([:heat, :d30, :d90, :d180], color_palette, 1:4)
        ax = Axis(f[i,1],
        xlabel = "Year",
        xlabelvisible = i == 4 ? true : false,
        xticklabelsvisible = i == 4 ? true : false,
        xticks=((df.Year[1]):5:(df.Year[end]),string.((df.Year[1]):5:(df.Year[end])))
        )
        single_event_bar!(ax, df, var, color; slope = true)
    end
    linkaxes!(f.content[1], f.content[2], f.content[3], f.content[4],)
    Label(f[1:4,0],
        "Percentage of annual days and land area", rotation = pi/2,
    )
    plots = [
        f.content[1].scene[1], f.content[2].scene[1], f.content[3].scene[1], f.content[4].scene[1],
        f.content[1].scene[2], f.content[2].scene[2], f.content[3].scene[2], f.content[4].scene[2],
    ]
    f[5,1] = Legend(f,
            plots,
            map(x -> x.label, plots), 
            position = :lb, 
            orientation = :horizontal,
            nbanks = 4,
            framevisible = false,
            )
    xlims!(df.Year[1]-1,df.Year[end]+1)
    # f.content[4].xticks=((df.Year[1]):5:(df.Year[end]),string.((df.Year[1]):5:(df.Year[end])))
    f
end
f = figtypsub(dfp1)
save("$path2v/fig/landArea_by_EventTypeSubplot.png",f)

# since 1970 
f = figtypsub(dfp2)
save("$path2v/fig/landArea_by_EventTypeSubplot_1970.png", f)


# single continent
function single_event_bar!(ax::Makie.Axis, df::DataFrame, var::String, color; slope=false)
    # df = df |> 
    #     (df -> subset(df, var)) |>
    #     (df -> groupby(df, :Year)) |> 
    #     (gdf -> combine(gdf, :Area_pc => sum))
    barplot!(ax, df.Year, df.Area_pc_sum,
        color = color, label = var)
        if slope
            reject_null_hypothesis, p_value, Tau, slope, intercept = mann_kendall(df.Year, df.Area_pc_sum)
            lines!(ax, df.Year,  (slope .* df.Year .+ intercept), color = color, label = "Theil-Sen estimator: $(round(slope; sigdigits = 2)) * Year + ($(round(intercept; sigdigits = 2))) \n Mann-Kendall test: p-value = $(round(p_value; sigdigits = 2))")
        end
    return ax
end

figcontsub = function(df)
    f = Figure(size = (800,660));
for (var, color, i) in zip(continents, cont_cols, 1:length(continents))
        ax = Axis(f[i,1],
        xlabel = "Year",
        xlabelvisible = i == length(continents) ? true : false,
        xticklabelsvisible = i == length(continents) ? true : false,
        xticks = ((df.Year[1]):5:(df.Year[end]),string.((df.Year[1]):5:(df.Year[end])))
        )
        single_event_bar!(ax, 
            df |> (df -> subset(df, :Continent => y -> y .== var[2])),
            var[2], color; slope = true
            )
    end
    linkaxes!(Tuple([f.content[i] for i in 1: length(f.content)])...)
    plots = vcat(([f.content[i].scene[1], f.content[i].scene[2]] for i in 1: length(f.content))...)
        # Tuple([(f.content[i].scene[1], f.content[i].scene[2]) for i in 1: length(f.content)])

    Label(f[1:length(continents),0],
        "Percentage of annual days and land area", rotation = pi/2,
    )
    f[1:length(continents), 2] = Legend(f,
            plots,
            map(x -> x.label, plots), 
            position = :lb, 
            # orientation = :horizontal,
            # nbanks = length(continents),
            framevisible = false,
            )
    xlims!(df.Year[1]-1,df.Year[end]+1)
    # f.content[length(continents)].xticks=((df.Year[1]):5:(df.Year[end]),string.((df.Year[1]):5:(df.Year[end])))
    f
end
f = figcontsub(df22 |> (df -> subset(df, :Year => x -> x.>= 1970) ))
save("$path2v/fig/landArea_by_ContSubplot_1970.png",f)
# function single_cont_bar(p::Tuple, df::DataFrame, var::String, color::Color; slope=false)
#     # df = df |> 
#     #     # (df -> subset(df, var)) |>
#     #     (df -> groupby(df, :Year)) |> 
#     #     (gdf -> combine(gdf, :Area_pc => sum))
#     b = @df df bar(:Year, :Area_pc_sum, 
#         label = var, 
#         legend=:outerright,
#         color = color,
#         alpha = 0.8,
#         ylim = (0.0,1.51),
#         lw=0,
#         size=(1000,800), dpi=300, 
#         # left_margin = (5, :mm), bottom_margin = (5, :mm),
#         )
#     if slope
#         m, b1 = theilsen(df[!, :Year], df[!,:Area_pc_sum])
#         plot!(b, df[!, :Year],  m .* df[!, :Year] .+ b1, colour = color, label = "Theil-Sen estimator: $(round(m; sigdigits = 2)) * Year + $(round(b1; sigdigits = 2))")
#     end
#     return p = (p..., b)
# end
# l = @layout [a;b;c;d;e;f;g;h]
# p = ();
# for (var, color) in zip(continents, cont_cols)
#     p = single_cont_bar(p, df22 |> (df -> subset(df, [:Year, :Continent] => (x,y) -> (x .>= 1970) .& (y .== var[2]))),
#         var[2], color; slope = true)
# end

# plot(p... , layout = l)
# png("$path2v/fig/landArea_by_ContSubplot_1970.png")

# by continent and event type
sorted_cont = stack(DataFrame(continents), [:1, :2, :3, :4, :5, :6, :7, :8]) |>
    (df -> select(df, :value)) |>
    (df -> sort(df))
dfpc = df |> 
    (df -> DataFrames.groupby(df, [:Year, :Continent])) |> 
    (df -> DataFrames.transform(df, :Area => (x -> x./sum(x) .*100) => :Area_pc)) |>
    (df -> subset(df, :Type => x-> x.>0 .&& x.<16)) |>
    # macro type
    (df -> DataFrames.transform(df, :Type => ByRow(x -> map(macrotype,x)) => :MacroType)) |>
    # group by year, continent and MacroType
    (df -> groupby(df, [:Year, :Continent, :MacroType])) |>
    # sum area
    (gdf -> combine(gdf, :Area_pc => sum)) |>
    # subset only dry
    (df -> filter(:MacroType => ==(3), df)) |>
    (df -> subset(df, :Year => x -> x .>=1970)) 
# data check
d0 = dfpc |> 
    (df -> groupby(df, :Continent)) |> 
    (df -> combine(df, :Area_pc_sum => mean)) |>
    (df -> rename(df, :Area_pc_sum_mean => :Years_1970_2022))
d1 = dfpc |> 
    (df -> filter(:Year => x -> x .< 2000, df)) |>
    (df -> groupby(df, :Continent)) |> 
    (df -> combine(df, :Area_pc_sum => mean))|>
    (df -> rename(df, :Area_pc_sum_mean => :Years_1970_1999))
d2 = dfpc |> 
    (df -> filter(:Year => x -> x .>= 2000, df)) |>
    (df -> groupby(df, :Continent)) |> 
    (df -> combine(df, :Area_pc_sum => mean))|>
    (df -> rename(df, :Area_pc_sum_mean => :Years_2000_2022))
d = leftjoin(d0,leftjoin(d1,d2, on = :Continent), on = :Continent)
# # data check
# d0 = dfpc |> 
#     (df -> groupby(df, :Continent)) |> 
#     (df -> combine(df, :Area_pc_sum => mean)) |>
#     (df -> rename(df, :Area_pc_sum_mean => :Years_1970_2022))
# d1 = dfpc |> 
#     (df -> filter(:Year => x -> x .< 1996, df)) |>
#     (df -> groupby(df, :Continent)) |> 
#     (df -> combine(df, :Area_pc_sum => mean))|>
#     (df -> rename(df, :Area_pc_sum_mean => :Years_1970_1995))
# d2 = dfpc |> 
#     (df -> filter(:Year => x -> x .>= 1996, df)) |>
#     (df -> groupby(df, :Continent)) |> 
#     (df -> combine(df, :Area_pc_sum => mean))|>
#     (df -> rename(df, :Area_pc_sum_mean => :Years_1996_2022))
# d = leftjoin(d0,leftjoin(d1,d2, on = :Continent), on = :Continent)

# global
d0 = dfpp |> 
    (df -> combine(df, :Area_pc_sum => mean)) |>
    (df -> rename(df, :Area_pc_sum_mean => :Years_1970_2022))
d1 = dfpp |> 
    (df -> filter(:Year => x -> x .< 2000, df)) |>
    (df -> combine(df, :Area_pc_sum => mean))|>
    (df -> rename(df, :Area_pc_sum_mean => :Years_1970_1999))
d2 = dfpp |> 
    (df -> filter(:Year => x -> x .>= 2000, df)) |>
    (df -> combine(df, :Area_pc_sum => mean))|>
    (df -> rename(df, :Area_pc_sum_mean => :Years_2000_2022))

show(stdout, MIME("text/latex"),vcat(d, hcat(DataFrame(Continent = "Global"), d0, d1, d2)))

# l = @layout [a;b;c;d;e;f;g;h]
# p = ();
# for (var, color) in zip(sorted_cont[!, :value], cont_cols)
#     p = single_cont_bar(p, dfpc |> (df -> subset(df, [:Year, :Continent] => (x,y) -> (x .>= 1970) .& (y .== var))),
#         var, color; slope = true)
# end
# plot(p... , layout = l)
# png("$path2v/fig/landArea_hotndry_by_ContSubplot_1970.png")

f = figcontsub(dfpc |> (df -> subset(df, :Year => x -> x.>= 1970) ))
ylims!(-0.05, 1.0)
save("$path2v/fig/landArea_hotndry_by_ContSubplot_1970.png",f)



# # trend on last 50 years (1972 to 2021)
# using RollingFunctions
# function single_event_smooth_plot!(p::Plots.Plot, df::DataFrame, var::Symbol, color::Color)
#     df = df |> 
#         (df -> subset(df, var)) |>
#         (df -> subset(df, :Year => x -> x .> "1991")) |>
#         # (df -> groupby(df, :Year)) |> 
#         (gdf -> combine(gdf, :Area_pc => sum)) |>
#         (df -> DataFrames.transform(df, :Area_pc_sum => x -> runmean(x, 10)))
#     # print(df)
#     @df df plot!(p, :Year, :Area_pc_sum_function, label = String(var), legend=:top, color=color)
#     return p
# end
# p = plot(title = "Running Mean (10)", xlabel = "Year", 
#     ylabel = "Percentage of annual days and land area affected",
#     left_margin = (5, :mm), dpi = 300);
# for (var, color) in zip([:heat, :d30, :d90, :d180], color_palette)
#     single_event_smooth_plot!(p, dfp1, var, color)
# end
# p
# savefig(p, "$path2v/fig/RM10_LandAreaDays_by_EventType.png")

# map plots of extreme years

# function mapplot(data;kwargs...)
#     data = data[[721:1440;1:720],end:-1:1]
#     heatmap(data';kwargs...)
# end

# lsmmask = lsm_notime.lsm[:,:];

# hwyear = ds.layer[time=1952:1952][:,:,:];
# m = sum(==(0x01),hwyear,dims=1);
# mapplot(m[1,:,:] .* lsmmask, title = "Hot days on land in 1952", dpi = 300)
# savefig("$path2v/fig/HotDays1952.png")

# hwyear2 = ds.layer[time=2020:2020][:,:,:];
# m2 = sum(==(0x01),hwyear2,dims=1);
# mapplot(m2[1,:,:].* lsmmask, title = "Hot days on land in 2020", dpi = 300)
# savefig("$path2v/fig/HotDays2020.png")

# dyear = ds.layer[time=1951:1951][:,:,:];
# md = sum(x-> x >=(0x01) && x < (0x10),dyear,dims=1);
# mapplot(md[1,:,:] .* lsmmask, title = "Dry days in 1951", dpi = 300)
# savefig("$path2v/fig/DryDays1951.png")

# dyear1 = ds.layer[time=2020:2020][:,:,:];
# md1 = sum(x-> x >=(0x01) && x < (0x10),dyear1,dims=1);
# mapplot(md1[1,:,:] .* lsmmask, title = "Dry days in 2020", dpi = 300)
# savefig("$path2v/fig/DryDays2020.png")

# dyear2 = ds.layer[time=2021:2021][:,:,:];
# md2 = sum(x-> x >=(0x01) && x < (0x10),dyear2,dims=1);
# mapplot(md2[1,:,:] .* lsmmask, title = "Dry days in 2021", dpi = 300)
# savefig("$path2v/fig/DryDays2021.png")


# ii = CartesianIndex(1035,399)
# zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
# era = open_dataset(zg)
# era.t2mmax
# ts1952 = era.t2mmax.data[1035,399,:]
# q1 = quantile(ts1952,0.99)

# n = size(ts1952)[1]
# sub = 1:3000;
# p = plot(era.time[sub],ts1952[sub], labels=false)
# hline!(p,[q1],labels=false)

# #savefig(p,"n_extremes.png")