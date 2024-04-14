# plot indicators annually by continent

import CSV
using DataFrames

path2v = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3"
# import preliminary back extension stats
back = "$path2v/indicators_annual_wstats_continents_backext.csv"
dfb = CSV.read(back, DataFrame);
dfb[!, "version"] .= "preliminary" ;

# import current version (v3)
currentv = "$path2v/indicators_annual_wstats_continents.csv"
dfc = CSV.read(currentv, DataFrame);
dfc[!, "version"] .=  "current" ;

# combine both data DataFrames
df = DataFrames.vcat(dfb, dfc)

# modify tp to get mm instead of M
df[df.variable .== "tp", :value] .*= 1e3;

# Do plots. Group by version

using StatsPlots

continents = Dict(
        1 => "Africa", 
        2 => "Asia", 
        3 => "Australia", 
        4 => "North America",
        5 => "Oceania", 
        6 => "South America", 
        7 => "Antarctica",  
        8 => "Europe",
        )

function subplot(df, varname, region, rf; p = plot(), kwargs...)
    plot!(; kwargs...)
    sdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => ==(varname), df))
        # @show sdf
    @df sdf |> 
            (df -> filter(:yr => <(1979), df)
            ) scatter!(
            :yr, :value, group = (:variable, :version),# smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \n $varname over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            )
    @df sdf |>
        (df -> filter(:yr => >=(1979), df)
            ) scatter!(
            :yr, :value, #group = :version,
            label = "$varname 1979-2022"
            )
    # add Theil-Sen estimator
    cdf = filter([:yr, :version] => (x, y) -> (x < 1979 && y == "current") || (x >= 1979), sdf)
    m, b = theilsen(cdf.yr,cdf.value)
    plot!(cdf.yr, map(x -> m*x+b, cdf.yr), label = "$varname Theil-Sen current")

    pdf = filter([:yr, :version] => (x, y) -> (x < 1979 && y == "preliminary") || (x >= 1979), sdf)
    m, b = theilsen(pdf.yr,pdf.value)
    plot!(pdf.yr, map(x -> m*x+b, pdf.yr), label = "$varname Theil-Sen preliminary")

end


include("mytheme.jl")
theme(:mytheme)
for (rf, rn) in continents
    println("$rn has number $rf")
    region = rn
    # t2mmax
    theme(:mytheme; seriesalpha = 0.66)
    p = subplot(df, "t2mmax", region, rf);
    # t2mmean
    theme(:mytheme)
    p = subplot(df, "t2m", region, rf; p);
    # t2mmin
    theme(:mytheme; seriesalpha = 0.33)
    p = subplot(df, "t2mmin", region, rf; p)
    # change ylabel
    p = plot!(ylabel = "2m temperature yearly \n continental average of over land",)
    # reset theme
    theme(:mytheme)
    savefig(p, "$path2v/fig/scatter_annual_$(rn)_t2m.png")

    # tp
    p = subplot(df, "tp", region, rf)
    p = plot!(ylabel = "Total precipitation yearly \n continental average of over land",)
    savefig(p, "$path2v/fig/scatter_annual_$(rn)_tp.png")
    
    # pet
    p = subplot(df, "pet", region, rf)
    savefig(p, "$path2v/fig/scatter_annual_$(rn)_pet.png")

    # pei 
    p = subplot(df, "pei_30", region, rf)
    theme(:mytheme; seriesalpha = 0.66)
    subplot(df, "pei_90", region, rf; p)
    theme(:mytheme; seriesalpha = 0.33)
    subplot(df, "pei_180", region, rf; p)
    plot!(p, ylabel = "PEI yearly continental \n average of over land")
    # reset theme
    theme(:mytheme)
    savefig(p, "$path2v/fig/scatter_annual_$(rn)_pei.png")
end

