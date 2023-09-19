# plot indicators annually by continent

import CSV
using DataFrames

# import preliminary back extension stats
back = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/indicators_annual_wstats_continents_backext.csv"
dfb = CSV.read(back, DataFrame);
dfb[!, "version"] .= "preliminary" ;

# import current version (v3)
currentv = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/indicators_annual_wstats_continents.csv"
dfc = CSV.read(currentv, DataFrame);
dfc[!, "version"] .=  "current" ;

# combine both data DataFrames
df = DataFrames.vcat(dfb, dfc)

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

function subplot(df, varname, region, rf)
    sdf = df |> 
        (df -> filter(:continent => ==(rf), df)) |>
        (df -> filter(:stat => ==("mean"), df)) |>
        (df -> filter(:variable => ==(varname), df))
        # @show sdf
    p = @df sdf |> (df -> filter(:yr => <(1979), df)) scatter(
            :yr, :value, group = :version, smooth = true,
            legend = :outerright, lw = 1,
            title = region,
            xlabel = "Year", 
            ylabel = "Yearly continental average of \n $varname over land",
            size=(900,500), dpi=300, left_margin = (5, :mm), bottom_margin = (5, :mm),
            )
end

for (rf, rn) in continents
    println("$rn has number $rf")
    region = rn
    # t2mmax
    display(subplot(df, "t2mmax", region, rf))
    # t2mmean
    display(subplot(df, "t2m", region, rf))
    # t2mmin
    display(subplot(df, "t2mmin", region, rf))
    # tp
    display(subplot(df, "tp", region, rf))

    # pet
    # tp+pet
    # pei

end



df