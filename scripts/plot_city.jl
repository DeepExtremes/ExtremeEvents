##################################
## Investigate single locations ##
##################################
using YAXArrays, OnlineStats, WeightedOnlineStats, Zarr
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates #DateFormats
import CSV
# import StatsBase
# using Measures
using CairoMakie, GeoMakie

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/v3"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

trial = "ranked_pot0.01_ne0.1"
etrial = "$(trial)_cmp_S1_T3"
startyears = 1970:10:2010 
intervals = map( y -> (y, y+12), startyears)
landonly = "landonly"

mutable struct City
    name::String
    title::String
    lon::Float64
    lat::Float64
    period::Any
end

# load also ERA5 and PEICube
zg = zopen("$(path)ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
# For whatever reason, the time axis of pet is 
# ↗ Time      Sampled{Int64} 1:26663 ForwardOrdered Regular Points
# instead of 
# ↗ Ti        Sampled{DateTime} [1950-01-01T00:00:00, …, 2022-12-31T00:00:00] ForwardOrdered Irregular Points
# solved manually by editing .zmetadata
tmax = era.t2mmax
rt = Cube("$(path)tmax_ranked.zarr")

peis = open_dataset(zopen("$(path)PEICube.zarr",consolidated=true, fill_as_missing = false))
rp = open_dataset(zopen("$(path)pei_ranks.zarr",consolidated=true, fill_as_missing = false))

eec = open_dataset(zopen("$(path)EventCube_$(trial).zarr",consolidated=true, fill_as_missing = false))

import Statistics
function plot_city(city::City)
    deo = eec.layer[time = city.period, latitude = At(city.lat, atol=0.25), longitude = At(city.lon, atol=0.25)];

    qtN = Statistics.quantile(tmax[lon = At(city.lon, atol=0.25), lat = At(city.lat, atol=0.25)][:], [0.99, 0.9])
    qtN .- 273.15
    qpe30N = Statistics.quantile(skipmissing(peis.pei_30[lon = At(city.lon, atol=0.25), lat = At(city.lat, atol=0.25)][:]), [0.01, 0.1])
    qpe90N = Statistics.quantile(skipmissing(peis.pei_90[lon = At(city.lon, atol=0.25), lat = At(city.lat, atol=0.25)][:]), [0.01, 0.1])
    qpe180N = Statistics.quantile(skipmissing(peis.pei_180[lon = At(city.lon, atol=0.25), lat = At(city.lat, atol=0.25)][:]), [0.01, 0.1])

    stp = era.tp[time = city.period, latitude = At(city.lat, atol=0.25), longitude = At(city.lon, atol=0.25)]
    spet = era.pet[time = city.period, latitude = At(city.lat, atol=0.25), longitude = At(city.lon, atol=0.25)]
    stmx = tmax[time = city.period, latitude = At(city.lat, atol=0.25), longitude = At(city.lon, atol=0.25)]
    spei = peis[time = city.period, latitude = At(city.lat, atol=0.25), longitude = At(city.lon, atol=0.25)]
    
    f = Figure(size=(600,800))
    
    tempo = lookup(stp, :Ti)
    ti = 1:length(tempo)

    ax1 = Axis(f[1,1:4],
        xlabel = "Time [day]",
        ylabel = "°C",)
    Label(f[1,1, Top()], 
        text = city.title,
        halign = :left
        )
    Label(f[1,1, TopLeft()], 
        text = "(a)",
        halign = :left
        )
    
    # Tmax
    tm = scatter!(ax1,
        # twinx(), 
        ti, stmx.data[:] .- 273.15, 
        # yaxis = "Temperature [Kelvin]", 
        label = "Max. temperature at 2m",
        markersize = 4,
        # markeralpha = 0.3,
        color = :grey50,
        );
    # thresholds
    t90 = hlines!(ax1, qtN[2]-273.15, label = "90th percentile", color = :grey50, linestyle = :dot) 
    t99 = hlines!(ax1, qtN[1]-273.15, label = "99th percentile", color = :grey50, linestyle = :dash) #; xmin = ti[1], xmax = ti[end], )
    
    ax2 = Axis(f[2,1:4],
        xlabel = "Time [day]",
        ylabel ="mm/day"
    )
    Label(f[2,1, TopLeft()], 
        text = "(b)",
        halign = :left
        )
    # tp
    btp = barplot!(ax2, ti, stp.data[:].*1e3, 
        label = "Total precipitation", 
        strokewidth=0, 
        # linealpha = 1.0,
        color = Colors.JULIA_LOGO_COLORS.blue,
        );
    
    # pet
    bpet = barplot!(ax2, ti, convert(Vector{Float64},  spet.data[:]), 
        label = "Ref. evapotranspiration",
        strokewidth=0,
        # linealpha = 1.0,
        color = Colors.JULIA_LOGO_COLORS.red,
        );

    ax3 = Axis(f[3,1:4],
        xlabel = "Time [day]",
        ylabel ="mm/day"
    )
    Label(f[3,1, TopLeft()], 
        text = "(c)",
        halign = :left
        )
    # PE
    pe30 = lines!(ax3, ti, spei["pei_30"].data[:], label = "pei_30",);
    pe90 = lines!(ax3, ti, spei["pei_90"].data[:], label = "pei_90",);
    pe180 = lines!(ax3, ti, spei["pei_180"].data[:], label = "pei_180",);
    # thresholds
    pe30_90 = hlines!(ax3, qpe30N[2], label = "pei_30 10th percentile", color = 1, colormap = :tab10, colorrange = (1, 10), linestyle = :dot) #; xmin = ti[1], xmax = ti[end],  )
    pe30_99 = hlines!(ax3, qpe30N[1], label = "pei_30 1st percentile", color = 1, colormap = :tab10, colorrange = (1, 10), linestyle = :dash) #; xmin = ti[1], xmax = ti[end], )
    pe90_90 = hlines!(ax3, qpe90N[2], label = "pei_90 10th percentile", color = 2, colormap = :tab10, colorrange = (1, 10), linestyle = :dot) #; xmin = ti[1], xmax = ti[end],  )
    pe90_99 = hlines!(ax3, qpe90N[1], label = "pei_90 1st percentile", color = 2, colormap = :tab10, colorrange = (1, 10), linestyle = :dash) #; xmin = ti[1], xmax = ti[end], )
    pe180_90 = hlines!(ax3, qpe180N[2], label = "pei_180 10th percentile", color = 3, colormap = :tab10, colorrange = (1, 10), linestyle = :dot) #; xmin = ti[1], xmax = ti[end],  )
    pe180_99 = hlines!(ax3, qpe180N[1], label = "pei_180 1st percentile", color = 3, colormap = :tab10, colorrange = (1, 10), linestyle = :dash) #; xmin = ti[1], xmax = ti[end], )
    
    ax4 = Axis(f[4,1:4],
        backgroundcolor = :transparent,
        xlabel = "Time [day]",
        ylabel = "DEO", #"Event \n type"
        )
    Label(f[4,1, TopLeft()], 
        text = "(d)",
        halign = :left
        )
    cols = [colorant"#FFFFFF",  
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
    b = barplot!(ax4, ti, repeat([1], length(deo)), 
        gap = 0, 
        color = map(x -> cols[x+1], deo.data[:]),
        )
    translate!(b, 0, 0, -1000)

    # Axes
    linkxaxes!(ax1, ax2, ax3, ax4)
    hidexdecorations!(ax1, grid = false)
    hidexdecorations!(ax2, grid = false)
    hidexdecorations!(ax3, grid = false)
    hidespines!(ax4)
    hideydecorations!(ax4, label = false)
    # xpos, ticks = time_ticks(tempo)
    # ax1.xticks = ax2.xticks = ax3.xticks = ax4.xticks = (xpos, ticks)
    ax1.xticks = ax2.xticks = ax3.xticks = ax4.xticks = setticks(tempo, 5, x -> string(Date(x)); ticksvalue=false)
    # ax4.xticklabelrotation = π / 4
    # ax4.xticklabelalign = (:right, :center)
    # arrange plots on
    rowsize!(f.layout,1,Relative(1/4))
    rowsize!(f.layout,2,Relative(1/4))
    rowsize!(f.layout,3,Relative(1/4))

    plots = [tm, t90, t99, btp, bpet, pe30, pe90, pe180, pe30_90, pe30_99]
    f[5,1:2] = Legend(f,
        plots,
        map(x -> x.label, plots), 
        position = :lb, 
        orientation = :horizontal,
        nbanks = 5,
        framevisible = false,
        )
    ecbar = Colorbar(f[5,3], 
            colormap = cgrad(cols[[17,1,2,9,5,3,10,6,4]], categorical=true),
            limits = (-0.5,8.5),
            halign = :left,
            spinewidth = 0,
            ticksvisible = false,
        )
    ecbar.ticks = (
        [0,1,2,5,8,],#0:4, 
        [
            "no extreme",
            "10th/90th percentile",
            "only hot",
            "only dry",
            "dry and hot",
        ],
    )
    ecbar1 = Colorbar(f[5,4], 
            colormap = cgrad([:white, :white], categorical=true),
            halign = :left,
            spinewidth = 0,
            labelvisible = false,
            ticksvisible = false,
            ticklabelsvisible = false,
        )
    
    save("$path/fig/City_$(city.name).png", f)
    return f
end

function time_ticks(dates; frac=8)
    tempo = string.(Date.(dates))
    lentime = length(tempo)
    slice_dates = range(1, lentime, step=lentime ÷ frac)
    return slice_dates, tempo[slice_dates]
end

# setticks
function setticks(x::Any, n::Int, f::Any; ticksvalue = true)
    ticks = [1] #; ticks[n] = x[end];
    nx = length(x)
    for i in 2:n
        push!(ticks, convert(Int, round(nx/(n-1) * (i-1))))
    end
    ticklabels = map(f, x[ticks])
    if ticksvalue
        return((x[ticks], ticklabels))
    else
        return((ticks, ticklabels))
    end
end
# Event 8: heatwave in British Columbia around 29 June 2021
Lytton  = (-121.5885 +360, 50.2260284,)
period = Date("2021-06-21") .. Date("2021-08-17")
city = City("Lytton", "Lytton, BC, Canada", -121.5885 +360, 50.2260284, Date("2021-06-21") .. Date("2021-08-20"))
f = plot_city(city)
# When the strong heatwave occurs, PEIs are still ok.
    # Temperature reaches 33.5 degrees on 29th-30th June, far from reported maximum of 40+
    # Tempertaures remain above 90th percentile for a while
    # but are not so high when the PEI30 drops below the 1% threshold.
    # Temperatures rise again reaching 
# Event 3: heatwve in Pakistan
f = plot_city(City("Karachi", "Karachi, Pakistan", 67.0011, 24.8607, Date("2018-05-12") .. Date("2018-05-22")))
# Event 14: heatwave in Northern India. Temperature in Delhi reaching 47 ?
f = plot_city(City("Delhi", "Delhi, India", 77.1025, 28.7041, Date("2018-05-12") .. Date("2018-06-10")))

# Event 25: heatwave in Tunisia
f = plot_city(City("Tunis", "Tunis, Tunisia", 10.181667, 36.806389, Date("2022-07-10") .. Date("2022-07-25")))

# Event 28: drought + heat in Texas, 2011
f = plot_city(City("SanAngelo", "San Angelo, TX, USA", -100.4370 +360, 31.4638, Date("2011-05-01") .. Date("2011-08-31")))

# Event 33: heatwave in France 2003. Lyon - Latitude : 45.750000. Longitude : 4.850000
f = plot_city(City("Lyon", "Lyon, France", 4.85, 45.75, Date("2003-07-10")..Date("2003-08-31")))
f = plot_city(City("Clermont", "Clermont-Ferrand, France",3.071508, 45.7871015, Date("2003-07-10") .. Date("2003-09-09")))
f = plot_city(City("Clermont1", "Clermont-Ferrand, France",3.071508, 45.7871015, Date("2003-07-10") .. Date("2003-09-30")))
# Carcassone: Latitude : 43.216667. Longitude : 2.350000
# Event 33 - France/Europe heatwave of 2003
f = plot_city(City("Carcassonne", "Carcassone, France", 2.35, 43.22, Date("2003-07-10") .. Date("2003-08-31")))
# Beauraing
f = plot_city(City("Beauraing_22", "Beauraing, Belgium", 4.9554, 50.1102, Date("2022-06-21") .. Date("2022-09-20")))
f = plot_city(City("Beauraing_21", "Beauraing, Belgium", 4.9554, 50.1102, Date("2021-06-21") .. Date("2021-09-20")))
f = plot_city(City("Beauraing_20", "Beauraing, Belgium", 4.9554, 50.1102, Date("2020-06-21") .. Date("2020-09-20")))
f = plot_city(City("Beauraing_19", "Beauraing, Belgium", 4.9554, 50.1102, Date("2019-06-21") .. Date("2019-09-20")))
f = plot_city(City("Beauraing_18", "Beauraing, Belgium", 4.9554, 50.1102, Date("2018-06-21") .. Date("2018-09-20")))

Niameylon = 2.1254; Niameylat = 13.5116
f = plot_city(City("Niamey_81_85", "Niamey, Niger", Niameylon, Niameylat, Date(1981) .. Date(1985,12,31)))
f = plot_city(City("Niamey_83", "Niamey, Niger", Niameylon, Niameylat, Date(1983) .. Date(1983,12,31)))

Jenalat = 50.92; Jenalon = 11.59
f = plot_city(City("Jena", "Jena, Germany", Jenalon, Jenalat, Date(2018) .. Date(2022,12,31)))
f = plot_city(City("Jena_20", "Jena, Germany", Jenalon, Jenalat, Date(2020) .. Date(2020,12,31)))
