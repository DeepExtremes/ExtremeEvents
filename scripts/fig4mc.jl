## figure for DeepExtremesCubes paper (Ji et al. 2024)
# time series for event obs_event=19. 7 days (25/07/2019 to 01/07/2019). 
# first row: t2mmax + legend
# second row: pei_30? + legend
# third row: EventCube + legend
# fourth row: labelcube + legend

# using ExtremeEvents
using YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using HTTP
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import StatsBase
using Measures
using CairoMakie, GeoMakie
import Random

obs_event = 19
trial = "ranked_pot0.01_ne0.1"
etrial = "$(trial)_cmp_2016_2021"

if haskey(ENV, "https_proxy") && occursin( "bgc-jena", ENV["https_proxy"])
    path = "/Net/Groups/BGI/work_1/scratch/s3/deepextremes/v2/"
else
    path = "https://s3.bgc-jena.mpg.de:9000/deepextremes/v2/"
end

labelpath = path * "labelcube_$etrial.zarr"
lzg = zopen(labelpath,consolidated=true, fill_as_missing = false)
labels = open_dataset(labelpath)

zg = zopen("$(path)ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
tmax = era.t2mmax
rt = Cube("$(path)tmax_ranked.zarr")

peis = open_dataset(zopen("$(path)PEICube.zarr",consolidated=true, fill_as_missing = false))
rp = open_dataset(zopen("$(path)pei_ranks.zarr",consolidated=true, fill_as_missing = false))

eec = open_dataset(zopen("$(path)EventCube_$(trial).zarr",consolidated=true, fill_as_missing = false))

if haskey(ENV, "https_proxy") && occursin( "bgc-jena", ENV["https_proxy"])
    df0 = CSV.read("$(path)SanityCheck_$etrial.csv", DataFrame, header=1)
else
    df0_http = HTTP.get("$(path)SanityCheck_$etrial.csv")
    df0 = CSV.read(df0_http.body, DataFrame, header=1)
end 
df = subset(df0, :obs_event => x -> x .== obs_event, :volume => x -> x .>= 10.0, :duration => x -> x .>= "2 days", :area => x -> x .>= 5.0)

# helper functions
include("../src/stats.jl")
include("../src/plots.jl")

lat = (34.0, 60);
lon = (-10,25);
period = (Date("2019-06-27"), Date("2019-07-01"))
# transform lon to match cube
if lon[1] < 0 
    if lon[2] <= 0
        # shift longitudes
        lon0 = lon
        lon = broadcast(x->x+360,lon)
    else # lon[2] > 0
        # split bbox into 2
        lon1 = (lon[1]+360,360)
        lon2 = (0,lon[2])
        lon0 = lon
        lon = [lon1,lon2]
    end
else lon0 = lon
end
# number of days in obs_event
nd = Dates.value(period[2]-period[1]) + 1

# subset label cube with maximum intersecting bounding box of labelled events
lato = (minimum(df[: ,:latitude_min]), maximum(df[: ,:latitude_max]))
lono = (minimum(df[:, :longitude_min]), maximum(df[:, :longitude_max]))
periodo = (minimum(df[:, :start_time]), maximum(df[:,:end_time])) #???

n = (periodo[2] - periodo[1]).value + 1
time_lapse = Day(1)

# bbox obs
xlims = typeof(lon) <: Vector ? expand(extrema(lon0)) : expand(extrema((lon..., lono...)))
xlims = xlims[1] >= 180 ? xlims.-360 : xlims
ylims = expand(extrema((lat...,lato...)))

# selected labels
lbls = df.label
ulbls = sort(unique(lbls))
nlb = length(ulbls)
labcols = cgrad(:tab20, 20, categorical=true)[1:nlb]
# EventCube colorscale - limit to 5 colours : 10th/90th percentile; only hot; only dry; dry & hot; no extreme
etcols = [colorant"#FFFFFF",  
    colorant"#FFB86F", # 1 Light Orange
    # colorant"#A6C5E8", # 2 Light Blue
    colorant"#4C7FB8", # 2 Medium Blue
    colorant"#8464A5", # 3 Medium Purple (Medium)
    # colorant"#A386BB", # 3 Medium Purple (Light)
    colorant"#4C7FB8", # 4 Medium Blue
    colorant"#8464A5", # 5 Medium Purple (Medium)
    colorant"#4C7FB8", # 6 Medium Blue
    colorant"#8464A5", # 7 Medium Purple (Medium)
    # colorant"#002D5A", # 8 Dark Blue
    colorant"#4C7FB8", # 8 Medium Blue
    colorant"#8464A5", # 9 Medium Purple (Medium)
    # colorant"#65498C", # 9 Medium Purple (Dark)
    # colorant"#002D5A", # 10 Dark Blue
    colorant"#4C7FB8", # 10 Medium Blue
    colorant"#8464A5", # 11 Medium Purple (Medium)
    # colorant"#65498C", # 11 Medium Purple (Dark)
    # colorant"#002D5A", # 12 Dark Blue
    colorant"#4C7FB8", # 12 Medium Blue
    colorant"#8464A5", # 13 Medium Purple (Medium)
    # colorant"#65498C", # 13 Medium Purple (Dark)
    colorant"#4C7FB8", # 14 Medium Blue
    # colorant"#002D5A", # 14 Dark Blue
    colorant"#8464A5", # 15 Medium Purple (Medium)
    # colorant"#65498C", # 15 Medium Purple (Dark)
    colorant"#BBBBBB", # 16
    ] 

# Threshold colours
tthcols = cgrad([colorant"#FFB86F", colorant"#BBBBBB", colorant"#FFFFFF",], [0.02,0.5], categorical = true)
pthcols = cgrad([colorant"#A6C5E8", colorant"#BBBBBB", colorant"#FFFFFF",], [0.02,0.5], categorical = true)

# ratio = diff([xlims[1],xlims[2]]) ./ diff([ylims[1], ylims[2]])
# fig = Figure(size = (round(1000 * ratio[1]), 50+250*F[end][1]));
n = nd
function myfig(;size = (2400, 1500), kwargs...)
fig = Figure(;size = size, kwargs...);
for t in 1:n
    periodt = (period[1] + (t - 1) * time_lapse, period[1] + t * time_lapse - Day(1))
    Label(fig[0, t], string(periodt[1]); padding=(2, 2, 2, 2))
    lvls = [0.01, 0.1, 0.9]
    for i in 1:4
        # aggregate over time by mode
        # plot bounding box
        axt = GeoAxis(fig[i,t],
            dest = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs",
        );
        limits!(axt, (xlims[1],xlims[2]+15), ylims,)
        cl=lines!(axt, 
            GeoMakie.coastlines(),
            color = :grey20, linewidth=0.5)
        translate!(cl, 0, 0, 1000)

        # indicators
        # Tmax
        if i == 1
            # heatmap of t2mmax
            ax, h = cubeplot!(axt, tmax, periodt, ylims, lon; colormap = Reverse(:lajolla), colorrange = (288.15, 318.15)) # :lajolla # vik 

            # contour of rt 0.01
            ax, c = cubeplot!(axt, rt, periodt, ylims, lon; plotfn = contour!, levels = lvls, colormap = tthcols, colorrange = (0.0, 1.0))

        end 

        # PEICube
        if i == 2
            # plot pei_30
            cubeplot!(axt, peis.pei_30, periodt, ylims, lon; colormap = Reverse(:berlin), colorrange = (-5,5)) # :bilbao # :managua
            # contour rp 0.01
            cubeplot!(axt, rp.pei_30, periodt, ylims, lon; plotfn = contour!, levels = lvls, colormap = pthcols, colorrange = (0.0, 1.0))

        end 

        # EventCube
        if i == 3
            cubeplot!(axt,eec.layer, periodt, ylims, lon;
                colormap = Makie.Categorical(etcols), 
                zlims = (0,16),
                colorrange = (0,16))
        end

        # labelcube
        if i == 4
            # skip timestep if no data
            ind = df.start_time .<= periodt[2] .&& df.end_time .>= periodt[1];  
            if any(ind)
                latt = (
                    minimum(df.latitude_min[ind]),
                    maximum(df.latitude_max[ind])
                    )
                lont = (minimum(df.longitude_min[ind]), maximum(df.longitude_max[ind])+.25)
                # labels in this time step
                lblt = sort(unique(df.label[ind]))
                axt, h = labelplot!(axt, labels, periodt, latt, lont, lblt; colormap = labcols[indexin(lblt, ulbls)])
            end
        end
        # remove ticks and grid
        axt.xticklabelsvisible = false;
        axt.xgridcolor[] = colorant"transparent";
        axt.ygridcolor[] = colorant"transparent";
        axt.xticklabelsvisible = false;
        axt.yticklabelsvisible = false;  

        if t==1 
            rowsize!(fig.layout,i,Relative(1/4))
        end
        if i == 1
            colsize!(fig.layout,t,Relative(1/6))
        end
    end
end
# colorbar
# t2mmax
Label(fig[1, n+1, Left()], 
    "Tmax (°C)",
    rotation = π/ 2, padding=(2, 2, 2, 2))
fg = fig[1,n+1] = GridLayout()
cbar1 = Colorbar(fg[1,1],
        colormap = Reverse(:lajolla), 
        colorrange = (288.15, 318.15),
        ticks = ([293.15, 303.15, 313.15], ["20", "30", "40"])
    )

lt = Legend(fg[1,2],
    [LineElement(color = tthcols[1], linestyle = nothing), 
    LineElement(color = tthcols[2], linestyle = nothing), 
    LineElement(color = tthcols[3], linestyle = nothing), ],
    ["0.01", "0.1", "0.9"],
    "Tmax Rank",
    titlefont = :regular,
    patchsize = (25, 25),
    backgroundcolor = RGB(1, 0.9978, 0.79425),
    framecolor = colorant"#FFFFFF",
    )
    
# pei
Label(fig[2, n+1, Left()], 
    "PE30 (mm day⁻¹)",
    rotation = π/ 2, padding = (4,4,4,4))

pcbar1 = Colorbar(fig[2,n+1][1,1],
        # label = L"\text{PE30 (mm day}^{-1})",
        colormap = Reverse(:berlin), #:bilbao, # :managua
        colorrange = (-5, 5),
        halign = :left,
    )

lp = Legend(fig[2,n+1][1,2],
    [LineElement(color = pthcols[1], linestyle = nothing), 
    LineElement(color = pthcols[2], linestyle = nothing), 
    LineElement(color = pthcols[3], linestyle = nothing), ],
    ["0.01", "0.1", "0.9"],
    "PEI30 Rank",
    titlefont = :regular,
    patchsize = (25, 25), rowgap = 10,
    backgroundcolor = RGB(0.99987, 0.68007, 0.67995),
    framecolor = colorant"#FFFFFF",
    )

# EventCube
Label(fig[3, n+1, Left()], "Event-Cube", rotation = π/ 2, padding = (4,4,4,4))

ecbar = Colorbar(fig[3,n+1], 
            colormap = cgrad(etcols[[1,2,3,4,17]], categorical=true),
            limits = (-0.5,4.5),
            halign = :left,
        )
ecbar.ticks = (
    # 0:16,
    0:4, 
    [
        "10th/90th percentile",
        "only hot",
        "only dry",
        "dry and hot",
        "no extreme",
    ],
    # [
    #     "no event : rank > 0.1 and rank < 0.9",
    #     "only hot",
    #     "only dry (30d)",
    #     "dry (30d) and hot",
    #     "only dry (90d)", 
    #     "dry (90d) and hot", 
    #     "dry (30d and 90d)", 
    #     "dry (30d and 90d) and hot", 
    #     "only dry (180d)", 
    #     "dry (180d) and hot", 
    #     "dry (30d and 180d)", 
    #     "dry (30d and 180d) and hot", 
    #     "dry (90d and 180d)", 
    #     "dry (90d and 180d) and hot", 
    #     "dry (30d, 90d and 180d)", 
    #     "dry (30d, 90d and 180d) and hot", 
    #     "no event : rank < 0.1 and rank > 0.9"]
    )

# labels
Label(fig[4, n+1, Left()], "Label-Cube", rotation = π/ 2, padding = (4,4,4,4))

lcbar = Colorbar(fig[4,n+1], 
        # label = L"Labelled CHD events lasting $> 2$ days)",
        colormap = cgrad(labcols[1:nlb], nlb, categorical=true),
        halign = :left,
    )
if length(ulbls) > 1 
    lcbar.limits = (1,nlb)
    lcbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), string.(ulbls))
else
    lcbar.ticks = ([0.5], string.(ulbls))
end

colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)
fig
end

fontsize_theme = Theme(fontsize = 32)
fig = with_theme(fontsize_theme) do
    fig = myfig(font = "DejaVu Sans")
end

if haskey(ENV, "https_proxy") && occursin( "bgc-jena", ENV["https_proxy"])
    save("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v2/fig/plot" * "_" * trial * "_Event_$obs_event" * "_GeoAxis.png", fig, dpi = 300) 
else
    save("plot" * "_" * trial * "_Event_$obs_event" * "_GeoAxis.png", fig, dpi = 300) 
end

#### without GeoAxis
function myfig(;size = (2400, 1500), kwargs...)
fig = Figure(;size = size, kwargs...);
for t in 1:n
    periodt = (period[1] + (t - 1) * time_lapse, period[1] + t * time_lapse - Day(1))
    Label(fig[0, t], string(periodt[1]); padding = (4,4,4,4))
    lvls = [0.01, 0.1, 0.9]
    for i in 1:4
        # aggregate over time by mode
        # plot bounding box
        axt = Axis(fig[i,t],
            # topspinevisible = i == 4 ? true : false,
            # leftspinevisible = i == 4 ? true : false,
            # bottomspinevisible = i == 4 ? true : false,
            # rightspinevisible = i == 4 ? true : false,
            # spinewidth = 0.5,
            topspinecolor = colorant"#aaaaaa",
            leftspinecolor = colorant"#aaaaaa",
            bottomspinecolor = colorant"#aaaaaa",
            rightspinecolor = colorant"#aaaaaa",
            aspect = DataAspect(),
        );
        limits!(axt, xlims, ylims,)
        cl=lines!(axt, 
            GeoMakie.coastlines(),
            color = :grey20, linewidth=0.5)
        translate!(cl, 0, 0, 1000)

        # indicators
        # Tmax
        if i == 1
            # heatmap of t2mmax
            ax, h = cubeplot!(axt, tmax, periodt, ylims, lon; colormap = Reverse(:lajolla), colorrange = (288.15, 318.15)) # :lajolla # vik 

            # contour of rt 0.01
            ax, c = cubeplot!(axt, rt, periodt, ylims, lon; plotfn = contour!, levels = lvls, colormap = tthcols, colorrange = (0.0, 1.0))

        end 

        # PEICube
        if i == 2
            # plot pei_30
            cubeplot!(axt, peis.pei_30, periodt, ylims, lon; colormap = Reverse(:berlin), colorrange = (-5,5)) # :bilbao # :managua
            # contour rp 0.01
            cubeplot!(axt, rp.pei_30, periodt, ylims, lon; plotfn = contour!, levels = lvls, colormap = pthcols, colorrange = (0.0, 1.0))

        end 

        # EventCube
        if i == 3
            cubeplot!(axt,eec.layer, periodt, ylims, lon;
                colormap = Makie.Categorical(etcols), 
                zlims = (0,16),
                colorrange = (0,16))
        end

        # labelcube
        if i == 4
            # skip timestep if no data
            ind = df.start_time .<= periodt[2] .&& df.end_time .>= periodt[1];  
            if any(ind)
                latt = (
                    minimum(df.latitude_min[ind]),
                    maximum(df.latitude_max[ind])
                    )
                lont = (minimum(df.longitude_min[ind]), maximum(df.longitude_max[ind])+.25)
                # labels in this time step
                lblt = sort(unique(df.label[ind]))
                axt, h = labelplot!(axt, labels, periodt, latt, lont, lblt; colormap = labcols[indexin(lblt, ulbls)])
            end
        end
        # remove grid and axes styles
        hidedecorations!(axt)  # Hide the default grid and decorations for customization

        if t==1 
            # rowsize!(fig.layout,i,Relative(1/4))
            rowsize!(fig.layout,i,Fixed(345))
        end
        if i == 1
            # colsize!(fig.layout,t,Relative(1/6))
            colsize!(fig.layout,t,Fixed(345))
        end
    end
end
# colorbar
# t2mmax
Label(fig[1, n+1, Left()], 
    "Tmax (°C)",
    rotation = π/ 2, padding = (4,4,4,4))
fg = fig[1,n+1] = GridLayout()
cbar1 = Colorbar(fg[1,1],
        colormap = Reverse(:lajolla), 
        colorrange = (288.15, 318.15),
        ticks = ([293.15, 303.15, 313.15], ["20", "30", "40"])
    )

lt = Legend(fg[1,2],
    [LineElement(color = tthcols[1], linestyle = nothing), 
    LineElement(color = tthcols[2], linestyle = nothing), 
    LineElement(color = tthcols[3], linestyle = nothing), ],
    ["0.01", "0.1", "0.9"],
    "Tmax Rank",
    titlefont = :regular,
    patchsize = (25, 25),
    backgroundcolor = RGB(1, 0.9978, 0.79425),
    framecolor = colorant"#FFFFFF",
    )
    
# pei
Label(fig[2, n+1, Left()], 
    "PE30 (mm day⁻¹)",
    rotation = π/ 2, padding = (4,4,4,4))

pcbar1 = Colorbar(fig[2,n+1][1,1],
        colormap = Reverse(:berlin), #:bilbao, # :managua
        colorrange = (-5, 5),
        halign = :left,
    )

lp = Legend(fig[2,n+1][1,2],
    [LineElement(color = pthcols[1], linestyle = nothing), 
    LineElement(color = pthcols[2], linestyle = nothing), 
    LineElement(color = pthcols[3], linestyle = nothing), ],
    ["0.01", "0.1", "0.9"],
    "PEI30 Rank",
    titlefont = :regular,
    patchsize = (25, 25), rowgap = 10,
    backgroundcolor = RGB(0.99987, 0.68007, 0.67995),
    framecolor = colorant"#FFFFFF",
    )

# EventCube
Label(fig[3, n+1, Left()], "Event-Cube", rotation = π/ 2, padding=(2, 2, 2, 2))

ecbar = Colorbar(fig[3,n+1], 
            colormap = cgrad(etcols[[1,2,3,4,17]], categorical=true),
            limits = (-0.5,4.5),
            halign = :left,
        )
ecbar.ticks = (
    # 0:16,
    0:4, 
    [
        "10th/90th percentile",
        "only hot",
        "only dry",
        "dry and hot",
        "no extreme",
    ],
    )

# labels
Label(fig[4, n+1, Left()], "Label-Cube", rotation = π/ 2, padding=(2, 2, 2, 2))

lcbar = Colorbar(fig[4,n+1], 
        colormap = cgrad(labcols[1:nlb], nlb, categorical=true),
        halign = :left,
    )
if length(ulbls) > 1 
    lcbar.limits = (1,nlb)
    lcbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), string.(ulbls))
else
    lcbar.ticks = ([0.5], string.(ulbls))
end

# Set gap between rows and columns
colgap!(fig.layout, 4)
rowgap!(fig.layout, 4)
return fig
end

fontsize_theme = Theme(fontsize = 40)
fig = with_theme(fontsize_theme) do
    fig = myfig(size = (2500, 1500), font = "DejaVu Sans")
end

if haskey(ENV, "https_proxy") && occursin( "bgc-jena", ENV["https_proxy"])
    save("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v2/fig/plot" * "_" * trial * "_Event_$obs_event" * ".png", fig, dpi = 300) 
else
    save("plot" * "_" * trial * "_Event_$obs_event" * ".png", fig, dpi = 300) 
end

