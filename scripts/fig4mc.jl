## figure for minicube paper:
# time series if event 22. 7 days (25/07/2019 to 01/07/2019). 
# first row: t2mmax + legend
# second row: pei_30? + legend
# third row: pei_90
# fourth row: pei_180
# fifth row: EventCube + legend
# sixth row: labelcube + legend

# using ExtremeEvents
using YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import StatsBase
using Measures
using CairoMakie, GeoMakie
import Random
# using PerceptualColourMaps

if haskey(ENV, "https_proxy") && occursin( "bgc-jena", ENV["https_proxy"])
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "https://s3.bgc-jena.mpg.de:9000/deepextremes/v2/"
end


etrial = "ranked_pot0.01_ne0.1"
trial = "$(etrial)_cmp_2016_2021" 
landonly = "" #"_landonly"
events = CSV.read("$(path)EventStats_$(trial)$(landonly).csv", DataFrame)

labelpath = path * "labelcube_$trial.zarr"
lzg = zopen(labelpath,consolidated=true, fill_as_missing = false)
labels = open_dataset(labelpath)

zg = zopen("$(path)ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
tmax = era.t2mmax
rt = Cube("$(path)tmax_ranked.zarr")

peis = open_dataset(zopen("$(path)PEICube.zarr",consolidated=true, fill_as_missing = false))
rp = open_dataset(zopen("$(path)pei_ranks.zarr",consolidated=true, fill_as_missing = false))

eec = open_dataset(zopen("$(path)EventCube_$(etrial).zarr",consolidated=true, fill_as_missing = false))

df0 = CSV.read("$(path)/SanityCheck_v2_$trial.csv", DataFrame, header=1)

obs_event = 22 #6 #22

df = subset(df0, :obs_event => x -> x .== obs_event, :volume => x -> x .>= 10.0, :duration => x -> x .>= "2 days", :area => x -> x .>= 5.0)
period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End]) ) #+Day(1)?
# lat = (obs[obs_event,:South], obs[obs_event,:North])
# lon = (obs[obs_event,:West], obs[obs_event,:East])
lat = (34.0, 60);
lon = (-10,25);
# period = (Date("2019-06-24"), Date("2019-07-01"))
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
# EventCube colorscale
etcols = [colorant"#FFFFFF",  
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

# Threshold colours
tthcols = cgrad([colorant"#FFB86F", colorant"#BBBBBB", colorant"#FFFFFF",], [0.02,0.5], categorical = true)
pthcols = cgrad([colorant"#A6C5E8", colorant"#BBBBBB", colorant"#FFFFFF",], [0.02,0.5], categorical = true)

# ratio = diff([xlims[1],xlims[2]]) ./ diff([ylims[1], ylims[2]])
# fig = Figure(size = (round(1000 * ratio[1]), 50+250*F[end][1]));
function myfig()
fig = Figure(size = (2400,1050));
for t in 1:n
    periodt = (periodo[1] + (t - 1) * time_lapse, periodo[1] + t * time_lapse - Day(1))
    Label(fig[1, t, Top()], string(periodt[1]); fontsize=18, padding=(2, 2, 2, 2))
    lvls = [0.01, 0.1, 0.9]
    for i in 1:4
        # aggregate over time by mode
        # plot bounding box
        axt = GeoAxis(fig[i,t],
            # title = i == 1 ? string(periodt[1]) : "", 
            # titlesize=10 ,; 
            dest = "+proj=longlat +datum=WGS84"
        );
        limits!(axt, xlims, ylims,)
        cl=lines!(axt, 
            GeoMakie.coastlines(),
            # x1,y1,
            color = :grey20, linewidth=0.5)
        translate!(cl, 0, 0, 1000)

        @show periodt

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
            # subcube1, axsp1 = getsubcube(eec.layer, periodt, lato, lon[1]);
            # subcube2, axsp2 = getsubcube(eec.layer, periodt, lato, lon[2]);
            # subcube = cat(subcube1, subcube2, dims = 2);
            # axs = (axsp1[1], Dim{:longitude}(vcat(axsp1[2].val, axsp2[2][1]:0.25:axsp2[2][end])), axsp1[3])    
            # # plot
            # h = heatmap!(axt, lookup(axs, :longitude), lookup(axs, :latitude), dropdims(subcube; dims = 1),
            #     colormap = Makie.Categorical(etcols), 
            #     zlims = (0,16),
            #     #  colorbar_ticks = (0:16, ["no event - ", "only hot", "only dry (30d)", "dry and hot", "only dry (90d)", "dry and hot", "dry", "dry and hot", "only dry (180d)", "dry and hot", "dry", "dry and hot", "dry", "dry and hot", "dry", "dry and hot", "no event - "])
            #     )
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
        # remove ticks
        axt.xticklabelsvisible = false;
        axt.xgridcolor[] = colorant"transparent";
        axt.ygridcolor[] = colorant"transparent";
        axt.xticklabelsvisible = false;
        axt.yticklabelsvisible = false;  
    end
    
end
# colorbar
# t2mmax
Label(fig[1, 1, Left()], L"\text{Tmax} (\degree \text{C})", rotation = π/ 2, padding=(2, 2, 2, 2), fontsize=18)
fg = fig[1,n+1] = GridLayout()
cbar1 = Colorbar(fg[1,1],
        # label = L"\text{Tmax} (\degree \text{C})",
        colormap = Reverse(:lajolla), 
        colorrange = (288.15, 318.15),
        ticks = ([293.15, 303.15, 313.15], ["20", "30", "40"])
    )

lt = Legend(fg[1,2],
    [LineElement(color = tthcols[1], linestyle = nothing), 
    LineElement(color = tthcols[2], linestyle = nothing), 
    LineElement(color = tthcols[3], linestyle = nothing), ],
    ["0.01", "0.1", "0.9"],
    L"\text{Tmax Rank}",
    patchsize = (25, 25),
    #  rowgap = 5,
    backgroundcolor = RGB(1, 0.9978, 0.79425),
    framecolor = colorant"#FFFFFF",
    # framevisible = false,
    )
    
# pei
Label(fig[2, 1, Left()], L"\text{PE30 (mm day}^{-1})", rotation = π/ 2, padding = (2, 2, 2, 2), fontsize = 18)

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
    L"\text{PEI30 Rank}",
    patchsize = (25, 25), rowgap = 10,
    backgroundcolor = RGB(0.99987, 0.68007, 0.67995),
    framecolor = colorant"#FFFFFF",
    # framevisible = false,
    )

# EventCube
Label(fig[3, 1, Left()], L"\text{Event-Cube}", rotation = π/ 2, padding=(2, 2, 2, 2), fontsize=18)

ecbar = Colorbar(fig[3,n+1], 
            # label = L"\text{Event type}",
            colormap = cgrad(etcols, categorical=true),
            # size = 12,
            limits = (-0.5,16.5),
            halign = :left,
            # ticklabelrotation = - π / 3,
            # vertical = false,
        )
ecbar.ticks = (0:16, [
        "no event : rank > 0.1 and rank < 0.9",
        "only hot",
        "only dry (30d)",
        "dry (30d) and hot",
        "only dry (90d)", 
        "dry (90d) and hot", 
        "dry (30d and 90d)", 
        "dry (30d and 90d) and hot", 
        "only dry (180d)", 
        "dry (180d) and hot", 
        "dry (30d and 180d)", 
        "dry (30d and 180d) and hot", 
        "dry (90d and 180d)", 
        "dry (90d and 180d) and hot", 
        "dry (30d, 90d and 180d)", 
        "dry (30d, 90d and 180d) and hot", 
        "no event : rank < 0.1 and rank > 0.9"])

# labels
Label(fig[4, 1, Left()], L"\text{Label-Cube}", rotation = π/ 2, padding=(2, 2, 2, 2), fontsize=18)

lcbar = Colorbar(fig[4,n+1], 
        label = L"Labelled CHD events lasting $> 2$ days)",
        colormap = cgrad(labcols[1:nlb], nlb, categorical=true),
        # size = 40,
        # limits = (1,nlb),
        halign = :left,
        labelrotation = 0,
        # flip_vertical_label = true
    )
if length(ulbls) > 1 
    lcbar.limits = (1,nlb)
    lcbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), string.(ulbls))
else
    lcbar.ticks = ([0.5], string.(ulbls))
end
# # cbar axis and label to the left
# lcbar.flipaxis = false

colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)
fig
end
fig = with_theme(theme_latexfonts()) do
    fig = myfig()
end
save(path2 * "v2/fig/plot" * "_" * trial * "_Event_$obs_event" * "_full.pdf", fig) 
