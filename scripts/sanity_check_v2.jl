# check that documented events get detected by analyse_events
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

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
    path2 = "$path" 
else
    # path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/"
    path = "/Net/Groups/BGI/work_1/scratch/s3/deepextremes/v2/"
    path2 = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/"
end


# load list of events
# obs = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart1_csv_Ltime.csv", DataFrame)
# remove spaces in columns names. DataConvenience.cleannames! from DataConvenience.jl
# import DataConvenience
# DataConvenience.cleannames!(obs)
# sort!(obs, :when_from, rev=true)
# # select drought and heatwave (drop floods)
# filter!(:Event_type=>!=("flood"), obs)

# for each event documented in the literature, subset eventcube according to time and spatial extent, 
# mask > 
# compute weighted "volume" of each type of event, (min, mean, max) for t2mmax, 

# run stats.jl first
include("../src/stats.jl")

# # load also ERA5 and PEICube
# peis = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")
# zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
# era = open_dataset(zg)
    
# for trial in ["ranked_pot0.01_ne0.1"]
#     eventspath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_$trial.zarr"
#     @time obs_check = sanity_check(obs, eventspath, peis, era);

#     # join data frames df and obs
#     obs.rowid = rownumber.(eachrow(obs));
#     obs_with_check = leftjoin(obs, obs_check, on=:rowid=>:label )
#     # export to csv
#     CSV.write("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart1_SanityCheck_$trial.csv", obs_with_check)

# end

# obs_with_check[:,[:Name,:Event_type, :Year, :heat, :drought30, :drought90, :drought180, :compound]]

### EventPart2
# obs = CSV.read("$(path)./EventPart2.csv", DataFrame; header=3)
# obs = CSV.read("$(path)RecentEvents.csv", DataFrame; header=3)
obs = CSV.read("$(path2)v3/EventPart2.csv", DataFrame; header=3)
include("../src/plots.jl")
# clean obs
obs = dropmissing!(obs)
obs.Start .= replace.(obs.Start, r"\." => "-");
obs.End .= replace.(obs.End, r"\." => "-");
# export to latex table
show(stdout, MIME("text/latex"),obs)

#### start helper fns
function labobs(df0::DataFrame, lon, lat, period, obs_event)
    #
    # plot labelled events flattened over time
    
    # retrive labels and count them
    tab = CubeTable(
        label    = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
    )
    labcount = fitalllabelcount(tab);
    sort!(labcount, by=i->i[end].c, rev=true);
    # toDF
    labcountdf = DataFrame(map(collectresults, labcount));
    # remove event with few voxels (arbitrarily I had 99, but maybe lower to 14) and empty lines
    labcountdf = labcountdf[map(>(0), labcountdf.count), :]
    # extract stats from events
    df = filter(:label => in(labcountdf.label), events)
    if isempty(df)
        return df0
    else
        df.start_time = Date.(df.start_time);
        df.end_time = Date.(df.end_time);
        # add column obs_event
        df.obs_event .= obs_event;
        # vcat df
        df0 = vcat(df0, df)
        return df0
    end
end

function labplot!(ax, lon, lat, period; reduced = :Ti, kwargs...)
    # sublabels = labels.layer[time=periodo[1]..periodo[2], latitude=lato[1]..lato[2], longitude=lono[1]..lono[2]]
    sublabels = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
    data = ( in(df.label).(sublabels.data))[:,:,:];
    if lon[1] >= 180
        # modify axes
        axs = modaxs(sublabels.axes)
    else
        axs = sublabels.axes
    end 
    if reduced == :longitude
        axs = (Ti((nd - Dates.value(Date(obs[obs_event,:End])-periodo[1])) : nd - Dates.value(Date(obs[obs_event,:End])-periodo[2])), axs[2], axs[3])
    end
    h = hm!(ax, data, axs; reduced, kwargs...)  
    # but this approach doesn't show if labelled events span outside the observed event bbox
    return h
end

function labplot3!(ax3, lon, lat, period; kwargs...)
    # sublabels = labels.layer[time=timlim[1]..timlim[2], latitude=latlim[1]..latlim[2], longitude=lon1[1]..lon1[2]]
    if typeof(lon) <: Vector
        sublabels1 = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1][1]..lon[1][2]]
        data1 = ( in(df.label).(sublabels1.data))[:,:,:];
        sublabels2 = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[2][1]..lon[2][2]]
        data2 = ( in(df.label).(sublabels2.data))[:,:,:];
        # concatenate over :longitude
        data = cat(data1, data2; dims = 1);

        x1 = lookup(sublabels1, :longitude) .- 360
        x2 = lookup(sublabels2, :longitude)
        x = cat(x1, x2, dims = 1)
        y = lookup(sublabels1, :latitude)
        tempo = lookup(sublabels1, :Ti)
    else
        sublabels = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
        data = ( in(df.label).(sublabels.data))[:,:,:];

        x = lookup(sublabels, :longitude)
        y = lookup(sublabels, :latitude)
        tempo = lookup(sublabels, :Ti)
    end

    δlon = abs(x[2] - x[1])
    δlat = abs(y[2] - y[1])
    # fix time
    z = 1 : length(tempo)

    points3d = [Point3f(ix, iz, iy) for ix in x, iz in z, iy in y];
    data_vec = [data[ix, iy, iz] for (ix, i) in enumerate(x), (iz, k) in enumerate(tempo), (iy, j) in enumerate(y)];

    nt = length(tempo)
    y_ticks_pos = nt > 5 ? (1:(nt÷5):nt) : (1:nt);
    ticks_time = string.(Date.(tempo[y_ticks_pos]));

    m = meshscatter!(ax3, points3d[:]; 
        color = data_vec[:], 
        colorrange = (0,1), 
        marker=Rect3f(Vec3f(-0.5), Vec3f(1)),
        markersize=Vec3f(δlon - 0.05*δlon, 0.95, δlat - 0.05δlat),
        transparency=true,
        # colormap,)
        kwargs...
    );
    ax3.yticks = (y_ticks_pos, ticks_time);
    ax3.xlabeloffset = 50;
    return ax3, m
end    

function getx(lon)
    x = [lon[1], lon[1], lon[2], lon[2], lon[1]];
    #### need to fix for Europe etc. e.g obs_event = 5
    x = lon[1] >= 180 ? x.-360 : x
end

function getsublabels(labels, period, lat, lon, lblt)
    sublabels = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
        # load to memory and set all other values to 0 (so that they will be set to NaN by prephm)
        # small events have been discarded from the plot
        sublabels1 = map(x -> x in lblt ? x : 0, (sublabels.data)[:,:,:]);
    
        if any(lon.>180)
            # modify axes
            axs = modaxs(sublabels.axes)
            # modify sublabels1: shift lon
            shifts = getshifts(axs)
            sublabels1 = circshift(sublabels1, shifts)
        else
            axs = sublabels.axes
        end
        return sublabels1, axs
    end

function labelplot!(ax, labels, period, lat, lon, lblt; kwargs...)
    # subset labelcube
    if typeof(lon) <: Vector
        # load 2 parts and join them
        sublabelsp1, axsp1 = getsublabels(labels, periodt, lato, lon[1], lblt)
        sublabelsp2, axsp2 = getsublabels(labels, periodt, lato, lon[2], lblt)
        sublabels1 = cat(sublabelsp1, sublabelsp2, dims = 2)
        axs = (axsp1[1], Dim{:longitude}(vcat(axsp1[2].val, axsp2[2][1]:0.25:axsp2[2][end])), axsp1[3])
    else
        sublabels1, axs = getsublabels(labels, period, lat, lon, lblt)
    end
    # cols = Makie.Categorical(:tab20)
    # fg = Figure(); ax1 = Axis(fg)
    x,y,z = prephm(sublabels1, axs, mode)
    res = Dict()
    foreach((x, y) -> push!(res, x => Float64(y)), lblt, eachindex(lblt))
    # res
    replace!(z, Tuple(res)...);
    # heatmap(x, y, z; colormap = cols[1:length(lblt)])
    h = heatmap!(ax, x, y, z; kwargs...)
    # h = hm!(ax, sublabels1, axs; fn = mode, kwargs...)
    return ax, h
end
function expand(x::Tuple{Int, Int})
    convert(Tuple{Float64, Float64}, x)
end
function expand(x::Tuple{Float64, Float64})
    x1 = round(x[1] - 1, RoundDown; digits = -1, base = 5)
    x2 = round(x[2] + 1, RoundUp; digits = -1, base = 5)
    return (x1, x2)
end
##### end helper fns

# for trial in ("ranked_pot0.01_ne0.1_cmp_2016_2021")#,"ranked_pot0.01_ne0.1_tcmp_2016_2021"# "ranked_pot0.01_ne0.1_tcmp_2016_2021","ranked_pot0.01_ne0.1_cmp_Sdiam3_T5_new_2016_2021")
trial = "ranked_pot0.01_ne0.1_cmp_2016_2021" # "ranked_pot0.01_ne0.1_cmp_2016_2021_land"#"ranked_pot0.01_ne0.1_tcmp_Sdiam3_T5_2016_2021" # "ranked_pot0.01_ne0.1_cmp_2016_2021" # 
# trial = "ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022"
landonly = "" #"_landonly"
# events_all = CSV.read(path * "EventStats_ranked_pot0.01_ne0.1_cmp_2016_2021.csv", DataFrame)
events = CSV.read("$(path)EventStats_$(trial)$(landonly).csv", DataFrame)
# look for intersection between spatial and temporal range of events from the table or directly in the labelcube
labelpath = path * "labelcube_$trial.zarr"
# labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_$trial.zarr"
labels = open_dataset(labelpath)# labels = Cube(labelpath)
# labels_all = open_dataset(path * "labelcube_ranked_pot0.01_ne0.1_cmp_2016_2021.zarr")

global df0 = DataFrame()
 # df0 = CSV.read("$(path)SanityCheck_$trial.csv", DataFrame, header=1)

# loop over observed events
for obs_event in 22#1: nrow(obs)
    if Date(obs[obs_event,:Start]) > labels.time[end]
        continue
    end
    # obs_event=10
    print(obs[obs_event,:])
    period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End])+Day(1))
    lat = (obs[obs_event,:South], obs[obs_event,:North])
    lon = (obs[obs_event,:West], obs[obs_event,:East])
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
    
    # get labels intersecting obs_event
    if typeof(lon) <: Vector{}
        # do everything twice...
        df1 = labobs(df0, lon[1], lat, period, obs_event);
        df1 = labobs(df1, lon[2], lat, period, obs_event);
    else
        # do only once
        df1 = labobs(df0, lon, lat, period, obs_event);
    end
    global df0 = df1;
end
unique!(df0)
# export to csv
CSV.write(path2*"SanityCheck_v2_$trial.csv", df0)
# df0 = CSV.read("$(path2)SanityCheck_$trial.csv", DataFrame, header=1)

for obs_event in 22 #1:nrow(obs)
    println(obs_event)
    # subset df0 with obs_event and minimum volume
    df = subset(df0, :obs_event => x -> x .== obs_event, :volume => x -> x .>= 10.0, :duration => x -> x .>= "2 days", :area => x -> x .>= 5.0);
    period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End]) ) #+Day(1)?
    lat = (obs[obs_event,:South], obs[obs_event,:North])
    lon = (obs[obs_event,:West], obs[obs_event,:East])
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

    # === lon-lat view ===
    figT = Figure();
    axT = Axis(figT[1,1], title = obs[obs_event,:Event] * " in " * obs[obs_event, :Area] * "\n from " * obs[obs_event, :Start] * " to " * obs[obs_event, :End],
            xlabel = "Longitude",
            ylabel = "Latitude");

    # fix the longitudes
    # event spatial bounding box
    if typeof(lon) <: Vector
        x = getx(lon0)
    else
        x = getx(lon)
    end
    y = [lat[1], lat[2], lat[2], lat[1], lat[1]];
    
    # plot bounding box
    bbT = lines!(axT, x, y , label = "");

    # === Time-lat view ===
    figL = Figure();
    axL = Axis(figL[1,1], title = obs[obs_event,:Event] * " in " * obs[obs_event, :Area] * "\n from " * obs[obs_event, :Start] * " to " * obs[obs_event, :End],
            xlabel = "Time",
            ylabel = "Latitude");
    # event time-lat bounding box
    t = [period[1], period[1], period[2], period[2], period[1]]
    nd = Dates.value(period[2]-period[1]) + 1
    # Makie doesn't support DateTime
    tplot = [1, 1, nd, nd, 1];
    # plot bounding box
    bbL = lines!(axL, tplot, y , label = "");# 
    ntticks = 4
    x_ticks_pos = nd > ntticks ? (1:(nd÷ntticks):nd) : (1:nd);
    ticks_time = string.(Date(period[1]): (nd > ntticks ? (Day(nd÷ntticks)) : Day(1)) : Date(period[2]));
    axL.xticks = (x_ticks_pos, ticks_time);
    axL.xticklabelrotation = π/4

    # # ==== lat-lon-time 3D view
    # # colormap
    # colormap = [(:white, 0.0), (:orange, 0.5)] 

    # fig3 = Figure();
    # ax3 = Axis3(fig3[1, 1], 
    #     perspectiveness = 0.5,
    #     azimuth = 6.64,
    #     elevation = 0.57, aspect = (1, 1, 1), #(1,2,1)
    #     xlabel = "Longitude", ylabel = "Time", zlabel = "Latitude");
    # # 
    # # === get labelled events maximum bounding box ===
    
    if !isempty(df)
        # subset label cube with maximum intersecting bounding box of labelled events
        lato = (minimum(df[: ,:latitude_min]), maximum(df[: ,:latitude_max]))
        lono = (minimum(df[:, :longitude_min]), maximum(df[:, :longitude_max]))
        periodo = (minimum(df[:, :start_time]), maximum(df[:,:end_time]) +Day(2))
        nl = step(labels.axes[:longitude])
        nL = (diff(collect(lon0))/nl)[1];

        latlim = (minimum([lat[1], lato[1]])-nl, maximum([lat[2], lato[2]])+nl)
        lonlim = typeof(lon) <: Vector ? (lon0[1]-nl, lon0[2]+nl) : (minimum([lon[1], lono[1]])-nl, maximum([lon[2], lono[2]])+nl)
        timlim = (minimum([period[1], periodo[1]])-Day(1), maximum([period[2], periodo[2]])+Day(1))
        timlimplot = ((nd - Dates.value(period[2]-timlim[1])), nd - Dates.value(period[2]-timlim[2]))

    #     # 3d plot limits and bbox
    #     ax3.limits = (lonlim..., timlimplot..., latlim...)
    #     hyperrect = Rect3f(Vec3f(lon0[1], 1, lat[1]), Vec3f(lon0[2]-lon0[1], nd-1, lat[2]-lat[1]))
    #     # ideally one should plot the back of the hyperrectangle first and then the front after plotting the events.
    #     bb3 = wireframe!(ax3, hyperrect, alpha = 0.5)
        
        # plot labelled events flattened over :Ti and :longitude
        # handle neg longitude
        if typeof(lon) <: Vector
            # do twice
            # > 180 ( or <0)
            # not possible to plot over longitude outside bbox

            hT = labplot!(axT, lon[2], lato, periodo; reduced = :Ti, colormap = cgrad(:inferno, nd, categorical = true), colorrange = (1,nd))
            hT1 = labplot!(axT, lon[1], lato, periodo; reduced = :Ti, colormap = cgrad(:inferno, nd, categorical = true), colorrange = (1,nd))
            hL = labplot!(axL, lon[2], lato, periodo; reduced = :longitude, colormap = cgrad(:inferno, Int(round(nL*nl)), categorical = true), colorrange = (1,Int(round(nL*nl))))
            hL1 = labplot!(axL, lon[1], lato, periodo; reduced = :longitude, colormap = cgrad(:inferno, Int(round(nL*nl)), categorical = true), colorrange = (1,Int(round(nL*nl))))
            # m, ax3 = labplot3!(ax3, lon, latlim, timlim; colormap)
            # m1, ax3 = labplot3(ax3, lon2, latlim, timlim; colormap)
        else
            # do once
            hT = labplot!(axT, lono, lato, periodo; reduced = :Ti, colormap = cgrad(:inferno, nd, categorical = true), colorrange = (1,nd)) # colormap = Makie.Categorical(:inferno))
            hL = labplot!(axL, lono, lato, periodo; reduced = :longitude, colormap = cgrad(:inferno, Int(round(nL*nl)), categorical = true), colorrange = (1,Int(round(nL*nl)))) # colormap = Makie.Categorical(:inferno))#
            # m, ax3 = labplot3!(ax3, lonlim, latlim, timlim; colormap)
        end
    
        axT.limits = (lonlim...,latlim...)
        axL.limits = (timlimplot...,latlim...)
        cbarT = Colorbar(figT[1,2], hT, label = "Number of days in labelled events")
        cbarL = Colorbar(figL[1,2], hL, label = "Number of longitudinal increments in labelled events")
        if nd > 5
            cbarT.ticks = ((1+(nd-1)/nd/2):((nd÷5) * (nd-1)/nd):(nd), string.(1:(nd÷5):nd)) 
        else
            cbarT.ticks = ((1+(nd-1)/nd/2):((nd-1)/nd):(nd), string.(1:nd)) 
        end
        if nL*nl > 5
            n = round(nL*nl)
            cbarL.ticks = ((1+(n-1)/n/2) : ((n-1)/n * (n÷5)) : n, string.(1:Int(n÷5):Int(n))) # not needed with Makie.Categorical
        else
            cbarL.ticks = ((1+(n-1)/n/2) : ((n-1)/n) : n, string.(1:Int(n))) # not needed with Makie.Categorical
        end
    
    end
    save(path2 * "v2/fig/plot_" * trial * "_HistEvent_$(obs_event)_LatLon.png", figT)
    save(path2 * "v2/fig/plot_" * trial * "_HistEvent_$(obs_event)_LatTime.png", figL)
    # save(path2 * "v2/fig/plot_" * trial * "event_$(obs_event)_LatLonTime.png", fig3)

    # or 
    # plot them simultaneously in different colours and plot times separately
    
    if !isempty(df) #any(df0.obs_event .== obs_event)
        # bbox obs
        xlims = typeof(lon) <: Vector ? expand(extrema(lon0)) : expand(extrema((lon..., lono...)))
        xlims = xlims[1] >= 180 ? xlims.-360 : xlims
        ylims = expand(extrema((lat...,lato...)))
        
        # selected labels
        lbls = df.label
        ulbls = sort(unique(lbls));# sort # if I do not sort, labels will be sorted by volume
        nlb = length(ulbls)
        # if nlb == 1
            # cols = cgrad(:viridis,1, categorical=true)
        # else
            # cols = cgrad(:viridis, length(ulbls), categorical=true)
            # cols = Makie.Categorical(:tab20)[1:nlb]
            cols = cgrad(:tab20, 20, categorical=true)[1:nlb]
        # end

        # maximum 14 graphs per obs_event (Fig 3x5, with Fig(1,1) = colorbar)
        n = minimum([14, (periodo[2] - periodo[1]).value])
        # period
        time_lapse = maximum(((periodo[2] - periodo[1] + Day(1)) ÷ n, Day(1)));
        # subfigures [1,2],[1,3],[1,4],[1,5],[2,1]
        F = Vector{Union{Nothing, Vector}}(nothing, n)
        i = 1; j = 2;
        for t in 1:n
            F[t]  = [i,j]
            (i, j) = j == 5 ? ((i+1, 1)) : (i, j+1)
        end
        # Set up Figure size as a function of xlims, ylims and number of subplots
        ratio = diff([xlims[1],xlims[2]]) ./ diff([ylims[1], ylims[2]])
        fig = Figure(size = (round(1000 * ratio[1]), 50+250*F[end][1]));
        # colorbar
        if length(ulbls) <= 9
            cbar = Colorbar(fig[1,1], 
                label = "Events labels \n" * obs[obs_event,:Event] * " in " * obs[obs_event, :Area],
                # colormap = cols,
                colormap = cgrad(cols[1:nlb], nlb, categorical=true)
                # size = 40,
                # limits = (1,nlb),
            )
            if length(ulbls) > 1 
                cbar.limits = (1,nlb)
                cbar.ticks = ((1+(nlb-1)/nlb/2):((nlb-1)/nlb):(nlb), string.(ulbls))
            else
                cbar.ticks = ([0.5], string.(ulbls))
            end
            # cbar axis and label to the left
            cbar.flipaxis = false
        else
            # split cbar into 2
            cbar1 = Colorbar(fig[1,1][1,1], 
                label = "Events labels \n" * obs[obs_event,:Event] * " in " * obs[obs_event, :Area],
                colormap = cgrad(cols[1:9], 9, categorical=true),
            )
            cbar1.limits = (1,9)
            cbar1.ticks = ((1+(9-1)/9/2):((9-1)/9):(9), string.(ulbls[1:9]))
            cbar1.flipaxis = false
            cbar2 = Colorbar(fig[1,1][1,2], 
                colormap = cgrad(cols[10:end], 9, categorical=true),
            )
            cbar2.limits = (1,nlb-9)
            cbar2.ticks = ((1+(nlb-9-1)/(nlb-9)/2):((nlb-10)/(nlb-9)):(nlb-9), string.(ulbls[10:end]))
            cbar2.flipaxis = false
        end
        # force column to be 1/5 of fig
        colsize!(fig.layout, 1, Relative(1/(F[end][2]))) # Fixed(150)
        
        for t in 1:n
            # aggregate over time by mode
            periodt = (periodo[1] + (t - 1) * time_lapse, periodo[1] + t * time_lapse -Day(1))
            # plot bounding box
            axt = GeoAxis(fig[F[t][1], F[t][2]],
                # xlabel = F[t][1] == F[end][1] ? "Longitude" : "", 
                # ylabel = F[t][2] == 1 ? "Latitude" : "",
                # aspect = AxisAspect(1), 
                title ="from " * string(periodt[1]) * "\nto " * string(periodt[2]), 
                titlesize=10 ,
                );
            limits!(axt, xlims, ylims,)
            cl=lines!(axt, 
                GeoMakie.coastlines(),
                # x1,y1,
                color = :grey80, linewidth=0.85)
            translate!(cl, 0, 0, 1000)
            # skip timestep if no data
            ind = df.obs_event .== obs_event .&& df.start_time .< periodt[2] .&& df.end_time .>= periodt[1];  
            if any(ind)
                # print("$t : $(any(ind))")
                # pl = DataFrame(x=x,y=y, order=1:length(x)) |> @vlplot(:line, x=:x, y=:y, order=:order) # something wrong with the order in which data are plotted
                latt = (
                    minimum(df.latitude_min[ind]),
                    maximum(df.latitude_max[ind])
                    )
                lont = (minimum(df.longitude_min[ind]), maximum(df.longitude_max[ind])+.25)
                # use lat, lon from obs_event to subset cube
                
                # labels in this time step
                lblt = sort(unique(df.label[ind]))
                # println(lblt)

                if typeof(lon) <: Vector
                    # for L in lon
                    # axt, h = labelplot!(axt, labels, periodt, lato, L, lblt; colormap = Makie.Categorical(cols[indexin(lblt, ulbls)]))
                    # end
                    axt, h = labelplot!(axt, labels, periodt, lato, lon, lblt; colormap = cols[indexin(lblt, ulbls)])
                else 
                    # if length(lblt) == 1
                    #     axt, h = labelplot!(axt, labels, periodt, lato, lono, lblt; color = cols[indexin(lblt, ulbls)])
                    # else
                        axt, h = labelplot!(axt, labels, periodt, lato, lono, lblt; colormap = Makie.Categorical(cols[indexin(lblt, ulbls)]))
                    # end
                end
            end
            # remove ticks
            if F[t][1] != F[end][1]
                axt.xticklabelsvisible = false;
            end
            
            # if F[t][2] > 1
            #     if F[t][2] != 5
            #         axt.yticklabelsvisible = false; 
            #     else
            #         axt.yaxisposition = :right;
            #     end
            # end
            axt.xgridcolor[] = colorant"transparent";
            axt.ygridcolor[] = colorant"transparent";
            axt.xticklabelsvisible = false;
            axt.yticklabelsvisible = false;  
        end
        
        save(path2 * "v2/fig/plot" * "_" * trial * "_Event_$obs_event" * "_labels.pdf", fig) 
    end
end


    

# CSV.write(path * "RecentEvents.csv", obs)

# 
# df0 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", DataFrame, header=1)

 