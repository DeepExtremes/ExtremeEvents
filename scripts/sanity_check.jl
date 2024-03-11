# check that documented events get detected by analyse_events
using YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import StatsBase
using Measures
using CairoMakie, GeoMakie
# using PerceptualColourMaps

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/v3"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
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
obs = CSV.read("$(path)../EventPart2.csv", DataFrame; header=3)
# obs = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart2.csv", DataFrame; header=3)
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
    end
    df.start_time = Date.(df.start_time);
    df.end_time = Date.(df.end_time);
    # add column obs_event
    df.obs_event .= obs_event;
    # vcat df
    df0 = vcat(df0, df)
    return df0
end

function labplot(ax, lon, lat, period; reduced = :Ti, kwargs...)
    sublabels = labels.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
    data = ( in(df.label).(sublabels.data))[:,:,:];
    if lon[1] >= 180
        # modify axes
        axs = modaxs(sublabels.axes)
    else
        axs = sublabels.axes
    end 
    if reduced == :longitude
        axs = (axs[1], axs[2], Ti((nd - Dates.value(period[2]-periodl[1])) : nd - Dates.value(period[2]-periodl[2])))
    end
    h = hm!(ax, data, axs; reduced, kwargs...)  
    # but this approach doesn't show if labelled events span outside the observed event bbox
    return h
end

# function labplotTi(ax, period, lat, lon)
#     sublabels = labels.layer[time=period[1]..period[2], latitude=lat[2]..lat[1], longitude=lon[1]..lon[2]]
#     data = ( in(df.label).(sublabels.data))[:,:,:];
#     if lon[1] >= 180
#         # modify axes
#         axs = modaxs(sublabels.axes)
#     else
#         axs = sublabels.axes
#     end 
#     h = hm!(ax, data, axs; colormap = Makie.Categorical(:inferno))  
#     # but this approach doesn't show if labelled events span outside the observed event bbox
#     return h
# end

# function  labplotLon(ax, periodl::Tuple{2}, latl::Tuple{2}, lonl::Tuple{2})
#     sublabels = labels.layer[time=periodl[1]..periodl[2], latitude=latl[1]..latl[2], longitude=lonl[1]..lonl[2]]
#     # load to memory and flag pixels equal to any of the labels
#     data = (sublabels.data .∈ Ref(df[:, :label]))[:,:,:];
#     if lonl[1] >= 180
#         # modify axes
#         axs = modaxs(sublabels.axes) 
#     else
#         axs = sublabels.axes 
#     end
#     h = hm!(ax, data, (axs[1], axs[2], Ti((nd - Dates.value(period[2]-periodl[1])) : nd - Dates.value(period[2]-periodl[2]))); 
#         reduced = :longitude,
#         colormap = cgrad(:inferno, categorical = true), colorrange = (1,length(sublabels.axes[1])))
#     return h
# end
function labplot3(ax3, lon, lat, period; kwargs...)
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

    y_ticks_pos = 1:5:length(tempo);
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

function labplot!(pl, labels, period, lat, lon, ulbls, lblt, cols)
    # subset labelcube
    sublabels = Cube(subsetcube(labels, time=period, latitude=lat, longitude=lon))
    # sublabels = Cube(subsetcube(labels, time=periodt, latitude=latt, longitude=lont))
    # load to memory and set all other values to 0 (so that they will be set to NaN by prephm)
    sublabels1 = (sublabels.data)[:,:,:];
    # @time 
    sublabels1 = map(x -> x in ulbls ? x : 0, sublabels1);
    # cols = ("#28828F", "#6E6E6E", "#9E9E9E", "#C8C8C8", "#366570", "#8C8C8C", "#57A9BA", "#FFD966", "#EAF1F3")
    if length(lblt) == 1
        if length(ulbls) == 1
            colst = palette(unique(vcat([RGBA{Float64}(.5,.5,.5,1)],cols)))
        else
            indc = map(x -> findfirst(y -> y == x, ulbls), lblt)
            colst = palette(unique(vcat([RGBA{Float64}(.5,.5,.5,1)],cols[indc])))
        end
    else
        indc = map(x -> findfirst(y -> y == x, ulbls), lblt)
        tickval = (lblt.-minimum(lblt))./(maximum(lblt) +1 - (minimum(lblt)-1))
        colst = cgrad(cols[indc], tickval, categorical = true)
    end
    if any(lon.>180)
        # modify axes
        axs = modaxs(sublabels.axes)
        # modify sublabels1: shift lon
        shifts = getshifts(axs)
        sublabels1 = circshift(sublabels1, shifts)
    else
        axs = sublabels.axes
    end
    # println(unique(sublabels1))

    pl = hm!(sublabels1, axs = axs, fn = mode, reduced = "tim", 
        c=colst, aspect_ratio=:equal, colorbar = :none)
    # pl = cf!(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", colours=colst, colorbar_ticks = (lblt, lblt),aspect_ratio=:equal)
    # pl = labelplot(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", 
        # title = "from " * string(periodt[1]) * " to " * string(periodt[2]))
    # ml = makielabel(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", 
    #     title = "from " * string(periodt[1]) * " to " * string(periodt[2]))
end

function expand(x::Tuple{Float64, Float64})
    x1 = round(x[1] - 1, RoundDown; digits = -1, base = 5)
    x2 = round(x[2] + 1, RoundUp; digits = -1, base = 5)
    return (x1, x2)
end
##### end helper fns

# for trial in ("ranked_pot0.01_ne0.1_cmp_2016_2021")#,"ranked_pot0.01_ne0.1_tcmp_2016_2021"# "ranked_pot0.01_ne0.1_tcmp_2016_2021","ranked_pot0.01_ne0.1_cmp_Sdiam3_T5_new_2016_2021")
# trial = "ranked_pot0.01_ne0.1_cmp_2016_2021" # "ranked_pot0.01_ne0.1_cmp_2016_2021_land"#"ranked_pot0.01_ne0.1_tcmp_Sdiam3_T5_2016_2021" # "ranked_pot0.01_ne0.1_cmp_2016_2021" # 
trial = "ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022"
landonly = "_landonly"
# events_all = CSV.read(path * "EventStats_ranked_pot0.01_ne0.1_cmp_2016_2021.csv", DataFrame)
    events = CSV.read("$(path)EventStats_$(trial)$(landonly).csv", DataFrame)
    # look for intersection between spatial and temporal range of events from the table or directly in the labelcube
    labelpath = path * "labelcube_$trial.zarr"
    # labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_$trial.zarr"
    labels = open_dataset(labelpath)# labels = Cube(labelpath)
    # labels_all = open_dataset(path * "labelcube_ranked_pot0.01_ne0.1_cmp_2016_2021.zarr")

    global df0 = DataFrame()
    # df0 = CSV.read("/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/SanityCheck_$trial.csv", DataFrame, header=1)
    # df0 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", DataFrame, header=1)

    # loop over observed events
    # for obs_event in 1 : nrow(obs)
    obs_event=19
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
    axL.xticks = [1,nd], string.([period[1],period[2]]);
    # ax.xticklabelrotation = π/4

    # ==== lat-lon-time 3D view
    # colormap
    colormap = [(:white, 0.0), (:orange, 0.5)] 

    fig3 = Figure();
    ax3 = Axis3(fig3[1, 1], 
        perspectiveness = 0.5,
        azimuth = 6.64,
        elevation = 0.57, aspect = (1, 1, 1), #(1,2,1)
        xlabel = "Longitude", ylabel = "Time", zlabel = "Latitude");
    # 
    # === get labelled events maximum bounding box ===
    df = subset(df0, :obs_event => x -> x .== obs_event);
    # if !isempty(df)
        # subset label cube with maximum intersecting bounding box 
        periodl = (minimum(df[: ,:start_time]), maximum(df[: ,:end_time]))
        latl = (minimum(df[: ,:latitude_max]), maximum(df[: ,:latitude_min]))
        lonl = (minimum(df[: ,:longitude_min]), maximum(df[: ,:longitude_max]))
        nl = step(labels.axes[:longitude])

        latlim = (minimum([lat[1], latl[1]])-nl, maximum([lat[2], latl[2]])+nl)
        lonlim = typeof(lon) <: Vector ? (lon0[1]-nl, lon0[2]+nl) : (minimum([lon[1], lonl[1]])-nl, maximum([lon[2], lonl[2]])+nl)
        timlim = (minimum([period[1], periodl[1]])-Day(1), maximum([period[2], periodl[2]])+Day(1))
        timlimplot = ((nd - Dates.value(period[2]-timlim[1])), nd - Dates.value(period[2]-timlim[2]))

        # 3d plot limits and bbox
        ax3.limits = (lonlim..., timlimplot..., latlim...)
        cube = Rect3f(Vec3f(lon0[1], 1, lat[1]), Vec3f(lon0[2]-lon0[1], nd-1, lat[2]-lat[1]))
        bb3 = wireframe!(ax3, cube, alpha = 0.5)
        
        # plot labelled events flattened over :Ti and :longitude
        # handle neg longitude
        if typeof(lon) <: Vector
            # do twice
            # > 180 ( or <0)
            # not possible to plot over longitude outside bbox

            hT = labplot(axT, lon[2], latl, periodl; reduced = :Ti, colormap = cgrad(:inferno, nd, categorical = true), colorrange = (1,nd))
            hT1 = labplot(axT, lon[1], latl, periodl; reduced = :Ti, colormap = cgrad(:inferno, nd, categorical = true), colorrange = (1,nd))
            hL = labplot(axL, lon[2], latl, periodl; reduced = :longitude, colormap = cgrad(:inferno, (diff(collect(lon0))/nl)[1], categorical = true), colorrange = (1,(diff(collect(lon0))/nl)[1]))
            hL1 = labplot(axL, lon[1], latl, periodl; reduced = :longitude, colormap = cgrad(:inferno, (diff(collect(lon0))/nl)[1], categorical = true), colorrange = (1,(diff(collect(lon0))/nl)[1]))
            m, ax3 = labplot3(ax3, lon, latlim, timlim; colormap)
            # m1, ax3 = labplot3(ax3, lon2, latlim, timlim; colormap)
        else
            # do once
            hT = labplot(axT, lonl, latl, periodl; reduced = :Ti, colormap = Makie.Categorical(:inferno))
            hL = labplot(axL, lonl, latl, periodl; reduced = :longitude, colormap = Makie.Categorical(:inferno))#colormap = cgrad(:inferno, categorical = true), colorrange = (1,(diff(collect(lonl))/nl)[1]))
            m, ax3 = labplot3(ax3, lonlim, latlim, timlim; colormap)
        end
    
    axT.limits = (lonlim...,latlim...)
    axL.limits = (timlimplot...,latlim...)
    cbarT = Colorbar(figT[1,2], hT, label = "Number of days in labelled events")
    cbarL = Colorbar(figL[1,2], hL, label = "Number of longitudinal increments in labelled events")
    if typeof(lon) <: Vector
        cbarT.ticks = ((1+(nd-1)/nd/2):((nd-1)/nd):(nd), string.(1:nd)) # not needed with Makie.Categorical
        nL = (diff(collect(lon0))/nl)[1];
        cbarL.ticks = ((1+(nL-1)/nL/2):((nL-1)/nL*5):(nL), string.(1:5:Int(nL))) # not needed with Makie.Categorical
    end
    figT
    figL
    fig3

    save(path * "fig/plot_" * trial * "_HistEvent_$obs_event" * "_LatLon.png", figT)
    save(path * "fig/plot_" * trial * "_HistEvent_$obs_event" * "_LatTime.png", figL)
    save("figs/plot_event_$(obs_event)_LatLonTime.png", fig3)

# end
    # end

    # or 
    # plot them simultaneously in different colours and plot times separately
    
    if typeof(lon) <: Vector
        x = getx(lon0)
    else
        x = getx(lon)
    end
    y = [lat[1], lat[2], lat[2], lat[1], lat[1]];
    
    if any(df0.obs_event .== obs_event)
    # bbox obs
    lato = (minimum(df0.latitude_min[df0.obs_event .== obs_event]), maximum(df0.latitude_max[df0.obs_event .== obs_event]))
    lono = (minimum(df0.longitude_min[df0.obs_event .== obs_event]), maximum(df0.longitude_max[df0.obs_event .== obs_event]))
    periodo = (minimum(df0.start_time[df0.obs_event .== obs_event]), maximum(df0.end_time[df0.obs_event .== obs_event])+Day(1))
    # selected labels
    lbls = df0.label[df0.obs_event .== obs_event]
    ulbls = sort(unique(lbls));
    if length(ulbls) == 1
        cols = cgrad(:viridis,2, categorical=true)[1]
    else
        cols = cgrad(:viridis, length(ulbls), categorical=true)
    end

    # maximum 10 graphs per obs_event
    n = 14
    # period
    time_lapse = maximum(((period[2] - period[1] + Day(1)) ÷ n, Day(1)));
    # l = @layout [a b c d e; f g h i j; k l m n o]
    # global p = ()
    fig = Figure();
    # colorbar
    pl = Plots.heatmap([1:length(ulbls)' 1:length(ulbls)'  ],
         yticks = (1:length(ulbls),ulbls),
         colorbar=:none,
         xticks=:none,
         title = "Events labels \n" * obs[obs_event,:Event] * " in " * obs[obs_event, :Area],
         titlefontsize=10,
         c=cols,
         left_margin = 5mm,
         )
    xlims = typeof(lon) <: Vector ? expand(extrema(lon0)) : expand(extrema((lon..., lono...)))
    xlims = xlims[1] >= 180 ? xlims.-360 : xlims
    ylims = expand(extrema((lat...,lato...)))
    p = (p..., pl)
    for t in 1:n
        # aggregate over time by mode
        periodt = (period[1] + (t - 1) * time_lapse, period[1] + t * time_lapse)
        # println(periodt)
        # plot bounding box
        pl = plot(x, y, 
            xlabel = t in (8,9) ? "Longitude" : "", 
            ylabel = t in (2,4,6,8) ? "Latitude" : "",
            aspect_ratio=:equal, 
            label = "", 
            title ="from " * string(periodt[1]) * "\nto " * string(periodt[2]-Day(1)), 
            titlefontsize=10 ,
            xlims = xlims,
            ylims = ylims,
            );
        # skip timestep if no data
        ind = df0.obs_event .== obs_event .&& df0.start_time .< periodt[2] .&& df0.end_time .>= periodt[1];   
        if !any(ind)
            global p = (p..., pl);
            continue
        else
            # pl = DataFrame(x=x,y=y, order=1:length(x)) |> @vlplot(:line, x=:x, y=:y, order=:order) # something wrong with the order in which data are plotted
            latt = (
                minimum(df0.latitude_min[ind]),
                maximum(df0.latitude_max[ind])
                )
            lont = (minimum(df0.longitude_min[ind]), maximum(df0.longitude_max[ind])+.25)
            # use lat, lon from obs_event to subset cube
            
            # labels in this time step
            lblt = sort(unique(df0.label[ind]))
            # println(lblt)

            if typeof(lon) <: Vector
                for L in lon
                pl = labplot!(pl, labels, periodt, lato, L, ulbls, lblt, cols)
                end
            else 
                pl = labplot!(pl, labels, periodt, lato, lono, ulbls, lblt, cols)
            end
            global p = (p..., pl)

        end
    end
    p1 = plot(p..., layout=l, size = (1000,800));
    # for pl in p
    #     display(pl)
    # end   
    Plots.savefig(p1, path * "fig/plot" * "_" * trial * "_HistEvent_$obs_event" * "_labels.png")
    # Plots.savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/plot" * "_" * trial * "_HistEvent_$obs_event" * "_labels.png")
 

    end
end


    unique!(df0)
    # export to csv
    CSV.write(path * "SanityCheck_$trial.csv", df0)

    # end
# CSV.write(path * "RecentEvents.csv", obs)

# 
# df0 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", DataFrame, header=1)

