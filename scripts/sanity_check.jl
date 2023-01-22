# check that documented events get detected by analyse_events
using Revise, YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, Zarr
using DataFrames, Dates
import CSV
import StatsBase


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
obs = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventPart2.csv", DataFrame; header=3)
include("../src/plots.jl")
# clean obs
obs = dropmissing!(obs)
obs.Start .= replace.(obs.Start, r"\." => "-");
obs.End .= replace.(obs.End, r"\." => "-");

# helper fn
function labobs(df0::DataFrame, lon, lat, period, obs_event; p=p::Plots.Plot{})
    #
    # plot labelled events flattened over time
    # or select all labels
    # fix the longitudes
    # event spatial bounding box
    x = [lon[1], lon[1], lon[2], lon[2], lon[1]];
    x = lon[1] >= 180 ? x.-360 : x
    y = [lat[1], lat[2], lat[2], lat[1], lat[1]];
    p = plot!(x, y , label = "Historic event");# 

    
    # retrive labels and count them
    tab = CubeTable(
        label    = labels.layer[time=period, latitude=lat, longitude=lon]
    )
    labcount = fitalllabelcount(tab);
    sort!(labcount, by=i->i[end].c, rev=true);
    # toDF
    labcountdf = DataFrame(map(collectresults, labcount));
    # remove event with few voxels and empty lines
    labcountdf = labcountdf[map(>(99), labcountdf.count), :]
    # extract stats from events
    df = filter(:label => in(labcountdf.label), events)
    if isempty(df)
        return p, df0
    end
    df.start_time = Date.(df.start_time);
    df.end_time = Date.(df.end_time);
    # add column obs_event
    df.obs_event .= obs_event;
    # vcat df
    df0 = vcat(df0, df)
    sublabels = Cube(subsetcube(labels, time=period, latitude=lat, longitude=lon))
    sublabels1 = ( in(df.label).(sublabels.data))[:,:,:];
    if lon[1] >= 180
        # modify axes
        axs = modaxs(sublabels.axes)
        p = hm!(sublabels1, axs = axs)
    else
        p = hm!(sublabels1, axs = sublabels.axes)
    end   
    # but this approach doesn't show if labelled events span outside the observed event bbox
    return p, df0
end

# for trial in ("ranked_pot0.01_ne0.1_cmp_2016_2021")#, "ranked_pot0.01_ne0.1_tcmp_2016_2021","ranked_pot0.01_ne0.1_cmp_Sdiam3_T5_new_2016_2021")
trial = "ranked_pot0.01_ne0.1_tcmp_Sdiam3_T5_2016_2021"
    events = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventStats_$trial.csv", DataFrame)
    # look for intersection between spatial and temporal range of events from the table or directly in the labelcube
    labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_$trial.zarr"
    labels = open_dataset(labelpath)# labels = Cube(labelpath)

    # df0 = DataFrame()
    df0 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", DataFrame, header=1)

    # loop over observed events
    for obs_event in 1 : nrow(obs)
    # obs[obs_event,:]
    period =( Date(obs[obs_event,:Start]), Date(obs[obs_event,:End])+Day(1))
    lat = (obs[obs_event,:North], obs[obs_event,:South])
    lon = (obs[obs_event,:West], obs[obs_event,:East])
    # transform lon to match cube
    if lon[1] < 0 
        if lon[2] <= 0
            # shift longitudes
            lon = broadcast(x->x+360,lon)
        else # lon[2] > 0
            # split bbox into 2
            lon1 = (lon[1]+360,360)
            lon2 = (0,lon[2])
            lon0 = lon
            lon = [lon1,lon2]
        end
    end
    # p = plot(title = obs[obs_event,:Event] * " in " * obs[obs_event, :Area] * "\n from " * obs[obs_event, :Start] * " to " * obs[obs_event, :End]);
    
    # if typeof(lon) <: Vector{}
    #     # do everything twice...
    #     p, df0 = labobs(df0, lon[1], lat, period, obs_event);
    #     p, df0 = labobs(df0, lon[2], lat, period, obs_event);
    # else
    #     # do only once
    #     p, df0 = labobs(df0, lon, lat, period, obs_event);
    # end

    # # print(size(df0))
    # # print("\n")

    # # display(p)
    # # but this approach doesn't show events that overlap in time.
    # Plots.savefig(p, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/plot" * "_" * trial * "_HistEvent_$obs_event" * "_LatLon.png")

    # # plot all labelled events separately with time dimension flattened + bbox of obs_event
    # p = plot(title = obs[obs_event,:Event] * " in " * obs[obs_event, :Area] * "\n from " * obs[obs_event, :Start] * " to " * obs[obs_event, :End]);
    # # event spatial bounding box
    # x = [period[1], period[1], period[2], period[2], period[1]]
    # y = [lat[1], lat[2], lat[2], lat[1], lat[1]];
    # plot!(p, x, y , label = "Historic event");
    # df = subset(df0, :obs_event => x -> x .== obs_event);
    # if !isempty(df)
    #     for label_row in 1:nrow(df)
    #             # subset labelcube and load to memory
    #             periodl =( df[label_row,:start_time], df[label_row,:end_time]+Day(1))
    #             #### !!! need to fix longitude to match 0-360 !!!
    #             latl = (df[label_row,:latitude_max], df[label_row,:latitude_min])
    #             lonl = (df[label_row,:longitude_min], df[label_row,:longitude_max])
    #             sublabels = Cube(subsetcube(labels, time=periodl, latitude=latl, longitude=lonl))
    #             # load to memory and flag pixels equal to label
    #             label = df[label_row, :label];
    #             sublabels1 = (sublabels.data .== label)[:,:,:];
    #             # view over time
    #             p = hm!(sublabels1, axs = sublabels.axes, reduced = "lon", xlab="Time", ylab = "Latitude")
    #     end
    # end
    # # display(p)

    # Plots.savefig(p, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/plot" * "_" * trial * "_HistEvent_$obs_event" * "_LatTime.png")

    # or 
    # plot them simultaneously in different colours and plot times separately
    x = [lon[1], lon[1], lon[2], lon[2], lon[1]];
    #### need to fix for Europe etc. e.g obs_event = 5
    x = lon[1] >= 180 ? x.-360 : x
    y = [lat[1], lat[2], lat[2], lat[1], lat[1]];
    
    if any(df0.obs_event .== obs_event)
    # bbox obs
    lato = (minimum(df0.latitude_min[df0.obs_event .== obs_event]), maximum(df0.latitude_max[df0.obs_event .== obs_event]))
    lono = (minimum(df0.longitude_min[df0.obs_event .== obs_event]), maximum(df0.longitude_max[df0.obs_event .== obs_event]))
    periodo = (minimum(df0.start_time[df0.obs_event .== obs_event]), maximum(df0.end_time[df0.obs_event .== obs_event])+Day(1))
    # selected labels
    lbls = df0.label[df0.obs_event .== obs_event]
    ulbls = sort(unique(lbls));
    cols = cgrad(:inferno, length(ulbls), categorical=true)

    # maximum 10 graphs per obs_event
    n = 9
    # period
    time_lapse = (period[2] - period[1] + Day(1)) ÷ n;
    l = @layout [a b ;c d ;e f ;g h ;i j]
    p = ()
    pl = heatmap(collect(1:length(ulbls))',
         xticks = (1:9,ulbls),
         colorbar=:none,
         yticks=:none,
         title = "Events labels \n" * obs[obs_event,:Event] * " in " * obs[obs_event, :Area],
         titlefontsize=10,
        #  xtickfontrotation = 90.0,
         )
    p = (p..., pl)
    for t in 1:n
        # aggregate over time by mode
        periodt = (period[1] + (t - 1) * time_lapse, period[1] + t * time_lapse - Day(1))
        println(periodt)
        # plot bounding box
        pl = plot(x, y, 
            xlabel = t in (8,9) ? "Longitude" : "", 
            ylabel = t in (2,4,6,8) ? "Latitude" : "",
            aspect_ratio=:equal, 
            label = "", 
            title ="from " * string(periodt[1]) * "\nto " * string(periodt[2]), 
            titlefontsize=10 );
        # skip timestep if no data
        ind = df0.obs_event .== obs_event .&& df0.start_time .< periodt[2] .&& df0.end_time .>= periodt[1];   
        if !any(ind)
            p = (p..., pl);
            continue
        else
            # pl = DataFrame(x=x,y=y, order=1:length(x)) |> @vlplot(:line, x=:x, y=:y, order=:order) # something wrong with the order in which data are plotted
            latt = (
                minimum(df0.latitude_min[ind]),
                maximum(df0.latitude_max[ind])
                )
            lont = (minimum(df0.longitude_min[ind]), maximum(df0.longitude_max[ind]))
            # subset labelcube
            sublabels = Cube(subsetcube(labels, time=periodt, latitude=lat, longitude=lon))
            # sublabels = Cube(subsetcube(labels, time=periodt, latitude=latt, longitude=lont))
            # load to memory and set all other values to 0 (so that they will be set to NaN by prephm)
            sublabels1 = (sublabels.data)[:,:,:];
            @time sublabels1 = map(x -> x in ulbls ? x : 0, sublabels1);
            # need to fix mode!!!
            lblt = sort(unique(df0.label[ind]))
            println(lblt)
            indc = map(x -> findfirst(y -> y == x, lbls), lblt)
            # cols = ("#28828F", "#6E6E6E", "#9E9E9E", "#C8C8C8", "#366570", "#8C8C8C", "#57A9BA", "#FFD966", "#EAF1F3")
            if length(lblt) == 1
                colst = palette(unique(vcat(cols[1],cols[indc],cols[end])))
                #,  (unique(vcat(lbls[1],lblt, lbls[end])).-minimum(lbls))./(maximum(lbls) - (minimum(lbls))), rev = true, categorical = true)
            else
                tickval = (lblt.-minimum(lblt))./(maximum(lblt) +1 - (minimum(lblt)-1))
                colst = cgrad(cols[indc], tickval, rev = true, categorical = true)
            end
            pl = hm!(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", 
                c=colst, aspect_ratio=:equal, colorbar = :none)
            # pl = cf!(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", colours=colst, colorbar_ticks = (lblt, lblt),aspect_ratio=:equal)
            # pl = labelplot(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", 
                # title = "from " * string(periodt[1]) * " to " * string(periodt[2]))
            # ml = makielabel(sublabels1, axs = sublabels.axes, fn = mode, reduced = "tim", xlab="Longitude", ylab = "Latitude", 
            #     title = "from " * string(periodt[1]) * " to " * string(periodt[2]))
            p = (p..., pl)

        end
    end
    p1 = plot(p..., layout=l, size = (450,1000));
    # for pl in p
    #     display(pl)
    # end   
    Plots.savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/plot" * "_" * trial * "_HistEvent_$obs_event" * "_labels.png")
 

    end
end


    unique!(df0)
    # export to csv
    # CSV.write("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", df0)

    # end
# CSV.write("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/RecentEvents.csv", obs)

# 
# df0 = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/SanityCheck_$trial.csv", DataFrame, header=1)
# for obs_event = 1:maximum(df0.obs_event)
#     # bbox obs
#     lato = (minimum(df0.latitude_min[df0.obs_event .== obs_event]), maximum(df0.latitude_max[df0.obs_event .== obs_event]))
#     lono = (minimum(df0.longitude_min[df0.obs_event .== obs_event]), maximum(df0.longitude_max[df0.obs_event .== obs_event]))
#     periodo = (minimum(df0.start_time[df0.obs_event .== obs_event]), maximum(df0.end_time[df0.obs_event .== obs_event])+Day(1))
#     # subset labelcube
#     sublabels = Cube(subsetcube(labels, time=periodo, latitude=lato, longitude=lono))


