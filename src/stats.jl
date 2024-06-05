using OnlineStats, WeightedOnlineStats
using DataFrames
import Dates

# create function to define the basic stats we want to run on the events as pairs of variable => statistic()
function defineallstats()
    (EventLabel(),
    :time => MinMaxTime(),
    :longitude => Extrema(),
    :latitude => Extrema(),
    :t2mmax => WeightedMean(),
    :t2mmax => Extrema(),
    :pei_30 => WeightedMean(),
    :pei_30 => Extrema(),
    :pei_90 => WeightedMean(),
    :pei_90 => Extrema(),
    :pei_180 => WeightedMean(),
    :pei_180 => Extrema(),
    # :gpp => WeightedMean(),
    # :gpp => Extrema(),
    # :nee => WeightedMean(),
    # :nee => Extrema(),
    # :ter => WeightedMean(),
    # :ter => Extrema(),
    EventType(),
    LandShare(),
    Volume(),
    )
end

function defineintensitystats()
    (EventLabel(),
    Intensity(),
    Volume(),
    )
end

function definetairstats()
    (EventLabel(),
    :time => MinMaxTime(),
    :longitude => Extrema(),
    :latitude => Extrema(),
    :t2mmax => WeightedMean(),
    :t2mmax => Extrema(),
    Volume(),
    )
end

function definefluxstats()
    (EventLabel(),
    :gpp => WeightedMean(),
    :gpp => Extrema(),
    :nee => WeightedMean(),
    :nee => Extrema(),
    :ter => WeightedMean(),
    :ter => Extrema(),
    )
end

function definelabelcount()
    (EventLabelCount(),)
end

# define statitsics as mutable structures (fields can be set after construction)

# Event label
# (better to define type in construct rather than in field)
mutable struct EventLabel{I<:Int64}
    el::I
end
EventLabel() = EventLabel(0)
function computestat(el::EventLabel,row)
    el.el = row.label
end
function newlabel(el::EventLabel, ievent)
    el.el = ievent
end

# Event Year
# (better to define type in construct rather than in field)
mutable struct EventYear{I<:Int64}
    ey::I
end
EventYear() = EventYear(Dates.year(Date(1900,1,1)))
function computestat(ey::EventYear,row)
    ey.ey = Dates.year(row.time)
end

# Event label count
mutable struct EventLabelCount{I<:Int64,II<:Int64}
    el::I
    c::II
end
EventLabelCount() = EventLabelCount(0,0)
function computestat(elc::EventLabelCount,row)
    elc.el = row.label
    elc.c = elc.c + 1
end


# Volume of event is proportional to grid cell area and to duration of event
mutable struct Volume{F<:Float64}
    v::F
end
# initialize Volume with value 0.0
Volume() = Volume(0.0)
# Define function computestat for objects of type Volume:
# Update volume by adding grid cells' area (approximated by cosine of latitude)
function computestat(v::Volume,row)
    v.v = v.v + cosd(row.latitude)
end

# Intensity of event is cumulative anomalies for all 4 indicators
mutable struct Intensity{F1<:Float64, F2<:Float64, F3<:Float64, F4<:Float64}
    h::F1
    d30::F2
    d90::F3
    d180::F4
end
Intensity() = Intensity(0.0, 0.0, 0.0, 0.0)
# Define function computestat for objects of type Intensity:
# Update Intensity by adding 1-rank weighted by grid cells' area (approximated by cosine of latitude)
function computestat(in::Intensity, row)
    in.h = in.h .+ (1 .- row.rt) .* cosd.(row.latitude)
    in.d30 = in.d30 .+ (1 .- row.rd30) .* cosd.(row.latitude)
    in.d90 = in.d90 .+ (1 .- row.rd90) .* cosd.(row.latitude)
    in.d180 = in.d180 .+ (1 .- row.rd180) .* cosd.(row.latitude)
end

# Define function computestat for objects of type Pair 
# compute stats as functions of Pair (variable, stat) and row
# use OnlineStats.fit! for running computation of stats, updating data along iterations
function computestat(st::Pair,row)
    variable,s = st
    # skip missing
    if !ismissing(row[variable])
        # skip NaN
        if isa(row[variable], Number) 
            if isnan(row[variable])
                return
            end
        end
        if s isa WeightedOnlineStats.WeightedOnlineStat
            fit!(s,row[variable],cosd(row.latitude))
        else
            fit!(s,row[variable])
        end
    end
end

# define MinMaxTime as mutable structure 
mutable struct MinMaxTime{D<:DateTime,T<:DateTime}
    min::D
    max::T
end
# initialize values of min and max
MinMaxTime() = MinMaxTime(DateTime(3000,1,1), DateTime(1600,1,1))
# update value for each iteration
function OnlineStats.fit!(m::MinMaxTime,v)
    if v < m.min
        m.min = v
    end
    if v > m.max
        m.max = v
    end
end

# stat on proportion of event that is maxt, spei30, spei90, spei_180 and combinations
# Weighted sum of pixels for each type of event. To be divided par Volume()
mutable struct EventType{F1<:Float64,F2<:Float64,F3<:Float64,F4<:Float64,F5<:Float64}
    tm::F1
    d30::F2
    d90::F3
    d180::F4
    td::F5
end
EventType() = EventType(0.0,0.0,0.0,0.0,0.0)
function computestat(et::EventType, row)
    E_NONE=0x00
    E_TAIR=0x01
    E_SPEI30=0x02
    E_SPEI90=0x04
    E_SPEI180=0x08
    # bitwise comparison # bitstring(3) # parse(Int, "00000011", base = 2)
    et_old = EventType(et.tm,et.d30,et.d90,et.d180,et.td)
    et.tm = et.tm + ((E_TAIR & row.event) == E_TAIR) * cosd(row.latitude)
    et.d30 = et.d30 + ((E_SPEI30 & row.event) == E_SPEI30) * cosd(row.latitude)
    et.d90 = et.d90 + ((E_SPEI90 & row.event) == E_SPEI90) * cosd(row.latitude)
    et.d180 = et.d180 + ((E_SPEI180 & row.event) == E_SPEI180) * cosd(row.latitude)
    et.td = et.td + (!iseven(row.event) && row.event > 1) * cosd(row.latitude)
    #@assert et.td <= et.tm && et.td <= et.d30 && et.td <= et.d90 && et.td <= et.d180
    if !(et.td <= et.tm && et.td <= et.d30+et.d90+et.d180)
        println(row)
        println(et)
        print(et_old)
        error()
    end
end

# Land share
mutable struct LandShare{F<:Float64}
    ls::F
end
LandShare() = LandShare(0.0)
function computestat(ls::LandShare,row)
    # println("LandShare old is " * string(ls.ls))
    # print(row)
    # println(row.landmask)
    try  
        ls.ls = ls.ls +  (row.landmask > 0.5 ? 1 : 0 ) * cosd(row.latitude) 
    catch
        # println("no data")
    end
    # println("new LandShare value is $ls")
end

# define function appendresults for all types present in res
appendresult(x::EventLabel,t) = (;t...,label=x.el)

function appendresult(x::Pair{Symbol,<:OnlineStats.Extrema},t)
    namemin = Symbol(x[1],"_min")
    namemax = Symbol(x[1],"_max")
    toappend = (namemin=>minimum(x[2]), namemax=>maximum(x[2]))
    (;t...,toappend...)
end
function appendresult(x::Pair{Symbol, <: MinMaxTime}, t)
    (;t...,start_time = x[2].min, end_time = x[2].max)
end
function appendresult(x::Pair{Symbol,<:WeightedOnlineStats.WeightedMean},t)
    namenew = Symbol(x[1],"_mean")
    toappend = (namenew=>OnlineStats.value(x[2]),)
    (;t...,toappend...)
end
function appendresult(x::Volume,t)
    (;t...,volume = x.v)
end
function appendresult(x::EventType, t)
    (;t...,heat = x.tm, drought30 = x.d30, drought90 = x.d90, drought180 = x.d180, compound = x.td)
end
function appendresult(x::LandShare, t)
    (;t...,land_share = x.ls)
end

function appendresult(x::EventLabelCount, t)
    (;t...,label = x.el, count = x.c)
end

appendresult(x::EventYear,t) = (;t...,year=x.ey)

function appendresult(x::Intensity, t)
    (;t..., inth = x.h, intd30 = x.d30, intd90 = x.d90, intd180 = x.d180,  )
end

"""
    This function is used to create a named tuple from a tuple, appending the elements in a loop
"""
function collectresults(x)
    r = NamedTuple{}()
    for i in x
        #@show typeof(i)
        r = appendresult(i,r)
    end
    r
end


function fitalldata(tab)
    allstats = [defineallstats() for i in 1:1e6];
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            # compute stats for labels > 0 and over land only
            if row.label > 0 && row.landmask > 0.5
                # compute stats for current label
                stat = allstats[row.label]
                map(stat) do st
                    # @show st
                    # skip computestat if stat can't be executed on row (e.g. :gpp => Extrema, for time>2021)
                    # try
                        computestat(st,row)
                    # catch
                        # do nothing
                    # end
                end
            end
        end
    end
    # return allstats
    allstats
end

function fitalldata1(tab)
    allstats = [defineintensitystats() for i in 1:1e6];
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            # compute stats for labels > 0 and over land only
            if row.label > 0 && row.landmask > 0.5
                # compute stats for current label
                stat = allstats[row.label]
                map(stat) do st
                    # @show st
                    # skip computestat if stat can't be executed on row (e.g. :gpp => Extrema, for time>2021)
                    # try
                        computestat(st,row)
                    # catch
                        # do nothing
                    # end
                end
            end
        end
    end
    # return allstats
    allstats
end

# apply collect results (convert to DataFrame) and add derived stats
function toDF(results)
    df = DataFrame(map(collectresults, results))
    # drop empty lines if any
    df = df[map(>(0), df.label), :]

    # add derived stats to df
    df.duration = map(Dates.Day, (df.end_time - df.start_time )) .+ Dates.Day(1);
    # average area (should be multiplied by (Rt=6.371 * 1e3)^2 * (2pi/360 * 0.25)^2) to get sqkm)
    df.area = map(((x,y) -> x/(Dates.value.(y))), df.volume, df.duration);
    # Weighted sums need to be divided by Volume
    df.heat = map((x,y) -> round(x/y*100, digits = 2), df.heat, df.volume);
    df.drought30 = map((x,y) -> round(x/y*100, digits = 2), df.drought30, df.volume);
    df.drought90 = map((x,y) -> round(x/y*100, digits = 2), df.drought90, df.volume);
    df.drought180 = map((x,y) -> round(x/y*100, digits = 2), df.drought180, df.volume);
    df.compound = map((x,y) -> round(x/y*100, digits = 2), df.compound, df.volume);
    try df.land_share = map((x,y) -> round(x/y*100, digits = 2), df.land_share, df.volume); catch; end
    # temperature in deg Celsius
    df.t2mmax_mean = map((x) -> x - 273.15, df.t2mmax_mean);
    df.t2mmax_min = map((x) -> x - 273.15, df.t2mmax_min);
    df.t2mmax_max = map((x) -> x - 273.15, df.t2mmax_max);

    return df

end

function toDFtair(results)
    df = DataFrame(map(collectresults, results))
    # drop empty lines if any
    df = df[map(>(0), df.label), :]

    # add derived stats to df
    df.duration = map(Dates.Day, (df.end_time - df.start_time + Dates.Day(1)));
    # average area (should be multiplied by (Rt=6.371 * 1e3)^2 * (2pi/360 * 0.25)^2) to get sqkm)
    df.area = map(((x,y) -> x/(Dates.value.(y) + 1)), df.volume, df.duration);
    # temperature in deg Celsius
    df.t2mmax_mean = map((x) -> x - 273.15, df.t2mmax_mean);
    df.t2mmax_min = map((x) -> x - 273.15, df.t2mmax_min);
    df.t2mmax_max = map((x) -> x - 273.15, df.t2mmax_max);

    return df

end


function fittairdata(tab)
    allstats = [definetairstats() for i in 1:1e6];
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            if row.label > 0
                # compute stats for current label
                stat = allstats[row.label]
                map(stat) do st
                    # print(st)
                    computestat(st,row)
                end
            end
        end
    end
    # return allstats
    allstats
end

function fitfluxdata(tab)
    allstats = [definefluxstats() for i in 1:1e6];
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            if row.label > 0
                # compute stats for current label
                stat = allstats[row.label]
                map(stat) do st
                    # print(st)
                    computestat(st,row)
                end
            end
        end
    end
    # return allstats
    allstats
end

function fitsanitycheck(tab, allstats, ievent)
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            if row.event > 0x00 && row.event < 0x10
                # compute stats for current label
                stat = allstats[ievent]
                map(stat) do st
                    try
                        computestat(st,row)
                    catch
                        # label doesn't exist, replace by ievent
                        newlabel(st,ievent)
                    end
                end
            end
        end
    end
    # return allstats
    return allstats
end

function sanity_check(obs, eventspath, peis, era)
    eventcube = open_dataset(eventspath)
    allstats = [defineallstats() for i in 1:(size(obs)[1])];
    for ievent in 1: (size(obs)[1])
        period = (Date(replace(obs.when_from[ievent], "." => "-")) , Date(replace(obs.when_until[ievent], "." => "-")))
        lat = (obs.where_NE[ievent], obs.where_NW[ievent])
        lon = (obs.where_SW[ievent],obs.where_SE[ievent])
        cube_table = CubeTable(
            event    = eventcube.layer[time=period, latitude=lat, longitude=lon],
            pei_30  = peis.pei_30[time=period, latitude=lat, longitude=lon], 
            pei_90  = peis.pei_90[time=period, latitude=lat, longitude=lon], #
            pei_180 = peis.pei_180[time=period, latitude=lat, longitude=lon], #
            t2mmax   = era.t2mmax[time=period, latitude=lat, longitude=lon], 
            t2m      = era.t2m[time=period, latitude=lat, longitude=lon],
            t2mmin   = era.t2mmin[time=period, latitude=lat, longitude=lon],
            tp       = era.tp[time=period, latitude=lat, longitude=lon],
            pet      = era.pet[time=period, latitude=lat, longitude=lon],
            landmask = lsmask.lsm[latitude=lat, longitude=lon],
        );
        # run some statitsics (follow stats_extremes?)
        allstats = fitsanitycheck(cube_table, allstats, ievent);
    end
    df = toDF(allstats);
    return df
end

function fitalllabelcount(tab)
    allstats = [definelabelcount() for i in 1:1e6];
    # iterate over CubeTable
    for t in tab
        # loop on rows
        for row in Tables.rows(t)
            if row.label > 0
                # compute stats for current label
                stat = allstats[row.label]
                map(stat) do st
                    # print(st)
                    computestat(st,row)
                end
            end
        end
    end
    # return allstats
    allstats
end

function toDFlab(results)
    df = DataFrame(map(collectresults, results))
    # drop empty lines if any
    df = df[map(>(0), df.label), :]

    # add derived stats to df
    df.duration = map(Dates.Day, (df.end_time - df.start_time));
    # average area (should be multiplied by (Rt=6.371 * 1e3)^2 * (2pi/360 * 0.25)^2) to get sqkm)
    df.area = map(((x,y) -> x/(Dates.value.(y) + 1)), df.volume, df.duration);
    # temperature in deg Celsius
    df.t2mmax_mean = map((x) -> x - 273.15, df.t2mmax_mean);
    df.t2mmax_min = map((x) -> x - 273.15, df.t2mmax_min);
    df.t2mmax_max = map((x) -> x - 273.15, df.t2mmax_max);

    return df

end

# mode function, with init and skipmissing
# init = 0
# build freq table
# look for unique values
# if only one value, set mode = value
# else 
# count them
# select the unique value with the largest count, discarding 0
function mode(itr::AbstractVector; init = missing)
    if isempty(itr)
        return init
    end
    vals = unique(itr)
    if length(vals) == 1
        return vals[1]
    end
    freq = map(vals) do v
        count(itr .== v)
    end
    sp = sortperm(freq, rev=true)
    if vals[sp[1]] != 0
        return vals[sp[1]]
    else
        return vals[sp[2]]
    end
end

function mode(itr::AbstractArray; dims=:, init=0)
    # apply mode, reducing by dim.
    if typeof(dims) <: Colon
        # transform to Vector
        return mode(itr[dims], init=init)
    else
    # elseif typeof(dims) <: Int
    #     # reduce over single dimension
    #     sz = collect(size(itr))
    #     sz[dims] = 1
    #     output = convert(Array{Int,length(sz)}, zeros(tuple(sz...)))
    #     # transform array into array of arrays
    #     ##### ????
        return mapslices(mode, itr, dims)
    end
    
end

# aggregate count by Year
# use EventType and Year
function definecountyear()
    (
    EventYear(),
    EventType(),
    )
end

function countyear(eventcube::Dataset)
    @show cube_table = CubeTable(
        event = eventcube.layer,
        )
    years = unique(map(Dates.year,eventcube.axes[:time]))
    allstats = [definecountyear() for iyear in 1:length(years)];
    # iterate over CubeTable
    for t in cube_table
        # loop on rows
        for row in Tables.rows(t)
            if row.event > 0
                # @show row
                y = Dates.year(row.time)
                iyear = findfirst(x-> x==y,years)
                # compute stats for current event
                stat = allstats[iyear]
                map(stat) do st
                    computestat(st,row)
                end
            end
        end
    end
    # toDF
    df = DataFrame(map(collectresults, allstats))
    # drop empty lines if any
    df = df[map(>(0), df.year), :]
    return df    
end

# sanity check helper functions
function expand(x::Tuple{Int, Int})
    convert(Tuple{Float64, Float64}, x)
end
function expand(x::Tuple{Float64, Float64})
    x1 = round(x[1] - 1, RoundDown; digits = -1, base = 5)
    x2 = round(x[2] + 1, RoundUp; digits = -1, base = 5)
    return (x1, x2)
end