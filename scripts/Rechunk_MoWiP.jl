using Pkg
Pkg.activate("$(@__DIR__)/..")

using Distributed
using SlurmClusterManager

using NetCDF
using Printf
using Dates

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere begin
using YAXArrays, Zarr, DiskArrays
using Dates, YAXArrayBase
using NetCDF
using DataStructures
using DiskArrayTools
using Printf: @sprintf
using Statistics: mean
startyear = 2020#1950
endyear = 2020#2023
end

@everywhere function era5varcube(varname;ext="an")
    onenc = Cube("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$varname/1950/$(varname).hh.$(ext).era5.01.1950.nc")
    lonax = onenc.longitude
    latax = onenc.latitude
    T_CF = typeof(onenc.data)
    T_NC = typeof(onenc.data.a)
    props = onenc.properties

    cfargs = (onenc.data.mv,onenc.data.add_offset,onenc.data.scale_factor)
    spatsize = (length(lonax),length(latax))

    settings = (;
        varname,lonax, latax, T_CF, T_NC, props, cfargs, spatsize
    )
    allars = [diskarfromyearmonth(settings,y,m) for m in 1:12, y in startyear:endyear]

    onear = ConcatDiskArray(reshape(allars,1,1,length(allars)))
    timax = Dim{:Ti}(DateTime(startyear,1,1):Hour(1):DateTime(endyear,12,31,23,1,1))
    YAXArray((lonax,latax,timax),onear,props)
end



@everywhere function diskarfromyearmonth(settings,y,m)
    filelist = readdir("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$(settings.varname)/$(y)/")
    m2 = @sprintf("%02d",m)
    fshort = filter(endswith("$m2.$y.nc"),filelist)
    if length(fshort)>1
        filter!(!contains("back_extension"),fshort)
    end
    @assert length(fshort)==1
    fnow = joinpath("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$(settings.varname)/$(y)/",fshort[1])
    ntsa = daysinmonth(y,m)*24
    # check that data array has the right size
    mfd = open_dataset(fnow)
    nts = length(mfd.axes[:time])
    nlons = length(mfd.axes[:longitude])
    nlats = length(mfd.axes[:latitude])
    @assert nts == ntsa "Problem with $fshort : time axis has length $nts instead of $ntsa"
    @assert nlons == 1440 "Problem with $fshort : longitude axis has length $(nlons) instead of 1440"
    @assert nlats == 721 "Problem with $fshort : latitude axis has length $(nlats) instead of 721"
    ## 
    a = settings.T_NC(fnow,settings.varname,(settings.spatsize...,nts))
    settings.T_CF(a,settings.cfargs...)
end

varname= "tp"; outvarname = "tp"; ext = "fc"; aggfun = maximum

tp = era5varcube(varname,ext=ext)
thrtp = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_eratp_1971_2000.zarr").layer
pei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr/"))
thrpei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_pei_1971_2000.zarr"))

@everywhere function rolling_sum(a, npre::Int, npost::Int; oob_value = 0.0) # inspired by https://discourse.julialang.org/t/rolling-sum/30793/4
    @assert 1<=(npre + 1 + npost)<=length(a)
    out = similar(a)
    out[1] = sum(a[1:npost + 1])
    for i in eachindex(out)[2:end]
        if i <= (npre + 1)
            out[i] = out[i-1] + a[i+npost]
        elseif i > (npre + 1) && i <= (length(a) - npost)
            out[i] = out[i-1] - a[i-npre-1] + a[i+npost]
        else 
            out[i] = out[i-1] - a[i-npre-1] 
        end
    end
    return out
end

@everywhere function rolling_sum!(out, a, npre::Int, npost::Int; oob_value = 0.0) # inspired by https://discourse.julialang.org/t/rolling-sum/30793/4
    @assert 1<=(npre + 1 + npost)<=length(a)
    # @assert typeof(out) == typeof(a) && size(out) == size(a)
    out[1] = sum(a[1:npost + 1])
    for i in eachindex(out)[2:end]
        if i <= (npre + 1)
            out[i] = out[i-1] + a[i+npost]
        elseif i > (npre + 1) && i <= (length(a) - npost)
            out[i] = out[i-1] - a[i-npre-1] + a[i+npost]
        else 
            out[i] = out[i-1] - a[i-npre-1] 
        end
    end
    return out[:]
end

# @everywhere function rolling_sum(a, npre::Int, npost::Int; oob_value = 0.0, pre = similar(a)) # inspired by https://discourse.julialang.org/t/rolling-sum/30793/4
#     @assert 1<=(npre + 1 + npost)<=length(a)
#     out = similar(a)
#     out[1] = sum(a[1:npost + 1])
#     for i in eachindex(out)[2:end]
#         if i <= (npre + 1)
#             out[i] = out[i-1] + a[i+npost]
#         elseif i > (npre + 1) && i <= (length(a) - npost)
#             out[i] = out[i-1] - a[i-npre-1] + a[i+npost]
#         else 
#             out[i] = out[i-1] - a[i-npre-1] 
#         end
#     end
#     return out
# end

function trolling_sum(a, n::Int, nthreads) # https://discourse.julialang.org/t/rolling-sum/30793/4
    @assert 1<=n<=length(a)
    nseg=nthreads
    if nseg*n >= length(a)
        return rolling_sum(a,n)
    else
        out = similar(a, length(a)-n+1)
        lseg = (length(out)-1)÷nseg+1
        segments = [(i*lseg+1, min(length(out),(i+1)*lseg)) for i in 0:nseg-1]
        Threads.@threads for (start, stop) in segments
            out[start] = sum(a[start:start+n-1])
            for i in start+1:stop
                out[i] = out[i-1]-a[i-1]+a[i+n-1]
            end
        end
        return out
    end
end

@everywhere function aggregate1d(a, aggfun::Function, aggsize::Int)
    # check that aggregation works: size(a) must be multiples of aggsizes
    # @assert all(broadcast(gcd, size(a), aggsizes) .>= aggsizes)
    nin = length(a)
    @assert gcd(nin, aggsize) >= aggsize
    # ranges = broadcast(x -> UnitRange(1,x), (size(a) ./ aggsizes))
    nout = fld(nin,  aggsize)
    rng = UnitRange(1, nout)
    b = similar(a, nout)
    for i in eachindex(rng)
        # aggsize * (i-1) + 1 : aggsize * i
        b[i] = aggfun(skipmissing(a[aggsize * (i-1) + 1 : aggsize * i]))
    end
    return b
end

@everywhere function rolling_all(a, istart::Int, istop::Int, fun::Function, args...; oob_value = false)
    @assert 1<=(istop - istart)<length(a)
    out = BitArray(undef, size(a))
    # apply function to a
    b = vcat(istart < 0 ? fill(oob_value, -istart) : [], fun(a, args...), istop > 0 ? fill(oob_value, istop) : [])
    # rolling all (shorcircuit and)
    for i in eachindex(a)
        out[i] = all(b[(istart < 0 ? i :  (i + (istop - istart))])
    end
    return out
end

# test on small subset
# brisbane south -27.72,152.73
lon = 152.5 .. 153.5; lat = -28.0 .. -27.0;
lon1 = 152.73; lat1 = -27.72;  
# c = era5varcube(varname,ext=ext)[longitude = lon, latitude = lat, Ti = 1:744]    
# tmp = tp[longitude = At(150), latitude = At(-34), Ti = DateTime(2020,1,1,) .. DateTime(2020,12,31,23,1,1)]
tmp = tp[longitude = lon, latitude = lat, Ti = DateTime(2020,1,1,) .. DateTime(2020,12,31,23,1,1)]
tmp1 = pei[longitude = lon, latitude = lat, Ti = DateTime(2020,1,1,) .. DateTime(2020,12,31,)]
@time a = mapCube(
        tmp, indims = InDims("Ti"), outdims = OutDims("Ti"),
    ) do xout, xin 
        xout[:] = rolling_sum(xin, 11,12) 
    end
# 10.48 M allocations for 2 months???? because first run? needs to open all netcdfs because data are stored as maps?
# 42.462680 seconds (221.87 k allocations: 16.149 MiB, 7.46% gc time, 17.55% compilation time) for 12 months and 1 grid cell
# 49.474261 seconds (4.43 M allocations: 310.067 MiB, 7.58% gc time, 88.78% compilation time) for 12 months and 5x5 grid cells

# now a is in memory and i should be able to test agg on it... only 68.62 KB
aggsize = 24; aggfun = maximum
@time b = mapCube(
        a, indims = InDims("Ti"), outdims = OutDims(tmp1.axes[1]),
    ) do xout, xin
        xout[:] = aggregate1d(xin, aggfun, aggsize)
end
# first time: 4.665948 seconds (1.13 M allocations: 82.109 MiB, 62.51% gc time, 357.35% compilation time)
# second time: 0.290074 seconds (203.06 k allocations: 17.854 MiB, 42.83% gc time, 1310.05% compilation time)
# why so many allocations??? 
# @time aggregate1d(rand(24000), aggfun, aggsize)
# frist time: 0.102384 seconds (36.86 k allocations: 2.858 MiB, 99.73% compilation time)
# second time: 0.000244 seconds (1.00 k allocations: 445.484 KiB)
thr = thrtp[quantiles = At(0.95), longitude = lon, latitude = lat,]
@time c = mapCube(
        (b, thr), indims = (InDims("Ti"), InDims()), outdims = OutDims("Ti"),
    ) do xout, xin...
        xout[:] = xin[1] .> xin[2]
end
# 5.372753 seconds (6.69 M allocations: 454.580 MiB, 6.37% gc time, 591.47% compilation time)

@time d = mapCube(
        (tmp1, thrpei[quantiles = At(0.05), longitude = lon, latitude = lat,]),
        indims = (InDims(:Ti, :Variable), InDims(:Variable)),
        outdims = OutDims(:Ti),
    ) do xout, xin...
    xout[:] = rolling_all(xin[1][:,3], -6, -1, (x,y) -> x .< y, xin[2][3])
end
# 1.966096 seconds (300.37 k allocations: 266.186 MiB, 15.80% gc time, 347.45% compilation time)
# 366×5×5 YAXArray{Union{Missing, Float64},3} 
# why is it not returning Bool?

# OutDims(Dim{:Ti}(Date("2020-01-01") : Day(1) : Date("2020-02-29")))
# DateTime(startyear,1,1):Hour(1):DateTime(endyear,12,31,23,1,1)

@everywhere function getpextremes!(xout, xin...)
    # xin is a tuple whose 1st element has size (Nh,) and 2nd element has size (Nd, 3,)
    p = xin[1]
    pe = xin[2]
    pthr = xin[3]
    pethr = xin[4]

    # Moving sum of 24
    sp = rolling_sum(p,11,12)
    # Max over 24 h
    msp = aggregate1d(sp, maximum, 24)
    # apply precalc tp threshold
    xp = broadcast((x,y) -> x .> y && x .> 0.00786, msp, pthr) # gthr global threshold weighted q0.9 over land
    # pei_30 should be xtr (5th percentile) for the past 5 days: soil surface shoud be dry
    xpe30 = rolling_all(pe[:,3], -6, -1, (x,y) -> x .< y, pethr[3])
    # pei_90 should be xtr (below 5th percentile)
    xpe90 = pe[:,2] .< pethr[2]
    # pei_180 too
    xpe180 = pe[:,1] .< pethr[1]
    # keep all info in an EventCube
    ints = map(x -> UInt8(2^x), 0:4)
    bytear = broadcast(*, hcat(xp, xpe30, xpe90, xpe180), ints[1:4]')
    # reduce along variable dimension
    jointar = reduce(|, bytear, dims=2, init = 0x00) 
    # xout is bitwise OR
    xout[:] = jointar[:]
    
    # 0.220003 seconds (195.92 k allocations: 12.965 MiB, 99.42% compilation time)
end

@time r1 = mapCube( 
        (tmp, tmp1, thrtp[quantiles = At(0.95), longitude = lon, latitude = lat,], thrpei[quantiles = At(0.05), longitude = lon, latitude = lat,]), 
        indims = (InDims(:Ti), InDims(:Ti, :Variable), InDims(), InDims(:Variable)), 
        outdims = OutDims(tmp1.axes[1]),
        ) do xout, xin...
        xout[:] = getpextremes!(xout, xin...)
    end
# 78.911865 seconds (1.85 M allocations: 316.386 MiB, 0.43% gc time, 64.37% compilation time: <1% of which was recompilation)


# other approach to limit allocations:
# Do everything for one day in one go
@everywhere function getpextremes1!(xout, xin...)
    # xout is of type SubArray{Union{Missing, Float64}, 1, Array{Union{Missing, Float64}, 3}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64, Int64}, true} and size (366,)
    # xin is a tuple whose 1st element has size (Nh,), 2nd element has size (Nd, 3,), 3rd has size (1,), 4th has size (1,3,)
    p = xin[1]
    # p is of type SubArray{Union{Missing, Float64}, 1, Array{Union{Missing, Float64}, 3}, Tuple{Int64, Int64, Base.Slice{Base.OneTo{Int64}}}, true} and size (8784,)
    pe = xin[2]
    # pe is of type SubArray{Union{Missing, Float64}, 2, Array{Union{Missing, Float64}, 4}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64, Int64, Base.Slice{Base.OneTo{Int64}}}, false} and size (366, 3)
    pthr = xin[3][1]
    # pthr is of type Float64 and size ()
    pethr = xin[4]
    # pethr is of type SubArray{Union{Missing, Float64}, 1, Array{Union{Missing, Float64}, 3}, Tuple{Int64, Int64, Base.Slice{Base.OneTo{Int64}}}, true} and size (3,)
    # 
    E_P=0x01
    E_PE30=0x02
    E_PE90=0x04
    E_PE180=0x08
    # can be out of function and passed as arg

    # initialise rolling sum    
    rs = similar(p) # 2 allocations # Array{Union{Missing, Float64}}(missing, 8784)
    rs = rolling_sum!(rs, p, 12, 11) # 3 allocations
    # loop over days
    for dday in eachindex(pe[:,1])
            drange = 1 + (dday-1) * 24 : dday * 24
            msp = maximum(rs[drange])
            # extreme Precipitation
            r = (msp > pthr && msp > 0.00786) * E_P  # gthr global threshold weighted q0.9 over land
            # pei_30 should be xtr (5th percentile) for the past 5 days: soil surface shoud be dry
            r = r | ((all(x -> x .< pethr[3], dday > 6 ? pe[dday-6 : dday-1] : pe[1:dday])) * E_PE30)
            ## needs to be made more efficient...
            
            # pei_90 should be xtr (below 5th percentile)
            r = r | ((pe[dday,2] < pethr[2]) * E_PE90)
            # pei_180 too
            r = r | ((pe[dday,1] < pethr[1]) * E_PE180)
            
            # xout is bitwise OR
            xout[dday] = r
    end
    # 0.041947 seconds (10.46 k allocations: 654.133 KiB, 97.71% compilation time)
    return xout[:]
end

@time r2 = mapCube( 
        (tmp, tmp1, thrtp[quantiles = At(0.95), longitude = lon, latitude = lat,], thrpei[quantiles = At(0.05), longitude = lon, latitude = lat,]), 
        indims = (InDims(:Ti), InDims(:Ti, :Variable), InDims(), InDims(:Variable)), 
        outdims = OutDims(tmp1.axes[1]),
        ) do xout, xin...
        xout[:] = getpextremes1!(xout, xin...)
    end
# 50.496656 seconds (606.37 k allocations: 231.477 MiB, 0.51% gc time, 22.03% compilation time)

# # compute next day's 24h rolling_sum
#         rslast = rs[end]
#         rs = Array{Union{Missing, Float64}}(missing, 24) # 1 allocation
#         for hhour in 1:24
#             print("Hour: $hour " )
#             println( (dday - 1) * 24 + hhour + 12 : (dday - 1) * 24 + hhour + 35)
#             if dday < ndays - 1
#                 rs[hhour] = rslast - p[(dday - 1) * 24 + hhour + 12] + p[(dday - 1) * 24 + hhour + 35]
#             elseif dday == ndays - 1
#                 rs[hhour] = rslast - p[(dday - 1) * 24 + hhour + 12]
#             else
#                 # do nothing
#                 break
#             end
#             # update rslast
#             rslast = rs[hhour]
#         end


    # still to do
do
    nt = (Symbol(outvarname)=>c2,)
    ds = Dataset(;nt...)
    ds = setchunks(ds,(lon=60,lat=60,time=5844))
    # data are Float64, so 1 chunk (for one variable) would be around 64*60*60*5844/8 = 168 MB

    # writefac: read/write speed factor. Estimate of time ratio between reading input and writing output.
    # Usually reading input is faster than writing output, hence default value is 4.0
    # Here reading is much slower than writing because we aggregate on the fly the hourly data to daily, therefore value is set to 0.1.
    savedataset(ds, path="/Net/Groups/BGI/scratch/mweynants/ARCEME/output/tp_max24.zarr",append=true, max_cache=3e8,writefac=0.1)
    nothing
end

# test on small subset
# Figure
# https://knowledge.aidr.org.au/resources/hailstorm-south-east-queensland-2020/
    using CairoMakie
    a1 = a[Ti = DateTime(2020,10,15) .. DateTime(2020,11,5,23,0,0), longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25)]
    zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
    era = open_dataset(zg)[Ti = DateTime(a1.Ti[1]) .. DateTime(a1.Ti[end-23]), longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25)]
    b1 = b[Ti = DateTime(a1.Ti[1]) .. DateTime(a1.Ti[end-23]), longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25)]
    pe30 = pei[Ti = DateTime(a1.Ti[1]) .. DateTime(a1.Ti[end-23]), longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25), Variable = At("pei_30")]
    
    fig = Figure(size = (1000,800))
    ax = Axis(fig[1,1])
    bp1 = barplot!(ax,13:24:size(a1)[1], era.tp .*1e3, color=colorant"#7bc5cf", 
        alpha = 0.9, 
        label = "TP");
    bp2 = barplot!(ax, 13:24:size(a1)[1], era.pet, color = colorant"#ffaaa8", 
        alpha = 0.9, 
        label = "ETref")
    
    l1 = lines!(ax,1:(size(a1)[1]), a1 .*1e3, label = "Moving sum 24h", color=colorant"#40acbe", linewidth=2);
    l2 = scatterlines!(ax,13:24:size(a1)[1], b1 .*1e3, label = "Daily Max", color = colorant"#0097a7");
        
    l3 = lines!(ax,13:24:size(a1)[1], pe30, label = "PE_30", color=colorant"#b02f2c", linewidth=2)

    # thresholds
    # q = lookup(thrtp, :quantiles)
    hl1 = hlines!(ax,thrtp[longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25)].data[8:end] .*1e3, color = colorant"#0097a7", linestyle = :dash, label = "TP q0.9, q0.95, q0.975, q0.99")
    hl2 = hlines!(ax,thrpei[longitude = At(lon1, atol=0.25), latitude = At(lat1, atol=0.25), Variable = At("pei_30")].data[1:4], color=colorant"#b02f2c", linestyle = :dash, label = "PE_30 q0.01, 0.025, q0.05, q0.1")
    
    # xticks
    function time_ticks(dates; frac=8, start = 1)
        tempo = string.(Date.(dates))
        lentime = length(tempo)
        slice_dates = range(start, lentime, step=lentime ÷ frac)
        return slice_dates, tempo[slice_dates]
    end 
    tempo = lookup(a1, :Ti)
    ax.xticks = time_ticks(tempo, frac=11, start = 13)
    ax.xticklabelrotation = π / 4
    ax.xticklabelalign = (:right, :center)

    plots = [l1,l2,bp1,bp2, hl1, hl2,l3];
    fig[2,1] = Legend(fig,
        plots,
        map(x -> x.label, plots), 
        position = :lb, 
        orientation = :horizontal,
        nbanks = 2,
        framevisible = false,
        )
    Label(fig[1,1, Top()], 
        text = "Brisbane, $(Date(a1.Ti[1])) to $(Date(a1.Ti[end]))",
        halign = :left
        )
    fig
save("/Net/Groups/BGI/scratch/mweynants/ARCEME/fig/tp_versus_mvdmx.png", fig)
