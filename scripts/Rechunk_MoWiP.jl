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
startyear = 1971
endyear = 2022
end

ref_period = "1971_2000"
# threshold (for PE, computed for tp as 1-q)
q = 0.05
tp = Cube(open_dataset("/Net/Groups/BGI/work_2/scratch/mweynants/ARCEME/era_0d25_hourly_tp_tschunked.zarr"))
thrtp = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_eratp_raindays_$(ref_period).zarr")
pei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr/"))
thrpei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_pei_$(ref_period).zarr"))

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
        b[i] = aggfun(skipmissing(view(a, aggsize * (i-1) + 1 : aggsize * i)))
    end
    return b
end

@everywhere function rolling_all(a, istart::Int, istop::Int, fun::Function, args...; oob_value = false)
    @assert 1<=(istop - istart)<length(a)
    out = BitArray(undef, size(a))
    # apply function to a
    b = fun(a, args...)
    # rolling all (shorcircuit AND)
    for i in axes(a,1)
        out[i] = all(view(b, max(1, min(i + istart, length(a))) :  max(1, min(i + istop, length(a)))))
    end
    return out
end

# test on small subset
# brisbane south -27.72,152.73
lon = 152.5 .. 153.5; lat = -28.0 .. -27.0;
lon1 = 152.73; lat1 = -27.72;  

tmp = tp[longitude = lon, latitude = lat, Ti = DateTime(2020,1,1,) .. DateTime(2020,12,31,23,1,1)]
tmp1 = pei[longitude = lon, latitude = lat, Ti = DateTime(2020,1,1,) .. DateTime(2020,12,31,)]
pthr = thrtp[quantiles = At(0.95), longitude = lon, latitude = lat,]
pethr = thrpei[quantiles = At(0.05), longitude = lon, latitude = lat,]

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
        (tmp, tmp1, thrtp[quantiles = At(q-1), longitude = lon, latitude = lat,], thrpei[quantiles = At(q), longitude = lon, latitude = lat,]), 
        indims = (InDims(:Ti), InDims(:Ti, :Variable), InDims(), InDims(:Variable)), 
        outdims = OutDims(tmp1.axes[1]),
        ) do xout, xin...
        xout[:] = getpextremes!(xout, xin...)
    end
# 78.911865 seconds (1.85 M allocations: 316.386 MiB, 0.43% gc time, 64.37% compilation time: <1% of which was recompilation)


# other approach to limit allocations:
# Do everything for one day in one go
@everywhere function getpextremes1!(xout, p, pe, pthr, pethr, rschannel)
    
    E_P=0x01
    E_PE30=0x02
    E_PE90=0x04
    E_PE180=0x08

    # hourly rolling sum    
    rs = take!(rschannel) 
    rs = rolling_sum!(rs, p, 12, 11)

    # loop over days
    for dday in axes(pe,1)
            drange = 1 + (dday-1) * 24 : dday * 24
            msp = maximum(view(rs,drange))
            # extreme Precipitation
            r = (msp > pthr[1] && msp > 0.00786) * E_P  # gthr global threshold weighted q0.9 over land
            # pei_30 should be xtr (5th percentile) for the past 5 days: soil surface shoud be dry
            r = r | ((all(x -> x .< pethr[3], view(pe, max(1, min(dday-5, length(pe))) :  max(1, min(dday-1, length(pe))) ))) * E_PE30)
                        
            # pei_90 should be xtr (below 5th percentile)
            r = r | ((pe[dday,2] < pethr[2]) * E_PE90)
            # pei_180 too
            r = r | ((pe[dday,1] < pethr[1]) * E_PE180)
            
            # xout is bitwise OR
            xout[dday] = r
    end
    put!(rschannel,rs)
    # 0.041947 seconds (10.46 k allocations: 654.133 KiB, 97.71% compilation time)
    return xout
end
rschannel = Channel{Vector{Float64}}(Threads.nthreads())
for _ in 1:Threads.nthreads()
    put!(rschannel,zeros(length(tmp1.Ti)))
end

# mapCube
# xin1: hourly precipitation
xin1 = tmp;
# xin2 : daily PE (pei_180, pei_90, pe30)
xin2 = tmp1;
# xin3: daily total precipitation q-quantile
xin3 = thrtp[quantiles = At(1-q), longitude = lon, latitude = lat,];
# daily PE_x q-quantiles 
xin4 = thrpei[quantiles = At(q), longitude = lon, latitude = lat,];

# set time chunks to 4 years
@time r2 = mapCube( 
        getpextremes1!,
        (xin1, xin2, xin3, xin4), 
        rschannel,
        indims = (InDims(:Ti), InDims(:Ti, :Variable), InDims(), InDims(:Variable)), 
        outdims = OutDims(tmp1.axes[1],
            propoerties = Dict(
                :Name => "hpd_events",
                :Description => "Heavy precipitation (>q$(1-q) of rain days, encoded as UInt8(1)) following dry conditions: last 5 days with 30 days evaporative stress (Precipitation - Ref. Evapotranspiration) <q$(q)",
                :DataSource => "Hersbach, H. et al. (2023): ERA5 hourly data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47 (Accessed on 16-04-2024)",
                :Attribution => "Contains modified Copernicus Climate Change Service information (2024)",
                :Authors => "M. Weynants and F. Gans, Max Planck Institute for Biogeochemistry.",
                :CreationDate => "$(Date(now()))",
                :TemporalExtent => "$((Date(startyear),Date(endyear)))",
                :SpatialExtent => sprint(show,(xin1.longitude[1], xin1.latitude[1], xin1.longitude[end], xin1.latitude[end]), context=:compact=>true),
                :CRS => "EPSG:4326",
                ),
            chunksize=Dict("Ti"=>1461),backend=:zarr,
            path="/Net/Groups/BGI/work_2/scratch/mweynants/ARCEME/hpd_eventcube_q$(q)_ref.zarr"),
        ) 
# 50.496656 seconds (606.37 k allocations: 231.477 MiB, 0.51% gc time, 22.03% compilation time)

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
