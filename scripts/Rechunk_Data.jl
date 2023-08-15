using SlurmClusterManager, Distributed

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
using EarthDataLab, YAXArrays, Zarr, DiskArrays
using Dates, YAXArrayBase
using DataStructures, DiskArrayTools
using Printf: @sprintf
using Statistics: mean
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
    allars = [diskarfromyearmonth(settings,y,m) for m in 1:12, y in 1950:2022]

    onear = ConcatDiskArray(reshape(allars,1,1,length(allars)))
    timax = RangeAxis("time",DateTime(1950,1,1):Hour(1):DateTime(2022,12,31,23,0,0))
    YAXArray([lonax,latax,timax],onear,props)
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


to_run = [
    ("t2m","t2m","an",mean),
    ("t2m","t2mmax","an",maximum),
    ("t2m","t2mmin","an",minimum),
    # ("tp","tp","fc",sum),
    # ("ssrd","ssrd","fc",sum),
]

pmap(to_run) do (varname,outvarname,ext,aggfun)
    c = era5varcube(varname,ext=ext)    
    a2 = DiskArrayTools.AggregatedDiskArray(c.data,(1,1,24),aggfun) # aggfun instead of mean

    # #Fix bug in 1950 January precip, just place in March
    # if varname in ("tp","ssrd")
    #     a2.a.parents[1]=a2.a.parents[3]
    # end

    # add aggregation function to properties
    c.properties["aggfun"] = string(aggfun)

    # set time axis 
    c2 = YAXArray([c.longitude,c.latitude,RangeAxis("time", Date(1950,1,1):Day(1):Date(2022,12,31))],a2,c.properties)

    nt = (Symbol(outvarname)=>c2,)
    ds = Dataset(;nt...)
    ds = setchunks(ds,(lon=60,lat=60,time=5844))
    # data are Float64, so 1 chunk (for one variable) would be around 64*60*60*5844/8 = 168 MB

    # writefac: read/write speed factor. Estimate of time ratio between reading input and writing output.
    # Usually reading input is faster than writing output, hence default value is 4.0
    # Here reading is much slower than writing because we aggregate on the fly the hourly data to daily, therefore value is set to 0.1.
    savedataset(ds, path="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",append=true, max_cache=3e8,writefac=0.1)
    nothing
end
