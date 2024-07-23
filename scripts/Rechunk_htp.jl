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
using YAXArrays, Zarr, DiskArrays
using Dates, YAXArrayBase
using NetCDF
using DataStructures
using DiskArrayTools
using Printf: @sprintf
startyear = 1971
endyear = 2023
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


varname = "tp"
outvarname = "tp"
ext = "fc"

# build data cube with netcdfs
c = era5varcube(varname,ext=ext) 
# set chunks
target_chunks = (longitude = 60, latitude = 7, Ti = length(c.Ti))
# add to properties
c.properties["processing"] = "Rechunked with YAXArrays"
# create dataset to get proper layer name
nt = (Symbol(outvarname)=>c,)
ds = Dataset(;nt...)
ds = setchunks(ds,target_chunks)
# data are Float64, so 1 chunk (for one variable) would be around (64*60*7*53*365+(53÷4))/8 *1e-6 = 65 MB

# writefac: read/write speed factor. Estimate of time ratio between reading input and writing output.
# Usually reading input is faster than writing output, hence default value is 4.0
# Here reading is much slower than writing because we aggregate on the fly the hourly data to daily, therefore value is set to 0.1.
savedataset(ds, path="/Net/Groups/BGI/work_2/scratch/mweynants/ARCEME/era_0d25_hourly_tp_tschunked.zarr", max_cache=1e9,writefac=0.1)