tp = era.tp
tp.chunks[24,13,1]
tp.data[tp.chunks[24,13,1]...]
tp.data[1,1,1:5845]
findfirst(map(ismissing,tp.data[1,1,1:5845]))
tp.axes[3][3191]

tpfix = Cube("/Net/Groups/data_BGC/era5/e1/0d25_hourly/tp/1958/tp.hh.fc.era5.09.1958.nc")


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
    allars = [diskarfromyearmonth(settings,y,m) for m in 9:12, y in 1958]

    onear = ConcatDiskArray(reshape(allars,1,1,length(allars)))
    timax = RangeAxis("time",DateTime(1958,9,1):Hour(1):DateTime(1958,12,31,23,0,0))
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
    nts = daysinmonth(y,m)*24
    a = settings.T_NC(fnow,settings.varname,(settings.spatsize...,nts))
    settings.T_CF(a,settings.cfargs...)
end



varname = "tp"; outvarname = "tp"; ext = "fc"; aggfun=sum;
    c = era5varcube(varname,ext=ext)    
    a2 = DiskArrayTools.AggregatedDiskArray(c.data,(1,1,24),aggfun) # aggfun instead of mean

    # #Fix bug in 1950 January precip, just place in March
    # if varname in ("tp","ssrd")
    #     a2.a.parents[1]=a2.a.parents[3]
    # end

    # add aggregation function to properties
    c.properties["aggfun"] = string(aggfun)

    # set time axis 
    c2 = YAXArray([c.longitude,c.latitude,RangeAxis("time", Date(1958,9,1):Day(1):Date(1958,12,31))],a2,c.properties)

    nt = (Symbol(outvarname)=>c2,)
    ds = Dataset(;nt...)
    ds = setchunks(ds,(lon=60,lat=60,time=5844))
    # data are Float64, so 1 chunk (for one variable) would be around 64*60*60*5844/8 = 168 MB

    ds.tp.data[:,:,:]

    # writefac: read/write speed factor. Estimate of time ratio between reading input and writing output.
    # Usually reading input is faster than writing output, hence default value is 4.0
    # Here reading is much slower than writing because we aggregate on the fly the hourly data to daily, therefore value is set to 0.1.
    
    # savedataset(ds, path="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",append=true, max_cache=3e8,writefac=0.1)
    # nothing


# check data size.
# loop over years and months
# check whether data array has expected size.
using YAXArray, NetCDF

tpfix = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$varname/1959/$(varname).hh.$(ext).era5.01.1959.nc")

function checkarrayfromyearmonth(y,m)
    d = "/Net/Groups/data_BGC/era5/e1/0d25_hourly/$(varname)/$(y)/"
    filelist = readdir(d)
    m2 = @sprintf("%02d",m)
    fshort = filter(endswith("$m2.$y.nc"),filelist)
    if length(fshort)>1
        filter!(!contains("back_extension"),fshort)
    end
    @assert length(fshort)==1
    fnow = joinpath(d,fshort[1])
    mfd = open_dataset(fnow)
    nts = length(mfd.axes[:time])
    nlons = length(mfd.axes[:longitude])
    nlats = length(mfd.axes[:latitude])
    ntsa = daysinmonth(y,m)*24
    if nts != ntsa println("Problem with $fshort : time axis has length $nts instead of $ntsa") end
    if nlons != 1440 println("Problem with $fshort : longitude axis has length $(nlons) instead of 1440") end
    if nlats != 721 println("Problem with $fshort : latitude axis has length $(nlats) instead of 721") end
    return nothing
end

[checkarrayfromyearmonth(y,m) for m in 1:12, y in 1950:2022];
tpfix.axes[:time]
