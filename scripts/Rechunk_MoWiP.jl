using Pkg
Pkg.activate("$(@__DIR__)/..")

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
using Statistics: mean
startyear = 1950
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


to_run = [
    # ("t2m","t2m","an",mean),
    # ("t2m","t2mmax","an",maximum),
    # ("t2m","t2mmin","an",minimum),
    ("tp","tp","fc",maximum),
    # ("ssrd","ssrd","fc",sum),
]

# parallel mapping over variables # here i have only one variable....
# p
map(to_run) do (varname ,outvarname, ext, aggfun)

varname= "tp"; outvarname = "tp"; ext = "fc"; aggfun = maximum

# test on small subset
# lon = 150 .. 151; lat = -34.0 .. -33.0;     
# c = era5varcube(varname,ext=ext)[longitude = lon, latitude = lat, Ti = 1:744]    

c = era5varcube(varname,ext=ext)

    # Moving sum of 24
    r1 = mapCube( 
        c, 
        indims = InDims(MovingWindow("Ti", 12, 11), window_oob_value = missing), 
        outdims = OutDims()
        ) do xout, xin
        xout[] = sum(view(skipmissing(xin)))
    end

    # Max over
    a2 = DiskArrayTools.AggregatedDiskArray(r1.data,(1,1,24), aggfun) # aggfun instead of mean

    # #Fix bug in 1950 January precip, just place in March
    # if varname in ("tp","ssrd")
    #     a2.a.parents[1]=a2.a.parents[3]
    # end

    # add aggregation function to properties
    c.properties["aggfun"] = string(aggfun)

    # set time axis 
    c2 = YAXArray((c.longitude,c.latitude,Dim{:Ti}(Date(c.Ti[1]):Day(1):Date(c.Ti[end]))),a2,c.properties)

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

# # test on small subset
#     # total precipitation daily aggregation
#     a3 = DiskArrayTools.AggregatedDiskArray(c.data,(1,1,24), sum)
#     c3 = YAXArray((c.longitude,c.latitude,Dim{:Ti}(Date(c.Ti[1]):Day(1):Date(c.Ti[end]))),a3,c.properties)
#
# # Figure
#     using CairoMakie
#     fig,ax,l1 = lines(1:(size(r1)[3]), r1[1,1,:] .*1e3, label = "Moving sum 24h", color=:blue, linewidth=2);
#     l2 = lines!(ax,13:24:size(r1)[3], c2[1,1,:] .*1e3, label = "Daily Max", color = :red, alpha = 0.9);
#     b1 = barplot!(ax,13:24:size(r1)[3], c3[1,1,:] .*1e3, color=:blue, 
#         alpha = 0.5, label = "TP")

#     zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
#     era = open_dataset(zg)
#     pet = era.pet[Ti = DateTime(c.Ti[1]) .. DateTime(c.Ti[end]), longitude = lon, latitude = lat, Variable = At("pei_30")]
#     b2 = barplot!(ax, 13:24:size(r1)[3], pet[1,1,:], color = :red, 
#         alpha = 0.5, label = "PET")
    
#     pei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr/"))
#     pe30 = pei[Ti = DateTime(c.Ti[1]) .. DateTime(c.Ti[end]), longitude = lon, latitude = lat, Variable = At("pei_30")]
#     l3 = lines!(ax,13:24:size(r1)[3], pe30[:,1,1], label = "PE_30", color=:black, linewidth=2)
    
#     plots = [l1,l2,b1,b2,l3];
#     fig[2,1] = Legend(fig,
#         plots,
#         map(x -> x.label, plots), 
#         position = :lb, 
#         orientation = :horizontal,
#         nbanks = 2,
#         framevisible = false,
#         )
#     Label(fig[1,1, Top()], 
#         text = "Sydney, $(Date(c.Ti[1])) to $(Date(c.Ti[end]))",
#         halign = :left
#         )
#     fig
# save("/Net/Groups/BGI/scratch/mweynants/ARCEME/fig/tp_versus_mvdmx.png", fig)
