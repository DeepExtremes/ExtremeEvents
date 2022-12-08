using Zarr, YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats
using DataFrames
import Tables
import CSV
# set threshold for extremes detection
pot = 0.005
ne = 0.1
# smoothed events?
smoothed = true
sm = smoothed ? "smoothed_" : "ranked_"

# restrict to small area for testing
# Germany
# "Europe"
# "World"
# region = "World"
# region doesn't work with longitude 0.0 to 360
region="Germany"

period = (Date("2019-07-20"), Date("2019-07-30"))#2016:2021 # 
aperiod = "_2016_2021"#"Summer2019"#

# Open all Data Cubes
labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * ".zarr"
labels = open_dataset(labelpath)

pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)

eventspath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * ".zarr"
eventcube = open_dataset(eventspath)

# LandSeaMask# LandSeaMask
# not in the same grid now !!!
# lsmask = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.720.static.nc")

# Fluxcom carbon fluxes (gC m^(-2) year^(-1))
# gross primary productivity
gpp = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_GPP_final.zarr"))
renameaxis!(gpp, "Time" => "time")
# net ecosystem exchange
nee = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_NEE_final.zarr"))
renameaxis!(nee, "Time" => "time")
# terrestrial ecosystem respiration
ter = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_TER_final.zarr"))
renameaxis!(ter, "Time" => "time")

# create iterable table with data cube label
# each chunk can be read in memory
tab = CubeTable(
    label    = labels.layer[region=region, time=period], # [region=region,time=period]
    pei_30  = pei.pei_30[region=region, time=period], 
    pei_90  = pei.pei_90[region=region, time=period], #
    pei_180 = pei.pei_180[region=region, time=period], #
    t2mmax   = era.t2mmax[region=region, time=period], 
    t2m      = era.t2m[region=region, time=period],
    t2mmin   = era.t2mmin[region=region, time=period],
    tp       = era.tp[region=region, time=period],
    pet      = era.pet[region=region, time=period],
    event    = eventcube.layer[region=region, time=period],
    gpp      = gpp[region=region, time=period],
    nee      = nee[region=region, time=period],
    ter      = ter[region=region, time=period],
    #landmask = lsmask.lsm[region=region],
    )

#tab[1]

# include functions to compute stats and create DataFrame
include("../src/stats.jl")

# compute all stats on CubeTable
res = fitalldata(tab);
# delete empty lines 
# sort results by decreasing Volume of events
sort!(res, by=i->i[end].v, rev=true);


# convert res to data frame
# res is a vector of unnamed tuples
# first convert tuple to named tuple


df = toDF(res);
# write DataFrame out to CSV file
CSV.write("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventStats_$sm$pot$ne$aperiod.csv", df)
# df = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventStats_$sm$pot$ne$aperiod.csv", DataFrame)
# convert back to NamedTuple
# Tables.rowtable(df)

