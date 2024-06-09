# Compute summary statistics on labelled extrene events
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using Zarr, YAXArrays, EarthDataLab, OnlineStats, WeightedOnlineStats, NetCDF
@everywhere using DataFrames
import Tables
import CSV

@everywhere begin
# set threshold for extremes detection
pot = 0.01
ne = 0.1
# smoothed events?
smoothed = false
sm = smoothed ? "smoothed_" : "ranked_"
# compound events
compound_events = true
cmp = compound_events ? "_cmp" : ""

filter_events = true
filter = filter_events ? "_S1_T3" : ""  #"_Sdiam3_T5_new" # "_Sdiam3_T5" #  "_Sdiam4_T4" #

filter_land = false
land = filter_land ? "_land" : ""

# stats run on land only
landonly = "_landonly"

# restrict to small area for testing
# Germany
# "Europe"
# "World"
# region = "World"
# region doesn't work with longitude 0.0 to 360
# region="Luxembourg"
end #begin

# loop over period
# split period into 13 years periods with 3 years overlap
# period = 1950:2022 #2016:2022 # (Date("2018-06-01"), Date("2018-09-16"))#could be added to file name "*"_"*replace("$period", ":" => "_")*"
# aperiod = "_1950_2022" #"_2016_2021"
y = 1950; p = [y:(y+12)]
while y < 2010
    global y += 10
    append!(p,[y:(y+12)])
end
# period = 1950:2021 # 2016:2021 # (Date("2019-07-20"), Date("2019-07-30"))#
# aperiod = "_1950_2021" # "_2016_2021"#"Summer2019"#

@everywhere begin

    # include functions to compute stats and create DataFrame
include("../src/stats.jl")

# Open all Data Cubes
pei = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr")

zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)

eventspath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/EventCube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * ".zarr"
eventcube = open_dataset(eventspath)

# LandSeaMask# LandSeaMask
# not in the same grid now !!!
# lsmask = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.720.static.nc")
lsmask = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
lsmask_notime = subsetcube(lsmask, time=DateTime("2019-01-01T13:00:00"))

# # Fluxcom carbon fluxes (gC m^(-2) year^(-1))
# # gross primary productivity
# gpp = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_GPP_final.zarr"))
# renameaxis!(gpp, "Time" => "time")
# # net ecosystem exchange
# nee = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_NEE_final.zarr"))
# renameaxis!(nee, "Time" => "time")
# # terrestrial ecosystem respiration
# ter = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_TER_final.zarr"))
# renameaxis!(ter, "Time" => "time")

end # begin

@sync @distributed for period in p
    @show aperiod = replace(string(period), ":" => "_", r"^" => "_")

labelpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/labelcube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * cmp * filter * aperiod * land * ".zarr"
labels = open_dataset(labelpath)


# create iterable table with data cube label
# each chunk can be read in memory
tab = CubeTable(
    label    = labels.layer[time=period], # [region=region,time=period] # [time=period]
    pei_30  = pei.pei_30[time=period], 
    pei_90  = pei.pei_90[time=period], #
    pei_180 = pei.pei_180[time=period], #
    t2mmax   = era.t2mmax[time=period], 
    t2m      = era.t2m[time=period],
    t2mmin   = era.t2mmin[time=period],
    tp       = era.tp[time=period],
    pet      = era.pet[time=period],
    event    = eventcube.layer[time=period],
    # gpp      = gpp[ time=period],
    # nee      = nee[ time=period],
    # ter      = ter[ time=period],
    landmask = lsmask_notime.lsm#[region=region],
    )

#tab[1]

# compute all stats on CubeTable
@time res = fitalldata(tab);
# delete empty lines 
# sort results by decreasing Volume of events
sort!(res, by=i->i[end].v, rev=true);

# convert res to data frame
# res is a vector of unnamed tuples
# first convert tuple to named tuple
df = toDF(res)
# write DataFrame out to CSV file
outname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/EventStats_" * sm * "pot" * string(pot) * "_ne" * string(ne)  * cmp * filter * aperiod * land * landonly * ".csv"
CSV.write(outname, df)
print(outname)
print("\n done!")
# df = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventStats_" * sm * "pot" * string(pot) * "_ne" * string(ne)  * cmp * filter * aperiod * ".csv", DataFrame)
# convert back to NamedTuple
# Tables.rowtable(df)

end # pmap

# close workers
rmprocs(workers())

# # combine tables
# # add colum with aperiod
# stats_df = DataFrame()
# for period in p
#     @show aperiod = replace(string(period), ":" => "_", r"^" => "_") 
#     df = CSV.read("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/EventStats_" * sm * "pot" * string(pot) * "_ne" * string(ne)  * cmp * filter * aperiod * land * landonly * ".csv", DataFrame)
#     df.aperiod = repeat(aperiod, size(df, 1))
#     stats_df = vcat(stats_df, df)
# end