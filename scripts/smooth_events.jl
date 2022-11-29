using Distributed

sleep(1)
addprocs(20)
sleep(1)

@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere begin
    using YAXArrays, EarthDataLab, Statistics, SphericalConvolutions
    include("detection.jl")
end

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
tair = era.t2mmax#[time=1980:2021]

pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")
peicube = Cube(pei)

# daily maximum temperature extremes of interest are in the upper range, so we multiply values by -1 before we normalize them
ranked_t = rescale(tair,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr", multiplier = -1)
# PEI values of interest are in the lower range (P-E << 0), so we don't transform them before we normalize them
ranked_pei = rescale(peicube,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")

smoothed = smooth(ranked_t, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr")
smoothed = smooth(ranked_pei, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_smoothed.zarr")


