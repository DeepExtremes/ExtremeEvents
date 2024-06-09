using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using YAXArrays, Statistics, Zarr#, DiskArrays
# 
@everywhere include("../src/detection.jl")

zg = zopen("$(path)ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
tmax = era.t2mmax

qdoy(
    tmax, 
    "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qdoy_tmax.zarr/";
    ref = (1971,2000), 
    w = 15, 
    q = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
    )

# close workers

println("done!")

# test plot
# tqdoy = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qdoy_test.zarr/")
# f = Figure();
# axs = ()
# lon = lookup(tqdoy,:longitude)
# lat = lookup(tqdoy,:latitude)
# for i in eachindex(lat)
#     for j in eachindex(lon)
#         ax = Axis(f[i,j]);
#         for iq in q
#             lines!(ax, lookup(tqdoy,:doy), 
#                 tqdoy[latitude = At(lat[i]), longitude = At(lon[j]), quantiles = At(iq)][:])
#         end
#         hidedecorations!(ax, grid = false)
#         axs = (axs..., ax)
#     end
# end
# linkxaxes!(axs...)
# linkyaxes!(axs...)
# # seems to be a jump at 29 Feb...