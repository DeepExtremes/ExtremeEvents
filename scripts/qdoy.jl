using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:ENV["SLURM_NTASKS"]
        addprocs(1)
        sleep(0.001)
    end
end


@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using YAXArrays, Statistics, Zarr, DimensionalData#, DiskArrays
# 
@everywhere include("../src/detection.jl")

path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"

zg = zopen("$(path)ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = Cube(open_dataset(zg))
# tmax = era.t2mmax

@time qdoy(
    era, #[latitude = 40.0 .. 42.0, longitude = 10.0 .. 12.0], #tmax
    "$(path)qdoy_era.zarr/";
    ref = (1971,2000), 
    w = 15, 
    q = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
    )
# Progress: 100%|█████████████████████████████████████████| Time: 6:59:44[K
# 25198.262774 seconds (28.88 M allocations: 1.909 GiB, 0.00% gc time, 0.07% compilation time)

println("done!")

# # test plot
# using CairoMakie
# tqdoy = Cube("$(path)qdoy_era.zarr/")
# q = lookup(tqdoy, :quantiles)
# lon = [20.0, 50.0, 80.0, 290.0]
# lat = [50.0, 30.0, 0.0, -30.0]
# doy = lookup(tqdoy,:doy)
# ndays = length(doy)
# f = Figure(size = (600, 600));
# axs = ()
# for i in eachindex(lat)
#     for j in eachindex(lon)
#         ax = Axis(f[i,j]);
#         for iq in q
#             lines!(ax, 1:ndays, 
#                 tqdoy[latitude = At(lat[i]), longitude = At(lon[j]), quantiles = At(iq)][:].-273.15,
#                 label = string(iq))
#         end
#         hidedecorations!(ax, grid = false)
#         # xticks
#         ax.xticks = (1:90:ndays, doy[1:90:ndays])
#         if i == length(lat)
#             ax.xticksvisible = true
#             ax.xticklabelsvisible = true
#             ax.xticklabelrotation = π/3
#         end
#         # yticks
#         if j == 1
#             ax.yticksvisible = true
#             ax.yticklabelsvisible = true
#         end
#         axs = (axs..., ax)
#     end
# end
# linkxaxes!(axs...)
# linkyaxes!(axs...)
# Label(f[:,0], "Maximum daily temperature at 2m",
#     rotation = π/2)
# for ilon in eachindex(lon)
#     Label(f[0, ilon], "Lon $(lon[ilon])",
#     tellwidth = false)
# end
# for ilat in eachindex(lat)
#     Label(f[ilat,5], "Lat $(lat[ilat])",
#     rotation = -π/2,
#     tellheight = false)
# end
# f[5,:] = Legend(f, axs[1], "Quantiles", nbanks = 5)
# f
# save("$(path)fig/tmax_qdoy.png")
