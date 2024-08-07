using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
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

pei = Cube(open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr/"))
@time qref(
    pei,
    "$(path)qref_pei_1971_2000.zarr/";
    ref = (1971,2000), 
    q = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
)
@time qref(
    pei,
    "$(path)qref_pei_1981_2010.zarr/";
    ref = (1981,2010), 
    q = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
)
@time qref(
    pei,
    "$(path)qref_pei_1991_2020.zarr/";
    ref = (1991,2020), 
    q = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975, 0.99],
    overwrite = true,
    backend = :zarr
)
println("done!")
