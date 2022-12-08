# detect extremes
using Distributed

sleep(1)
addprocs(20)
sleep(1)

pot = 0.005;
ne = 0.1;
# smoothed events?
smoothed = true
sm = smoothed ? "smoothed_" : "ranked_"

@everywhere begin
    using Pkg
    Pkg.activate(".")
end
@everywhere begin
    using YAXArrays, EarthDataLab
    include("../src/detection.jl")
end

inpath_t = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr"
inpath_pei = "/Net/Groups/BGI/scratch/fgans/DeepExtremes/smoothed_pei_ranks.zarr"
r_t = open_dataset(inpath_t)
r_pei = open_dataset(inpath_pei)
inputs = (r_t.layer,
        r_pei.pei_30, 
        r_pei.pei_90,
        r_pei.pei_180)


outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * ".zarr"
tmp = compute_extremes(inputs, pot, outpath; tresne = ne);

print("done.")

# check plot
# using Plots
# tmp = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_ranked_pot0.01_ne0.1.zarr")
tmp = Cube(outpath)
include("../src/plots.jl")
simpleplot(tmp,182, 2019, 4)

