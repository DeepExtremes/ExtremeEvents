# detect extremes
using Distributed

# sleep(1)
# addprocs(20)
# sleep(1)

pot = 0.01;
ne = 0.1;
# smoothed events?
smoothed = false
sm = smoothed ? "smoothed_" : "ranked_"
# start year?
start_year = "_2016"#"" #

@everywhere begin
    using Pkg
    Pkg.activate("..")
end
@everywhere begin
    using YAXArrays, EarthDataLab
    include("../src/detection.jl")
end

inpath_t = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr"#tmax_smoothed_40"*start_year*".zarr"
inpath_pei = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr"#"/Net/Groups/BGI/scratch/fgans/DeepExtremes/smoothed_pei_ranks.zarr"
r_t = open_dataset(inpath_t)
r_t = r_t[Time = (2016,2022)]
# r_t_old = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr")
r_pei = open_dataset(inpath_pei)
using Dates
inputs = (r_t.layer,
        r_pei.pei_30[Time = (r_t.time[1], r_t.time[end]+Day(1))], 
        r_pei.pei_90[Time = (r_t.time[1], r_t.time[end]+Day(1))],
        r_pei.pei_180[time = (r_t.time[1], r_t.time[end]+Day(1))])


outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_" * sm * "pot" * string(pot) * "_ne" * string(ne) * start_year *".zarr"
tmp = compute_extremes(inputs, pot, outpath; tresne = ne);

print("done.")

# # check plot
# # using Plots
# # tmp = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_ranked_pot0.01_ne0.1.zarr")
# tmp = Cube(outpath)
# include("../src/plots.jl")
# simpleplot(tmp,182, 2019, 4)

