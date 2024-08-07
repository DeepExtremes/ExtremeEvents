using Distributed, SlurmClusterManager

if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
        addprocs(1)
        sleep(0.001)
    end
end

@everywhere using YAXArrays, Zarr, NetCDF

lsm = Cube("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")[
    time = At(DateTime("2019-01-01T13:00:00")),
    # region = region,
    ]

labelpath = "/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/MergeLabelCubes/mergedlabels.zarr"
labels = Cube(labelpath )# 
largest = [42561  51252  24092  55983  55755  25632  18958  44770  36790  53015]
longest = [31767  31866  49981  44843  18998  24109  42561  28223  50340  54071]

# test on small subset (42561 (Russian heatwave of 2010) + 18958)
lon1 = 50 .. 51; lat1 = 54 .. 55;  

# function getllabels!(xout,lbls,largest,longest)
#     # more events in recent years, hence reverse lbls
#     lbls1 = reverse(lbls)
#     islargest = in.(largest, Ref(lbls1))
#     xout[1] = all(islargest .== 0) ? 0 : first(largest[islargest])
#     islongest = in.(longest, Ref(lbls1))
#     xout[2] = all(islongest .== 0) ? 0 : first(longest[islongest])
#     return xout
# end
# # without reverse, with last
# # 25.532521 seconds (2.42 M allocations: 39.085 GiB, 6.54% gc time, 2.92% compilation time)
# # with reverse, with first
# # 23.977383 seconds (3.71 M allocations: 38.867 GiB, 7.15% gc time, 6.13% compilation time)
# # 2-element Vector{Int64}:
# #  18958
# #  42561

@everywhere function getllabelsfaster!(xout,lbls,ls,largest,longest)
    if ls[1] < 0.5
        return xout[:] .= 0
    end
    # more events in recent years, hence reverse lbls
    lbls1 = reverse(lbls)
    # i want the smallest events among the largest. Hence reverse as input
    islargest = 0; islongest = 0
    for i in largest
        if in(i, lbls1)
            islargest = i
            break
        end
    end
    xout[1] = islargest
    for i in longest
        if in(i, lbls1)
            islongest = i
            break
        end
    end
    xout[2] = islongest
    return xout
end
# 8.996266 seconds (758.49 k allocations: 15.124 GiB, 7.74% gc time, 0.18% compilation time)
# 2-element Vector{Int64}:
#  18958
#  42561

# xout = [0,0];
# lbls = labels[Ti = DateTime(1970,1,1,) .. DateTime(2022,12,31,23,1,1), longitude = At(lon1, atol = 0.25), latitude = At(lat1, atol = 0.25),]
# lbls = labels[Ti = DateTime(1970,1,1,) .. DateTime(2022,12,31,23,1,1), longitude = lon1, latitude = lat1,]
# # lablels[Ti = DateTime(1970,1,1,) .. DateTime(2022,12,31,23,1,1)]
# ls = lsm[longitude = At(lon1, atol = 0.25), latitude = At(lat1, atol = 0.25)]
# ls = lsm[longitude = lon1, latitude = lat1]
# getllabelsfaster!(xout,lbls,ls,reverse(largest),reverse(longest))

@time llabels = mapCube(
    getllabelsfaster!,
    (labels, lsm),
    reverse(largest),
    reverse(longest);
    indims = (InDims(:Ti), InDims()),
    outdims = OutDims(
        Dim{:Variable}(["largest", "longest"]),
        outtype = Int64,
        path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/largest_longest_labels.zarr",
        backend = :zarr,
        overwrite = true,
        chunksize=Dict("Variable"=>1),
        )
)

# up to 30 seconds per grid cell (when no event)
# calc * lon * lat / seconds / minutes / hours/ nthreads / nprocs
# 30 * 1440 * 721 / 60 / 60 / 24 /8 /8 # 5.6 days!
# over land only: would reduce time by half

# completely wrong calculations!
# 224.836979 seconds (25.62 M allocations: 1.699 GiB, 0.49% gc time, 6.65% compilation time)
