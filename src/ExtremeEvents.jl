module ExtremeEvents

using Reexport: @reexport

# # using Pkg;Pkg.add(url="https://github.com/meggart/SPhericalConvolutions.jl.git")

# # load all needed packages
# using YAXArrays, Statistics, Zarr, DiskArrays, DiskArrayTools
# # distributed computing
# using Distributed, SlurmClusterManager
# # label connected components
# using ImageMorphology, ImageFiltering
# # stats
# using OnlineStats, WeightedOnlineStats
# using DataFrames
# import Dates

include("preprocessing.jl")
include("detection.jl")
include("stats.jl")
include("plots.jl")

@reexport using .detection

end
