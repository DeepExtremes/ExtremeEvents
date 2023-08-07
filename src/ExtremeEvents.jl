module ExtremeEvents
using Revise
# load all needed packages
using YAXArrays, EarthDataLab, Statistics, Zarr, DiskArrays, DiskArrayTools
# distributed computing
using Distributed, SlurmClusterManager
# spatial filter
using SphericalConvolutions
# label connected components
using ImageMorphology, ImageFiltering
# stats
using OnlineStats, WeightedOnlineStats
using DataFrames
import Dates

"""
    splat(f)
Defined as
```julia
    splat(f) = args->f(args...)
```
i.e. given a function returns a new function that takes one argument and splats
its argument into the original function. This is useful as an adaptor to pass
a multi-argument function in a context that expects a single argument, but
passes a tuple as that single argument.
# Example usage:
```jldoctest
julia> map(Base.splat(+), zip(1:3,4:6))
3-element Vector{Int64}:
 5
 7
 9
```
source: https://github.com/JuliaLang/julia/blob/76952a88ea650ae8c6b6b1d010ef695ed9c8244d/base/operators.jl#L1150
"""
# splat(f) = args ->f(args...)

include("preprocessing.jl")
include("detection.jl")
include("stats.jl")
include("plots.jl")
#include("sanity_check.jl")
end
