using YAXArrays, OnlineStats, WeightedOnlineStats, Zarr, NetCDF
using DimensionalData
using DimensionalData.LookupArrays
using DataFrames, Dates
import CSV
import Statistics #StatsBase
# import Plots
# using Measures
using CairoMakie, GeoMakie

include("../src/plots.jl")

# load absolute percentiles computed over 30 years w/o seasonality
qref = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_era_1971_2000.zarr")
qreftp = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/qref_eratp_1971_2000.zarr")

## plot 90 percentile of total daily precipitation

# get tp as a YAXArray
cube = Cube(qreftp)
# get dimensions
lon = lookup(cube, :longitude);
lat = lookup(cube, :latitude)[end:-1:1];
axs = modaxs(cube.axes)
lon_dim =  dimnum(axs, :longitude);
nlon = collect(axs[lon_dim]);
# drop quantiles dimension from axs
qdim = dimnum(axs, :quantiles)
vaxs = [axs...]
deleteat!(vaxs, qdim)
axs = Tuple(x for x in vaxs)
# shift longitude
shifts = getshifts(axs)

# map given percentile
function pfig(x; colorrange=(0.001,0.02), highclip = :yellow, lowclip = nothing, colormap = :viridis)
f = Figure();
ax = GeoAxis(f[1,1])
s = surface!(ax, nlon, lat, circshift(cube[quantiles = At(x)].data[:,end:-1:1], shifts),
        # colorscale = Makie.pseudolog10,
        colormap = colormap,
        colorrange = colorrange,
        highclip = highclip,
        lowclip = lowclip,
        # kwargs...
        );
cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :grey50, linewidth=0.85,
        )
translate!(cl, 0, 0, 1000);

Colorbar(f[2,1], s, vertical = false, label = "Daily precipitation $(round(x*100))th percentile [m/day]")
return f
end

f90 = pfig(0.9; colorrange=(0.001,0.05))
f50 = pfig(0.5; colorrange=(0.001,0.01))

# get median of 50th percentile over land
lsm = Cube(open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc"))[time = At(DateTime("2019-01-01T13:00:00"))]
lsmask = ifelse.(lsm .> 0.5, 1, 0);
# @time maskedarray = ifelse.(lsmask, qreftp, 0)

# weighted online mean
t = CubeTable(values = qreftp.layer,
        lsm = lsmask)
c_tbl = DataFrame(t[1]);
first(c_tbl, 5)

fitcube = cubefittable(t, WeightedMean, :values, weight=(i->abs(cosd(i.latitude))), by = (:lsm, :quantiles))
# global threshold: area weighted mean of 0.9 quantile over land
gth = fitcube[lsm = At(1), quantiles = At(0.9)].data[]
# 0.00786, i.e. 8 mm/day

# see how many pixels have their 0.99 quantile > gth
cmap1 = vcat(colorant"grey99", resample_cmap(:berlin,50))
fg = pfig(0.99, colorrange = (gth, 0.075), colormap = Reverse(:lapaz), lowclip = colorant"grey99", highclip = :black)
save("/Net/Groups/BGI/scratch/mweynants/ARCEME/fig/era_ref_tp_99.png", fg, dpi = 300)