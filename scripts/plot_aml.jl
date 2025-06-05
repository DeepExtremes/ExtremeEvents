# plot annual max label
using YAXArrays, Zarr
using CairoMakie, GeoMakie

path = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"
aml = open_dataset(joinpath(path,"AnnualMaxLabel_ranked_pot0.01_ne0.1_cmp_S1_T3.zarr")).layer

f = Figure();
# gax = GeoAxis(f[1,1], 
#         # limits=(extrema(lon), extrema(lat)), 
#         source="+proj=latlong +datum=WGS84", # src CRS
#         dest="+proj=eqearth", # destination CRS, in which you want to plot
#         # coastlines = true # plot coastlines from Natural Earth, as a reference.
#     )
ax = Axis(f[1,1])
data = circshift(convert(Array{Float64},aml[year=Near(2020)]), (180/0.25, 0));
replace!(data, 0 => NaN);
lon = circshift(map(x -> x >= 180 ? x-360 : x, lookup(aml, :longitude)), (180/0.25));
lat = lookup(aml, :latitude);
h = heatmap!(ax,lon,lat,data);
# coastlines
cl=lines!(ax, 
    GeoMakie.coastlines(),
    color = :black, linewidth=0.85)
translate!(cl, 0, 0, 1000)

cb = Colorbar(f[2,1], h, vertical=false)
f
