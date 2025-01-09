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
lon = circshift(map(x -> x .> 180 ? x-360 : x, lookup(aml, :longitude)), (180/0.25));
lat = lookup(aml, :latitude);
h = heatmap!(ax,lon,lat,data);
# coastlines
cl=lines!(ax, 
    GeoMakie.coastlines(),
    color = :black, linewidth=0.85)
translate!(cl, 0, 0, 1000)

cb = Colorbar(f[2,1], h, vertical=false)
f

# why are there artefacts that don't appear when plotting without circshift, from 0 to 360?
# if I plot the map in two steps, OK:
f1 = Figure();
# gax = GeoAxis(f[1,1], 
#         # limits=(extrema(lon), extrema(lat)), 
#         source="+proj=latlong +datum=WGS84", # src CRS
#         dest="+proj=eqearth", # destination CRS, in which you want to plot
#         # coastlines = true # plot coastlines from Natural Earth, as a reference.
#     )
ax1 = Axis(f1[1,1])
data = convert(Array{Float64},aml[year=Near(2020)]);
replace!(data, 0 => NaN);
lon = lookup(aml, :longitude);
lat = lookup(aml, :latitude);
h1 = heatmap!(ax1,lon[1:Int(180/0.25)],lat,data[1:Int(180/0.25),:]);
h2 = heatmap!(ax1,lon[Int(180/0.25+1):end].-360, lat, data[Int(180/0.25+1):end,:]);
# coastlines
cl=lines!(ax1, 
    GeoMakie.coastlines(),
    color = :black, linewidth=0.85)
translate!(cl, 0, 0, 1000)

cb = Colorbar(f1[2,1], h1, vertical=false)
f1