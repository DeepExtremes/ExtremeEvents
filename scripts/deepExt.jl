using YAXArrays, Zarr, EarthDataLab
include("../src/plots.jl")

# make plots for DeepExtremes review meeting
# Pass over Threshold
x = sub_tmax[lon=10,lat=51][:];
m = (markershape = :circle,
    markersize = 2,
    markeralpha = 0.6,
    #markercolor = :green,
    #markerstrokewidth = 0,
    #markerstrokealpha = 0.2,
    #markerstrokecolor = :black,
    markerstrokestyle = :none
    )
p = plot(sub_tmax.time,x, seriestype = :scatter, color = :grey, markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
    xlabel = "Days", ylabel = "Tmax", label = "Events",
    legend = :bottomright)
savefig(p,"../tmax_series.png")
q = quantile(x, 0.99)
ind = x .>= q;
sum(ind)
plot!(p,sub_tmax.time[[1,end]], [q,q], color= Colors.colorant"#28828F", lw=3, label = "Threshold")
savefig(p,"../tmax_tres.png")
plot!(p,sub_tmax.time[ind], x[ind], seriestype = :scatter, 
    color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
    label = "Extremes")
savefig(p,"../tmax_extremes.png")


## connected components
# eec = open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/EventCube_0.01.zarr", consolidated = true))
eec = open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr", consolidated = true))

# eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[:,:,:] .> 0;
# Plots.heatmap(eec_bin[182,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), title = "Boolean layer")
# r = label_components(eec_bin);
#  # should be 70 labels
# findmax(r)
# Plots.heatmap(r[197,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [1, maximum(r)], title = "Connected components")
# # all same event over 1 year...

# # plot events
# p0 = simpleplot(subsetcube(eec, time = 2019:2019, region="Germany"), 182, 2019, 4)
# savefig(p0, "../events.png")
# # if selection on 1 day
# eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[182,:,:] .> 0;
# p1 = Plots.heatmap(eec_bin'[end:-1:1,:], c = cgrad(:gist_gray, categorical = true), title = "Boolean layer")
# savefig(p1,"../boolean_layer.png")
# r = label_components(eec_bin);
# p2 = Plots.heatmap(r'[end:-1:1,:], c = cgrad(:gist_earth, categorical = true), title = "Connected components")
# savefig(p2,"../connected_comp.png")


# plot events
cols = [colorant"#FFFFFF",  
    colorant"#FFB86F", # 1 Light Orange
    colorant"#A6C5E8", # 2 Light Blue
    colorant"#A386BB", # 3 Medium Purple (Light)
    colorant"#4C7FB8", # 4 Medium Blue
    colorant"#8464A5", # 5 Medium Purple (Medium)
    colorant"#4C7FB8", # 6 Medium Blue
    colorant"#8464A5", # 7 Medium Purple (Medium)
    colorant"#002D5A", # 8 Dark Blue
    colorant"#65498C", # 9 Medium Purple (Dark)
    colorant"#002D5A", # 10 Dark Blue
    colorant"#65498C", # 11 Medium Purple (Dark)
    colorant"#002D5A", # 12 Dark Blue
    colorant"#65498C", # 13 Medium Purple (Dark)
    colorant"#002D5A", # 14 Dark Blue
    colorant"#65498C", # 15 Medium Purple (Dark)
    colorant"#BBBBBB", # 16
    ] 
year = 2019; d = 180;

function getsubcube(cube, period, lat, lon)
    subcube = cube.layer[time=period[1]..period[2], latitude=lat[1]..lat[2], longitude=lon[1]..lon[2]]
    subcubeData = subcube.data[:,:,:]    
    if any(lon.>180)
        # modify axes
        axs = modaxs(subcube.axes)
        # modify sublabels1: shift lon
        shifts = getshifts(axs)
        subcubeData = circshift(subcubeData, shifts)
    else
        axs = subcube.axes
    end
    return subcubeData, axs
end

# bbox
lat = (34.0, 60);
lon = [(350.0, 360.0)
 (0, 25.0)];
 date = Date("2019-06-30")
 period = (date, date);#Date("2019-06-30")

subcube1, axsp1 = getsubcube(eec, period, lat, lon[1])
subcube2, axsp2 = getsubcube(eec, period, lat, lon[2])
subcube = cat(subcube1, subcube2, dims = 2)
axs = (axsp1[1], Dim{:longitude}(vcat(axsp1[2].val, axsp2[2][1]:0.25:axsp2[2][end])), axsp1[3])    

# eec_1 = subsetcube(eec, time = year:year, lon = (0,25), lat = (35, 60));#region="Europe");
# eec_data = eec_1.data[d,:,:]'[end:-1:1,:];
# p0 = Plots.heatmap(eec_1.axes[2], eec_1.axes[3][end:-1:1], eec_data,
# Makie plot
fig = Figure(size = (800, 500));

axt = GeoAxis(fig[1,1],
                # xlabel = "Longitude", 
                # ylabel =  "Latitude",
                # aspect = AxisAspect(1), 
                title ="Subset of EventCube on " * string(Date("2019-06-30")) , 
                titlesize=14 ,
                );
limits!(axt, extrema(lookup(axs, :longitude)), lat,)
cl=lines!(axt, 
    GeoMakie.coastlines(),
    color = :grey20, linewidth=0.85)
translate!(cl, 0, 0, 1000)
h = heatmap!(axt, lookup(axs, :longitude), lookup(axs, :latitude), dropdims(subcube; dims = 1),
     colormap = Makie.Categorical(cols), 
     zlims = (0,16),
    #  colorbar_ticks = (0:16, ["no event - ", "only hot", "only dry (30d)", "dry and hot", "only dry (90d)", "dry and hot", "dry", "dry and hot", "only dry (180d)", "dry and hot", "dry", "dry and hot", "dry", "dry and hot", "dry", "dry and hot", "no event - "])
     )
axt.xgridcolor[] = colorant"transparent";
axt.ygridcolor[] = colorant"transparent";
axt.xticklabelsvisible = false;
axt.yticklabelsvisible = false;
# colorbar
cbar = Colorbar(fig[1,2], 
                label = "Event type",
                colormap = cgrad(cols, categorical=true),
                size = 14,
                limits = (-0.5,16.5),
            )
cbar.ticks = (0:16, [
    "no event : rank > 0.1 and rank < 0.9",
    "only hot",
    "only dry (30d)",
    "dry (30d) and hot",
    "only dry (90d)", 
    "dry (90d) and hot", 
    "dry (30d and 90d)", 
    "dry (30d and 90d) and hot", 
    "only dry (180d)", 
    "dry (180d) and hot", 
    "dry (30d and 180d)", 
    "dry (30d and 180d) and hot", 
    "dry (90d and 180d)", 
    "dry (90d and 180d) and hot", 
    "dry (30d, 90d and 180d)", 
    "dry (30d, 90d and 180d) and hot", 
    "no event : rank < 0.1 and rank > 0.9"])

save("../DeepExtremesOutput/v2/fig/Europe_events_" * string(Date("$year") + Day(d)) * ".png", fig)

# # if selection on 1 day
# eec_bin = eec_1.data[d,:,:]'[end:-1:1,:] .> 0 & ;
# p1 = Plots.heatmap(eec_1.axes[2], eec_1.axes[3][end:-1:1], eec_bin, c = cgrad(:gist_gray, categorical = true), title = "Boolean layer")
# savefig(p1,"../boolean_layer.png")
# r = label_components(eec_bin);
# p2 = Plots.heatmap(r'[end:-1:1,:], c = cgrad(:gist_earth, categorical = true), title = "Connected components")
# savefig(p2,"../connected_comp.png")


