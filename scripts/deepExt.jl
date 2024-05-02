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
eec = Cube(open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/deepextremes/v2/EventCube_ranked_pot0.01_ne0.1.zarr/", consolidated = true)))
eec = Cube(open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr"))

eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[:,:,:] .> 0;
Plots.heatmap(eec_bin[182,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), title = "Boolean layer")
r = label_components(eec_bin);
 # should be 70 labels
findmax(r)
Plots.heatmap(r[197,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [1, maximum(r)], title = "Connected components")
# all same event over 1 year...

# plot events
cols = [:white,  
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
eec_1 = subsetcube(eec, time = year:year, region="Europe");
p0 = Plots.heatmap(eec_1.axes[2], eec_1.axes[3][end:-1:1], eec_1.data[d,:,:]'[end:-1:1,:],
     c = cgrad(cols, categorical = true), 
     title = "Subset of Event Cube on " * string(Date("$year") + Day(d)), 
     ylabel="latitude", xlabel="longitude")
# p0 = simpleplot(subsetcube(eec, time = 2019:2019, region="Germany"), 182, 2019, 4)
savefig(p0, "../Europe_events_" * string(Date("$year") + Day(d)) * ".png")
# if selection on 1 day
eec_bin = eec_1 .> 0;
p1 = Plots.heatmap(eec_bin'[end:-1:1,:], c = cgrad(:gist_gray, categorical = true), title = "Boolean layer")
savefig(p1,"../boolean_layer.png")
r = label_components(eec_bin);
p2 = Plots.heatmap(r'[end:-1:1,:], c = cgrad(:gist_earth, categorical = true), title = "Connected components")
savefig(p2,"../connected_comp.png")


