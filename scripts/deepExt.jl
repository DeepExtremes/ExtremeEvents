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
eec = Cube(open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/EventCube_0.01.zarr", consolidated = true)))

eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[:,:,:] .> 0;
Plots.heatmap(eec_bin[182,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), title = "Boolean layer")
r = label_components(eec_bin);
 # should be 70 labels
findmax(r)
Plots.heatmap(r[197,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), zlims = [1, maximum(r)], title = "Connected components")
# all same event over 1 year...

# plot events
p0 = simpleplot(subsetcube(eec, time = 2019:2019, region="Germany"), 182, 2019, 4)
savefig(p0, "../events.png")
# if selection on 1 day
eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[182,:,:] .> 0;
p1 = Plots.heatmap(eec_bin'[end:-1:1,:], c = cgrad(:gist_gray, categorical = true), title = "Boolean layer")
savefig(p1,"../boolean_layer.png")
r = label_components(eec_bin);
p2 = Plots.heatmap(r'[end:-1:1,:], c = cgrad(:gist_earth, categorical = true), title = "Connected components")
savefig(p2,"../connected_comp.png")


