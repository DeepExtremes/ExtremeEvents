using YAXArrays, Zarr
using Plots
import Statistics

pei = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr")
rpei = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/pei_ranks.zarr")

Jenalat = 50.92; Jenalon = 11.59

spei = pei[lon = Jenalon, lat = Jenalat]
srpei = rpei[lon = Jenalon, lat = Jenalat]

x = spei.data;

color_palette = [ colorant"#A6C5E8", colorant"#4C7FB8", colorant"#002D5A"]

# setticks
function setticks(x::Any, n::Int, f::Any)
    ticks = [x[1]]#; ticks[n] = x[end];
    for i in 2:n
        push!(ticks, x[convert(Int, round(end/(n-1) * (i-1)))])
    end
    ticklabels = map(f, ticks)
    return((ticks, ticklabels))
end
# Pass over Threshold
ttt = setticks(pei.time,5,Date)
p = plot(pei.time, x, seriestype = :scatter, 
    markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
    xlabel = "Time", ylabel = "Daily P-E", #label = "Occurrences",
    # xlims = (st.time[1], st.time[end]),
    xticks = ttt,# ([ st.time[1],  st.time[end]], [ Date(st.time[1]),  Date(st.time[end])]),
    # xtickfonthalign = :right,
    legend = :bottomright)
# savefig(p,"../v3/fig/pei_series.png")
# q = Statistics.quantile(x, 0.99)
# ind = x .>= q;
# sum(ind)
# plot!(p,st.time[[1,end]], [q,q], color= Colors.colorant"#28828F", lw=3, label = "Threshold")
# savefig(p,"../v3/fig/tmax_tres.png")
# plot!(p,st.time[ind], x[ind], seriestype = :scatter, 
#     color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
#     label = "Extremes")
# savefig(p,"../v3/fig/tmax_extremes.png")

# distribution
# p1 = Plots.scatter(srpei.data, x,
#     label = srpei.Variable,
#     xlabel = "ranked",
#     ylabel = "values",
#     color = color_palette, 
#     markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
#     )

# do something like this:
p1 = plot(
    xlabel = "ranked",
    ylabel = "values",
    );
# for (var, color) in zip(["pei_30", "pei_90", "pei_180"], color_palette)
#     plot!(srpei[var].data, rpei[var].data, label = var, c = color,
#     markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,)
# end
# p1

savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/scatter_pei_vs_ranked.png")

# plot!(p1,[0.01,0.01], Statistics.quantile(x,[0,1]), color= Colors.colorant"#28828F", lw=1, label = "Threshold")
# savefig(p,"../v3/fig/tmax_ranked_tres.png")

# plot!(p1,srt[ind], x[ind], seriestype = :scatter, 
#     color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
#     label = "Extremes")
# savefig(p,"../v3/fig/tmax_extremes.png")



