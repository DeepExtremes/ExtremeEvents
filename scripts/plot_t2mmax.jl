using YAXArrays, Zarr
using Plots
import Statistics

zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
t = era.t2mmax
rt = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/tmax_ranked.zarr")

Jenalat = 50.92; Jenalon = 11.59

st = t[lon = Jenalon, lat = Jenalat]
srt = rt[lon = Jenalon, lat = Jenalat]

x = st[:];

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
ttt = setticks(st.time,5,Date)
p = plot(st.time, x, seriestype = :scatter, 
    color = :grey, markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
    xlabel = "Time", ylabel = "Maximum daily temperature", label = "Occurrences",
    # xlims = (st.time[1], st.time[end]),
    xticks = ttt,# ([ st.time[1],  st.time[end]], [ Date(st.time[1]),  Date(st.time[end])]),
    # xtickfonthalign = :right,
    legend = :bottomright)
savefig(p,"../v3/fig/tmax_series.png")
q = Statistics.quantile(x, 0.99)
ind = x .>= q;
sum(ind)
plot!(p,st.time[[1,end]], [q,q], color= Colors.colorant"#28828F", lw=3, label = "Threshold")
savefig(p,"../v3/fig/tmax_tres.png")
plot!(p,st.time[ind], x[ind], seriestype = :scatter, 
    color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
    label = "Extremes")
savefig(p,"../v3/fig/tmax_extremes.png")

# distribution
p1 = Plots.scatter(srt[:], x,
    label = "t2mmax",
    xlabel = "ranked",
    ylabel = "values",
    color = :grey, markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
    )
savefig(p1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig/scatter_tmax_vs_ranked.png")

plot!(p1,[0.01,0.01], Statistics.quantile(x,[0,1]), color= Colors.colorant"#28828F", lw=1, label = "Threshold")
savefig(p,"../v3/fig/tmax_ranked_tres.png")

plot!(p1,srt[ind], x[ind], seriestype = :scatter, 
    color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
    label = "Extremes")
savefig(p,"../v3/fig/tmax_extremes.png")



