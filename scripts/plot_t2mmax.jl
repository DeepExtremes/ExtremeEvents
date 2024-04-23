using YAXArrays, Zarr
using Plots
import Statistics

zg = zopen("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
t = era.t2mmax
tp = era.tp
pet = era.pet 
rt = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/tmax_ranked.zarr")
pei = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr")

Jenalat = 50.92; Jenalon = 11.59

path2fig = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/fig"

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
savefig(p,"$(path2fig)/tmax_series.png")
q = Statistics.quantile(x, 0.99)
ind = x .>= q;
sum(ind)
plot!(p,st.time[[1,end]], [q,q], color= Colors.colorant"#28828F", lw=3, label = "Threshold")
savefig(p,"$(path2fig)/tmax_tres.png")
plot!(p,st.time[ind], x[ind], seriestype = :scatter, 
    color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
    label = "Extremes")
savefig(p,"$(path2fig)/tmax_extremes.png")

# distribution
p1 = Plots.scatter(srt[:], x,
    label = "t2mmax",
    xlabel = "ranked",
    ylabel = "values",
    color = :grey, markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
    )
savefig(p1, "$(path2fig)/scatter_tmax_vs_ranked.png")

plot!(p1,[0.01,0.01], Statistics.quantile(x,[0,1]), color= Colors.colorant"#28828F", lw=1, label = "Threshold")
savefig(p,"$(path2fig)/tmax_ranked_tres.png")

plot!(p1,srt[ind], x[ind], seriestype = :scatter, 
    color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
    label = "Extremes")
savefig(p,"$(path2fig)/tmax_extremes.png")

# plot pet and tp for 7 years
startyear = 2012; stopyear = 2023
stp = tp[lon = Jenalon, lat = Jenalat, time = (startyear,stopyear)]
spet = pet[lon = Jenalon, lat = Jenalat, time = (startyear,stopyear)]
stmx = era.t2mmax[lon = Jenalon, lat = Jenalat, time = (startyear,stopyear)]
spei = pei[lon = Jenalon, lat = Jenalat, time = (startyear,stopyear)]
qt = Statistics.quantile(era.t2mmax[lon = Jenalon, lat = Jenalat][:], 0.99)
qt - 273.15
qpei = mapslices(x -> Statistics.quantile(skipmissing(x), 0.01), pei[lon = Jenalon, lat = Jenalat][:,:], dims = 1)
# ind = x .>= q;

include("mytheme.jl");
theme(:mytheme)
function tspl(stp, spet, stmx, spei)
p2 = bar(stp.time, stp[:].*1e3, 
    label = "Total precipitation", 
    linewidth=0, 
    linealpha = 1.0,
    );
bar!(spet.time, spet[:], label = "Ref. evapotranspiration",
    xlabel = "Time [day]",
    ylabel = "Water flux [mm/day] \n Temperature [°C]",
    linewidth=0,
    linealpha = 1.0,
    );
scatter!(
    # twinx(), 
    stmx.time, stmx[:] .- 273.15, #yguideposition = :right, 
    # yaxis = "Temperature [Kelvin]", 
    label = "Max. temperature at 2m",
    markersize = 2, markeralpha = 0.3);
plot!(spei.time, spei.data[:,3], label = "pei_30",);
plot!(spei.time, spei.data[:,2], label = "pei_90",);
plot!(spei.time, spei.data[:,1], label = "pei_180",);
plot!(
    dpi = 300, 
    size = (800, 600),
    legend = :outerbottom,
    legend_column = 3,
    xlims = (spei.time[1], spei.time[end]),
    );
return p2
end
p2 = tspl(stp, spet, stmx, spei)
savefig(p2,"$(path2fig)/timeseries_$(Jenalon)_$(Jenalat)_$(startyear)_$(stopyear).png", )

# Niamey 13.5116° N, 2.1254° E
Niameylon = 2.1254; Niameylat = 13.5116
qtN = Statistics.quantile(era.t2mmax[lon = Niameylon, lat = Niameylat][:], 0.99)
qtN - 273.15
qpeiN = mapslices(x -> Statistics.quantile(skipmissing(x), 0.01), pei[lon = Niameylon, lat = Niameylat][:,:], dims = 1)
stp = tp[lon = Niameylon, lat = Niameylat, time = (startyear,stopyear)]
spet = pet[lon = Niameylon, lat = Niameylat, time = (startyear,stopyear)]
stmx = era.t2mmax[lon = Niameylon, lat = Niameylat, time = (startyear,stopyear)]
spei = pei[lon = Niameylon, lat = Niameylat, time = (startyear,stopyear)]
p3 = tspl(stp, spet, stmx, spei)
savefig(p3,"$(path2fig)/timeseries_$(Niameylon)_$(Niameylat)_$(startyear)_$(stopyear).png", )
