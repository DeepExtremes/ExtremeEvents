using YAXArrays, Zarr
# using Plots
using CairoMakie
import Statistics
import CairoMakie.Colors.RGB

eraold = open_dataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr")
pet = eraold.pet 
told = eraold.t2mmax
eranew = open_dataset("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/ERA5Cube.zarr/")
t = eranew.t2mmax
tp = eranew.tp
rtold = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/tmax_ranked.zarr")
rtnew = Cube("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/tmax_ranked.zarr/")
peiold = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PEICube.zarr")
peinew = Cube("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/PEICube.zarr/")
rpeiold = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/pei_ranks.zarr")
rpeinew = Cube("/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/pei_ranks.zarr/")


# Jenalat = 50.92; Jenalon = 11.59
# Niameylon = 2.1254; Niameylat = 13.5116
Jenalat =  13.5; Jenalon = 2.1254

path2fig = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/fig"

st = t[lon = Near(Jenalon), lat = Near(Jenalat)]
stold = told[lon = Near(Jenalon), lat = Near(Jenalat)]
srtnew = rtnew[lon = Near(Jenalon), lat = Near(Jenalat)]
srtold = rtold[lon = Near(Jenalon), lat = Near(Jenalat)]

# setticks
function setticks(x::Any, n::Int, f::Any; ticksvalue = true)
    ticks = [1] #; ticks[n] = x[end];
    nx = length(x)
    for i in 2:n
        push!(ticks, convert(Int, round(nx/(n-1) * (i-1))))
    end
    ticklabels = map(f, x[ticks])
    if ticksvalue
        return((x[ticks], ticklabels))
    else
        return((ticks, ticklabels))
    end
end
# # Pass over Threshold
# ttt = setticks(st.Ti,5,Date)
# p = plot(st.time, x, seriestype = :scatter, 
#     color = :grey, markershape = :circle, markersize = 2, markeralpha = 0.6, markerstrokewidth = 0,
#     xlabel = "Time", ylabel = "Maximum daily temperature", label = "Occurrences",
#     # xlims = (st.time[1], st.time[end]),
#     xticks = ttt,# ([ st.time[1],  st.time[end]], [ Date(st.time[1]),  Date(st.time[end])]),
#     # xtickfonthalign = :right,
#     legend = :bottomright)
# savefig(p,"$(path2fig)/tmax_series.png")
# q = Statistics.quantile(x, 0.99)
# ind = x .>= q;
# sum(ind)
# plot!(p,st.time[[1,end]], [q,q], color= Colors.colorant"#28828F", lw=3, label = "Threshold")
# savefig(p,"$(path2fig)/tmax_tres.png")
# plot!(p,st.time[ind], x[ind], seriestype = :scatter, 
#     color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
#     label = "Extremes")
# savefig(p,"$(path2fig)/tmax_extremes.png")

# distribution
f = Figure(size=(800,800));
ax = Axis(f[1,1], 
    xlabel = "ranked",
    ylabel = "values",
    )
scatter!(ax,srtold.data[:], stold.data[1:length(srtold)],
    label = "old",
    color = :grey, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax,srtnew.data[:], st.data[:],
    label = "new",
    color = :blue, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
axislegend(ax, framevisible = true)
ax1 = Axis(f[2,1], 
    xlabel = "ranked",
    ylabel = "values",)
ind = map(x -> x.<=0.01, srtold.data[:])
scatter!(ax1,srtnew.data[1:length(srtold)][ind], st.data[1:length(srtold)][ind],
    label = "new as old",
    color = :red, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
indn = map(x -> x.<=0.01, srtnew.data[:])
scatter!(ax1,srtold.data[indn[1:length(srtold)]], stold.data[1:length(srtold)][indn[1:length(srtold)]],
    label = "old as new",
    color = :orange, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax1, srtold.data[ind], stold.data[1:length(srtold)][ind],
    label = "old",
    color = :grey, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax1,srtnew.data[indn], st.data[indn],
    label = "new",
    color = :blue, 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
axislegend(ax1, framevisible = true)
f
save("$(path2fig)/scatter_tmax_vs_ranked.png",f)

# scatter old versus new
f = Figure();
ax = Axis(f[1,1], xlabel = "old", ylabel = "new")
scatter!(ax, stold.data[:], st.data[:], markersize = 2)
f
save("$(path2fig)/scatter_tmax_old_vs_new.png",f)

# do the same with PEIs...
speinew = peinew[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)]
speiold = peiold[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)]
srpeinew = rpeinew[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)]
srpeiold = rpeiold[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)]
# scatter old versus new
f = Figure();
for (i,variable) in zip(1:3,lookup(speinew, :Variable,))
    ax = Axis(f[1,i], title = variable, xlabel = "old", ylabel = "new")
    scatter!(ax, speiold[Variable = At(variable)].data[:], speinew[Variable = At(variable)].data[:], markersize = 2)
end
f
save("$(path2fig)/scatter_pei_old_vs_new.png",f)

cols = Makie.wong_colors()
f = Figure(size=(800,600));
for (i,variable) in zip(1:3,lookup(speinew, :Variable,))
ax = Axis(f[1,i], 
    xlabel = "ranked",
    ylabel = "values",
    )
scatter!(ax,srpeiold[Variable = At(variable)].data[:], speiold[Variable = At(variable)].data[1:size(srpeiold,1)],
    label = "old",
    color = cols[1], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax,srpeinew[Variable = At(variable)].data[:], speinew[Variable = At(variable)].data[:],
    label = "new",
    color = cols[2], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
# axislegend(ax, framevisible = true)
ax1 = Axis(f[2,i], 
    xlabel = "ranked",
    ylabel = "values",
    title = variable)
ind = map(x -> x.<=0.01, srpeiold[Variable = At(variable)].data[:])
indn = map(x -> x.<=0.01, srpeinew[Variable = At(variable)].data[:])
scatter!(ax1,srpeinew.data[1:size(srpeiold,1)][ind], speinew[Variable = At(variable)].data[1:size(srpeiold,1)][ind],
    label = "new as old",
    color = cols[3], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax1,srpeiold[Variable = At(variable)].data[indn[1:size(srpeiold,1)]], speiold[Variable = At(variable)].data[1:size(srpeiold,1)][indn[1:size(srpeiold,1)]],
    label = "old as new",
    color = cols[4], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax1, srpeiold[Variable = At(variable)].data[ind], speiold[Variable = At(variable)].data[1:size(srtold,1)][ind],
    label = "old",
    color = cols[1], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
scatter!(ax1,srpeinew[Variable = At(variable)].data[indn], speinew[Variable = At(variable)].data[indn],
    label = "new",
    color = cols[2], 
    marker = :circle, markersize = 2, alpha = 0.6, strokewidth = 0,
    )
# axislegend(ax1, framevisible = true)
end
Legend(f[3,:], 
    [
        MarkerElement(color = cols[1], marker = :circle, markersize = 5,),
        MarkerElement(color = cols[2], marker = :circle, markersize = 5,),
        MarkerElement(color = cols[3], marker = :circle, markersize = 5,),
        MarkerElement(color = cols[4], marker = :circle, markersize = 5,),
    ],
    ["old", "new", "new as old", "old as new"],
    orientation=:horizontal,
)
f
save("$(path2fig)/scatter_pei_vs_ranked.png",f)
###########

# plot!(p1,[0.01,0.01], Statistics.quantile(x,[0,1]), color= Colors.colorant"#28828F", lw=1, label = "Threshold")
# savefig(p,"$(path2fig)/tmax_ranked_tres.png")

# plot!(p1,srt[ind], x[ind], seriestype = :scatter, 
#     color = Colors.colorant"#366570", markersize = 3, markeralpha = 0.6, markerstrokewidth = 0,
#     label = "Extremes")
# savefig(p,"$(path2fig)/tmax_extremes.png")


# modified from wong_palette
my_palette = [
    RGB(([230, 159,   0] / 255)...), # orange
    RGB(([ 86, 180, 233] / 255)...), # sky blue
    RGB(([  0, 158, 115] / 255)...), # blueish green
    RGB(([213,  94,   0] / 255)...), # vermillion
    RGB(([  0, 114, 178] / 255)...), # blue
    # RGB(([204, 121, 167] / 255)...), # reddish purple
    ];

# Makie 
function mtspl!(fig, stp, spet, stmx, spei; i = 1)
ax1 = Axis(fig[1,i],
    xlabel = "Time [day]",
    ylabel = "Water flux [mm/day]",
    ); #\nTemperature [°C]
ax2 = Axis(fig[1,i], yaxisposition = :right,
    ylabel = "\nTemperature [°C]", #"PEI [mm/day]", 
    flip_ylabel = true,
    )
hidespines!(ax2)
hidexdecorations!(ax2)
linkxaxes!(ax1, ax2)
tempo = lookup(stp, :Ti); nt = length(tempo)
btp = barplot!(ax1, 1:nt, stp.data[:] .*1e3, 
    label = "Total precipitation", 
    # strokewidth=0, 
    # linealpha = 1.0,
    color = my_palette[1],
    );
# pet is Float32 and tp is Float64 # is this why it won't be plotted?
bpe = barplot!(ax1, 1:nt, map(x-> convert(Float64, x), spet.data[:]), 
    label = "Ref. evapotranspiration",
    # strokewidth=0,
    # linealpha = 1.0,
    color = my_palette[2],
    );
st = scatter!(ax2,
    1:nt, stmx.data[:] .- 273.15, #yguideposition = :right, 
    label = "Max. temperature at 2m",
    markersize = 2, alpha = 0.3,
    color = my_palette[3]);
#  with pei, also problem of conversion to Float64
lpe30 = lines!(ax1, 1:nt, map(x-> convert(Float64, x), spei.data[:,3]), linewidth = 0.5, label = "pei_30", color = my_palette[5],);
lpe90 = lines!(ax1, 1:nt, map(x-> convert(Float64, x), spei.data[:,2]), linewidth = 0.5, label = "pei_90", color = my_palette[4],);
lpe180 = lines!(ax1, 1:nt, map(x-> convert(Float64, x), spei.data[:,1]), linewidth = 0.5, label = "pei_180", color = my_palette[3],);

# ticks
ax1.xticks = ax2.xticks = setticks(tempo, 8, x -> string(Date(x)); ticksvalue=false)
ax1.xticklabelrotation = π / 4

plots = [btp, bpe, st, lpe30, lpe90, lpe180]
return fig, plots
end

# plot pet and tp for 10 years
startyear = 2013; stopyear = 2022
jstp = tp[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
jspet = pet[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
jstmx = era.t2mmax[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
jspei = pei[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
jqt = Statistics.quantile(era.t2mmax[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)][:], 0.99)
jqt - 273.15
jqpei = map(x -> Statistics.quantile(skipmissing(x), 0.01), eachslice(pei[lon = At(Jenalon, atol=0.25), lat = At(Jenalat, atol=0.25)], dims = 2))
# ind = x .>= q;

fig = Figure(size = (1000,500));
fig, plots = mtspl!(fig,jstp, jspet, jstmx, jspei);

# Niamey 13.5116° N, 2.1254° E
Niameylon = 2.1254; Niameylat = 13.5116
nqt = Statistics.quantile(era.t2mmax[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25)][:], 0.99)
nqt - 273.15
nqpei = map(x -> Statistics.quantile(skipmissing(x), 0.01), eachslice(pei[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25)], dims = 2))
nstp = tp[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
nspet = pet[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
nstmx = era.t2mmax[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]
nspei = pei[lon = At(Niameylon, atol=0.25), lat = At(Niameylat, atol=0.25), time = Date(startyear) .. Date(stopyear, 12, 31)]

fig, plots = mtspl!(fig,nstp, nspet, nstmx, nspei; i=2);

# add legend
fig[2,1:2] = Legend(fig,
        plots,
        map(x -> x.label, plots), 
        position = :lb, 
        orientation = :horizontal,
        nbanks = 3,
        framevisible = false,
        )
linkyaxes!(fig.content[1], fig.content[3])
linkyaxes!(fig.content[2], fig.content[4])
fig.content[2].ylabelvisible = false
fig.content[3].ylabelvisible = false
# ylims!(fig.content[1], (-5,20))
Label(fig[1,1, Top()], 
        text = "(a) Jena, Germany",
        halign = :left
        )
Label(fig[1,2, Top()], 
        text = "(b) Niamey, Niger",
        halign = :left
        )
fig
save("$(path2fig)/timeseries_$(startyear)_$(stopyear)_1.png", fig)

# to do: split figure into 3 subfig, with one variable in each
# split plot into subplots
function mpl!(fig, stp, spet, stmx, spei, qt, qpei; i = 1)
ax1 = Axis(fig[1,i], yaxisposition = i==1 ? :left : :right,
    ylabel = "Temperature [°C]", #"PEI [mm/day]", 
    flip_ylabel = i==1 ? false : true,
    )
ax2 = Axis(fig[2,i],
    ylabel = "Water flux [ l m^{-2} d^{-1}]",
    yaxisposition = i==1 ? :left : :right,
    flip_ylabel = i==1 ? false : true,
    ); #\nTemperature [°C]
ax3 = Axis(fig[3,i],
    xlabel = "Time [d]",
    ylabel = L"PEI [l m^{-2} d^{-1}]",
    yaxisposition = i==1 ? :left : :right,
    flip_ylabel = i==1 ? false : true,
)
# hidespines!(ax2)
# hidexdecorations!(ax2)
linkxaxes!(ax1, ax2, ax3)
tempo = lookup(stp, :Ti); nt = length(tempo)
btp = barplot!(ax2, 1:nt, stp.data[:] .*1e3, 
    label = "Total precipitation", 
    # strokewidth=0, 
    # linealpha = 1.0,
    color = my_palette[1],
    );
# pet is Float32 and tp is Float64 # is this why it won't be plotted?
bpe = barplot!(ax2, 1:nt, map(x-> convert(Float64, x), spet.data[:]), 
    label = "Ref. evapotranspiration",
    # strokewidth=0,
    # linealpha = 1.0,
    color = my_palette[2],
    );
st = scatter!(ax1,
    1:nt, stmx.data[:] .- 273.15, #yguideposition = :right, 
    label = "Max. temperature at 2m",
    markersize = 2, alpha = 0.3,
    color = my_palette[3]);
hlt = hlines!(ax1, qt-273.15, linewidth = 1, linestyle = :dash, label = "threshold T\${2m,max}\$", color = my_palette[3],);
#  with pei, also problem of conversion to Float64
lpe30 = lines!(ax3, 1:nt, map(x-> convert(Float64, x), spei.data[:,3]), linewidth = 1, label = "pei_30", color = my_palette[5],);
hlpe30 = hlines!(ax3, qpei.data[3], linewidth = 1, linestyle = :dash, label = "threshold PEI_30", color = my_palette[5],);
lpe90 = lines!(ax3, 1:nt, map(x-> convert(Float64, x), spei.data[:,2]), linewidth = 1, label = "pei_90", color = my_palette[4],);
hlpe90 = hlines!(ax3, qpei.data[2], linewidth = 1, linestyle = :dash, label = "threshold PEI_90", color = my_palette[4],);
lpe180 = lines!(ax3, 1:nt, map(x-> convert(Float64, x), spei.data[:,1]), linewidth = 1, label = "pei_180", color = my_palette[3],);
hlpe180 = hlines!(ax3, qpei.data[1], linewidth = 1, linestyle = :dash, label = "threshold PEI_180", color = my_palette[3],);

# ticks
ax1.xticks = ax2.xticks = ax3.xticks = setticks(tempo, 5, x -> string(Date(x)); ticksvalue=false)
# ax1.xticklabelrotation = π / 4

plots = [btp, bpe, st, hlt, lpe30, hlpe30, lpe90, hlpe90, lpe180, hlpe180]
return fig, plots
end

f = with_theme(theme_latexfonts()) do
    f = Figure();
end
f, plots = mpl!(f, jstp, jspet, jstmx, jspei, jqt, jqpei)