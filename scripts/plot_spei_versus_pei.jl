# Script to plot comparisons of PEI and SPEI on non-deseasonalized daily data
# params of the SPEI have been precomputed by Fabian.
# quantiles of PEI_30 have been precomputed by Fabian.
using CairoMakie, YAXArrays, Zarr, NetCDF
using FHist
import GeoMakie
thresholds = open_dataset("/Net/Groups/BGI/work_2/scratch/fgans/pei_30_quantiles.zarr")
params = open_dataset("/Net/Groups/BGI/work_2/scratch/fgans/SPEI_30_params.zarr/")
path = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"
data = open_dataset(joinpath(path, "PEICube.zarr")).pei_30
lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
isoc = lsm.lsm.data[:,:,1] .< 0.5
tres = readcubedata(thresholds.layer[quantile=1]);

tres.data[repeat(vec(isoc), inner=2)] .= NaN;

# when γ is smaller than x, can't plot pdf (NaN)?
par = readcubedata(params.layer[param=1:3]);
par.data[repeat(vec(isoc), inner=3)] .= NaN;

function hmap(cube; f = Figure(), i=1, title="", clabel="", kwargs...)
    ax = GeoMakie.GeoAxis(f[i,1], title=title) 
    # GeoMakie.GeoAxis brings artefacts when doing circ shift, hence remove lon=180
    lon = circshift(map(x -> x >= 180 ? x-360 : x, YAXArrays.lookup(cube, :longitude)), (180/0.25));
    lat = YAXArrays.lookup(cube, :latitude);
    h1 = heatmap!(ax, lon[2:720], lat, cube[longitude=180.1..360].data; kwargs...)
    h2 = heatmap!(ax, lon[721:end], lat, cube[longitude=0..179.9].data; kwargs...)
    
    cb = Colorbar(f[i,2], h1, label=clabel)
    # remove gridlines
    ax.xgridcolor[] = colorant"transparent";
    ax.ygridcolor[] = colorant"transparent";
    ax.xticklabelsvisible = false;
    ax.yticklabelsvisible = false;

    f
end
function hmap!(f::Figure, cube; i=1, kwargs...)
    hmap(cube; f = f, i=i, kwargs...)

    Label(f[i,1, TopLeft()], 
        text = "($('`'+i))",
        halign = :left
        )
    
    return f
end

f = Figure(size=(600,850));
# Heatmap for a) fitted distribution; b) empirical
for i in 1:2
    hmap!(f, tres[method=i], i=i, title="$(tres.method[i]), q = $(thresholds.quantile[1])", colorrange=(-6.0,2.0), clabel = "PEI30 (mm/day)")
end
# Heatmap of differences
hmap!(f, tres[method=1].-tres[method=2], i=3, title = "Difference: $(tres.method[1]) - $(tres.method[2])", colorrange=(-3,3),colormap=:bam100, clabel = "(mm/day)")
save(joinpath(path,"fig/pei_vs_spei.png"), f)

f = Figure(size=(600,850));
# Heatmap for parameters of fitted distribution
for i in 1:3
    hmap!(f, par[param=i], i=i, colorrange=(-50.0,50.0), clabel = "$(par.param[i])", colormap=:bam)
end
f
save(joinpath(path,"fig/spei_30_param.png"), f)

function make_hist_plot!(fig,lons,lats,params,data,thresholds;quantile=1,i=1,j=1, title="")
    p = params.layer[lon=Near(lons),lat=Near(lats)]
    @show sp = SPEI(p[1:3]...)
    d = data[lon=Near(lons),lat=Near(lats)][Ti=Date(1984)..Date(2024)].data[:]
    t_dist, t_quant = thresholds.layer[lon=Near(lons),lat=Near(lats), quantile=quantile].data
    mi,ma = extrema(d)
    mi = min(mi,t_dist-3.0)
    binedges = range(mi,ma,length=100)
    h = Hist1D(d; binedges) |> normalize
    
    xr = bincenters(h)
    dens = pdf.(xr,(sp,));
    i_dist,i_quant = searchsortedfirst(xr,t_dist),searchsortedfirst(xr,t_quant)
    ax = Axis(fig[i,j], title = title, titlealign = :left)
    lines!(ax,xr,h.bincounts,color=:blue,linewidth=2);
    lines!(ax,xr,dens,color=:orange,linewidth=2)
    band!(ax,xr[1:i_quant],zeros(i_quant),h.bincounts[1:i_quant],color=(:blue,0.2))
    band!(ax,xr[1:i_dist],zeros(i_dist),dens[1:i_dist],color=(:orange,0.2))
    lines!(ax,[xr[i_quant],xr[i_quant]],[0,h.bincounts[i_quant]],color=:blue,linewidth=2)
    lines!(ax,[xr[i_dist],xr[i_dist]],[0,dens[i_dist]],color=:orange,linewidth=2)
    # vlines!(ax,t_quant)
    # vlines!(ax,t_dist)

    Label(f[i,j, TopLeft()], 
        text = "($('`'+((i-1)*2+j)))",
        halign = :left
        )
    fig
end

struct SPEI
    α::Float64
    β::Float64
    γ::Float64
end
function (spei::SPEI)(x)
    x < spei.γ && return NaN
    spei.β < 0.0 && return NaN
    spei.α < 0.0 && return NaN
    any(isnan,(spei.α, spei.β, spei.γ)) && return NaN
    F = inv(1 + (spei.α / (x-spei.γ))^spei.β)
    P = 1-F
    C0 = 2.515517; C1 = 0.802853; C2 = 0.010328; d1 = 1.432788; d2 = 0.189269; d3 = 0.001308
    if P<=0.5
        W = sqrt(-2*log(P))
        return W - (C0 + C1*W + C2*W*W)/(1 + d1*W + d2*W*W + d3*W*W*W)
    else
        W = sqrt(-2*log(1-P))
        return -W  + (C0 + C1*W + C2*W*W)/(1 + d1*W + d2*W*W + d3*W*W*W)
    end
end

function pdf(x,sp::SPEI)
    # x < sp.γ && return NaN
    # sp.β < 0.0 && return NaN
    # sp.α < 0.0 && return NaN
    # any(isnan,(sp.α, sp.β, sp.γ)) && return NaN

    # some filter should be adapted to be able to plot. The one above removes too much data.

    try
        sp.β / sp.α * ((x-sp.γ)/sp.α)^(sp.β-1) * (1 + ((x-sp.γ)/sp.α)^sp.β)^(-2)
    catch
        @show x
        NaN
    end
end


# The following 4 examples could be shown in a plot, or even more
f = Figure();
# Example western Europe
Jenalat = 50.92; Jenalon = 11.59
make_hist_plot!(f, Jenalon, Jenalat ,params,data,thresholds, i=1, j=1, title = "Jena, Germany ($Jenalon, $Jenalat)");
# Example dryland
Niameylon = 2.12; Niameylat = 13.51
make_hist_plot!(f, Niameylon, Niameylat ,params,data,thresholds, i=1, j=2, title = "Niamey, Niger ($Niameylon, $Niameylat)");
# Example Australia 
make_hist_plot!(f,118.21,-31.47,params,data,thresholds, i=2, j=2, title = "Merredin, Australia (118.21,-31.47)");
# Example: tropics where threshold for dist is much higher than from quantile
make_hist_plot!(f,360-69,5,params,data,thresholds, i=2, j=1, title="Casanare, Colombia (-69.0,5.0)");
Label(f[3,1:2], text="PEI_30 (mm/d)", halign=:center);
Label(f[1:2,0], text="Density", rotation=π/2);
f
save(joinpath(path, "fig/pei_vs_spei_density.png"), f)