# annual number of events by type and by grid cell (lat, lon)
using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    # addprocs(SlurmManager())
    # delay addprocs
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
        addprocs(1)
        sleep(0.001)
    end
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere begin
    using YAXArrays, OnlineStats, WeightedOnlineStats, Zarr, NetCDF
    using DimensionalData
    using DimensionalData.LookupArrays
    using DataFrames, Dates
    import CSV
    import Statistics #StatsBase
    # import Plots
    # using Measures
    using CairoMakie, GeoMakie

    include("../src/plots.jl")
end

if occursin("/Users", pwd())
    path = "/Users/mweynants/BGI/DeepExtremes/DeepExtremesOutput/"
else
    path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/"
end

trial = "ranked_pot0.01_ne0.1_cmp_S1_T3_2010_2022"
aec = Cube("$(path)AnnualEventCount_ranked_pot0.01_ne0.1.zarr")
aml= Cube("$(path)AnnualMaxLabel_cmp_S1_T3.zarr")

years = lookup(aec, :year, )

cmap = vcat(colorant"grey99", resample_cmap(cgrad(:viridis, 100, categorical = true), 100))
cmap1 = vcat(colorant"grey99", resample_cmap(Reverse(:acton50),50))#cgrad(:acton50, 50, categorical = true))

lon = lookup(aec, :longitude);
lat = lookup(aec, :latitude)[end:-1:1];
axs = modaxs(aec.axes)
ydim = dimnum(axs, :year)
tdim = dimnum(axs, :net)
vaxs = [axs...]
deleteat!(vaxs, [ydim,tdim])
axs = Tuple(x for x in vaxs)
lon_dim =  dimnum(axs, :longitude);
nlon = collect(axs[lon_dim]);
# shift longitude
shifts = getshifts(axs)
# axsl = modaxs(aml.axes)
# shifts_aml = getshifts(axsl)

plotorder = ["d30", "d90", "d180", "h", "dh"]
plottitles = ["dry (30 days)", "dry (90 days)", "dry (180 days)", "hot", "dry and hot"]

pmap(years) do y
f0 = Figure(;size=(1500,600), fontsize = 24);
for i in eachindex(plotorder)
    row = Int((i + 2.5) ÷ 3)
    col = i - 3 * (i ÷ 4)
    # @show (row, col)
    ax = #GeoAxis
        Axis(f0[row,col], title = plottitles[i]);
    s = surface!(ax, nlon, lat, circshift(aec[year = At(y), net = At(plotorder[i])].data[:,end:-1:1], shifts),
        colormap = cmap1,
        colorrange=(0,100),
        highclip=cmap1[end],
        );
    cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :grey50, linewidth=0.85,
        )
    translate!(cl, 0, 0, 1000);
    # # remove gridlines
    # ax.xgridcolor[] = colorant"transparent";
    # ax.ygridcolor[] = colorant"transparent";
    # ax.xticklabelsvisible = false;
    # ax.yticklabelsvisible = false;
    # remove decorations
    hidedecorations!(ax)
end
ax = #GeoAxis
    Axis(f0[2,3], title = "labels")
yml = circshift(aml[year = At(y)].data[:,end:-1:1], shifts)
# replace 0 by NaN
replace!(yml, 0 => missing)
l = surface!(ax, nlon, lat, yml,
        colormap = Makie.Categorical(:glasbey_bw_minc_20_minl_30_n256),
        );
    cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :grey50, linewidth=0.85,
        )
    translate!(cl, 0, 0, 1000);
    # # remove gridlines
    # ax.xgridcolor[] = colorant"transparent";
    # ax.ygridcolor[] = colorant"transparent";
    # ax.xticklabelsvisible = false;
    # ax.yticklabelsvisible = false;
    # remove decorations
    hidedecorations!(ax)


cbar = Colorbar(f0[1,4], colormap = cmap1,
        colorrange=(0,100),
        highclip=cmap1[end],
        label = "Number of dry and/or hot days")
cbarl = Colorbar(f0[2,4], l, label = "Labelled events", 
    ticklabelsvisible = false, ticksvisible = false)
Label(f0[0,:], "Year $y")
# f0
@time save("$(path)fig/annual/ndays_$(trial)_year$(y).png", f0, dpi = 300)
end # for year