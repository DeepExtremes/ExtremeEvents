using YAXArrays, Zarr
using Dates, DataFrames

path = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"

startyear = 1966
endyear = 2023

deo = Cube(joinpath(path,"EventCube_ranked_pot0.01_ne0.1.zarr"))[Ti=Date(startyear)..Date(endyear,12,31)]
tempo = year.(lookup(deo, :Ti))
# aggregate by decade? only 5 points...
tempo = Int.(round.(tempo, digits=-1))

include("../src/stats.jl")

function get_trend!(xout,xin,tempo::Vector{Int64}; rule::Function=x -> (x .>1) .& (iseven.(x.+1)))
    df = DataFrame(Year=tempo, Type=xin) |>
    (df -> DataFrames.subset(df, :Type => rule)) |>
    (df -> DataFrames.groupby(df, :Year)) |>
    (df -> DataFrames.combine(df, nrow)) |>
    (df -> leftjoin(DataFrame(Year=unique(tempo)), df, on=[:Year => :Year], order=:left)) |>
    (df -> DataFrames.select!(df, :Year, :nrow => x -> replace(x, missing=>0))) |>
    (df -> rename!(df, :nrow_function => "ndays"))

    m,_ = theilsen(df.Year, df.ndays)
    return xout[] = Float32(m)
end

# xout=[Float32(0.0)];
# xin=deo.data[:,10,550];
# dnh = get_trend!(xout,xin,tempo)
# hot = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(1)) .== Int8(1))
# dry30 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(2)) .== Int8(2))
# dry90 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(4)) .== Int8(4))
# dry180 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(8)) .== Int8(8))
# anyex = get_trend!(xout,xin,tempo,rule = x -> (x .> 0 .&& x .< 16))
# noex= get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(16)) .== Int8(16))

function trendmap(cube; title=cube.properties["name"], kwargs...)
    f = Figure()
    ax = Axis(f[1,1], title=title)
    data = circshift(convert(Array{Float32},cube), (180/0.25, 0));
    # replace!(data, 0 => NaN);
    lon = circshift(map(x -> x >= 180 ? x-360 : x, lookup(cube, :longitude)), (180/0.25));
    lat = lookup(cube, :latitude);
    h = heatmap!(ax,lon,lat,data; kwargs...);
    # coastlines
    cl=lines!(ax, 
        GeoMakie.coastlines(),
        color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 1000)

    cb = Colorbar(f[2,1], h, vertical=false)
    f
end


# # probably I can get a clear trend if I put all extremes together, but not the compound
@time trend_dh = mapCube(get_trend!,deo,tempo, indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_dh_decade_1966_2023.zarr"), chunksize=:max, overwrite=true, layername="dh"))
# # 2066.314741 seconds (3.65 G allocations: 480.355 GiB, 5.89% gc time, 0.01% compilation time)
# 190.135836 seconds (815.77 M allocations: 317.858 GiB, 12.77% gc time, 195.58% compilation time: <1% of which was recompilation)
trend_dh=Cube(joinpath(path,"trendmap_dh_decade_1966_2023.zarr"))
f = trendmap(trend_dh, title="Theil-Sen trend in decadal number of extremely dry and hot days (1970-2023)", colorrange=(-1,1), colormap=:vik)
save(joinpath(path,"fig/trendmap_dh_decadal_1966_2023.png"), f)

# @time trend_any = mapCube(get_trend!,deo,tempo,rule = x -> (x .> 0 .&& x .< 16), indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_any_1971_2023.zarr"), chunksize=:max, overwrite=true, layername="any"))
# # 2195.409061 seconds (3.67 G allocations: 498.512 GiB, 2.60% gc time, 0.01% compilation time)
# f = trendmap(trend_any, title="Theil-Sen trend in annual number of extremely dry or hot days (1971-2023)", colorrange=(-.5,.5), colormap=:vik)
# save(joinpath(path,"fig/trendmap_any_1971_2023.png"), f)

# @time mapCube(get_trend!,deo,tempo,rule = x -> (x .& Int8(1)) .== Int8(1), indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_hot_1971_2023.zarr"), chunksize=:max, overwrite=true, layername="hot"))
# # 2072.418393 seconds (3.65 G allocations: 487.068 GiB, 2.73% gc time, 0.01% compilation time)
# @time mapCube(get_trend!,deo,tempo,rule = x -> (x .& Int8(2)) .== Int8(2), indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_dry30_1971_2023.zarr"), chunksize=:max, overwrite=true, layername="dry30"))
# # 2016.432447 seconds (3.65 G allocations: 485.673 GiB, 2.75% gc time, 0.01% compilation time)
@time mapCube(get_trend!,deo,tempo,rule = x -> (x .& Int8(4)) .== Int8(4), indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_dry90_1971_2023.zarr"), chunksize=:max, overwrite=true, layername="dry90"))
# 1976.989004 seconds (3.65 G allocations: 484.997 GiB, 2.88% gc time, 0.01% compilation time)
@time mapCube(get_trend!,deo,tempo,rule = x -> (x .& Int8(8)) .== Int8(8), indims=InDims(:Ti), outdims=OutDims(outtype=Float32,path=joinpath(path,"trendmap_dry180_1971_2023.zarr"), chunksize=:max, overwrite=true, layername="dry180"))
# 1956.747982 seconds (3.65 G allocations: 484.381 GiB, 2.83% gc time, 0.01% compilation time)

f = trendmap(Cube(joinpath(path,"trendmap_hot_1971_2023.zarr")), title="Theil-Sen trend in annual number of extremely hot days (1971 - 2023)", colorrange=(-0.5,0.5), colormap=:vik)
save(joinpath(path,"fig/trendmap_hot_1971_2023.png"), f)
f = trendmap(Cube(joinpath(path,"trendmap_dry30_1971_2023.zarr")), title="Theil-Sen trend in annual number of extremely dry (last 30) days (1971 - 2023)", colorrange=(-0.5,0.5), colormap=:vik)
save(joinpath(path,"fig/trendmap_dry30_1971_2023.png"), f)
f = trendmap(Cube(joinpath(path,"trendmap_dry90_1971_2023.zarr")), title="Theil-Sen trend in annual number of extremely dry (last 90) days (1971 - 2023)", colorrange=(-0.1,0.1), colormap=:vik)
save(joinpath(path,"fig/trendmap_dry90_1971_2023.png"), f)
f = trendmap(Cube(joinpath(path,"trendmap_dry180_1971_2023.zarr")), title="Theil-Sen trend in annual number of extremely dry (last 180) days (1971 - 2023)", colorrange=(-0.1,0.1), colormap=:vik)
save(joinpath(path,"fig/trendmap_dry180_1971_2023.png"), f)

print("done!")