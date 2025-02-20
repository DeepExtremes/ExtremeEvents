using YAXArrays, Zarr
using Dates, DataFrames
using CairoMakie
import GeoMakie

path = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"

# pixelwise, there is mostly no significant trend bc some/many years have no extreme days
# Average number of annual extremely dry and hot days by spatial grid cell
# Compare two time periods:
# - one scale for number of days in 1st period (intensity?)
# - one scale for trend

startyear = 1970
endyear = 2023

deo = Cube(joinpath(path,"EventCube_ranked_pot0.01_ne0.1.zarr"))[Ti=Date(startyear)..Date(endyear,12,31)]
tempo = year.(lookup(deo, :Ti))
# aggregate by decade? only 5 points...

include("../src/stats.jl")

function get_diff!(xout,xin,tempo::Vector{Int64}; rule::Function=x -> (x .>1) .& (iseven.(x.+1)))
    df = DataFrame(Year=tempo, Type=xin) |>
    (df -> DataFrames.subset(df, :Type => rule)) |>
    (df -> DataFrames.groupby(df, :Year)) |>
    (df -> DataFrames.combine(df, nrow)) |>
    (df -> leftjoin(DataFrame(Year=unique(tempo)), df, on=[:Year => :Year], order=:left)) |>
    (df -> DataFrames.select!(df, :Year, :nrow => x -> replace(x, missing=>0))) |>
    (df -> rename!(df, :nrow_function => "ndays")) |>
    (df -> sort(df, order(:Year)))

    p1 = mean(df.ndays[1:30])
    p2 = mean(df.ndays[31:end])
    diff = (p2 - p1)#/p1 *100
    return xout .= [Float32(p1), Float32(p2), Float32(diff)]
end

xout=[Float32(0.0), Float32(0.0), Float32(0.0)];
xin=deo.data[:,10,550];
dnh = get_diff!(xout,xin,tempo)
# hot = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(1)) .== Int8(1))
# dry30 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(2)) .== Int8(2))
# dry90 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(4)) .== Int8(4))
# dry180 = get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(8)) .== Int8(8))
# anyex = get_trend!(xout,xin,tempo,rule = x -> (x .> 0 .&& x .< 16))
# noex= get_trend!(xout,xin,tempo,rule = x -> (x .& Int8(16)) .== Int8(16))

# to get bivariate palette, need to 1st bin the bivariate space
# scatter a subsample?


function diffmap(ds; title=ds.diff.properties["name"], kwargs...)
    f = Figure(size=(800,1200))
    ax1 = Axis(f[1,1], title=title, aspect=DataAspect())
    # replace!(data, 0 => NaN);
    lon = circshift(map(x -> x >= 180 ? x-360 : x, lookup(ds.p1, :longitude)), (180/0.25));
    lat = lookup(ds.p1, :latitude);
    data1 = circshift(convert(Array{Float32},ds.p1), (180/0.25, 0));
    h1 = heatmap!(ax1,lon,lat,data1; kwargs...);
    # coastlines
    cl=lines!(ax1, 
        GeoMakie.coastlines(),
        color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 1000)

    cb = Colorbar(f[1,2], h1)#, vertical=false)

    # 
    ax2 = Axis(f[2,1], aspect=DataAspect())
    data2 = circshift(convert(Array{Float32},ds.p2), (180/0.25, 0));
    h2 = heatmap!(ax2,lon,lat,data2; kwargs...);
    # coastlines
    cl=lines!(ax2, 
        GeoMakie.coastlines(),
        color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 1000)

    cb = Colorbar(f[2,2], h2)#, vertical=false)
    #
    ax3 = Axis(f[3,1], aspect=DataAspect())
    data3 = circshift(convert(Array{Float32},ds.diff), (180/0.25, 0));
    h3 = heatmap!(ax3,lon,lat,data3; kwargs...);
    # coastlines
    cl=lines!(ax3, 
        GeoMakie.coastlines(),
        color = :black, linewidth=0.85)
    translate!(cl, 0, 0, 1000)

    cb = Colorbar(f[3,2], h3)#, vertical=false)
    #    
    f
end


@time diff_dh = mapCube(get_diff!,deo,tempo, indims=InDims(:Ti), outdims=OutDims(Dim{:Variable}(["p1","p2","diff"]),outtype=Float32,path=joinpath(path,"diffmap_dh_1970to1999_2000to2023.zarr"), chunksize=:max, overwrite=true, layername="dh"))
# 225.314949 seconds (805.31 M allocations: 325.946 GiB, 9.45% gc time, 35.86% compilation time)
diff_dh=open_dataset(joinpath(path,"diffmap_dh_1970to1999_2000to2023.zarr"))
f = diffmap(diff_dh, title="Average annual number of extremely dry and hot days in 1970-1999 versus 2000-2023)")
save(joinpath(path,"fig/diffmap_dh_1970to1999_2000to2023.png"), f)

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