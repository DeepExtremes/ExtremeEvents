# compare PEI_x at ICOS sites with SPEI_x from Pohl (EOBS) and with SPEI_x from Liu (ERA5)
# I don't expect good alignment...
using  YAXArrays,Zarr,NetCDF
import CSV, Statistics
using DataFrames
using CairoMakie, Colors
using Dates, DimensionalData

path2Dheed = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"
peis = open_dataset(joinpath(path2Dheed,"PEICube.zarr/"))
deo =  open_dataset(joinpath(path2Dheed,"EventCube_ranked_pot0.01_ne0.1.zarr/"))

# download SPEI from Pohl https://zenodo.org/records/7509504
path2spei_pohl = "/Net/Groups/BGI/work_2/scratch/mweynants/Pohl_long-term_2023_SPEI/"
# download icos stations from https://www.icos-cp.eu/observations/station-network

# icos = CSV.read(joinpath(path2spei_pohl,"icos_stations.csv"), DataFrame)
# download again from  https://zenodo.org/records/7473637
icos = CSV.read(joinpath(path2spei_pohl,"icos_list_lon_lat.csv"), DataFrame)

# download SPEI from Liu https://zenodo.org/records/8060268
# ;cd /Net/Groups/BGI/work_2/scratch/mweynants/SPEI_Liu
# ;wget https://zenodo.org/records/8060268/files/30days.zip
# ;unzip 30days.zip -d 30days
# wget https://zenodo.org/records/8060268/files/90days.zip
# wget https://zenodo.org/records/8060268/files/180days.zip
path2spei_liu = "/Net/Groups/BGI/work_2/scratch/mweynants/SPEI_Liu"

cols = Makie.wong_colors()

elem_1 = LineElement(color = cols[1], linestyle = nothing,
        points = Point2f[(0.5, 0), (0.5, 1)], label = "SPEI ≤ -2 : Extremely dry")
elem_2 = LineElement(color = cols[2], linestyle = nothing,
        points = Point2f[(0.5, 0), (0.5, 1)], label = "SPEI ≤ -2 : Extremely dry")
elem_3 = LineElement(color = cols[4], linestyle = nothing,
        points = Point2f[(0.5, 0), (0.5, 1)], label = "PEI ≤ 1st percentile")

agree = DataFrame(site = "icos.site",yr=1, d=30, pohl=tuple([1]), liu=tuple([1]), dheed=tuple([1]))
# mkdir(joinpath(path2Dheed,"fig/icos"))
for yr in 2000:2021
lspies30 = open_dataset(joinpath(path2spei_liu,"30days/Daily_SPEI_$(yr)_30Day.nc"),)
lspies90 = open_dataset(joinpath(path2spei_liu,"90days/Daily_SPEI_$(yr)_90Day.nc"),)
lspies180 = open_dataset(joinpath(path2spei_liu,"180days/Daily_SPEI_$(yr)_180Day.nc"),)

# # chunking of this thing is terrible to look at time series, even for a year
# # would be better to extract all locations in one read operation
# slspies30 = lspies30[lon = Near(icos.lon), lat = Near(icos.lat)].spei
# data = slspies30.data[:,:,:];

peis20 = peis[Ti = Date(yr,1,1)..Date(yr,12,31)]

for i = eachindex(icos.ID)
ilon = icos.lon[i] < 0 ? icos.lon[i] +360 : icos.lon[i]

pspeis = CSV.read(joinpath(path2spei_pohl,"SPEI_$(icos.site[i]).csv"), DataFrame) |>
    (df -> subset(df, :Date => x -> (x .>= Date(yr,1,1)) .&& (x .<= Date(yr,12,31))))

qpe30N = Statistics.quantile(skipmissing(peis.pei_30[lon = Near(ilon), lat = Near(icos.lat[i])][:]), [0.01, 0.05, 0.1])
qpe90N = Statistics.quantile(skipmissing(peis.pei_90[lon = Near(ilon), lat = Near(icos.lat[i])][:]), [0.01, 0.05, 0.1])
qpe180N = Statistics.quantile(skipmissing(peis.pei_180[lon = Near(ilon), lat = Near(icos.lat[i])][:]), [0.01, 0.05, 0.1])

fig = Figure(size=(1000,600), );
Label(fig[1,1, Top()], 
        text = icos.site[i] * " Year $(yr)",
        halign = :left
        )
ax30 = Axis(fig[1,1], ylabel = "SPEI_30 \n PEI_30 (mm/d)")
Label(fig[1,1, Right()], 
        text = "30 days",
        rotation = π/2, 
        fontsize = 20,
        )
psp30 = lines!(ax30, parse.(Float64,pspeis.SPEI_30), label = "SPEI Pohl", color = cols[1])
lsp30 = lines!(ax30, lspies30.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:], label="SPEI Liu", color = cols[2])
# sp_2 = hlines!(ax30,-2, label="SPEI ≤ -2 : Extremely dry")
# sp_15 = hlines!(ax,-1.5, label="Severly dry")
qp30 = lines!(ax30, peis20.pei_30[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:], label="PEI", color=cols[4])
# qp30_01 = hlines!(ax30, qpe30N[1])
# plot threshold with transparent vertical lines
append!(agree, DataFrame(site = icos.site[i],yr=yr, d=30, pohl=tuple(), liu=tuple(), dheed=tuple()), promote=true)
x = parse.(Float64,pspeis.SPEI_30) .≤ -2
if any(x)
agree[end,:pohl] = tuple(findall(x))
vlines!(ax30, 
        # ti, repeat([1], length(ti)), 
        agree[end,:pohl][1],
        # gap = 0, 
        # color = map(x -> cols[x+1], parse.(Float64,pspeis.SPEI_30) .≤ -2),
        alpha = 0.5,
        label="SPEI ≤ -2 : Extremely dry", 
        color=cols[1],
        )
end
x = lspies30.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:] .≤ -2
if any(x)
agree[end,:liu] = tuple(findall(x))
vlines!(ax30, 
        agree[end,:liu][1],
        alpha = 0.5,
        color=cols[2],
        )
end
x = peis20.pei_30[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:] .≤ qpe30N[1]
if any(x)
agree[end,:dheed] = tuple(findall(x))
vlines!(ax30, 
        agree[end,:dheed][1],
        alpha = 0.5,
        color=cols[4],
        )
end
ax90 = Axis(fig[2,1], ylabel = "SPEI_90 \n PEI_90 (mm/d)")
Label(fig[2,1, Right()], 
        text = "90 days",
        rotation = π/2, 
        fontsize = 20,
        )
lines!(ax90, parse.(Float64,pspeis.SPEI_90), color=cols[1],)
lines!(ax90, lspies90.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:], color=cols[2],)
lines!(ax90, peis20.pei_90[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:], color=cols[4],)
append!(agree, DataFrame(site = icos.site[i],yr=yr, d=90, pohl=tuple(), liu=tuple(), dheed=tuple()), promote=true)
x = parse.(Float64,pspeis.SPEI_90) .≤ -2
if any(x)
agree[end,:pohl] = tuple(findall(x))
vlines!(ax90, 
        agree[end,:pohl][1],
        alpha = 0.5,
        label="SPEI ≤ -2 : Extremely dry",
        color=cols[1],
        )
end
x = lspies90.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:] .≤ -2
if any(x)
agree[end,:liu] = tuple(findall(x))
vlines!(ax90, 
        agree[end,:liu][1],
        alpha = 0.5,
        color=cols[2],
        )
end
x = peis20.pei_90[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:] .≤ qpe90N[1]
agree[end,:dheed] = tuple(findall(x))
if any(x)
vlines!(ax90, 
        agree[end,:dheed][1],
        alpha = 0.5,
        color=cols[4],
        )
end
ax180 = Axis(fig[3,1], ylabel = "SPEI_180 \n PEI_180 (mm/d)")
Label(fig[3,1, Right()], 
        text = "180 days",
        rotation = π/2, 
        fontsize = 20,
        )
lines!(ax180, parse.(Float64,pspeis.SPEI_180), color=cols[1],)
lines!(ax180, lspies180.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:], color=cols[2],)
lines!(ax180, peis20.pei_180[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:], color=cols[4],)
append!(agree, DataFrame(site = icos.site[i],yr=yr, d=180, pohl=tuple(), liu=tuple(), dheed=tuple()), promote=true)
x = parse.(Float64,pspeis.SPEI_180) .≤ -2
if any(x)
agree[end,:pohl] = tuple(findall(x))
vlines!(ax180, 
        # ti, repeat([1], length(ti)), 
        agree[end,:pohl][1],
        # gap = 0, 
        # color = map(x -> cols[x+1], parse.(Float64,pspeis.SPEI_30) .≤ -2),
        alpha = 0.5,
        label="SPEI ≤ -2 : Extremely dry",
        color = cols[1],
        )
end
x = lspies180.spei[lon = Near(icos.lon[i]), lat  = Near(icos.lat[i])].data[:] .≤ -2
if any(x)
agree[end,:liu] = tuple(findall(x))
vlines!(ax180, 
        # ti, repeat([1], length(ti)), 
        agree[end,:liu][1],
        # gap = 0, 
        # color = map(x -> cols[x+1], parse.(Float64,pspeis.SPEI_30) .≤ -2),
        alpha = 0.5,
        label="SPEI ≤ -2 : Extremely dry",
        color = cols[2],
        )
end
# translate!(b, 0, 0, 10)
x = peis20.pei_180[longitude = Near(ilon), latitude  = Near(icos.lat[i])].data[:] .≤ qpe180N[1]
if any(x)
agree[end,:dheed] = tuple(findall(x))
vlines!(ax180, 
        # ti, repeat([1], length(ti)), 
        agree[end,:dheed][1],
        # gap = 0, 
        # color = map(x -> cols[x+1], parse.(Float64,pspeis.SPEI_30) .≤ -2),
        alpha = 0.5,
        label="PEI ≤ 1st percentile",
        color = cols[4],
        )
end
# translate!(b, 0, 0, 10)

plots = [psp30,lsp30,qp30,elem_1,elem_2,elem_3]
fig[4,1] = Legend(fig,
        plots,
        map(x -> x.label, plots), 
        position = :lb, 
        orientation = :horizontal,
        nbanks = 3,
        framevisible = false,
        )

# fig
save(joinpath(path2Dheed,"fig/icos/$(icos.site[i])_$(yr).png"), fig)
end
end

CSV.write(joinpath(path2Dheed,"fig/icos/agree.csv"),agree[2:end,:],)
# 
# modify csv file offline
# 
agree = CSV.read(joinpath(path2Dheed,"fig/icos/agree.csv"), DataFrame, )

# count number of matches and non matches, for each year and each station 
str2vec = function (str)
        @assert str[1] == '[' && str[end] == ']'
        map(x->parse(Int16,x),split(str[2:end-1], ", "))
end

transform!(agree, :pohl => ByRow(x -> ismissing(x) ? missing : str2vec(x)), renamecols=false)
transform!(agree, :liu => ByRow(x -> ismissing(x) ? missing : str2vec(x)), renamecols=false)
transform!(agree, :dheed => ByRow(x -> ismissing(x) ? missing : str2vec(x)), renamecols=false)

isamatch = function (x,y,z)
        if !ismissing(y) && length(y) >= 365
                y = missing
        end
        if all(ismissing.((x,y,z)))
                return (dheed_in_pohl_and_liu, dheed_in_pohl, dheed_in_liu, dheed_nin, pohl_in_liu, pohl_nin , liu_nin) = zeros(7)
        end
        if !ismissing(z)
                if !ismissing(x) && !ismissing(y)
                        dheed_in_pohl_and_liu = sum(in.(z,Ref(x)) .& in.(z,Ref(y)))
                        dheed_in_pohl = sum(in.(z,Ref(x))) - dheed_in_pohl_and_liu
                        dheed_in_liu = sum(in.(z,Ref(y))) - dheed_in_pohl_and_liu
                        dheed_nin = length(z) - dheed_in_pohl_and_liu - dheed_in_pohl - dheed_in_liu
                        pohl_in_liu = sum(in.(x, Ref(y))) - dheed_in_pohl_and_liu
                        pohl_nin = length(x) - dheed_in_pohl_and_liu - dheed_in_pohl - pohl_in_liu
                        liu_nin = length(y) - dheed_in_pohl_and_liu - dheed_in_liu - pohl_in_liu
                elseif !ismissing(x) && ismissing(y)
                        (dheed_in_pohl_and_liu, dheed_in_liu, pohl_in_liu, liu_nin) = zeros(4)
                        dheed_in_pohl = sum(in.(z,Ref(x))) 
                        dheed_nin = length(z)  - dheed_in_pohl 
                        pohl_nin = length(x)  - dheed_in_pohl 
                elseif ismissing(x) & !ismissing(y)
                        (dheed_in_pohl_and_liu, dheed_in_pohl, pohl_in_liu, pohl_nin ) = zeros(4)
                        dheed_in_liu = sum(in.(z,Ref(y))) 
                        dheed_nin = length(z) - dheed_in_liu
                        liu_nin = length(y) - dheed_in_liu
                else # both missing
                        dheed_nin = length(z)
                        (dheed_in_pohl_and_liu, dheed_in_pohl, dheed_in_liu, pohl_in_liu, pohl_nin , liu_nin) = zeros(6)

                end
        else
                (dheed_in_pohl_and_liu, dheed_in_pohl, dheed_in_liu, dheed_nin) = zeros(4)
                if !ismissing(x) && !ismissing(y)
                        pohl_in_liu = sum(in.(x, Ref(y)))
                        pohl_nin = length(x)  - pohl_in_liu
                        liu_nin = length(y) - pohl_in_liu
                elseif !ismissing(x) && ismissing(y)
                        (pohl_in_liu, liu_nin) = zeros(2)
                        pohl_nin = length(x) - pohl_in_liu
                else#if ismissing(x) & !ismissing(y)
                        (pohl_in_liu, pohl_nin ) = zeros(2)
                        liu_nin = length(y) - pohl_in_liu
                end
        end
        return (dheed_in_pohl_and_liu, dheed_in_pohl, dheed_in_liu, dheed_nin, pohl_in_liu, pohl_nin , liu_nin)     
end

transform!(agree, [:pohl, :liu, :dheed] => ByRow((x, y, z) -> isamatch(x,y,z)) => [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin])
transform!(agree, [:pohl, :liu, :dheed] .=> ByRow(x -> isa(x, Vector)) .=> [:pohl_any, :liu_any, :dheed_any])

# df1 = groupby(agree, [:pohl_any, :liu_any, :dheed_any]) |>
#         (gdf -> combine(gdf, [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin] .=> sum , renamecols=false)) 


# merge all stations
df1 = agree |> #groupby(agree, :d, ) |>
        (gdf -> combine(gdf, [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin] .=> sum , renamecols=false)) |>
        (df -> stack(df, [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin]));

df1[!,:dheed] = vcat(ones(4), zeros(3));
df1[!,:pohl,] = [1,1,0,0,1,1,0];
df1[!,:liu] = [1,0,1,0,1,0,1];
transform!(df1, [:dheed, :pohl, :liu] => ByRow((x,y,z) -> RGB(x.*0.9,y.*0.9,z.*0.9)) => :color);
df1[!, :label] = [
        "Dheed ∩ Pohl ∩ Liu",
        "Dheed ∩ Pohl",
        "Dheed ∩ Liu",
        "Dheed only",
        "Pohl ∩ Liu",
        "Pohl only",
        "Liu only"]
# colorant"hsl(0, 100%, 25%)"
h=10; l = .7; s = 1
cols = [HSL(h,0,.75),
        HSL(h+60, s,l),
        HSL(h+300, s,l),#colorant"hsl(300, 100%, 75%)",
        HSL(h, s,l),#colorant"hsl(0, 100%, 75%)", 
        HSL(h+180,s,l),#colorant"hsl(180, 100%, 75%)",
        HSL(h+120, s,l),#colorant"hsl(120, 100%, 75%)", 
        HSL(h+240, s,l),#colorant"hsl(240, 100%, 75%)", 
        ]
# Colors.protanopic.(cols)
# Colors.deuteranopic.(cols)
# Colors.tritanopic.(cols)
# [colorant"hsl(0, 100%, 75%)", colorant"hsl(60, 100%, 75%)", colorant"hsl(120, 100%, 75%)", colorant"hsl(180, 100%, 75%)", colorant"hsl(240, 100%, 75%)"]

# range(HSL(colorant"red"), stop=HSL(colorant"green"), length=3)
# range(colorant"hsl(0, 80%, 25%)", stop = colorant"hsl(120, 80%, 25%)", length=3)


f,ax,plt = pie(df1.value, color=cols, #Makie.wong_colors(0.8)[1:7],
        radius = 4,
        inner_radius = 2,
        strokecolor = :white,
        strokewidth = 2,
        axis = (autolimitaspect = 1, )
        )
hidedecorations!(ax)
c = Colorbar(f[1,2], colormap=cgrad(cols, categorical=true), 
        limits = (0.5,7.5), ticks = (1:7, df1.label))                
f
save(joinpath(path2Dheed, "fig", "spei_validation.png"), f)
# barplot(gdf)
transform!(df1, :value => x -> round.(100 .* x ./ sum(x), digits=2) )
show(stdout, MIME("text/latex"), select(df1, :label, :value_function))

tmp = select(df1, :label, :value_function) |>
        (df -> filter(:label => x -> occursin("Dheed", x), df)) |>
        (df -> transform(df, :value_function => x -> round.(100 .* x ./ sum(x), digits=2) ))
show(stdout, MIME("text/latex"), tmp)

df2 = groupby(agree, :d, ) |>
        (gdf -> combine(gdf, [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin] .=> sum , renamecols=false)) |>
        (df -> stack(df, [:dheed_in_pohl_and_liu, :dheed_in_pohl, :dheed_in_liu, :dheed_nin, :pohl_in_liu, :pohl_nin , :liu_nin])) |>
        (df -> unstack(df, :variable, :d, :value)) |>
        (df -> select(df, [2,3,4] .=> x -> round.(100 .* x ./ sum(x), digits=2) ))

show(stdout, MIME("text/latex"), hcat(tmp, df2[1:4,:]))

show(stdout, MIME("text/latex"), hcat(select(df1, :label, :value_function)[5:7,:], df2[5:7,:]))
      