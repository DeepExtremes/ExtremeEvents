## example of workflow

# point to ERA5 daily data
# zarr group
# zg = zopen("https://s3.bgc-jena.mpg.de:9000/xaida/ERA5Data.zarr",consolidated=true)
zg = zopen("https://s3.bgc-jena.mpg.de:9000/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)

# YAXArray Dataset
ds = open_dataset(zg)
# YAXArray
era = Cube(ds)
# run preprocessing
## calc PET
## compute PEI = P - E
# run compute_pei.jl
pei = Cube(open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/SPEICube.zarr", consolidated = true)))

# subsetcube
sub_tmax = subsetcube(era,
    variable = ["t2mmax"],
    # region = "Germany",
    lon=(10,11.5),
    lat=(51,50),
    time=(Date("2000-01-01"),Date("2001-12-31")),
    )
sub_pei = subsetcube(pei,
    region = "Germany",
    # lon=(10,11.5),
    # lat=(51,50),
    #time=(Date("2000-01-01"),Date("2001-12-31")),
    )

# ranked: normalize data between 0 and 1, with extremes of interest in the lower range

# daily maximum temperature extremes of interest are in the upper range, so we multiply values by -1 before we normalize them
ranked_t = rescale(sub_tmax,"/Users/mweynants/BGI/temp/tmax_ranked.zarr", multiplier = -1)
# PEI values of interest are in the lower range (P-E << 0), so we don't transform them before we normalize them
ranked_pei = rescale(sub_pei,"/Users/mweynants/BGI/temp/pei_ranked.zarr")

# # with different data sources?
# using YAXArrays.Datasets: merge_datasets
# ds2 = open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/SPEICube.zarr", consolidated = true, fill_as_missing = true))
# mds = merge_datasets((ds,ds2))
# ranked1 = rescale(Cube(mds),"/Users/mweynants/BGI/temp/ranked1.zarr")


# smoothed: apply spatial filter on ranked data
## DO NOT RUN ON LOCAL MACHINE
# smoothed = smooth(ranked, "/Users/mweynants/BGI/temp/tmax_smoothed.zarr")
##

ranked_ds = open_dataset("/Users/mweynants/BGI/temp/ranked.zarr")[time=2000:2001, lon=(10,11.5), lat=(51,50)]
tmx_ds = open_dataset("/Users/mweynants/BGI/temp/tmax_ranked.zarr")
pei_ds = open_dataset("/Users/mweynants/BGI/temp/pei_ranked.zarr")
inputs = (tmx_ds.t2mmax[time=2019:2019],
        pei_ds.spei_30[time=2019:2019], pei_ds.spei_90[time=2019:2019], pei_ds.spei_180[time=2019:2019])

ranked_dc = Cube(ranked_ds)
# ranked_sc = subsetcube(ranked_dc,
#     #variable=["air_temperature_2m"], # if no [], drops the variable dimension.
#     lon=(10,11.5),
#     lat=(51,50),
#     time=(Date("2000-01-01"),Date("2001-12-31"))
#     )

extr = compute_extremes(ranked_dc, 0.01, "/Users/mweynants/BGI/temp/extr.zarr")
extr1 = compute_extremes(ranked_dc, 0.01, "/Users/mweynants/BGI/temp/extr1.zarr"; tresne = 0.1)

extr = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_ranked_pot0.01_ne0.1.zarr")
# extr = Cube("/Users/mweynants/BGI/temp/extr.zarr")

# simpleplot(dc::YAXArray, d::Int, year::Int, nlayer::Int; variable = nothing)
simpleplot(extr, 45, 2001, 4)
simpleplot(extr1, 40, 2001, 4)

#import Plots
Plots.plot(reshape(extr[:,:,:], (92,24)),legend = :outertopright)
Plots.plot(reshape(extr1[:,:,:], (92,24)),legend = :outertopright)

# other metheod for compute_extremes, with tuple as input 
tmp = compute_extremes(inputs,0.01, "/Users/mweynants/BGI/temp/tmp.zarr"; tresne = 0.1)
tmp1 = compute_extremes(inputs,0.01, "/Users/mweynants/BGI/temp/tmp1.zarr")
simpleplot(tmp, 182, 2019, 4)
simpleplot(tmp1, 182, 2019, 4)

eec = Cube(open_dataset(zopen("https://s3.bgc-jena.mpg.de:9000/xaida/EventCube_0.01.zarr", consolidated = true)))

## connected components
eec_bin = subsetcube(eec, time = 2019:2019, region="Germany").data[:,:,:] .> 0;
Plots.heatmap(eec_bin[182,:,:]'[end:-1:1,:], c = cgrad(:thermal, categorical = true), title = "Boolean layer")
# use ImageMorphology.label_components to label the connected blobs of events
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


ranked_pei = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")
c1 = ranked_pei[time=(Date(2020,6,15),Date(2020,6,17)), variable = "pei_30"]
tmp = smooth(c1, "smooth_pei.zarr")

tmp1 = tmp[:,:,:]
p = heatmap(tmp1[:,:,1]'[end:-1:1,:])
p = heatmap(tmp1[:,:,2]'[end:-1:1,:])

