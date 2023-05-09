using Revise
using YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates

# addprocs(8)
# @everywhere begin
#     using Pkg
#     Pkg.activate("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/LargeScaleEvents/")
# end
# @everywhere using ParallelUtilities, YAXArrays, Zarr, WeightedOnlineStats, OnlineStats, DataFrames, Dates
# @everywhere mergefun(h1,h2) = Dict(k=>merge!(h1[k],h2[k]) for k in keys(h1))
# @everywhere function fit1(df)
#     df.y = year.(df.time)
#     dfg = groupby(df,:y)
#     allhists = Dict(i=>WeightedHist(-0.5:1.0:32.5) for i in 1950:2021)
#     for k in keys(dfg)
#         dfs = dfg[k]
#         fit!(allhists[k[1]],dfs.event,cosd.(dfs.latitude))
#     end
#     @show myid()
#     allhists
# end

# ds = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr/")



# t = CubeTable(event = ds.layer)

# annualstats = pmapreduce(mergefun,t) do tab
#     fit1(DataFrame(tab))
# end
# rmprocs(workers())


# allx = [value(annualstats[yr]).x[i] for i in 1:17, yr in 1950:2021]
# allres = [value(annualstats[yr]).y[i] for i in 1:17, yr in 1950:2021]
# using Plots
# allres_norm = allres./sum(allres,dims=1)
# labels = map([UInt8(i) for _ in 1:1, i in 1:16]) do i
#     n = ""
#     ((i & 0x01) > 0) && (n = join((n,"HW"),"_"))
#     ((i & 0x02) > 0) && (n = join((n,"D30"),"_"))
#     ((i & 0x04) > 0) && (n = join((n,"D90"),"_"))
#     ((i & 0x08) > 0) && (n = join((n,"D180"),"_"))
#     lstrip(n,'_')
# end
# p = plot(1950:2021,allres_norm[1:end-1,:]',labels=labels,lw=1,size=(800,400),dpi=300)

# savefig(p,"n_extremes.png")


using DiskArrayEngine, StatsBase, DiskArrays, FillArrays, NetCDF
using DiskArrayEngine: ProductArray, InputArray, LoopWindows, create_userfunction, NoFilter, GMDWop
ds = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/EventCube_ranked_pot0.01_ne0.1.zarr/")

lsm = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")

a = ds.layer.data

t = ds.time
years, nts = rle(year.(t))
cums = [0;cumsum(nts)]

stepvectime = [cums[i]+1:cums[i+1] for i in 1:length(nts)]
rp = ProductArray((stepvectime,1:size(a,2),1:size(a,3)))

rplsm = ProductArray((1:size(a,2),1:size(a,3)))


rplat = ProductArray((1:size(a,3),))
lat = ds.latitude.values

# rangeproduct[3]
inars = (InputArray(a,LoopWindows(rp,Val((1,2,3)))),InputArray(lat,LoopWindows(rplat,Val((3,)))),InputArray(lsm.lsm.data[:,:,1],LoopWindows(rplsm,Val((2,3)))))

outrp = ProductArray((1:72,[1:17],[1:2]))
outwindows = (LoopWindows(outrp,Val((1,4,5))),)

function fit_online!(xout,x,lat,lsm)
  cdl = cosd(lat) 
  for ix in x
    xout[Int(ix)+1,(lsm > 0.5)+1] += cdl
  end
end

init = 0.0
filters = (NoFilter(),)
fin_onine(x) = nobs(x) == 0 ? missing : OnlineStats.value(x)
f = create_userfunction(
    fit_online!,
    typeof(init),
    is_mutating = true,
    red = (h1,h2)->merge!.(h1,h2), 
    init = init, 
    finalize=identity,
    buftype = typeof(init),  
)
  
  optotal = GMDWop(inars, outwindows, f)
  
optotal.windowsize

da = DiskArrayEngine.results_as_diskarrays(optotal)

da[1]

da[1][72,:,:]

  loopranges = ProductArray((DiskArrays.RegularChunks(10,0,72),eachchunk(a).chunks[2:3]...,[1:1],[1:1]))
  b = zeros(72,17,2);
  outars = (InputArray(b,outwindows[1]),)
  DiskArrayEngine.run_loop(optotal,loopranges,outars)
  

allres = b
allres_norm = allres./sum(allres,dims=2)
labels = map([UInt8(i) for _ in 1:1, i in 1:16]) do i
    n = ""
    ((i & 0x01) > 0) && (n = join((n,"HW"),"_"))
    ((i & 0x02) > 0) && (n = join((n,"D30"),"_"))
    ((i & 0x04) > 0) && (n = join((n,"D90"),"_"))
    ((i & 0x08) > 0) && (n = join((n,"D180"),"_"))
    lstrip(n,'_')
end
using Plots
p = plot(1950:2021,
  allres_norm[:,2:end-1,2],
  labels=labels,
  size=(800,400),
  dpi=300,
  xlabel = "year", 
  legend=:outerright
)

savefig(p,"n_extremes.png")



function mapplot(data)
    data = data[[721:1440;1:720],end:-1:1]
    heatmap(data')
end

hwyear = ds.layer[time=1952:1952][:,:,:]
m = sum(==(0x01),hwyear,dims=1)
mapplot(m[1,:,:])

hwyear2 = ds.layer[time=2020:2020][:,:,:]
m2 = sum(==(0x01),hwyear2,dims=1)
mapplot(m2[1,:,:])

ii = CartesianIndex(1035,399)
zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
era.t2mmax
ts1952 = era.t2mmax.data[1035,399,:]
q1 = quantile(ts1952,0.99)

sub = 1:3000
p = plot(era.time[sub],ts1952[sub], labels=false)
hline!(p,[q1],labels=false)

savefig(p,"n_extremes.png")


