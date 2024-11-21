## connected components

using SlurmClusterManager, Distributed

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    for iproc in 1:parse(Int,ENV["SLURM_NTASKS"])
        addprocs(1)
        sleep(0.001)
    end
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using YAXArrays, Statistics, DiskArrays, ImageMorphology, Dates, Zarr, NetCDF
@everywhere import ImageFiltering

# 
@everywhere include("../src/detection.jl")

# use ImageMorphology.label_components to label the connected blobs of events
# use Fabian's version to wrap around the globe
@everywhere include("../src/label.jl")


# split period into 13 years periods with 3 years overlap
# period = 1950:2022 #2016:2022 # (Date("2018-06-01"), Date("2018-09-16"))#could be added to file name "*"_"*replace("$period", ":" => "_")*"
# aperiod = "_1950_2022" #"_2016_2021"
## 16 years, 2 days overlap to get filtering correctly
p=1950:16:2020
# period = p
# aperiod = map( x -> replace(string(x), ":" => "_", r"^" => "_"),p)
@everywhere begin
    pot=0.01
    ne=0.1
    
    compound_events = true
    cmp = compound_events ? "cmp" : ""
    filter_events = true
    filtern = filter_events ? "S1_T3" : "" #"_Sdiam3_T5"
    filter_land = false
    land = filter_land ? "_land" : ""
    region = ""; # "Germany"

    if filter_events
        # filter maskarray to remove small events
        # diamond in space 100% + at least 3 contiguous times
        # set window (time window: Tx: x*2-1)
        window = (1,1,5)#(1,1,7) #(3,3,9)
        # # myfilter function needs following global variables (avoid computing them multiple times)
        # central get_diamond_indices
        Nh = map(x->Int((x+1)/2), window);
        # # only time filter
        # Nh = Int((window[3] + 1) /2 )
        # indices of 2D diamond
        diamondindices = get_diamond_indices(window[1])
        t = length(diamondindices); # 0.6 * length(diamondindices);
    end

    # set window and t accordingly
    #outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_ranked_pot" * string(pot) * "_ne" * string(ne) * "_sreld335.zarr"
    # outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/labelcube_ranked_pot$(pot)_ne$ne$cmp$filtern$aperiod$land.zarr"
    #outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/labelcube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"

    path2v4 = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4"
    # path2v3 = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3"
    inpath = "$path2v4/EventCube_ranked_pot$(pot)_ne$(ne).zarr"
    # inpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/EventCube_smoothed_pot" * string(pot) * "_ne" * string(ne) * ".zarr"
    c = Cube(inpath)
end

@sync @distributed for period in p
    # period=2014
    startyear=period;endyear=min(period+15, year(c.Ti[end]))
    @show aperiod = "$(startyear)_$(endyear)"
    outpath = "$path2v4/labelcube_ranked_pot$(pot)_ne$(ne)_$(cmp)_$(filtern)_$(aperiod)$(land)$(region).zarr"

clastyears = c[time=DateTime(startyear-1, 12, 30) .. DateTime(endyear+1, 1,2)]
# test:
# clastyears = c[time=period, latitude=45 ..50, longitude=10 .. 15]
# create binary array flagging all events (regardless of type of event) and load to memory
# keep out non-extremes, i.e. where event is < tresne or > 1-tresne
if compound_events
    # compound events are events where both temperature and any of the drought indicators are extremes
    # i.e. where the first bit (coding for temp extreme) is 1 and any other bit is also 1
    # function indsel(x)
    #     x > 0x01 && isodd(x)
    # end
    # @time maskarray = permutedims(
    #     map(indsel, clastyears.data)[:,:,:],
    #     (2,3,1)
        # );
    @time maskarray = permutedims(
        (broadcast(x -> isodd(x) && x > 0x01, view(clastyears.data, :,:,:))),
        (2,3,1)
        );
    # 136.463942 seconds (14.02 M allocations: 28.601 GiB, 0.68% gc time, 2.64% compilation time)
    # 1950-2021: 2084.692204 seconds (16.20 M allocations: 317.017 GiB, 0.35% gc time, 0.21% compilation time)
    # 1950-2022: 1784.978296 seconds (15.26 M allocations: 326.052 GiB, 0.09% gc time, 0.21% compilation time)
    # 2000-2022: 199.699789 seconds (215.27 k allocations: 72.378 GiB, 0.11% gc time)
    # test: 4.046847 seconds (5.23 M allocations: 399.573 MiB, 4.67% gc time, 90.28% compilation time)
    # versus with view 1.139852 seconds (272.27 k allocations: 56.643 MiB, 1.71% gc time, 37.92% compilation time)

else
    @time maskarray = permutedims(
        (broadcast(x -> x > 0x00 && x .< 0x10, view(clastyears.data,:,:,:))),
        (2,3,1)
        );
end

# import Plots
# red = reshape(reduce(+,maskarray, dims=3), size(maskarray)[1:2]);
# size(red)
# Plots.heatmap(red)

if filter_land
    lsm = 
    Cube(open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc"))[
    time=Near(DateTime("2019-01-01T13:00:00")),
    # latitude=45 ..50, longitude=10 .. 15
    # region = region,
    ]
    lsmask = lsm.data .> 0.5;
    @time maskedarray = maskarray .* lsmask
    # @time for t in 1: size(maskarray)[3]
    #     maskarray[:,:,t] = ifelse.(lsmask, maskarray[:,:,t], 0)
    # end
    else maskedarray = maskarray
end
# Plots.heatmap(lsmask)
# red1 = reshape(reduce(+,maskarray, dims=3), size(maskarray)[1:2]);
# Plots.heatmap(red1)

timeind = (ifelse(startyear == 1950, 1, 3), ifelse(endyear > 2020, 0, -2))

println("Filter events")
if filter_events
    # use ImageFiltering.mapwindow to apply filter
    # pass the loaded data to mapwindow and keep only the needed days
    @time filteredarray = ImageFiltering.mapwindow(
        myfilter, 
        maskedarray[:,:,:], 
        window; 
        border = ImageFiltering.Fill(false),
        )[:,:,timeind[1]:end+timeind[2]];
    ########
    # # 1916-2021 !!! 32027.975212 seconds (98.20 G allocations: 3.299 TiB, 1.51% gc time, 0.00% compilation time)
    # 2000-2022 ! 36920.514462 seconds (98.73 G allocations: 3.383 TiB, 0.82% gc time) (up to 44000 seconds)
    # @time filteredarray = ImageFiltering.mapwindow(mytimefilter, maskarray, window);

else
    @time filteredarray = maskedarray[:,:,timeind[1]:end+timeind[2]];
    # 77.158438 seconds (3.29 M allocations: 25.292 GiB, 0.71% gc time, 2.67% compilation time)
end
# axs = clastyears.axes[[2,3,1]]
# include("../src/plots.jl");
# hmf = hm(filteredarray, axs=axs)

# set d to longitude dimension 
# size(filteredarray)
# (lon, lat, time)
d = 1;
@time r = label_components(filteredarray, wrapdims=(d,));
# 1916-2021 8.853066 seconds (159.63 k allocations: 16.965 GiB, 1.89% gc time, 1.12% compilation time)
# 1950-2021 166.081424 seconds (162.74 k allocations: 203.452 GiB, 0.01% gc time, 0.11% compilation time)
# 1950-2022 24.839224 seconds (60 allocations: 36.733 GiB, 1.19% gc time
 
# # look at results
# Plots.heatmap(r[:,:,1]'[end:-1:1,:])
# Plots.heatmap(r[:,:,10]'[end:-1:1,:])
# Plots.heatmap(r[:,:,20]'[end:-1:1,:])
# Plots.heatmap(r[:,:,30]'[end:-1:1,:])
# Plots.heatmap(r[:,:,40]'[end:-1:1,:])
# Plots.heatmap(r[:,:,50]'[end:-1:1,:])
# Plots.heatmap(r[:,:,60]'[end:-1:1,:])
# Plots.heatmap(r[:,:,70]'[end:-1:1,:])
# Plots.heatmap(r[:,:,80]'[end:-1:1,:])
# Plots.heatmap(r[:,:,90]'[end:-1:1,:])
# Plots.heatmap(r[:,:,100]'[end:-1:1,:])

# define YAXArray
@time labelcube = YAXArray((c.longitude,c.latitude,Dim{:time}( Date(startyear,1,1):Day(1):Date(endyear,12,31))),r, chunks = DiskArrays.GridChunks(r,(120,120,90)))
nt = (:labels=>labelcube,)
ds = Dataset(;nt...)
savedataset(ds, path=outpath, overwrite=true)
nothing

print("done!")
# check
# using Plots
# labelcube = Cube(outpath)
# include("../src/plots.jl")
# simpleplot(labelcube,Date("2019-07-25"), 4, colours = cgrad(:darkterrain, 50, categorical = true), replacement = 0=>5e3)

end # pmap