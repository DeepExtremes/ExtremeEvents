using SlurmClusterManager, Distributed

varname = "TER"
orig_names = ("TER", "TER_HB")
unit = "gC/m^2/d"

#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate(@__DIR__)
end
@everywhere using EarthDataLab, YAXArrays, Zarr, DiskArrays, Statistics, Dates, Interpolations


@everywhere function ensmean(oc,ics...)
    v = Union{Missing,Float32}[]
    for i in eachindex(first(ics))
        empty!(v)
        foreach(ics) do ic
            push!(v, ic[i])
        end
        s = zero(Float32)
        n = 0
        for vv in v
            if !ismissing(vv) && (vv>0) && (vv<40)
                s+=Float32(vv)
                n+=1
            end
        end
        if n > 5
            oc[i] = s/n
        else
            oc[i] = missing
        end
    end
end

@everywhere function agg(xout,xin;fac=3)
    for ilat in axes(xout,2)
        for ilon in axes(xout,1)
            co = view(xin,(fac*(ilon-1)+1):ilon*fac,(fac*(ilat-1)+1):ilat*fac)
            xout[ilon,ilat] = mean(skipmissing(co))
        end
    end
end

function compute_fluxcom_ensemble()
    allcubes = []
    basep = "/Net/Groups/BGI/work_3/FluxcomDataStructure/internal/QOVbZtjf6sSIZowqIc99/CarbonFluxes/RS_V006/raw/8daily/"
    for method in readdir(basep)
        for v in orig_names
            p  = joinpath(basep, method,v,"*.nc")
            ds = YAXArrays.Datasets.open_mfdataset(p)
            ar = ds.cubes |> values |> first
            push!(allcubes,ar)
        end
    end

    id = [InDims("Lon","Lat") for _ in 1:length(allcubes)]
    od = OutDims("Lon","Lat",path = "/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_$(varname)_ens.zarr",backend=:zarr,overwrite=true)


    mapCube(ensmean,(allcubes...,),indims = id, outdims=od,max_cache=5e9)
end


#Load another dataset to get the chunk sizes

function aggregate_to_quarterdeg()
    dstest = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/ERA5Data.zarr/").t2m
    dstest = dstest[time=2001:2020]

    c = Cube("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_$(varname)_ens.zarr")

    #Aggregate in space first    
    
    mapCube(agg,c,
        indims=InDims("Lon","Lat"),
        outdims = OutDims(
            dstest.longitude,
            dstest.latitude,
            path = "/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_$(varname)_ens_agg.zarr",
            overwrite=true,
        )
    )
end



function interpolate_time()
    dstest = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/ERA5Data.zarr/").t2m
    dstest = dstest[time=2001:2020]

    c_quarter = Cube("/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_$(varname)_ens_agg.zarr")

    #The time axis always only stores the beginning of a time period. In order to make interpolation work, we need to shift this. 
    newtimes = copy(c_quarter.time.values)
    for i in 1:length(newtimes)-1
        newtimes[i] = newtimes[i] + (newtimes[i+1]-newtimes[i])/2
    end
    lastyear = year(newtimes[end])
    newtimes[end] = newtimes[end] + (DateTime(lastyear+1) - newtimes[end])/2

    renameaxis!(c_quarter,:time=>RangeAxis("Time",newtimes))


    icc = interpolatecube(c_quarter,Dict("Time"=>dstest.time.values),bc = Dict("Time"=>missing))

    dssym = Symbol(varname)
    dslist = (dssym=>icc,)

    ds = Dataset(;dslist...)
    ds.properties["Description"] = "Ensemble of Fluxcom RSv006 product $(varname) over $(orig_names) as well as all available machine learning methods. "
    ds[dssym].properties["name"] = varname
    ds[dssym].properties["units"] = unit
    dschunked = setchunks(ds,(lon=60,lat=60,time=dstest.chunks.chunks[3]))

    savedataset(dschunked,path = "/Net/Groups/BGI/work_1/scratch/fgans/DeepExtremes/Fluxcom_$(varname)_final.zarr",max_cache=5e9,overwrite=true)
end


#compute_fluxcom_ensemble()

#aggregate_to_quarterdeg()

interpolate_time()
