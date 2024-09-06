@show years = eval(Meta.parse(ARGS[1])) # years as String

using SlurmClusterManager, Distributed
#Quick check if we are in a slurm job
if haskey(ENV,"SLURM_CPUS_PER_TASK")
    addprocs(SlurmManager())
end

@everywhere begin
    using Pkg
    Pkg.activate("$(@__DIR__)/..")
end

@everywhere using YAXArrays, Dates, YAXArrayBase, Zarr, NetCDF
# Variables: 
# t2m temperature at 2 meters [K]
# u10 10 m u-wind component [m s^(-1)]
# v10 10 m v-wind component [m s^(-1)]
# sp surface pressure [Pa]
# snr surface net radiation [J m^(-2)] accumulated over 1 hour (ssr + str ?)
# vpd_cf = swvp-vp Saturation water vapour pressure - vapour pressure = vapour pressure deficit [hPa] cf ?
@everywhere begin 
    varlist = (:t2m, :u10, :v10, :sp, :snr, :vpd_cf)
    # load land sea mask (land>=1;sea==0)
    lsmask = open_dataset("/Net/Groups/data_BGC/era5/e1/0d25_static/lsm.1440.721.static.nc")
    # filter for year 2019
    lsm = lsmask.lsm[Ti = At(DateTime("2019-01-01T13:00:00"))]
    # load data to memory
    lsmdata = lsm[:,:]

"""
    This function is used to convert the units of the input variables to approprate units.
"""
function unit_conv(t2m_K, u10m_m_s,v10m_m_s,sp_Pa,snr_J_m2,vpd_hPa,lsm)
    
    
    # Convert units, etc of hourly data.
    t2m_C = t2m_K - 273.15
    
    # Wind speed at 2 m (use wind profile to scale from 10 m)
    ws10m_m_s = sqrt(u10m_m_s^2 + v10m_m_s^2)
    # u(z2) = u(z1) * log((z2-d)/z0) / log((z1-d)/z0)
    # ???
    ws2m_m_s = ws10m_m_s*(4.87/(log(67.8*10-5.42)))
    
    # Net downward radiation.
    snr_MJ_m2 = snr_J_m2 / 1e6
    
    # Surface pressure Pa to kPa.
    sp_kPa = sp_Pa / 1000.
 
    # Vapour pressure deficit hPa to kPa
    vpd_kPa = vpd_hPa / 10

    # soil heat flux (condition, day, night) 
    G_MJ_m2 = snr_MJ_m2 < 0.0 ? snr_MJ_m2*0.1 : snr_MJ_m2*0.5 
    # mask out sea
    G_MJ_m2 = lsm > 0.5 ? G_MJ_m2 : 0.0

    # add Cd denominator constant for time step (condition, day, night) 
    # For daily PET, Cd = 0.34
    # effect of changing Cd between day and night is negligeable when integrated over 24 hours.
    Cd = snr_MJ_m2 < 0.0 ? 0.24 : 0.96
    
    return t2m_C, ws2m_m_s, snr_MJ_m2, G_MJ_m2, sp_kPa, vpd_kPa, Cd
end


function calculate_pet(
        t2m_C,         # Daily mean temperature at 2 m
        ws2m_m_s,      # Windspeed at 2 m
        snr_MJ_m2,     # Total net downward radiation MJ/m2/day
        G_J_m2,        # Ground Heat Flux
        sp_kPa,        # Surface Pressure
        vpd_kPa,       # Vapour pressure deficit kPa
        Cd             # Cd denominator constant for time step
        )     
    """
    This is the function that calculate the PET based on the PM method.
    """
    # Constants.
    lmbda = 2.45  # Latent heat of vaporization [MJ kg -1] (simplification in the FAO PenMon (latent heat of about 20°C)
    cp = 1.013e-3 # Specific heat at constant pressure [MJ kg-1 °C-1]
    eps = 0.622   # Ratio molecular weight of water vapour/dry air
   
#     # Atmospheric pressure [kPa] eq 7 in FAO.
     P_kPa = sp_kPa #101.3*((293.0-0.0065*height_m) / 293.0)**5.26

#     # Psychrometric constant (gamma symbol in FAO) eq 8 in FAO.
     psychometric_kPa_c = cp*P_kPa / (eps*lmbda)

#     # Saturation vapour pressure, eq 11 in FAO.
    svp_kPa = 0.6108*exp((17.27*t2m_C) / (t2m_C+237.3))

#     # Delta (slope of saturation vapour pressure curve) eq 13 in FAO.
     delta_kPa_C = 4098.0*svp_kPa / (t2m_C+237.3)^2

#     # Actual vapour pressure, eq 14 in FAO.
#     avp_kPa = 0.6108*np.exp((17.27*dewpoint2m_C) / (dewpoint2m_C+237.3))

#     # Saturation vapour pressure deficit.
#     svpdeficit_kPa = svp_kPa - avp_kPa

    
    # Calculate ET0, equation 6 in FAO Crop ET - Chap. 2 (https://www.fao.org/3/x0490e/x0490e06.htm#equation)
    # also Walter et al. 2001 The ASCE Standardized Reference Evapotranspiration Equation 
    numerator = 0.408*delta_kPa_C*(snr_MJ_m2 - G_J_m2) + 
            psychometric_kPa_c*(37/(t2m_C+273))*ws2m_m_s*vpd_kPa

    denominator = delta_kPa_C + psychometric_kPa_c*(1 + Cd * ws2m_m_s)
    # according to ASCE, denominator should vary between day and night: 
    # coefficient 0.34 (which is valid for daily PM) 
    # should be replaced by
    # 0.24 during daytime
    # 0.96 during nighttime
    
    ET0_mm_hr = numerator / denominator
    
    # express ET0 as a negative (upward) flux
    ET0_mm_hr = -max(zero(ET0_mm_hr),ET0_mm_hr)
    return ET0_mm_hr
end    

function pet_with_units(t2m_K, u10m_m_s,v10m_m_s,sp_Pa,snr_J_m2,vpd_hPa,lsm)
    """
    This is the function that converts units and calculates PET.
    """
    any(ismissing,(t2m_K, u10m_m_s,v10m_m_s,sp_Pa,snr_J_m2,vpd_hPa,lsm)) && return missing
    o = unit_conv(t2m_K, u10m_m_s,v10m_m_s,sp_Pa,snr_J_m2,vpd_hPa,lsm)
    calculate_pet(o...)
end

end # begin

# run function pet_with_units on era5 0d25_hourly data
# for each year

pmap(years) do yr # 2022 # 1950:2022 #[1950:1952; 1958; 1960; 1966:2022] # [1953:1957; 1959; 1960:1965] # [1959; 1961; 1962; 1963; 1964; 1965;]
    # get all variables from varlist into new dataset
    allvars = map(varlist) do vn
        @show vn
        # select only files that are NOT "era5_backextention"
        filelist = readdir("/Net/Groups/data_BGC/era5/e1/0d25_hourly/$vn/$yr/", join = true)
        filter!(!contains("back_extension"),filelist)
        filter!(endswith(".nc"),filelist)
        ds = YAXArrays.Datasets.open_mfdataset(filelist)
        # replaces vn by first values
        vn=>first(values(ds.cubes))
    end
    ds = Dataset(;allvars...)
    # time resolution: day
    tr = Date(yr):Day(1):Date(yr,12,31)
    # define dimensions' axes
    outaxes = [ds.longitude, ds.latitude, Dim{:time}(tr)] # previously :Time and not "time"!
    # create empty dataset
    outds = YAXArrays.Datasets.createdataset(
        YAXArrayBase.ZarrDataset,
        outaxes,
        layername = "pet",
        path = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/$yr.zarr",
        persist = true,
        T=Union{Missing,Float32},
        overwrite=true,
        chunksize = (1440,721,1),
    )[1]

    # loop around days and agregate hourly data to daily
    for d in tr
        loadeddata = map(allvars) do v
            subcube = v[2][time=d..(d+Day(1))]
            if size(subcube,3) >= 24
                subcube[:,:,1:24]
            else
                missing
            end
        end
        outar = outds[time = At(d)].data
        if any(ismissing,loadeddata)
            outar[:,:] .= missing
        else
            resu = broadcast(pet_with_units, loadeddata...,lsmdata);
            outar[:,:] = sum(resu,dims=3)[:,:,1]
        end
    end
end
println("Done!")

zarrlist = readdir("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET")

map(1950:2023) do yr
    !( "$yr.zarr" in zarrlist) ? println("PET for $yr is missing") : missing
end

# using Plots
# p = heatmap(outar[:,:])
# figname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/PET_20210101_mela.png"
# Plots.savefig(p, figname)
# # compare with old
# pet_2021 = open_dataset("/Net/Groups/BGI/work_1/scratch/fgans/PET/2021.zarr")
# tmp = subsetcube(pet_2021, time=d);
# p = heatmap(tmp.layer[:,:])
# figname = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/fig/PET_20210101_fabian.png"
# Plots.savefig(p, figname)
# # diff
# p = heatmap(outar[:,:] - tmp.layer[:,:], title = "Effect of day/night coefficient")
# # negligible (<1e-6) so no need to implement it!