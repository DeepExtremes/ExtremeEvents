using Distributed

# sleep(1)
# addprocs(20)
# sleep(1)

@everywhere begin
    using Pkg
    Pkg.activate("..")
end
@everywhere begin
    using Zarr, YAXArrays, EarthDataLab, Statistics, SphericalConvolutions
    include("../src/detection.jl")
end

zg = zopen("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/ERA5Data.zarr",consolidated=true, fill_as_missing = false)
era = open_dataset(zg)
tair = era.t2mmax#[time=1980:2021]

# pei = open_dataset("/Net/Groups/BGI/work_1/scratch/s3/xaida/v2/PEICube.zarr")
# peicube = Cube(pei)

# daily maximum temperature extremes of interest are in the upper range, so we multiply values by -1 before we normalize them
ranked_t = if ispath("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr")
    Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr")    
else
    rescale(tair,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_ranked.zarr", multiplier = -1)
end
# # # PEI values of interest are in the lower range (P-E << 0), so we don't transform them before we normalize them
# ranked_pei = if ispath("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")
#     Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")    
# else
#     rescale(peicube,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_ranks.zarr")
# end

# tmp = subsetcube(ranked_t, time = (Date(2020,6,15), Date(2020,6,17)))
# tmp1 = smooth(tmp, "test.zarr"; lbord=40, width=4)
smoothed_t = smooth(subsetcube(ranked_t, time=(2016,2022)), "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed_80_2016.zarr"; lbord=80, width=4)

# # close workers
# t = rmprocs(workers(), waitfor=0)
# wait(t)
# workers()

# smoothed_t = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed.zarr")
# smoothed_t40 = Cube("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/tmax_smoothed_40.zarr")
# # # my smooth function seems to make Julia crash...
# # smooth_pei =  = smooth(ranked_pei[variable="pei_30"], "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei_smoothed.zarr")
# # # smoothed_pei30 = smooth(ranked_pei[variable="pei_30"], "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei30_smoothed.zarr")
# # # smoothed_pei90 = smooth(ranked_pei[variable="pei_90"], "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei90_smoothed.zarr")
# # # smoothed_pei180 = smooth(ranked_pei[variable="pei_180"], "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/pei180_smoothed.zarr")

# # print("done! \n")
# # # # => try with old method.
# # # apply low pass filter to spatial data
# # function applylowpass(xout,xin;lbord = 20, width=2)
# #     im2 = [xin[j,i] for i in size(xin,2):-1:2,j in 1:(size(xin,1)-1)]
# #     imfiltered = SphericalConvolutions.lowpass(Float64.(im2),lbord = lbord,width=width)
# #     xout[1:end-1,end:-1:2] = permutedims(imfiltered)
# #     xout[end,:] = xout[end-1,:]
# #     xout[:,1] = xout[:,end]
# #     xout
# # end
# # indims = InDims("Lon","Lat")
# # outdims = OutDims("Lon","Lat",
# #     path="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/smoothed_pei_ranks.zarr",
# #     #path="smooth_pei.zarr",
# #     backend=:zarr,overwrite=true)

# # mapCube(applylowpass, ranked_pei, indims=indims, outdims=outdims, max_cache=5e8)


# # check
# ## plot rank, spectrum, filter ##
# t0 = readcubedata(tair[time=Date(2000,6,15)]);
# t1 = readcubedata(ranked_t[time=Date(2000,6,15)]);
# t2 = readcubedata(smoothed_t[time=Date(2020,6,15)]);
# t3 = readcubedata(smoothed_t40[time=Date(2020,6,15)]);

# # c1 = ranked_pei[time=Date(2020,6,15),variable="pei_30"]
# # c1 = readcubedata(c1)

# # c2 = smoothed_pei[time=Date(2020,6,15),variable="pei_30"]
# # c2 = readcubedata(c2)

# using Plots
# xin0 = t0.data[:,:];
# xin1 = t1.data[:,:];
# xin2 = t2.data[:,:];
# xin3 = t3.data[:,:];

# # xin1 = c1.data[:,:];
# # xin2 = c2.data[:,:];
# # # # reshuffle data
# # # im0 = [xin0[j,i] for i in size(xin0,2):-1:2,j in 1:(size(xin0,1)-1)]
# im1 = [xin1[j,i] for i in size(xin1,2):-1:2,j in 1:(size(xin1,1)-1)]
# im2 = [xin2[j,i] for i in size(xin2,2):-1:2,j in 1:(size(xin2,1)-1)]
# im3 = [xin3[j,i] for i in size(xin3,2):-1:2,j in 1:(size(xin3,1)-1)]
# lats = reverse(t1.latitude)[2:end]
# lons = t1.longitude[1:end-1]

# # # p0 = heatmap(lons,lats,im0,title="Daily max t2m")
# # # savefig(p0,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_t.png")
# p1 = Plots.heatmap(lons,lats,im1,title="Daily max t2m normalised")
# # savefig(p1,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_t_normalised.png")
# p2 = Plots.heatmap(lons,lats,im2,title="Daily max t2m smoothed")
# # savefig(p2,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_t_smoothed20.png")
# p3 = Plots.heatmap(lons,lats,im3,title="Daily max t2m smoothed 40")
# # savefig(p3,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_t_smoothed40.png")

# import PlotlyJS
# h1 = PlotlyJS.plot(PlotlyJS.histogram2d(x=reshape(im1,:),y=reshape(im2, :)), PlotlyJS.Layout(xaxis_title = "Normalised t2mmax", yaxis_title = "Smoothed t2mmax lbord 20", title= " Time " * string(Date(2000,6,15))))
# # PlotlyJS.savefig(h1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/ranked_vs_smooth20.png")
# h2 = PlotlyJS.plot(PlotlyJS.histogram2d(x=reshape(im1,:),y=reshape(im3, :)), PlotlyJS.Layout(xaxis_title = "Normalised t2mmax", yaxis_title = "Smoothed t2mmax lbord 40", title= " Time " * string(Date(2000,6,15))))
# # PlotlyJS.savefig(h1, "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/ranked_vs_smooth40.png")

# # p1 = heatmap(lons,lats,im1,title="PEI_30 normalised")
# # savefig(p1,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_pei_normalised.png")
# # p2 = heatmap(lons,lats,im2,title="PEI_30 smoothed")
# # savefig(p2,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/map_pei_smoothed.png")

# im_trans = do_transform(Float64.(im1))


# scalef(l;lbord=20,width=2) = tanh((lbord-l)/width)/2.0+0.5

# function getSpectrum(im)
#     im_trans = do_transform(Float64.(im))
#     map(0:SphericalConvolutions.lmax(im_trans)) do l
#         mm = mmax(im_trans,l)
#         sum(m->abs2(im_trans[l,m]),-mm:mm)
#     end
# end

# spec = getSpectrum(im1)

# p = Plots.plot(spec[2:end],yaxis=:log,xlabel="Wavenumber l",ylabel="Power",yrange=[1e-6,0.0],label="Unaltered Spectrum")
# Plots.plot!(p,clamp.(scalef.(1:720),1e-6+1e-7,1),yaxis=:log,label="Filter function")
# Plots.plot!(p,clamp.(scalef.(1:720, lbord=40, width=4),1e-6+1e-7,1),yaxis=:log,label="Filter function lbord 40")
# Plots.plot!(p,spec[2:end].*scalef.(1:719),label="Filtered Spectrum",title="Spherical Spectrum")
# Plots.plot!(p,spec[2:end].*scalef.(1:719, lbord=40, width=4),label="Filtered Spectrum",title="Spherical Spectrum")
# # savefig(p,"/Net/Groups/BGI/scratch/mweynants/DeepExtremes/spei_spectrum.png")

# imfiltered2 = SphericalConvolutions.lowpass(Float64.(im1),lbord = 20,width=2)
# imfiltered4 = SphericalConvolutions.lowpass(Float64.(im1),lbord = 40,width=4)
# h = PlotlyJS.plot(PlotlyJS.histogram2d(x=reshape(imfiltered2,:),y=reshape(imfiltered4, :)), Layout(xaxis_title = "Smoothed t2mmax lbord 20", yaxis_title = "Smoothed t2mmax lbord 40", title= "bbox Lon " * string(lon) * " Lat " * string(lat) * " Time " * string(period[1]) * " to " * string(period[2])))


# # p3 = heatmap(lons,lats,imfiltered,title = "Filtered PEI rank")
# # savefig(p3,"spei_filtered.png")