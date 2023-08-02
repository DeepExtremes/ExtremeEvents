using EarthDataLab, YAXArrays, Zarr, DiskArrays

using YAXArrayBase

outpath = "/scratch/fgans/XAIDAERACUBE.zarr"

ds = YAXArrays.Datasets.open_mfdataset("/Net/Groups/BGI/work_1/scratch/fgans/PET/*.zarr")
ds = setchunks(ds,(lon=60,lat=60,time=5844))

savedataset(ds, path=outpath,append=true,max_cache=4e9)