using EarthDataLab, YAXArrays, Zarr, DiskArrays

using YAXArrayBase

outpath = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr"

ds = YAXArrays.Datasets.open_mfdataset("/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/PET/*.zarr")
ds = setchunks(ds,(lon=60,lat=60,time=5844))

savedataset(ds, path=outpath,append=true,max_cache=4e9)