outpath = "/Net/Groups/BGI/work_2/scratch/mweynants/Dheed_v4/ERA5Cube.zarr"
using CondaPkg; CondaPkg.add("xarray"); CondaPkg.add("zarr")
using  PythonCall
zr = pyimport("zarr")
g = zr.open_group(outpath)
# rename dimensions
# lon to longitude
# lat to latitude
# zr.storage.rename(g, 'lat', 'latitude')
# zr.storage.rename(g, 'lon', 'longitude')
# consolidate_metadata
g.consolidate_metadata()


