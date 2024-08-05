import xarray as xr

def add2cube(oldcube, newcube):
    ds0 = xr.open_zarr(oldcube)
    ds1 = xr.open_zarr(newcube)
    ds1.to_zarr(
        oldcube,
        mode = "a",
        append_dim = "time",
        consolidated = True,
    )

    