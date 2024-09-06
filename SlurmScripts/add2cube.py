import xarray as xr

def add2cube(oldcube, newcube):
    ds = xr.open_zarr(newcube)
    ds.to_zarr(
        oldcube,
        mode = "a",
        append_dim = "time",
        consolidated = True,
    )

    