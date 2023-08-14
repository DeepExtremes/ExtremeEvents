path2cube = "/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr"
python -c "import zarr; g = zarr.open_group(${path2cube}); zarr.consolidate_metadata(g.store)"