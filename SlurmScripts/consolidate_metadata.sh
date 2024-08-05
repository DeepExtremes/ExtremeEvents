#!/bin/bash 
PATH2CUBE="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube_2023.zarr"
python -c "import zarr; g = zarr.open_group(\"${PATH2CUBE}\"); zarr.consolidate_metadata(g.store)"

# # manually edit .zmetadata: modify per dimensions :Time => :time
# # find
#  ".zattrs": {
#             "_ARRAY_DIMENSIONS": [
#                 "Time",
# # replace by
# ".zattrs": {
#             "_ARRAY_DIMENSIONS": [
#                 "time",
# # in ./time/.zattrs, replace Time  => time
# # in ./pet/.zattrs
# rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/Time
# rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/layer
# # consolidate metadata
# python -c "import zarr; g = zarr.open_group(\"${PATH2CUBE}\"); zarr.consolidate_metadata(g.store)"
