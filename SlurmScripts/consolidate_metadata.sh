#!/bin/bash 
PATH2CUBE="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr"
python -c "import zarr; g = zarr.open_group(\"${PATH2CUBE}\"); zarr.consolidate_metadata(g.store)"

# manually edit .zmetadata: modify per dimensions :Time => :time
# find
 "pet/.zattrs": {
            "_ARRAY_DIMENSIONS": [
                "Time",
# replace by
"pet/.zattrs": {
            "_ARRAY_DIMENSIONS": [
                "time",
rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/Time
rm -rf /Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr/layer
# consolidate metadata
python -c "import zarr; g = zarr.open_group(\"${PATH2CUBE}\"); zarr.consolidate_metadata(g.store)"
