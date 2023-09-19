#!/bin/bash 
PATH2CUBE="/Net/Groups/BGI/scratch/mweynants/DeepExtremes/v3/ERA5Cube.zarr"
python -c "import zarr; g = zarr.open_group(\"${PATH2CUBE}\"); zarr.consolidate_metadata(g.store)"