# rasterize Chena HUC for quick-and-dirty comparison
import os, shutil, rasterio
from rasterio.features import rasterize
import geopandas as gpd
import numpy as np

template_fn = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_testing_extent/wrf_extent_raw_TEST.tif'
shp_fn = '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_19040506.shp'

shp = gpd.read_file( shp_fn )

with rasterio.open( template_fn ) as tmp:
	meta = tmp.meta.copy()
	meta.update( compress='lzw' )
	shape = tmp.read(1).shape

shp = shp.to_crs( meta['crs'] )
mask = rasterize( list(shp.geometry), out_shape=shape, fill=0, transform=meta['transform'], all_touched=True, default_value=1, dtype=np.float32 )

output_filename = '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_19040506.tif'
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write( mask, 1 )
