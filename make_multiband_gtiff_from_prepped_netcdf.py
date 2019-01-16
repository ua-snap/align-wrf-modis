# restructure the formatted 20km data to GTiff and reproject to 3338

def make_gtiff( arr, meta, out_fn ):
	shape = arr.shape
	if len(shape) == 2:
		count = 1
		arr = arr[np.newaxis,...]
	else:
		count = shape[0]
	meta.update( count=count ) 
	with rasterio.open(out_fn,'w',**meta) as out:
		out.write( arr )
	return out_fn

if __name__ == '__main__':
	import os, rasterio, glob, subprocess
	import numpy as np
	import pandas as pd
	from affine import Affine
	import xarray as xr
	from rasterio.warp import reproject, Resampling

	wrf_prepped_dir = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
	modis_prepped_dir = '/workspace/Shared/Users/malindgren/MODIS_DATA/regrid_20km_wrf_netcdf'

	template_fn = '/workspace/Shared/Users/malindgren/MODIS_DATA/regrid_20km_wrf/MOD11A2_A2011137_InteriorAK_006_2016055185500_01_wrf20km.tif'

	with rasterio.open( template_fn ) as tmp:
		meta = tmp.meta.copy()
		meta.update(compress='lzw')
		# add in the counts of the number of files here as well

	# open the netcdf we want
	variable = 'tsk'
	ds = xr.open_dataset(os.path.join( wrf_prepped_dir,'tsk_8Day_daytime_wrf_ERA-Interim_historical_MOD11A2.nc' ))
	da = ds[variable]
	proj4 = '+a=6370000 +b=6370000 +k=1 +lat_0=90 +lat_ts=64 +lon_0=-152 +no_defs +proj=stere +units=m +x_0=0 +y_0=0'

	arr = da.values

	out_fn = '/workspace/Shared/Users/malindgren/MODIS_DATA/GTiff_multiband/tsk_8Day_daytime_wrf_ERA-Interim_historical_MOD11A2.tif'
	out_fn_proj = '/workspace/Shared/Users/malindgren/MODIS_DATA/GTiff_multiband/tsk_8Day_daytime_wrf_ERA-Interim_historical_MOD11A2_EPSG3338.tif'

	_ = make_gtiff( arr[:3,...], meta, out_fn2 )
	_ = subprocess.call([ 'gdalwarp', '-overwrite', '-multi', '-s_srs', proj4,'-t_srs','EPSG:3338', '-co', 'COMPRESS=LZW', '-r', 'cubicspline', out_fn, out_fn_proj ])


	


