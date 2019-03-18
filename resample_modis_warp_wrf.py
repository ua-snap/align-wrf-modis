# [Fish & Fire] make the data more accessible for R-workflow used for this project.
# ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ ------ 

# # # # # TEMPORARY NOTES:
# # how to run this properly...
# 1. run `resample_wrf_tsk_dayhours_to_8day_mean_min_max.py`
# 	- for both directories...  GFDL-CM3/ ERA-Interim AND NCAR-CCSM4
# 2. run the top script from this file
# 3. run the bottom script (separately) from this file


# # resample to MODIS prepped (EPSG:3338) to 20km resolution.
# def resample_prepped_modis_to_20km( fn, out_fn ):
# 	print(fn)
# 	_ = subprocess.call(['gdalwarp', '-q','-multi', '-overwrite', '-tr', '20000', '20000', '-t_srs', 'EPSG:3338', 
# 						'-overwrite', '-srcnodata', '0', '-dstnodata', '-9999', '-r', 'lanczos', fn, out_fn ])
# 	return out_fn

# def wrapper( x ):
# 	return resample_prepped_modis_to_20km( *x )

# if __name__ == '__main__':
# 	import os, glob, subprocess
# 	import pandas as pd
# 	import multiprocessing as mp

# 	modis_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/rescaled'
# 	out_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos'
# 	for sensor in ['MOD11A2','MYD11A2']:
# 		# list and sort the files chronologically
# 		modis_files = glob.glob( os.path.join( modis_path, '{}*_01.tif'.format(sensor) ) )
# 		dates = [ int(os.path.basename(fn).split('.tif')[0].split('_')[1].replace('A','')) for fn in modis_files ]
# 		modis_files = pd.DataFrame({'fn':modis_files,'date':dates}).fn.tolist()
# 		out_files = [ os.path.join(out_path, os.path.basename(fn).replace('.tif','_20km_3338.tif')) for fn in modis_files ]
		
# 		pool = mp.Pool( 32 )
# 		pool.map( wrapper, zip(modis_files, out_files ) )
# 		pool.close()
# 		pool.join()

# --------------------------------------------------------------------------------------------------------------------------------
# reproject the WRF TSK resampled to 8-day means, to the prepped 20km MODIS extent
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

def reproject_wrf_to_modis_20km( fn, out_fn ):
	print(fn)
	_ = subprocess.call(['gdalwarp', '-q', '-r', 'lanczos', fn, out_fn ])
	with rasterio.open( out_fn, mode='r+' ) as out:
		arr = out.read( 1 )
		arr[arr == 0] = -9999
		out.write(arr, 1)
	return out_fn

def run_reproject( arr, meta, out_fn, template_fn ):
	import copy
	fn = make_gtiff( arr, meta, out_fn )
	out_fn = copy.copy(fn).replace('.tif', '_3338.tif')
	with rasterio.open( template_fn ) as tmp:
		arr = np.empty_like(tmp.read())
		tmp_meta = tmp.meta
		tmp_meta.update( compress='lzw' )
		with rasterio.open( out_fn, 'w', **tmp_meta ) as rst:
			rst.write( arr )
		_ = reproject_wrf_to_modis_20km( fn, out_fn )
		os.unlink( fn )
	return out_fn

def wrapper( x ):
	return run_reproject( *x )

if __name__ == '__main__':
	import os, glob, rasterio, subprocess, itertools
	import xarray as xr
	import numpy as np
	import multiprocessing as mp

	variable = 't2'
	base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
	wrf_path = os.path.join(base_path,'WRF_DATA','WRF_Day_hours', variable)
	out_path = os.path.join(base_path, 'WRF_DATA','WRF_Day_hours','{}_20km_3338_lanczos'.format(variable))
	template_fn = '/workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos/MOD11A2_A2005257_InteriorAK_006_2015261161547_01_20km_3338.tif'
	groups = [ 'ERA-Interim_historical','GFDL-CM3_historical','GFDL-CM3_rcp85','NCAR-CCSM4_historical','NCAR-CCSM4_rcp85',]
	sensors = ['MOD11A2'] #'MYD11A2'] # # ,'MYD11A2'
	metrics = ['min','max','mean']

	for group, sensor, metric in itertools.product(groups, sensors, metrics):
		wrf_fn, = glob.glob( os.path.join( wrf_path, '*{}_{}_{}.nc'.format(group,metric,sensor) ) )
		ds = xr.open_dataset( wrf_fn ).load()
		
		da = ds[variable]
		bands,rows,cols = da.shape
		res = (20000.0,20000.0)
		transform = rasterio.transform.from_origin( ds.xc.data.min()-(res[0]/2.), ds.yc.data.max()+(res[0]/2.), res[0], res[1] )
		proj4 = '+a=6370000 +b=6370000 +k=1 +lat_0=90 +lat_ts=64 +lon_0=-152 +no_defs +proj=stere +units=m +x_0=0 +y_0=0'
		meta = {'count': 1,
				'crs': rasterio.crs.CRS.from_string(proj4),
				'driver': 'GTiff',
				'dtype': 'float32',
				'height': rows,
				'nodata': -9999,
				'transform': transform,
				'width': cols}

		# make output_filenames for wrf to gtiff
		times = ds.time.to_index()
		jdates = list(times.map(lambda x: x.strftime('%Y%j')))
		out_files = [os.path.join( out_path,'{}_8Day_daytime_wrf_{}_{}_{}_{}.tif'.format(variable,group,metric,sensor,jdate) ) for jdate in jdates]

		# make args
		args = [(a, meta, out_fn, template_fn) for a, out_fn in zip(list(da.values),out_files) ]

		# serial processing
		out = [ wrapper(x) for x in args ]
		
		# # multicore process is giving errors, but this is how it is deployed.
		# pool = mp.Pool( 25 )
		# out = pool.map( wrapper, args )
		# pool.close()
		# pool.join()

		# now put these new GTiffs into a large multiband file...
		def open_raster( fn, band=1 ):
			with rasterio.open( fn ) as rst:
				arr = rst.read( band ).copy()
			rst = None
			return arr
		
		newfiles = sorted( glob.glob( os.path.join( out_path, '*.tif' ) ) )
		new_arr = np.array([ open_raster(fn) for fn in newfiles ])

		with rasterio.open( template_fn ) as tmp:
			new_meta = tmp.meta.copy()
			new_meta.update(compress='lzw', count=new_arr.shape[0])

		# os.mkdir(out_path+'_multiband')
		with rasterio.open( os.path.join( out_path+'_multiband','{}_8Day_daytime_wrf_{}_{}_{}_3338_multiband.tif'.format(variable, group, metric, sensor)), 'w', **new_meta ) as out:
			out.write( new_arr )

		tmp = None
		out = None
		del new_arr, tmp, out


# # zip em
# cd /workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos
# zip -rj5 /workspace/Shared/Users/malindgren/MODIS_DATA/zips_for_stephen/lst_20km_3338_lanczos.zip *.tif

# cd /workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_20km_3338_lanczos
# zip -rj5 /workspace/Shared/Users/malindgren/MODIS_DATA/zips_for_stephen/tsk_20km_3338_lanczos.zip *.tif

