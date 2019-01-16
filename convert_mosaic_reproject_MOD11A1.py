# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# CONVERSION STRATEGY --> MODIS. 
# --> using GDAL CLI tools
# Michael Lindgren (malindgren@alaska.edu) - September 2018
# *PRELIMINARY VERSION*
# LICENSE: MIT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def run( fn ):
	print( fn )
	# translate to GTiff
	return subprocess.call([ 'gdal_translate', '-q', '-of', 'GTiff', '-sds', '-co', 'COMPRESS=LZW', fn, os.path.basename(fn).replace('.','_').replace('_hdf','.tif') ])


if __name__ == '__main__':
	import subprocess, glob, os, shutil
	import pandas as pd
	import numpy as np
	import rasterio
	import multiprocessing as mp

	# list all files
	base_dir = '/workspace/Shared/Users/malindgren/MODIS_DATA'
	files = [ os.path.join(r,fn) for r,s,files in os.walk( os.path.join(base_dir,'raw') ) \
					for fn in files if fn.endswith('.hdf') if 'MOD11A2' not in fn and 'MYD11A2' not in fn ]

	output_path = os.path.join(base_dir, 'MOD11A1', 'converted')
	if not os.path.exists( output_path ):
		_ = os.makedirs( output_path )

	# make a tempdir to dump MOD11A2 bands
	temp_dir = os.path.join(base_dir, 'MOD11A1', 'TEMP')
	if not os.path.exists( temp_dir ):
		_ = os.makedirs( temp_dir )

	os.chdir( temp_dir )

	# parallel it
	pool = mp.Pool(16)
	out = pool.map( run, files )
	pool.close()
	pool.join()

	name_lookup = {'01':'LST_Day_1km','02':'QC_Day','03':'Day_view_time','04':'Day_view_angl','05':'LST_Night_1km','06':'QC_Night',
					'07':'Night_view_time','08':'Night_view_angl','09':'Emis_31','10':'Emis_32','11':'Clear_sky_days','12':'Clear_sky_nights'}
	
	# # delete the unwanted bands -- this is just easier since there is no the granularity in the -sds flag to deal with single bands.
	# remove = [ fn for fn in glob.glob( '*.tif' ) if '_01.tif' not in fn and '_03.tif' not in fn ]
	# _ = [ os.unlink(i) for i in remove ]
	# [ os.rename(fn, fn.replace('_'+fn.split('.')[0].split('_')[-1], '_'+name_lookup[fn.split('.')[0].split('_')[-1]])) for fn in glob.glob( '*.tif*' ) ]
	
	# move the bands we want and dont want around
	_ = [ shutil.move( fn, os.path.join( output_path, os.path.basename(fn) ) ) for fn in glob.glob( '*.tif' ) if '_01.tif' in fn ] # and '_03.tif' in fn ]
	raw_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/MOD11A1/raw_extra_bands_extracted'
	if not os.path.exists(raw_path):
		_ = os.makedirs(raw_path)

	os.system( 'mv *.tif {}'.format(raw_path) )
	
	# # # # MOSAIC -- using gdal_merge.py
	files = glob.glob( os.path.join( base_dir, 'converted', '*.tif' ) )
	colnames = ['product', 'date', 'tile', 'version', 'production', 'band']
	df = pd.DataFrame([os.path.basename(fn).split('.')[0].split('_') for fn in files ], columns=colnames )
	df['fn'] = files

	mosaicked_dir = os.path.join( base_dir, 'mosaicked' )
	if not os.path.exists( mosaicked_dir ):
		_ = os.makedirs( mosaicked_dir )

	for group_name, sub_df in df.groupby([ 'product', 'date', 'band' ]):
		if sub_df.shape[0] == 2:
			fn1,fn2 = sorted(sub_df.fn.tolist())
			out_fn = os.path.join( mosaicked_dir, '_'.join(sub_df[['product', 'date','version', 'production', 'band']].iloc[0].tolist()) ) + '.tif'
			out_fn = out_fn.replace('_006_','_InteriorAK_006_')
			out_vrt_fn = out_fn.replace('.tif', '.vrt')

			# _ = subprocess.call(['gdalbuildvrt', out_vrt_fn, fn1, fn2 ])
			_ = subprocess.call([ 'gdal_merge.py','-n','0','-a_nodata','0', '-o', out_fn, fn1, fn2 ])
		

	# # # # WARP TO EPSG:3338 AK Albers 
	files = glob.glob( os.path.join( base_dir, 'mosaicked', '*.tif' ) )
	colnames = ['product', 'date', 'tile', 'version', 'production', 'band']

	warped_dir = os.path.join( base_dir, 'warped' )
	if not os.path.exists( warped_dir ):
		_ = os.makedirs( warped_dir )

	for fn in files:
		out_fn = os.path.join( warped_dir, os.path.basename(fn) ) 
		_ = subprocess.call([ 'gdalwarp','-overwrite','-multi','-t_srs', 'EPSG:3338','-co', 'COMPRESS=LZW', fn, out_fn ])


	# # # RESCALE --- UPDATES NEEDED HERE TO DO PROPER SCALING OF VALUES FOR DIFFERENT BANDS
	# https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
	files = glob.glob( os.path.join( base_dir, 'warped', '*.tif' ) )
	rescaled_dir = os.path.join( base_dir, 'rescaled' )
	if not os.path.exists( rescaled_dir ):
		_ = os.makedirs( rescaled_dir )

	for fn in files:
		with rasterio.open( fn ) as rst:
			meta = rst.meta.copy()
			if 'affine' in meta.keys():
				meta.pop('affine')
			meta.update( compress='lzw', dtype='float32', nodata=0 )
			arr = rst.read(1)
			arr_out = np.copy(arr).astype(np.float32)
			ind = np.where( arr != meta['nodata'] )

			# scale it:
			if '_01.tif' in fn or '_05.tif' in fn:
				arr_out[ind] = arr[ind]*0.02
			
			elif '_09.tif' in fn or '_10.tif' in fn:
				arr_out[ind] = (arr[ind]*0.0020) + 0.49
			
			elif '_03.tif' in fn:
				arr_out[ind] = arr[ind]*0.1

			else:
				raise BaseException('wrong bands')

			out_fn = fn.replace( 'warped', 'rescaled' )
			with rasterio.open( out_fn, 'w', **meta ) as out:
				out.write( arr_out.astype( np.float32 ), 1 )

