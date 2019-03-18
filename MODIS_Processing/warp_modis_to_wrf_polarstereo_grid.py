# convert MODIS LST to WRF domain/res for comparison
# Author: Michael Lindgren (malindgren@alaska.edu)
#
# [ NOT-USED ]: due to difficulties with working in WRF ref system,
# 				but this does indeed work.

def reproject( template_fn, in_fn, out_fn ):
	try:
		dirname = os.path.dirname( out_fn )
		if not os.path.exists( dirname ):
			_ = os.makedirs( dirname )
	except:
		pass
	
	_ = shutil.copy( template_fn, out_fn )
	with rasterio.open( out_fn, mode='r+' ) as tmp:
		arr = tmp.read(1)
		arr[:] = -9999
		tmp.write( arr, 1 )

	return subprocess.call(['gdalwarp','-q','-multi', '-r', 'cubicspline', '-srcnodata', '0', in_fn, out_fn ])

def run( x ):
	return reproject( *x )

if __name__ == '__main__':
	import rasterio
	import os, glob, subprocess, shutil
	import multiprocessing as mp

	base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/MODIS_DATA'
	wrf_fn = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/WRF_DATA/wrf_extent/wrf_extent_raw.tif'
	data_path = os.path.join(base_path, 'rescaled')
	out_path = os.path.join(base_path, 'regrid_to_wrf')
	
	wildcard = '_01.tif' # we only want band01 which is daily_lst

	# list the files from the akcan dir
	files = [ os.path.join(r, fn) for r,s,files in os.walk( data_path ) for fn in files if fn.endswith('.tif') and wildcard in fn ]

	# make output_filenames and paths
	output_filenames = [ fn.replace(data_path, out_path).replace('.tif','_wrf20km.tif') for fn in files ]

	args = list( zip([wrf_fn for fn in files], files, output_filenames ) )

	pool = mp.Pool( 14 )
	out = pool.map( run, args )
	pool.close()
	pool.join()
