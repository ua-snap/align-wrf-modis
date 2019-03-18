# resample to SNAP-prepped MODIS (EPSG:3338) to 20km resolution.
# Author: Michael Lindgren (malindgren@alaska.edu)

def resample_prepped_modis_to_20km( fn, out_fn ):
	print(fn)
	_ = subprocess.call(['gdalwarp', '-q','-multi', '-overwrite', '-tr', '20000', '20000', '-t_srs', 'EPSG:3338', 
						'-overwrite', '-srcnodata', '0', '-dstnodata', '-9999', '-r', 'lanczos', fn, out_fn ])
	return out_fn

def wrapper( x ):
	return resample_prepped_modis_to_20km( *x )

if __name__ == '__main__':
	import os, glob, subprocess
	import pandas as pd
	import multiprocessing as mp

	base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/MODIS_DATA'
	modis_path = os.path.join(base_path, 'rescaled')
	out_path = os.path.join( base_path, 'modis_regrid_20km_3338')

	if not os.path.exists( out_path ):
		_ = os.makedirs( out_path )

	for sensor in ['MOD11A2','MYD11A2']:
		# list and sort the files chronologically
		modis_files = glob.glob( os.path.join( modis_path, '{}*_01.tif'.format(sensor) ) )
		dates = [ int(os.path.basename(fn).split('.tif')[0].split('_')[1].replace('A','')) for fn in modis_files ]
		modis_files = pd.DataFrame({'fn':modis_files,'date':dates}).fn.tolist()
		out_files = [ os.path.join(out_path, os.path.basename(fn).replace('.tif','_20km_3338.tif')) for fn in modis_files ]
		
		pool = mp.Pool( 32 )
		pool.map( wrapper, zip(modis_files, out_files ) )
		pool.close()
		pool.join()