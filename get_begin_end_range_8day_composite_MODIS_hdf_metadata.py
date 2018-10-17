# pull range begin/end dates/times from the raw HDF5 files...
def get_ranges( fn ):
	elems = ['product', 'date', 'location', 'version', 'production_date']
	fn_dict = dict(zip(elems, os.path.basename(fn).split('.hdf')[0].split('.')))
	fn_dict['date'] = fn_dict['date'].replace('A','')
	fn_dict['fn'] = fn
	f = gdal.Open( fn )
	meta = f.GetMetadata_Dict()
	out = {'RANGEBEGINNINGDATE':meta['RANGEBEGINNINGDATE'],
	'RANGEBEGINNINGTIME':meta['RANGEBEGINNINGTIME'],
	'RANGEENDINGDATE':meta['RANGEENDINGDATE'],
	'RANGEENDINGTIME':meta['RANGEENDINGTIME']}
	out.update( fn_dict )
	f = None
	return out

if __name__ == '__main__':
	import os
	import gdal
	import multiprocessing as mp
	import pandas as pd
	import numpy as np

	# args
	base_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/raw'
	files = [ os.path.join(r,fn) for r,s,files in os.walk(base_path) for fn in files if fn.endswith('.hdf') and '11A2' in fn ]

	# run it
	pool = mp.Pool( 15 )
	out = pool.map( get_ranges, files )
	pool.close()
	pool.join()

	# make a table out of it and dump to disk...
	df = pd.DataFrame( out )
	df.to_csv( '/workspace/Shared/Users/malindgren/MODIS_DATA/ancillary/MODIS_LST_8dayComposite_begin_end_range_metadata.csv' )

	# subset it since we have posited that they use the same dates in each tile agg
	new_df = df.groupby(['product','date']).apply(lambda x: x.iloc[0].drop('location'))
	new_df.to_csv( '/workspace/Shared/Users/malindgren/MODIS_DATA/ancillary/MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv' )


	