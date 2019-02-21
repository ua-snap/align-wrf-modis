
# combine stations by variable for comparison with WRF and LST
if __name__ == '__main__':
	import os, glob
	import pandas as pd
	import numpy as np

	path = '/workspace/UA/malindgren/repos/modis_lst/ancillary/stations'
	out_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted'
	files = glob.glob( os.path.join( path, '*.csv') )

	for variable in ['tmin','tavg','tmax']:
		df_list = []
		for fn in files:
			station = os.path.basename(fn).split('_')[0]
			df = pd.read_csv(fn, index_col=0)
			sub_df = df[variable]
			sub_df.name = station
			df_list = df_list + [sub_df]

		# concat
		full_df = pd.concat(df_list, axis=1)
		full_df = full_df.replace('M', np.nan) # convert to np.nan
		out_fn = os.path.join( out_path,'{}_acis_chenariver_allstations.csv'.format(variable) )
		full_df.to_csv( out_fn )
		
