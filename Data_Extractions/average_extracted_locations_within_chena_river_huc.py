# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# get avg of all 'good' stations within the CHENA RIVER HUC
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_good_stations( fn ):
	df = pd.read_csv(fn, index_col=0)
	rows,cols = df.shape
	has_75_pct = ((rows - df[~df.isna()].count()) < rows*.25)
	return df[has_75_pct[has_75_pct == True].index]

def avg_stations( df ):
	mean = df.mean(axis=1)
	return mean.to_frame(name='Chena River HUC')

def write_to_disk( df, out_fn ):
	return df.to_csv( out_fn )

def make_station_averages( fn, out_fn ):
	df = avg_stations(get_good_stations(fn))
	return write_to_disk(df, out_fn)


if __name__ == '__main__':
	import os, glob
	import pandas as pd

	# input pathing
	group = 'snotel'
	base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/data_extractions/chena_river_huc_stations_extracted'
	data_path = os.path.join( base_path, group )
	output_path = os.path.join( base_path, group, 'chena_river_huc_station_avg')

	if not os.path.exists(output_path):
		_ = os.makedirs(output_path)

	# list the csvs we want
	files = glob.glob(os.path.join(data_path, '*.csv'))
	out_files = [ os.path.join(output_path, os.path.basename(fn).replace('.csv', '_stationavg.csv')) for fn in files ]
	
	# run the averaging
	out = [ make_station_averages( fn, out_fn ) for fn,out_fn in zip(files, out_files)]

