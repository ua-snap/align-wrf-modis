
if __name__ == '__main__':
	import os
	import pandas as pd
	import geopandas as gpd
	from shapely.geometry import Point
	import numpy as np

	fn = '/Users/malindgren/Documents/repos/modis_lst/ancillary/SNOTEL_points_all.csv'
	# output_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted/snotel'
	output_path = '/Users/malindgren/Documents/TEMP'

	df = pd.read_csv( fn )
	df['geometry'] = df.apply( lambda x: Point(x.Longitude, x.Latitude), axis=1 )
	df['ID'] = df.ID.apply(lambda x: x.replace('\t', '').replace(',', ''))
	df['Name'] = [''.join([j for j in i if j.isalnum()]) for i in df.Name ]

	# subset to remove some junk
	df = df[['Name','State' ,'ID', 'Longitude', 'Latitude', 'geometry']]

	pts = gpd.GeoDataFrame(df, geometry='geometry')
	pts.to_file('/Users/malindgren/Documents/repos/modis_lst/ancillary/SNOTEL_points_all.shp')
	
	huc_fn = '/Users/malindgren/Documents/repos/modis_lst/ancillary/chena_river_huc_19040506.shp'
	huc = gpd.read_file( huc_fn )
	huc = huc.to_crs( epsg=4326 )

	# now find which ones intersect
	ak_pts = pts[pts.geometry.apply(lambda x: x.intersects(huc.geometry[0]))]
	ak_pts.to_file('/Users/malindgren/Documents/repos/modis_lst/ancillary/SNOTEL_points_chena_river_huc.shp')

	# get the data
	import requests
	station_files = []
	for station_id in ak_pts.ID:
		# get the station name from the ID
		station_name = ak_pts[ak_pts.ID == station_id].Name
		station_name = ''.join([i for i in station_name if i.isalnum()])

		# download the data from the server
		url = 'https://wcc.sc.egov.usda.gov/reportGenerator/view_csv/customGroupByMonthReport/daily/{0}:AK:SNTL%7Cid=%22%22%7Cname/POR_BEGIN,POR_END/TMIN::value,TMAX::value,TAVG::value'.format(station_id)
		response = requests.get(url)
		out_fn = os.path.join(output_path, 'SNOTEL', 'raw', '{}.csv'.format(str(station_id)))
		with open(out_fn, 'w') as f:
			f.write( response.text )
		station_files = station_files + [out_fn]

	# NOW DO STUFF TO THE DOWNLOADED DATA TO MAKE IT CLEANER
	# separate variables and restructure for processing similar to the ACIS data
	# initial column-fu to deal with the nonsensical data format they are passing.
	colnames = ['Year','Day','Oct-TMIN','Oct-TMAX','Oct-TAVG','Nov-TMIN','Nov-TMAX','Nov-TAVG','Dec-TMIN','Dec-TMAX','Dec-TAVG','Jan-TMIN','Jan-TMAX',
				'Jan-TAVG','Feb-TMIN','Feb-TMAX','Feb-TAVG','Mar-TMIN','Mar-TMAX','Mar-TAVG','Apr-TMIN','Apr-TMAX','Apr-TAVG','May-TMIN','May-TMAX',
				'May-TAVG','Jun-TMIN','Jun-TMAX','Jun-TAVG','Jul-TMIN','Jul-TMAX','Jul-TAVG','Aug-TMIN','Aug-TMAX','Aug-TAVG','Sep-TMIN','Sep-TMAX','Sep-TAVG']

	for metric in ['TMAX', 'TMIN', 'TAVG']:
		out_data = []
		for station_fn in station_files:
			station_id = os.path.basename(station_fn).split('.')[0]
			station_name = ak_pts[ak_pts.ID == station_id].Name
			station_name = ''.join([i for i in station_name if i.isalnum()])
			print(station_name)
			df = pd.read_csv(station_fn, skiprows=56, header=None, names=colnames )
			
			cols = ['Year','Day'] + [ col for col in colnames if metric in col ]
			sub_df = df[ cols ]
			sub_df.columns = [i.split('-')[0] for i in sub_df.columns]
			
			# melt the dataframe to long-form
			melted = sub_df.melt(id_vars=['Year','Day'])
			melted.columns = ['Year', 'Day', 'Month', metric]
			
			# deal with month-y stuff
			month_days = {'Oct':31, 'Nov':30, 'Dec':31, 'Jan':31, 'Feb':28, 'Mar':31, 'Apr':30, 'May':31, 'Jun':30,'Jul':31, 'Aug':31, 'Sep':30}
			month_ids = {'Oct':10, 'Nov':11, 'Dec':12, 'Jan':1, 'Feb':2, 'Mar':3, 'Apr':4, 'May':5, 'Jun':6,'Jul':7, 'Aug':8, 'Sep':9}
			month_names = {v:k for k,v in month_ids.items()} # flip for inverse lookup

			# alter the month name to the month chronological id (1-12)
			months = melted.Month
			new_months = []
			for month in months:
				month_id = month_ids[month]
				new_months = new_months + [month_id]

			# update to numerical month
			melted.Month = new_months
			
			output_frames = []
			for idx, sub_df in list(melted.groupby(['Year', 'Month'])):
				# make a timestamp
				def make_timestamp(x):
					try:
						return pd.Timestamp('{}-{}-{}'.format(int(x.Month), int(x.Day), int(x.Year)))
					except:
						pass

				new_index = pd.DatetimeIndex( sub_df.apply( lambda x:make_timestamp(x) , axis=1 ) )
				# new_index = new_index[~pd.isna(new_index)]
				# drops errant extra days due to the asinine way they stored the data
				sub_df = sub_df[ metric ].to_frame(name=station_name)
				sub_df.index = new_index
				sub_df = sub_df[~pd.isna(sub_df.index)]
				output_frames = output_frames + [sub_df]
			station_df = pd.concat( output_frames )
			out_data = out_data + [station_df]
		# write to csv
		out_df = pd.concat(out_data, axis=1)
		# rescale the temps to celcius from fahrenheit
		out_df = ((out_df - 32.0) * 5.0 / 9.0).round(1) # make celcius
		out_fn = os.path.join(output_path, 'SNOTEL', 'prepped','SNOTEL_points_chena_river_huc_{}.csv'.format(metric.lower()) )
		out_df.to_csv( out_fn )
		

