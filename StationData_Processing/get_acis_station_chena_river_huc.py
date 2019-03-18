# # # # # # # # # # # # # # # # # # # # 
# GET ACIS Station Data:
# # # # # # # # # # # # # # # # # # # # 

def open_csv_change_colname_to_station( fn, metric ):
	station_name = os.path.basename(fn).split('.')[0].split('_')[0]
	df = pd.read_csv(fn, index_col=0)[metric]
	df.name = station_name
	return df

if __name__ == '__main__':
	import json, os, requests
	import pandas as pd
	import geopandas as gpd
	import numpy as np
	from shapely.geometry import Point

	# # figure out the station names that we need. -- overlap the Chena River HUC
	# stations = pd.read_csv('/workspace/UA/malindgren/repos/modis_lst/ancillary/ghcnd-stations_fixed.csv')
	# stations['geometry'] = stations.apply( lambda x: Point(x.lat,x.lon), axis=1 )
	# shp = gpd.GeoDataFrame(stations, crs={'init':'epsg:4326'}, geometry='geometry')

	# huc = gpd.read_file( '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_19040506.shp' )
	# huc = huc.to_crs( shp.crs )
	# over = shp.intersects( huc )
	base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
	acis_path = os.path.join(base_path,'STATION_DATA','acis')
	stations = gpd.read_file(os.path.join(acis_path, 'chena_river_huc_station_ids.shp'))
	output_path = os.path.join(acis_path, 'raw')

	if not os.path.exists(output_path):
		_ = os.makedirs(output_path)

	ids = [ int(''.join([ j for j in i if j.isnumeric()])) for i in stations.StationID ]

	# lookup key/vals for REST-like access
	out_names = []
	for idx,sid in enumerate(ids):
		sdate = '2002-01-01'
		edate = '2019-03-14' # pd.Timestamp.today().strftime('%Y-%m-%d')
		elems = 'mint,avgt,maxt'
		output_type = 'json'

		url = 'http://data.rcc-acis.org/StnData?sid={}&sdate={}&edate={}&elems={}&output={}'.format(sid, sdate, edate, elems, output_type )

		# get the station data json
		response = requests.get(url)
		json_data = json.loads(response.text)
		
		if list(json_data.keys())[0] != 'error':
			try:
				# make a dataframe from the returned json 
				df = pd.DataFrame( json_data['data'] )
				df.columns = ['index', 'tmin', 'tavg', 'tmax']
				df.index = df['index']
				df = df.drop('index', axis=1 )
				name = ''.join([i for i in stations.iloc[idx]['name'] if i.isalnum()])

				# fill missing vals... and make trace really small
				df = df.replace('M', np.nan).replace('T', '.01').astype(np.float32)
				
				# make celcius
				df = ((df - 32.0) * 5.0 / 9.0).round(2)

				out_fn = os.path.join( output_path, name+'_acis_stationdata.csv' )
				df.to_csv( out_fn, index=True )
				out_names = out_names + [out_fn]
			except:
				print(sid)
				pass


	# now lets make a single file out of each of the downloaded metrics across all stations
	stacked_path = os.path.join(acis_path, 'prepped')
	if not os.path.exists( stacked_path ):
		_ = os.makedirs( stacked_path )

	for metric in ['tmin', 'tavg', 'tmax']:
		final = pd.concat([ open_csv_change_colname_to_station( fn, metric ) for fn in out_names ], axis=1)
		final = final.round(2)
		final.to_csv( os.path.join(stacked_path, '{}_acis_stationdata.csv'.format(metric)) )


