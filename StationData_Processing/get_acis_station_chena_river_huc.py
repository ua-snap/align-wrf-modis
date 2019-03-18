# # # # # # # # # # # # # # # # # # # # 
# GET ACIS Station Data:
# # # # # # # # # # # # # # # # # # # # 


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


# lookup key/vals for REST Service
out = {}
for idx,sid in enumerate(ids):
	sdate = '2002-01-01'
	edate = '2018-12-31'
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
			
			out_fn = os.path.join( output_path, name+'_acis_stationdata.csv' )
			df.to_csv( out_fn, index=True )
		except:
			print(sid)
			pass

