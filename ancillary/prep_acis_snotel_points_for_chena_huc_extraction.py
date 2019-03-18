# prep the points layers for extractions from MODIS / WRF data
import geopandas as gpd

out_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/shapefiles'

# prep the ACIS points in the HUC we are interested in:
fn = '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_station_ids.shp'

df = gpd.read_file( fn )
df = df[['lat', 'lon', 'name', 'geometry']]

df.to_file(os.path.join(out_path, 'acis_points_chena_river_huc.shp'))

# prep the SNOTEL points in the HUC we are interested in:
fn = '/workspace/UA/malindgren/repos/modis_lst/ancillary/SNOTEL_points_chena_river_huc.shp'
df = gpd.read_file( fn )
df = df[['Latitude', 'Longitude', 'Name', 'geometry']]

df.columns = ['lat', 'lon', 'name', 'geometry']
df.to_file(os.path.join(out_path, 'snotel_points_chena_river_huc.shp'))
