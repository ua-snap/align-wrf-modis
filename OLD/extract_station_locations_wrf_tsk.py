# extract locations from WRF TSK
def open_raster( fn, band=1 ):
    with rasterio.open(fn) as rst:
        arr = rst.read(band)
    return arr


if __name__ == '__main__':
    import affine
    import os, glob
    import pandas as pd
    import numpy as np
    import geopandas as gpd
    from shapely.geometry import Point
    import rasterio
    from datetime import datetime as dt

    group = 'snotel' # input locations should be previously prepped with ancillary prep script
    path = '/workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos'
    out_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted/{}'.format(group)

    if not os.path.exists(out_path):
        _ = os.path.makedirs(out_path)

    # path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_20km_3338_lanczos'
    # out_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted/snotel'
    variables = ['min','mean','max']
    var_lookup = {'min':'tmin','mean':'tavg','max':'tmax'}

    # get the stations locations
    stations = gpd.read_file('/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/shapefiles/{}_points_chena_river_huc.shp'.format(group))
    points = {''.join([j for j in i if j.isalnum()]):{'geometry':Point(lat,lon)} 
                for i,lat,lon in zip(stations.name,stations.lon, stations.lat)}
    shp = gpd.GeoDataFrame(points, crs={'init':'epsg:4326'}).T
    shp = shp.to_crs(epsg=3338)
    points = shp.geometry.apply( lambda x: (x.x, x.y) ).to_dict()

    sensors = ['MOD11A2','MYD11A2' ]
    model_groups = ['NCAR-CCSM4_rcp85','GFDL-CM3_rcp85', 'ERA-Interim_historical']

    for sensor in sensors:
        for variable in variables:
            for model_group in model_groups:
                # list the files:
                files = sorted(glob.glob(os.path.join( path, '*{}*{}*{}*.tif'.format(model_group,variable,sensor,))))

                # fun with dates!
                dates = [ os.path.basename(fn).split('_')[-2] for fn in files ]
                datetimes = [ dt.strptime(dates[0], '%Y%j') for d in dates ]

                # stack the files
                arr = np.array([ open_raster(fn,band=1) for fn in files ])

                # get metadata
                with rasterio.open( files[0] ) as tmp:
                    a = tmp.transform

                out_vals = {}
                for station in points:
                    col,row = list(np.array(~a*points[station]).round(0).astype(int))
                    vals = arr[:,row,col]
                    out_vals[station] = dict(zip(dates,vals.tolist()))

                # write to dataframe
                out_fn = os.path.join( out_path, '{}_{}_{}_wrf_tsk_allstations.csv'.format(model_group, sensor, var_lookup[variable]) )
                df = pd.DataFrame( out_vals ) - 273.15
                df.to_csv( out_fn )
