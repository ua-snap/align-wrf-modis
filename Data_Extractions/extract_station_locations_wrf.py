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

    variable = 't2'
    group = 'snotel' # input locations should be previously prepped with ancillary prep script
    base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
    path = os.path.join( base_path, 'WRF_DATA','WRF_Day_hours','{}_20km_3338_lanczos'.format(variable) )
    out_path = os.path.join(base_path, 'data_extractions','chena_river_huc_stations_extracted',group)

    if not os.path.exists(out_path):
        _ = os.makedirs(out_path)

    metrics = ['min','mean','max']
    var_lookup = {'min':'tmin','mean':'tavg','max':'tmax'}

    # get the stations locations
    stations = gpd.read_file(os.path.join( base_path, 'STATION_DATA', group, '{}_points_chena_river_huc.shp'.format(group)))
    points = {''.join([j for j in i if j.isalnum()]):{'geometry':Point(lat,lon)} 
                for i,lat,lon in zip(stations.name,stations.lon,stations.lat)}
    shp = gpd.GeoDataFrame(points, crs={'init':'epsg:4326'}).T
    shp = shp.to_crs(epsg=3338)
    points = shp.geometry.apply( lambda x: (x.x, x.y) ).to_dict()

    sensors = ['MOD11A2','MYD11A2' ]
    model_groups = ['NCAR-CCSM4_rcp85','GFDL-CM3_rcp85', 'ERA-Interim_historical']

    for sensor in sensors:
        for metric in metrics:
            for model_group in model_groups:
                # list the files:
                files = sorted(glob.glob(os.path.join( path, '*{}*{}*{}*.tif'.format(model_group,metric,sensor))))
                
                # fun with dates!
                dates = [ os.path.basename(fn).split('_')[-2] for fn in files ]
                datetimes = pd.DatetimeIndex([ dt.strptime(d, '%Y%j') for d in dates ])

                # stack the files
                arr = np.array([ open_raster(fn, band=1) for fn in files ])

                # get metadata
                with rasterio.open( files[0] ) as tmp:
                    a = tmp.transform

                out_vals = {}
                for station in points:
                    col,row = list(np.array(~a*points[station]).round(0).astype(int))
                    vals = arr[:,row,col]
                    out_vals[station] = dict(zip(dates,vals.tolist()))

                # write to dataframe
                out_fn = os.path.join( out_path, '{}_{}_{}_{}_wrf_{}_chena_river_huc_stations.csv'.format(variable, model_group, sensor, metric, group) )
                df = pd.DataFrame( out_vals ) - 273.15
                df.index = datetimes
                df.to_csv( out_fn )
