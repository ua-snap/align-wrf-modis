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

    path = '/workspace/Shared/Users/malindgren/MODIS_DATA/modis_prepped_20km_3338_lanczos'
    out_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted'

    # get the stations locations
    stations = gpd.read_file('/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_station_ids.shp')
    points = {''.join([j for j in i if j.isalnum()]):{'geometry':Point(lat,lon)} for i,lat,lon in zip(stations.name,stations.lon, stations.lat)}
    shp = gpd.GeoDataFrame(points, crs={'init':'epsg:4326'}).T
    shp = shp.to_crs(epsg=3338)
    points = shp.geometry.apply( lambda x: (x.x, x.y) ).to_dict()

    sensors = ['MOD11A2']
    variable = 'lst'

    for sensor in sensors:
        # list the files:
        files = sorted(glob.glob(os.path.join( path, '{}*.tif'.format(sensor,))))

        # fun with dates!
        dates = [ os.path.basename(fn).split('_')[1][1:] for fn in files ]
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
        out_fn = os.path.join( out_path, '{}_{}_allstations.csv'.format(sensor, variable) )
        df = pd.DataFrame( out_vals ) - 273.15
        df = df.replace(-10272.150000, np.nan) # this is an odd oob value but I am removing it anyway.
        df.to_csv( out_fn )
