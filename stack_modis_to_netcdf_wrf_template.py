def open_raster( fn, band=1 ):
    with rasterio.open( fn ) as rst:
        arr = rst.read( band )
    return arr

if __name__ == '__main__':
    import os, glob, datetime
    import numpy as np
    import pandas as pd
    import xarray as xr
    import rasterio
    import multiprocessing as mp

    # get the files
    template_nc = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/daily/tsk/tsk_daily_wrf_ERA-Interim_historical_2000.nc'
    tmp_ds = xr.open_dataset( template_nc )
    path = '/workspace/Shared/Users/malindgren/MODIS_DATA/regrid_20km_wrf'
    files = glob.glob( os.path.join( path,'*.tif' ) )
    ncpus = 64
    
    # make a dataframe and group / sort
    elem_names =['sensor','date','region','version','process_date','band']
    df = pd.DataFrame([ dict(zip(elem_names, os.path.basename(fn).split('.')[0].split('_'))) for fn in files ])
    df['fn'] = files
    df['date'] = df.date.apply(lambda x: x.replace('A','')).astype(int)
    df = df.groupby('sensor').apply(lambda x: x.sort_values('date'))
    
    for idx,sensor in enumerate(df.sensor.unique()):
        fn_list = list(df.loc[sensor]['fn'])

        # stack to 3D array
        pool = mp.Pool( ncpus )
        arr = np.array( pool.map( open_raster, fn_list ) )
        pool.close()
        pool.join()

        # get the coordinates -- from the template WRF file since it is the same.
        xc, yc = tmp_ds.xc, tmp_ds.yc
        lon, lat = tmp_ds.lon, tmp_ds.lat

        # make a time coordinate
        dates = pd.DatetimeIndex(df.loc[sensor].date.astype(str).apply( lambda x: datetime.datetime.strptime(x, '%Y%j') ).tolist() )

        # make the new LST NetCDF dataset
        new_ds = xr.Dataset({'lst': (['time','yc', 'xc'], arr)},
                            coords={'xc': ('xc', xc.values),
                                    'yc': ('yc', yc.values),
                                    'lon':(['yc','xc'], lon.values ),
                                    'lat':(['yc','xc'], lat.values ),
                                    'time':dates })

        new_ds.lst.attrs = {'_FillValue':-9999}

        # add compression encoding
        encoding = new_ds[ 'lst' ].encoding
        encoding.update( zlib=True, complevel=5, contiguous=False, chunksizes=None, dtype='float32' )
        new_ds[ 'lst' ].encoding = encoding
        
        # output_filenaming-fu
        elems = df.loc[sensor].iloc[0]
        begin_date = str(elems['date'])
        end_date = str(df.loc[sensor].iloc[-1].date)
        elems = elems.loc[['sensor','region','version','band']].tolist()
        elems = ['lst'] + elems + ['{}-{}'.format(begin_date, end_date)]
        output_filename = os.path.join('/workspace/Shared/Users/malindgren/MODIS_DATA/regrid_20km_wrf_netcdf', '_'.join(elems)+'.nc')
        
        # write to disk
        new_ds.to_netcdf( output_filename )