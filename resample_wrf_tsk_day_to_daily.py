# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# resample the WRF skin-temp hourlies to daily day-times and then to 8-Day averages for comparison with MODIS
# Author: Michael Lindgren (malindgren@alaska.edu)
# License: MIT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_day_hours(month):
    return (month >= 9) & (month <= 17)

def modisify_wrf( files, output_path, sensor ):
    ''' make the WRF data temporally look like MODIS LST 8-Day composite Daytime'''
    ds = xr.open_mfdataset( files )
    fn = '_'.join(files[0].split('_')[:-1])+'_{}.nc'.format(sensor)
    output_filename = os.path.join( output_path, os.path.basename(fn).replace('hourly','8Day_daytime') )
    print( output_filename )

    with xr.open_mfdataset( files ) as ds:
        # get the hours we are interested in
        selected_hours = ds.sel( time=get_day_hours(ds['time.hour']) )
        
        # resample to daily the selected hours...
        new_ds = selected_hours.resample( time='D' ).mean()

        # make an 8-day average for comparison with MODIS
        new_mean = new_ds.resample( time='8D' ).mean()
    
    dirname = os.path.dirname( output_filename )
    if not os.path.exists( dirname ):
        _ = os.makedirs( dirname )

    # write to disk and cleanup
    new_mean.to_netcdf( output_filename )
    del new_mean, new_ds
    return output_filename


if __name__ == '__main__':
    import os, glob
    import numpy as np
    import xarray as xr
    import dask
    import pandas as pd
    import datetime
    import multiprocessing as mp
    from functools import partial

    path = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix'
    output_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
    variable = 'tsk'

    # mess with the begin / end dates for the MODIS data for comparison
    sensor_dates = {'MOD11A2': (2000049,2018233), 'MYD11A2':(2002185,2018233)}
    sensors = ['MOD11A2','MYD11A2']
    
    for sensor in sensors:
        # more date-fu
        begin,end = sensor_dates[ sensor ]
        begin_dt = datetime.datetime( int(str(begin)[:4]),1,1 ) + datetime.timedelta( int(str(begin)[4:]) )
        end_dt = datetime.datetime( int(str(end)[:4]),1,1 ) + datetime.timedelta( int(str(end)[4:]) )
        begin = begin_dt.strftime('%Y-%m-%d')
        end = end_dt.strftime('%Y-%m-%d')

        # list files we want within the range we want...
        files = [ fn for fn in glob.glob( os.path.join( path, variable, '*.nc' ) ) 
                    if int(fn.split('_')[-1].split('.')[0]) in list(range(begin_dt.year, end_dt.year+1))]

        # make a dataframe of the files
        elems = [ 'variable', 'agg', 'wrf', 'model', 'scenario', 'year' ]
        df = pd.DataFrame([ dict(zip(elems, os.path.splitext(os.path.basename(fn))[0].split('_'))) for fn in files ])
        df['fn'] = files

        # group and sort the files 
        grouped_files = df.sort_values('year').groupby(['model','scenario']).apply(lambda x: x.fn.tolist()).tolist()
        
        # run processing
        f = partial( modisify_wrf, output_path=output_path, sensor=sensor )
        pool = mp.Pool( 3 )
        done = pool.map( f, grouped_files )
        pool.close()
        pool.join()
