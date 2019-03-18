# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# resample the WRF skin-temp hourlies to daily day-times and then to 8-Day averages for comparison with MODIS
# Author: Michael Lindgren (malindgren@alaska.edu)
# License: MIT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_day_hours( month ):
    return (month >= 9) & (month <= 17)

def get_year_from_fn( fn ):
    return int(os.path.basename(fn).split('.nc')[0].split('_')[-1])

def fill_leap( fn, out_fn, variable ):
    with xr.open_dataset( fn ) as ds:
        ds_filled = ds.resample(time='1H').interpolate('slinear')
        ds_filled[variable].encoding = ds[variable].encoding
        ds_filled.to_netcdf(out_fn)
        ds_filled.close(); ds_filled = None
    return out_fn

def interpolate_leap_years( fn, year, variable, tmp_path='/atlas_scratch/malindgren/TMP' ):
    # make an output path to dump the temporary intermediate file
    out_fn = os.path.join( tmp_path, os.path.basename(fn) )
    if os.path.exists( out_fn ):
        _ = os.unlink( out_fn )
    if calendar.isleap( year ):
        # fill it with a spline
        _ = fill_leap( fn, out_fn, variable )
    else:
        # copy the file to the temp location for speed the symlinks were very slow.
        _ = shutil.copy( fn, out_fn )
    return out_fn

def run_interpolate_leap_years(x):
    return interpolate_leap_years(*x)

def make_wrf_like_modis( files, meta_df, metric, product, variable ):
    ''' resample (temporal) WRF hourlies to the MODIS ranges used in 8-day compositing. '''
    with xr.open_mfdataset( files ) as ds:
        ds_sel = ds.sel( time=get_day_hours(ds['time.hour']) ).compute()
    
    # sort the values -- probably unnecessary, but is here.
    meta_df = meta_df.sort_values('date')

    # make the slices we need from the metadata file
    slices = meta_df.apply( lambda row:slice(row['RANGEBEGINNINGDATE']+'T'+row['RANGEBEGINNINGTIME'], \
                                    row['RANGEENDINGDATE']+'T'+row['RANGEENDINGTIME']), axis=1 ).tolist()

    if metric == 'mean':
        out_arr = np.array([ ds_sel.sel(time=sl).mean('time')[variable].values for sl in slices ])
    elif metric == 'min':
        out_arr = np.array([ ds_sel.sel(time=sl).min('time')[variable].values for sl in slices ])
    elif metric == 'max':
        out_arr = np.array([ ds_sel.sel(time=sl).max('time')[variable].values for sl in slices ])

    new_times = pd.DatetimeIndex(meta_df['date'].astype(str).apply(lambda x: datetime.datetime.strptime(str(x), '%Y%j')).tolist())

    # make nc
    new_ds = xr.Dataset({variable.lower(): (['time','yc', 'xc'], out_arr)},
                        coords={'xc': ('xc', ds_sel.xc.values),
                                'yc': ('yc', ds_sel.yc.values),
                                'time':new_times })

    new_ds[variable].encoding = ds_sel[variable].encoding
    return new_ds

def modisify_wrf( files, output_path, sensor, meta_df, metric, variable, ncpus=7 ):
    ''' make the WRF data temporally look like MODIS LST 8-Day composite Daytime'''
    fn = '_'.join(files[0].split('_')[:-1])+'_{}_{}.nc'.format(metric,sensor)
    output_filename = os.path.join( output_path, os.path.basename(fn).replace('hourly','8Day_daytime') )
    print( output_filename )

    # this whole piece is a way for us to interpolate a leapday for comparison...
    # # WRF has NOLEAP calendar and there is no way to properly aggregate using xarray semantics without adding it.
    years = [ get_year_from_fn(fn) for fn in files ]
    args = [ (fn, year, variable) for fn,year in zip(files,years) ]
    pool = mp.Pool( ncpus )
    fixed_files = pool.map( run_interpolate_leap_years, args )
    pool.close()
    pool.join()

    # get the proper begin and end dates... this allows for proper subsetting and ~8-day computation like MODIS
    with xr.open_mfdataset( fixed_files ) as ds:
        times = ds.time.to_index()
        mod_times = [pd.Timestamp.strptime(i, '%Y-%m-%d') for i in times.strftime('%Y-%m-%d')]
        # cleanup
        ds = None
        del ds

    # resample the data following begin / end times from the MODIS file metadata for 8-day composites
    ds_res = make_wrf_like_modis( fixed_files, meta_df, metric, product=sensor, variable=variable )

    # make sure the output path exists
    dirname = os.path.dirname( output_filename )
    if not os.path.exists( dirname ):
        _ = os.makedirs( dirname )

    # write to disk and cleanup
    print('writing 8-Day netcdf to disk')
    ds_res.to_netcdf( output_filename, mode='w', format='NETCDF4' )
    ds_res.close()
    del ds_res

    # clear out the temp files we have made
    _ = [ os.unlink(fn) for fn in fixed_files ]
    return output_filename

def run( x ):
    return modisify_wrf( **x )

if __name__ == '__main__':
    import os, glob, calendar, shutil, itertools
    import numpy as np
    import xarray as xr
    import dask
    import pandas as pd
    import datetime
    import multiprocessing as mp

    # data_path = '/storage01/malindgren/wrf_ccsm4/hourly_fix' #'/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix' #,]
    data_path = '/rcs/project_data/wrf_data/hourly_fix'
    variable = 'tsk'
    base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
    output_path = os.path.join(base_path, 'WRF_DATA','WRF_Day_hours', variable)
    metrics = ['mean','max','min']
    sensors = ['MOD11A2','MYD11A2']
    
    # read in a metadata file built from metadata entries within each MOD11A2/MYD11A2 HDF file 
    #    --> see "MODIS_Processing/get_begin_end_range_8day_composite_MODIS_hdf_metadata.py" to build this file
    meta_df = pd.read_csv(os.path.join(base_path,'MODIS_DATA','ancillary','MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv'))
    meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )
    model_groups = ['ERA-Interim_historical', 'GFDL-CM3_historical', 'GFDL-CM3_rcp85', 'NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']
    # model_groups = ['NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']

    begin_year, end_year = 2000, 2019
    for sensor, metric, group in itertools.product(sensors, metrics, model_groups):
        print(sensor,metric,group)

        # subset the metadata DataFrame
        meta_df_sub = meta_df[meta_df['product'] == sensor].sort_values('date')
        
        # convert the RANGE TIMES to Pandas date objects for easier querying.
        meta_times_end = meta_df_sub['RANGEENDINGDATE'].apply(lambda x: pd.Timestamp.strptime(x+' 23:00:00','%Y-%m-%d %H:%M:%S') )
        meta_times_begin = meta_df_sub['RANGEBEGINNINGDATE'].apply(lambda x: pd.Timestamp.strptime(x+' 23:00:00','%Y-%m-%d %H:%M:%S') )

        # list files we want within the range we want...
        files = sorted( glob.glob( os.path.join( data_path, variable, '*{}*.nc'.format(group) ) ) )
        # filter to logical years to avoid over-processing...
        files = [fn for fn in files if fn.endswith(tuple('{}.nc'.format(i) for i in np.arange(begin_year, end_year+1)))]

        # pull the time objects from the file list
        with xr.open_mfdataset(files) as ds:
            times = ds.time.to_index()
        
        # slice meta_df_sub to the overlapping times
        meta_df_sub = meta_df_sub[(meta_times_begin > times[0])&(meta_times_end < times[-1])].copy(deep=True)
        print(meta_df_sub.iloc[0].date,meta_df_sub.iloc[-1].date)

        # run
        names = ['files', 'output_path', 'sensor','metric','meta_df','variable',]
        args = dict(zip(names,[files, output_path, sensor, metric, meta_df_sub, variable,])) 
        out = run( args )
        
