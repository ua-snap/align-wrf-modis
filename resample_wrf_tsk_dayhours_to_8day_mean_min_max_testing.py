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
    ds = xr.open_mfdataset( files )
    ds_sel = ds.sel( time=get_day_hours(ds['time.hour']) )
    
    # now subset to the times that match
    ds_sel = ds_sel.sel( time=slice(begin, ds_sel.time.to_index()[-1].strftime('%Y-%m-%d')) ).load()
    ds.close()
    ds = None

    # # make time objects for begin and end times
    # b = pd.Timestamp( begin )
    # e = pd.Timestamp( end + ' 23:59:59' )

    # read in and slice to the proper dates in the pre-built metadata file
    group = meta_df[meta_df['product'] == product].copy(deep=True).sort_values('date')
    # group = group[(group.date >= int(b.strftime('%Y%j'))) & (group.date <= int(e.strftime('%Y%j')))]

    # # handle range ending date for the 8-day interval and see if the ds includes that date
    # modrangeend = pd.Timestamp(group.iloc[-1]['RANGEENDINGDATE'])
    # dsend = ds_sel.time.to_index()[-1]
    # if np.abs((modrangeend - dsend)) < pd.Timedelta(6,'D'):
    #     group = group.iloc[:-1]

    # make the slices we need from the metadata file
    slices = group.apply( lambda row:slice(row['RANGEBEGINNINGDATE']+'T'+row['RANGEBEGINNINGTIME'], \
                                    row['RANGEENDINGDATE']+'T'+row['RANGEENDINGTIME']), axis=1 )

    if metric == 'mean':
        out_arr = np.array([ ds_sel.sel(time=sl).mean('time')[variable].values for sl in slices ])
    elif metric == 'min':
        out_arr = np.array([ ds_sel.sel(time=sl).min('time')[variable].values for sl in slices ])
    elif metric == 'max':
        out_arr = np.array([ ds_sel.sel(time=sl).max('time')[variable].values for sl in slices ])

    new_times = pd.DatetimeIndex(group['date'].astype(str).apply(lambda x: datetime.datetime.strptime(str(x), '%Y%j')).tolist())

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

    path = '/storage01/malindgren/wrf_ccsm4/hourly_fix' #'/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix' #,]
    variable = 't2'
    base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
    output_path = os.path.join(base_path, 'MODIS_DATA','WRF_Day_hours', variable)
    metrics = ['mean','max','min']
    sensors = ['MOD11A2','MYD11A2']

    # mess with the begin / end dates for the MODIS data for comparison [this is now only partially used...  Harwire was more sane.]
    # sensor_dates = {'MOD11A2': (2000049,2018233), 'MYD11A2':(2002185,2018233)}
    
    # read in a metadata file built from metadata entries within each MOD11A2/MYD11A2 HDF file
    meta_df = pd.read_csv(os.path.join(base_path,'MODIS_DATA','ancillary','MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv'))
    meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )
    # groups_lists = [['ERA-Interim_historical', 'GFDL-CM3_historical', 'GFDL-CM3_rcp85',],['NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']]
    # model_groups = ['ERA-Interim_historical', 'GFDL-CM3_historical', 'GFDL-CM3_rcp85', 'NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']
    model_groups = ['NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']

    begin_year, end_year = 2000, 2019
    # for path, groups in zip(paths, groups_lists): # WHEN ALL DATA ARE MOVED TO RCS THIS NEEDS UPDATING
    for sensor, metric, group in itertools.product(sensors, metrics, model_groups):
        # for sensor, metric, group in itertools.product(sensors, metrics, groups):
        print(sensor,metric,group)

        # subset the metadata DataFrame
        meta_df_sub = meta_df[meta_df['product'] == sensor].sort_values('date')
        
        # convert the RANGE TIMES to Pandas date objects for easier querying.
        meta_times_end = meta_df_sub['RANGEENDINGDATE'].apply(lambda x: pd.Timestamp.strptime(x+' 23:00:00','%Y-%m-%d %H:%M:%S') )
        meta_times_begin = meta_df_sub['RANGEBEGINNINGDATE'].apply(lambda x: pd.Timestamp.strptime(x+' 23:00:00','%Y-%m-%d %H:%M:%S') )


        # # # # NEW BELOW
        # # find the most appropriate overlapping dates
        # sensor_begin, sensor_end = sensor_dates[ sensor ]
        # sensor_begin = datetime.datetime( int(str(sensor_begin)[:4]),1,1 ) + datetime.timedelta( int(str(sensor_begin)[4:]) )
        # sensor_end = datetime.datetime( int(str(sensor_end)[:4]),1,1 ) + datetime.timedelta( int(str(sensor_end)[4:]) )

        # list files we want within the range we want...
        files = sorted( glob.glob( os.path.join( path, variable, '*{}*.nc'.format(group) ) ) )
        # filter to logical years to avoid over-processing...
        files = [fn for fn in files if fn.endswith(tuple('{}.nc'.format(i) for i in np.arange(begin_year, end_year+1)))]

        # pull the time objects from the file list
        with xr.open_mfdataset(files) as ds:
            times = ds.time.to_index()
        ds = None; del ds # cleanup

        # slice meta_df_sub to the overlapping times
        meta_df_sub = meta_df_sub[(meta_times_begin > times[0])&(meta_times_end < times[-1])].copy(deep=True)

        print(meta_df_sub.iloc[0].date,meta_df_sub.iloc[-1].date)

        # # make a dataframe of the files
        # elems = [ 'variable', 'agg', 'wrf', 'model', 'scenario', 'year' ]
        # df = pd.DataFrame([ dict(zip(elems, os.path.splitext(os.path.basename(fn))[0].split('_'))) for fn in files ])
        # df['fn'] = files

        # # get begin/end times from the dataset(s) being processed
        # begin_fn, = df[df.year == df.year.min()].fn.tolist()
        # with xr.open_dataset( begin_fn ) as tmp:
        #     ds_begin = tmp.time.to_pandas()[0]

        # end_fn, = df[df.year == df.year.max()].fn.tolist()
        # with xr.open_dataset( end_fn ) as tmp:
        #     ds_end = tmp.time.to_pandas()[-1]

        # # get the start date based on the meta_df
        # ds_begin_ordinal = int(ds_begin.strftime('%Y%j'))
        # start_idx = (np.abs(np.array(meta_df['date']) - ds_begin_ordinal)).argmin()
        # new_begin = pd.Timestamp.strptime(meta_df.loc[start_idx, 'date'].astype(str), '%Y%j')


        # ds_end_ordinal = int(ds_end.strftime('%Y%j'))
        # end_idx = (np.abs(np.array(meta_df['date']) - ds_end_ordinal)).argmin()
        # new_end = pd.Timestamp.strptime(meta_df.loc[end_idx, 'date'].astype(str), '%Y%j')
        # ds_end = pd.Timestamp.strptime(str(ds_end_ordinal),'%Y%j')
        # ndays = (new_end - ds_end).days

        # if (np.sign(ndays) == 1) and (ndays )
        # if new_end > ds_end:
        #     pd.Timestamp.strptime(meta_df.loc[end_idx-1, 'date'].astype(str), '%Y%j')




        # # now its time to play with times
        # if sensor_begin < ds_begin:
        #     begin = ds_begin
        # elif sensor_begin > ds_begin:
        #     begin = sensor_begin
        # else:
        #     begin = sensor_begin

        # if sensor_end < ds_end:
        #     end = sensor_end
        # elif sensor_end > ds_end:
        #     end = ds_end
        # else:
        #     end = sensor_end

        # # grab 8-days beyond -- for aggregation? THIS NEEDS SOME TESTING!
        # end = end + pd.Timedelta(8, 'D')

        # # now make the years an integer and return the files that match the begin/end period
        # df.year = df.year.astype(int)
        # files = df[(df.year >= begin.year ) & (df.year <= end.year)].fn.tolist()

        # # # # NEW ABOVE
        # run
        # names = ['files', 'output_path', 'sensor','metric','meta_df','begin','end','variable',]
        # args = dict(zip(names,[files, output_path, sensor, metric,meta_df,begin.strftime('%Y-%m-%d'),end.strftime('%Y-%m-%d'),variable,]))
        names = ['files', 'output_path', 'sensor','metric','meta_df','variable',]
        args = dict(zip(names,[files, output_path, sensor, metric, meta_df_sub, variable,])) 
        out = run( args )
        break
