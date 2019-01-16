# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# resample the WRF skin-temp hourlies to daily day-times and then to 8-Day averages for comparison with MODIS
# Author: Michael Lindgren (malindgren@alaska.edu)
# License: MIT
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

def get_day_hours( month ):
    return (month >= 9) & (month <= 17)

def get_year_from_fn( fn ):
    return int(os.path.basename(fn).split('.nc')[0].split('_')[-1])

def fill_leap( fn, out_fn ):
    with xr.open_dataset( fn ) as ds:
        ds_filled = ds.resample(time='1H').interpolate('slinear')
        ds_filled.tsk.encoding = ds.tsk.encoding
        ds_filled.to_netcdf(out_fn)
        ds_filled.close(); ds_filled = None
    return out_fn

def interpolate_leap_years( fn, year, tmp_path='/atlas_scratch/malindgren/TMP' ):
    # make an output path to dump the temporary intermediate file
    out_fn = os.path.join( tmp_path, os.path.basename(fn) )
    if os.path.exists( out_fn ):
        _ = os.unlink( out_fn )
    if calendar.isleap( year ):
        # fill it with a spline
        _ = fill_leap( fn, out_fn )
    else:
        # copy the file to the temp location for speed the symlinks were very slow.
        _ = shutil.copy( fn, out_fn )
    return out_fn

def run_interpolate_leap_years(x):
    return interpolate_leap_years(*x)

def make_wrf_like_modis( files, meta_df, metric, product, begin, end ):
    ds = xr.open_mfdataset( files )
    ds_sel = ds.sel( time=get_day_hours(ds['time.hour']) )
    
    # now subset to the times that match
    ds_sel = ds_sel.sel( time=slice(begin, ds_sel.time.to_index()[-1].strftime('%Y-%m-%d')) ).load()
    ds.close()
    ds = None

    # make time objects for begin and end times
    b = pd.Timestamp( begin )
    e = pd.Timestamp( end + ' 23:59:59' )

    # read in and slice to the proper dates in the pre-built metadata file
    group = meta_df[meta_df['product'] == product].copy(deep=True).sort_values('date')
    group = group[(group.date >= int(b.strftime('%Y%j'))) & (group.date <= int(e.strftime('%Y%j')))]

    # handle range ending date for the 8-day interval and see if the ds includes that date
    modrangeend = pd.Timestamp(group.iloc[-1]['RANGEENDINGDATE'])
    dsend = ds_sel.time.to_index()[-1]
    if np.abs((modrangeend - dsend)) < pd.Timedelta(6,'D'):
        group = group.iloc[:-1]

    # make the slices we need from the metadata file
    slices = group.apply( lambda row:slice(row['RANGEBEGINNINGDATE']+'T'+row['RANGEBEGINNINGTIME'], \
                                    row['RANGEENDINGDATE']+'T'+row['RANGEENDINGTIME']), axis=1 )

    if metric == 'mean':
        out_arr = np.array([ ds_sel.sel(time=sl).mean('time').tsk.values for sl in slices ])
    elif metric == 'min':
        out_arr = np.array([ ds_sel.sel(time=sl).min('time').tsk.values for sl in slices ])
    elif metric == 'max':
        out_arr = np.array([ ds_sel.sel(time=sl).max('time').tsk.values for sl in slices ])

    new_times = pd.DatetimeIndex(group['date'].astype(str).apply(lambda x: datetime.datetime.strptime(str(x), '%Y%j')).tolist())

    # make nc
    new_ds = xr.Dataset({variable.lower(): (['time','yc', 'xc'], out_arr)},
                        coords={'xc': ('xc', ds_sel.xc.values),
                                'yc': ('yc', ds_sel.yc.values),
                                'time':new_times })

    new_ds.tsk.encoding = ds_sel.tsk.encoding
    return new_ds

def modisify_wrf( files, output_path, sensor, meta_df, metric, begin, end, ncpus=7 ):
    ''' make the WRF data temporally look like MODIS LST 8-Day composite Daytime'''
    fn = '_'.join(files[0].split('_')[:-1])+'_{}_{}.nc'.format(metric,sensor)
    output_filename = os.path.join( output_path, os.path.basename(fn).replace('hourly','8Day_daytime') )
    print( output_filename )

    # this whole piece is a way for us to interpolate a leapday for comparison...
    # # WRF has NOLEAP calendar and there is no way to properly aggregate using xarray semantics without adding it.
    years = [ get_year_from_fn(fn) for fn in files ]
    args = [ (fn, year) for fn,year in zip(files,years) ]
    pool = mp.Pool( ncpus )
    fixed_files = pool.map( run_interpolate_leap_years, args )
    pool.close()
    pool.join()

    # resample the data following begin / end times from the MODIS file metadata for 8-day composites
    ds_res = make_wrf_like_modis( fixed_files, meta_df, metric, product=sensor, begin=begin, end=end )

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

    paths = ['/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix','/storage01/malindgren/wrf_ccsm4/hourly_fix']
    output_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
    variable = 'tsk'
    metrics = ['mean','max','min']

    # mess with the begin / end dates for the MODIS data for comparison [this is now only partially used...  Harwire was more sane.]
    sensor_dates = {'MOD11A2': (2000049,2018233), 'MYD11A2':(2002185,2018233)}
    sensors = ['MOD11A2','MYD11A2']
    
    # read in a metadata file built from metadata entries within each MOD11A2/MYD11A2 HDF file
    meta_df = pd.read_csv('/workspace/Shared/Users/malindgren/MODIS_DATA/ancillary/MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv')
    meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )
    groups_lists = [['ERA-Interim_historical', 'GFDL-CM3_historical', 'GFDL-CM3_rcp85',],['NCAR-CCSM4_historical', 'NCAR-CCSM4_rcp85']]

    for path, groups in zip(paths, groups_lists):
        for sensor, metric, group in itertools.product(sensors, metrics, groups):
            print(sensor,metric,group)
            
            # # more date-fun [from older parsing method. this could be rebuilt but needs work.]
            # begin,end = sensor_dates[ sensor ]
            # begin_dt = datetime.datetime( int(str(begin)[:4]),1,1 ) + datetime.timedelta( int(str(begin)[4:]) )
            # end_dt = datetime.datetime( int(str(end)[:4]),1,1 ) + datetime.timedelta( int(str(end)[4:]) )
            # begin = begin_dt.strftime( '%Y-%m-%d' )
            # end = end_dt.strftime( '%Y-%m-%d' )

            # print(begin, end)

            # if (begin_dt.year < 2006) and ('rcp85' in group):
            #     begin_year = 2006
            # elif ('historical' in group) and (('NCAR-CCSM4' in group) or ('GFDL-CM3' in group)):
            #     begin_year = begin_dt.year
            # elif ('ERA-Interim' in group):
            #     begin_year = begin_dt.year
            # else:
            #     BaseException('wrong group')

            # hardwired date groups since it was confusing and there are only so man groups to assess.
            date_groups = {
                            'NCAR-CCSM4_rcp85_MOD11A2':('2006-01-02', '2018-09-01'),
                            'NCAR-CCSM4_rcp85_MYD11A2':('2006-01-02', '2018-09-01'),
                            'NCAR-CCSM4_historical_MOD11A2':('2000-02-18', '2005-12-30'),
                            'NCAR-CCSM4_historical_MYD11A2':('2002-07-04', '2005-12-30'),
                            'GFDL-CM3_rcp85_MOD11A2':('2006-01-02', '2018-09-01'),
                            'GFDL-CM3_rcp85_MYD11A2':('2006-01-02', '2018-09-01'),
                            'GFDL-CM3_historical_MOD11A2':('2000-02-18', '2005-12-30'),
                            'GFDL-CM3_historical_MYD11A2':('2002-07-04', '2005-12-30'),
                            'ERA-Interim_historical_MOD11A2':('2000-02-18', '2015-10-29'),
                            'ERA-Interim_historical_MYD11A2':('2002-07-04','2015-10-29')
                        }

            begin, end = date_groups['_'.join([group,sensor])]
            begin_dt, end_dt = pd.Timestamp(begin), pd.Timestamp(end)

            # list files we want within the range we want...
            files = [ fn for fn in glob.glob( os.path.join( path, variable, '*{}*.nc'.format(group) ) ) 
                                    if int(fn.split('_')[-1].split('.')[0]) in list(range(begin_dt.year, end_dt.year+1))]

            # make a dataframe of the files
            elems = [ 'variable', 'agg', 'wrf', 'model', 'scenario', 'year' ]
            df = pd.DataFrame([ dict(zip(elems, os.path.splitext(os.path.basename(fn))[0].split('_'))) for fn in files ])
            df['fn'] = files

            # group and sort the files 
            files = df.sort_values('year').fn.tolist()
            names = ['files', 'output_path', 'sensor','metric','meta_df','begin','end',]
            args = dict(zip(names,[files, output_path, sensor, metric,meta_df,begin,end]))
            out = run( args )
            