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

def make_wrf_like_modis( files, meta_df, product ):
    ds = xr.open_mfdataset( files )
    ds_sel = ds.sel( time=get_day_hours(ds['time.hour']) ).load()
    ds.close()
    ds = None

    # read in an slice to the proper dates the pre-built metadata file...
    group = meta_df[meta_df['product'] == product ].copy(deep=True).sort_values('date')

    # what are the best begin / end dates given the file dates
    modbegin = pd.Timestamp(datetime.datetime.strptime(str(group.iloc[0]['date']), '%Y%j'))
    modend = pd.Timestamp(datetime.datetime.strptime(str(group.iloc[-1]['date']), '%Y%j'))
    dsbegin = pd.Timestamp( ds_sel.isel(time=0).time.values )
    dsend = pd.Timestamp( ds_sel.isel(time=-1).time.values )

    if modbegin <= dsbegin:
        begin = modbegin
    else:
        begin = dsbegin

    if modend <= dsend:
        end = modend
    else:
        end = dsend

    begin = int(begin.strftime('%Y%j'))
    end = int(end.strftime('%Y%j'))

    # slice it up and make the '8-day' means...
    group = group[(group.date >= begin) & (group.date <= end)]
    slices = group.apply( lambda row:slice(row['RANGEBEGINNINGDATE']+'T'+row['RANGEBEGINNINGTIME'], \
                                    row['RANGEENDINGDATE']+'T'+row['RANGEENDINGTIME']), axis=1 )

    mean_arr = np.array([ ds_sel.sel(time=sl).mean('time').tsk.values for sl in slices ])
    new_times = pd.DatetimeIndex(group['date'].astype(str).apply(lambda x: datetime.datetime.strptime(str(x), '%Y%j')).tolist())

    # make nc
    new_ds = xr.Dataset({variable.lower(): (['time','yc', 'xc'], mean_arr)},
                        coords={'xc': ('xc', ds_sel.xc.values),
                                'yc': ('yc', ds_sel.yc.values),
                                'time':new_times })

    new_ds.tsk.encoding = ds_sel.tsk.encoding
    return new_ds

def run_make_wrf_like_modis(x):
    return make_wrf_like_modis(**x)

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

def modisify_wrf( files, output_path, sensor, meta_df, ncpus=5 ):
    ''' make the WRF data temporally look like MODIS LST 8-Day composite Daytime'''
    fn = '_'.join(files[0].split('_')[:-1])+'_{}.nc'.format(sensor)
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
    ds_res = make_wrf_like_modis( fixed_files, meta_df, product=sensor )

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
    import os, glob, calendar, shutil
    import numpy as np
    import xarray as xr
    import dask
    import pandas as pd
    import datetime
    import multiprocessing as mp

    path = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix'
    output_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
    variable = 'tsk'

    # mess with the begin / end dates for the MODIS data for comparison
    sensor_dates = {'MOD11A2': (2000049,2018233), 'MYD11A2':(2002185,2018233)}
    sensors = ['MOD11A2','MYD11A2']
    
    # read in a metadata file built from metadata entries within each MOD11A2/MYD11A2 HDF file
    meta_df = pd.read_csv('/workspace/Shared/Users/malindgren/MODIS_DATA/ancillary/MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv')
    meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )

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
        names = ['files', 'output_path', 'sensor']
        args = [ dict(zip(names,[group, output_path, sensor])) for group in grouped_files ]
        _ = [a.update(meta_df=meta_df) for a in args]
        break
        out = [run(a) for a in args]
