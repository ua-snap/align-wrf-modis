# how to resample the hourlies to daily day-times
def get_day_hours(month):
    return (month >= 9) & (month <= 17)

if __name__ == '__main__':
    import os
    import numpy as np
    import xarray as xr 
    import datetime

    path = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix'
    output_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
    variable = 'tsk'

    # mess with the begin / end dates for the MODIS data for comparison
    begin,end = (2000049,2018233)
    begin_dt = datetime.datetime( int(str(begin)[:4]),1,1 ) + datetime.timedelta( int(str(begin)[4:]) )
    end_dt = datetime.datetime( int(str(end)[:4]),1,1 ) + datetime.timedelta( int(str(end)[4:]) )
    begin = begin_dt.strftime('%Y-%m-%d')
    end = end_dt.strftime('%Y-%m-%d')

    # list files we want within the range we want...
    files = [ fn for fn in glob.glob( os.path.join( path, variable, '*.nc' ) ) 
                if int(fn.split('_')[-1].split('.')[0]) in list(range(begin_dt.year, end_dt.year+1))]

    # make a dataframe of the files
    elems = ['variable', 'agg', 'wrf', 'model', 'scenario', 'year']
    df = pd.DataFrame([ dict(zip(elems, os.path.splitext(os.path.basename(fn))[0].split('_'))) for fn in files ])
    df['fn'] = files

    grouped_files = df.sort_values('year').groupby(['model', 'scenario']).apply(lambda x: x.fn.tolist() ).tolist()

    for group_files in grouped_files:
        ds = xr.open_mfdataset( group_files )
        fn = '_'.join(group_files[0].split('_')[:-1])+'.nc' # fn prep
        output_filename = os.path.join( output_path, os.path.basename(fn).replace('hourly','8Day_daytime') )

        with xr.open_mfdataset( group_files ) as ds:

            # get the hours we are interested in
            selected_hours = ds.sel( time=get_day_hours(ds['time.hour']) )
            
            # resample to daily the selected hours...
            new_ds = selected_hours.resample( time='D' ).mean()

            # make an 8-day average for comparison with MODIS
            new_mean = new_ds.resample( time='8D' ).mean()
        
        dirname = os.path.dirname( output_filename )
        if not os.path.exists( dirname ):
            _ = os.makedirs( dirname )

        new_mean.to_netcdf( output_filename )
        del new_mean, new_ds

