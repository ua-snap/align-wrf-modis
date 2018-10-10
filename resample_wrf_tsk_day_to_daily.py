# how to resample the hourlies to daily day-times
def get_day_hours(month):
    return (month >= 9) & (month <= 17)

if __name__ == '__main__':
    import os
    import numpy as np
    import xarray as xr 

    path = '/workspace/Shared/Tech_Projects/wrf_data/project_data/wrf_data/hourly_fix'
    variable = 'tsk'

    fn = os.path.join( path, variable, 'tsk_hourly_wrf_GFDL-CM3_rcp85_2016.nc' )
    ds = xr.open_dataset( fn )

    # get the hours we are interested in
    selected_hours = ds.sel( time=get_day_hours(ds['time.hour']) )
    
    # resample to daily the selected hours...
    new_ds = selected_hours.resample( time='D' ).mean()

    # dump to disk
    new_ds.to_netcdf( '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_day_dayhours_wrf_GFDL-CM3_rcp85_2016.nc' )

    # make an 8-day average for comparison with MODIS
    # new_ds = xr.open_dataset('/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_day_dayhours_wrf_GFDL-CM3_rcp85_2016.nc')
    new_mean = new_ds.resample( time='8D' ).mean()
    new_mean.to_netcdf( '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk_8day_dayhours_wrf_GFDL-CM3_rcp85_2016.nc' )