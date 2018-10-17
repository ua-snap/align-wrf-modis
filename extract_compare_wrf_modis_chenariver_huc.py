# extract masked regions for comparison between MODIS-LST and WRF-TSK
if __name__ == '__main__':
    import matplotlib
    matplotlib.use('agg')
    from matplotlib import pyplot as plt
    import os, glob
    import rasterio
    import xarray as xr
    import numpy as np
    import pandas as pd

    modis_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/regrid_20km_wrf_netcdf'
    wrf_path = '/workspace/Shared/Users/malindgren/MODIS_DATA/WRF_Day_hours/tsk'
    mask_fn = '/workspace/UA/malindgren/repos/modis_lst/ancillary/chena_river_huc_19040506.tif'

    modis_files = glob.glob( os.path.join( modis_path, '*.nc' ) )
    wrf_files = [ fn for fn in glob.glob( os.path.join( wrf_path, '*.nc' ) ) if 'MOD11A2' in fn or ('ERA-Interim' in fn and 'MOD11A2' in fn)]
    all_files = modis_files + wrf_files
    
    with rasterio.open( mask_fn ) as mask:
        mask_arr = mask.read(1)
        rows, cols = np.where( mask_arr == 1 )
    
    d = {}
    for fn in all_files:
        if 'tsk_' in fn:
            variable = 'tsk'
            elems = os.path.basename(fn).split('_')
            if elems[5] == 'historical':
                elem = 'hist'
            else:
                elem = elems[5]

            if 'GFDL' in fn:
                model = 'GFDL'
            elif 'ERA-Interim' in fn:
                model = 'ERA'
            else:
                AttributeError('wrong model')

            colname = '-'.join([ model, elem ])
        else:
            variable = 'lst'
            colname = os.path.basename(fn).split('_')[1]

        with xr.open_dataset( fn ) as ds:
            data = ds[ variable ][:,rows,cols].reduce( np.nanmean, axis=(1,2))
            d[ colname ] = data.to_series().to_frame( name=colname )
    
    # # make a dataframe of the values -- this is ugly
    keys = [i for i in d.keys() if 'lst_MOD11A2' not in i]
    df = d['MOD11A2'].copy()
    for key in sorted(keys):
        if key != 'MOD11A2':
            df = df.join(d[key], on='time')

    df = df - 273.15
    df.to_csv( '/workspace/Shared/Users/malindgren/MODIS_DATA/comparison_plots/chena_river_huc_19040506_spatialmean_extraction.csv' )

    ax = df[sorted(keys)].plot( kind='line', linewidth=0.5, title='Chena River HUC Spatial Average\n~8-Day Time-step' )
    plt.ylabel('$^\circ$C')
    plt.xlabel('time')
    # Shrink current axis's height by 10% on the bottom
    box = ax.get_position()

    # Shrink current axis a bit
    ax.set_position([box.x0, box.y0, box.width * 0.83, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.savefig( '/workspace/Shared/Users/malindgren/MODIS_DATA/comparison_plots/chena_river_huc_19040506_wrf_modis_compare.png', dpi=300, figsize=(25,9) )
    plt.close()
