# agg ACIS like MODIS

import pandas as pd
import os, glob

path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data/chena_river_huc_stations_extracted'

# list the acis files that have been output:
files = glob.glob(os.path.join(path, '*acis*.csv'))
sensor = 'MOD11A2'
meta_df = pd.read_csv('/workspace/Shared/Users/malindgren/MODIS_DATA/ancillary/MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv')
meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )
meta_df = meta_df[meta_df['product'] == sensor]
avg_groups = [ list(row[['RANGEBEGINNINGDATE','RANGEENDINGDATE']]) for idx,row in meta_df.iterrows() ]

for fn in files:
    df = pd.read_csv( fn, index_col=0 )
    df = (df - 32.0) * 5.0 / 9.0 # make celcius
    dat = {}
    for mod_begin, mod_end in avg_groups:
        test = df.loc[slice(mod_begin, mod_end)]
        if not test.empty:
            dat[mod_begin] = test.mean()
    out_fn = fn.replace('.csv', '_8day.csv')
    out_df = pd.DataFrame(dat).T
    out_df.to_csv(out_fn)
