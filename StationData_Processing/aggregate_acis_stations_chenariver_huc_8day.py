# agg ACIS like MODIS
import os, glob
import pandas as pd
import numpy as np

base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'
acis_path = os.path.join(base_path,'STATION_DATA','acis', 'raw')
output_path = os.path.join(base_path,'STATION_DATA','acis', '8day')

if not os.path.exists( output_path ):
    _ = os.makedirs( output_path )

# open the aggregation ranges harvested from the HDF raw MODIS Tiles
meta_df = pd.read_csv(os.path.join(base_path,'MODIS_DATA','ancillary','MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv'))

# list the acis files that have been output:
files = glob.glob(os.path.join(acis_path, '*.csv'))
sensor = 'MOD11A2'
meta_df['year'] = meta_df['date'].apply( lambda x: int(str(x)[:4]) )
meta_df = meta_df[meta_df['product'] == sensor]
avg_groups = [ list(row[['RANGEBEGINNINGDATE','RANGEENDINGDATE']]) for idx,row in meta_df.iterrows() ]

for fn in files:
    df = pd.read_csv( fn, index_col=0 )
    df = df.replace('M', np.nan)
    df = df.replace('T', '.001')
    df = df.astype(float)

    df = (df - 32.0) * 5.0 / 9.0 # make celcius
    dat = {}
    for mod_begin, mod_end in avg_groups:
        test = df.loc[slice(mod_begin, mod_end)]
        if not test.empty:
            dat[mod_begin] = test.mean()
    
    out_fn = os.path.join(output_path, os.path.basename(fn).replace('.csv', '_8day.csv'))
    out_df = pd.DataFrame(dat).T
    out_df.to_csv(out_fn)
