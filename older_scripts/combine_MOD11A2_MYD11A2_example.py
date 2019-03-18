# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# [UNTESTED -- YMMV] combine MOD11A2 and MYD11A2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
import rasterio
import numpy as np

mod_fn = '/workspace/Shared/Users/malindgren/MODIS_DATA/TEMP/MOD11A2_A2018233_h11v02_006_2018242041022_01.tif'
myd_fn = '/workspace/Shared/Users/malindgren/MODIS_DATA/TEMP/MYD11A2_A2018233_h11v02_006_2018242042802_01.tif'
output_filename = '/workspace/Shared/Users/malindgren/MODIS_DATA/TEST/combined_A2018233_h11v02_006_01.tif'

mod = rasterio.open(mod_fn)
myd = rasterio.open(myd_fn)

mod_arr = mod.read(1)
myd_arr = myd.read(1)

mean_ind = np.where((mod_arr != 0) & (myd_arr != 0))
mod_fill_ind = np.where((mod_arr != 0) & (myd_arr == 0))
myd_fill_ind = np.where((mod_arr == 0) & (myd_arr != 0))

# fill or mean...
out_arr = np.zeros_like( mod_arr )

out_arr[mean_ind] = (mod_arr[mean_ind] + myd_arr[mean_ind])/2
out_arr[mod_fill_ind] = mod_arr[mod_fill_ind]
out_arr[myd_fill_ind] = myd_arr[myd_fill_ind]

out_arr = out_arr * 0.02

meta = mod.meta.copy()
meta.update( compress='lzw', dtype='float32' )
with rasterio.open( output_filename, 'w', **meta ) as out:
	out.write( out_arr.astype(np.float32), 1 )
