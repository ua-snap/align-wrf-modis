# SEE USER GUIDE -- page 10 
# https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf

import rasterio
import numpy as np

with rasterio.open('/workspace/Shared/Users/malindgren/MOD11A1/prepped/MOD11A1_A2018247_h11v02_006_2018248084624_Emis_31.tif') as rst:
	meta = rst.meta.copy()
	meta.update( compress='lzw', dtype='float32' )
	arr = rst.read(1)
	arr_out = np.copy(arr).astype(np.float32)
	ind = np.where( arr != meta['nodata'] )
	arr_out[ind] = (arr[ind]*0.0020) + 0.49

	with rasterio.open('/workspace/Shared/Users/malindgren/MOD11A1/rescaled/MOD11A1_A2018247_h11v02_006_2018248084624_Emis_31.tif', 'w', **meta ) as out:
		out.write( arr.astype(np.float32), 1 )



