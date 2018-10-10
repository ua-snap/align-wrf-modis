# figure out what hours it typically grabs the data from MOD11A2
import os, rasterio, itertools, glob
import numpy as np

os.chdir('/workspace/Shared/Users/malindgren/MODIS_DATA/TEMP')

all_hours = [ list(np.unique(rasterio.open(fn).read(1))) for fn in glob.glob('./*.tif') if '_03.tif' in fn ]
common_hours = np.array(list(itertools.chain.from_iterable( all_hours )))
scale_factor = 0.1 # from the MOD11A2 metadata
common_hours_scaled = common_hours[ common_hours != 255 ] * scale_factor
unique_counts = np.unique(common_hours, return_counts=True)
