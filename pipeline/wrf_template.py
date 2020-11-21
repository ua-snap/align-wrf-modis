"""Create a WRF template GeoTIFF in epsg:3338 for clipping MODIS"""

import glob, os, subprocess
import numpy as np
import rasterio as rio
from helpers import check_env


def reproject_band(arr, meta, band_fp):
    """run the reprojection for a template WRF band"""
    proj_fp = band_fp.replace("4326", "3338")
    with rio.open(band_fp, "w", **meta) as src:
        src.write(arr, 1)
    _ = subprocess.call(
        ["gdalwarp", "-t_srs", "epsg:3338", "-q", "-overwrite", band_fp, proj_fp]
    )
    os.unlink(band_fp)
    return proj_fp


if __name__ == "__main__":
    _ = check_env()
    # setup dirs
    wrf_var = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    wrf_dir = os.path.join(scratch_dir, "WRF", "resampled", wrf_var)
    # template WRF filepath
    wrf_fp = glob.glob(os.path.join(wrf_dir, f"*gfdl_2007-2017.nc"))[0]
    # with xr.open_dataset(wrf_fp) as ds:
    with rio.open(f"netcdf:{wrf_fp}:Band1") as src:
        band_arr = src.read(1)
        meta = src.meta
    meta.update(
        {"driver": "GTiff", "nodata": -9999,}
    )
    band_fp = os.path.join(scratch_dir, "ancillary", "wrf_4326_template.tif")
    proj_fp = reproject_band(band_arr, meta, band_fp)
    print(f"Template WRF GeoTIFF saved to {proj_fp}")
