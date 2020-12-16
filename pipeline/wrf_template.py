"""Create a WRF template GeoTIFF in epsg:3338 for clipping MODIS"""

import glob, os, subprocess
import numpy as np
import rasterio as rio
from helpers import check_env


def crop_nodata(fp):
    """Crop GeoTIFF to extent of valid data"""
    with rio.open(fp) as src:
        arr = src.read(1)
        meta = src.meta
        psize = src.xy(0,1)[0] - src.xy(0,0)[0]
    valid = arr != 0
    valid_idx = np.where(valid)
    rmin, rmax = min(valid_idx[0]), max(valid_idx[0]) + 1
    cmin, cmax = min(valid_idx[1]), max(valid_idx[1]) + 1
    # new data cropped to extent of valid data
    arr = arr[rmin:rmax,cmin:cmax]
    # create new metadata for cropped data and overwrite
    west, north = rio.transform.xy(meta["transform"], rmin, cmin)
    meta["transform"] = rio.transform.from_origin(west, north, psize, psize)
    meta["height"], meta["width"] = arr.shape
    os.unlink(fp)
    with rio.open(fp, "w", **meta) as src:
        src.write(arr, 1)
    return None


def reproject_band(arr, meta, band_fp):
    """run the reprojection for a template WRF band"""
    with rio.open(band_fp, "w", **meta) as src:
        src.write(arr, 1)
    proj_fp = band_fp.replace("4326", "3338")
    _ = subprocess.call(
        ["gdalwarp", "-t_srs", "epsg:3338", "-q", "-overwrite", band_fp, proj_fp]
    )
    os.unlink(band_fp)
    crop_nodata(proj_fp)
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
        {"driver": "GTiff", "nodata": 0,}
    )
    band_fp = os.path.join(scratch_dir, "ancillary", "wrf_4326_template.tif")
    proj_fp = reproject_band(band_arr, meta, band_fp)
    print(f"Template WRF GeoTIFF saved to {proj_fp}")
