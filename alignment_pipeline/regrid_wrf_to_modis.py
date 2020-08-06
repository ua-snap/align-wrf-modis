"""
Re-grid WRF to the clipped MODIS, which comes at a lower resolution (downsampling)
"""

import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
import os, glob, subprocess, itertools, datetime
from helpers import check_env


def warp_wrf(temp_arr, temp_meta, fp, out_fp):
    with rio.open(out_fp, "w", **temp_meta) as rst:
        rst.write(temp_arr)
    _ = subprocess.call(
        ["gdalwarp", "-s_srs", "epsg:3338", "-q", "-r", "near", fp, out_fp]
    )
    with rio.open(out_fp, mode="r+") as out:
        arr = out.read(1)
        # values close to zero form from reasmpling
        arr[np.isclose(arr, 0)] = -9999
        out.write(arr, 1)


def open_raster(fn, band=1):
    with rio.open(fn) as rst:
        arr = rst.read(band).copy()
    rst = None
    return arr


def get_dates(fps):
    jd_lst = [os.path.basename(fp).split("_")[-2] for fp in fps]
    dt = [datetime.datetime.strptime(jd_str, "%Y%j") for jd_str in jd_lst]
    dt_arr = np.array(dt, dtype="datetime64")
    return dt_arr


if __name__ == "__main__":
    # check environment
    wrf_env_var = check_env()
    if not wrf_env_var:
        exit("Environment variables incorrectly setup, check README for requirements")
    wrf_var = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    out_base_dir = os.getenv("OUTPUT_DIR")
    template_fp = glob.glob(os.path.join(out_base_dir, "MODIS", "lst_1km_3338", "*"))[0]
    out_dir = os.path.join(out_base_dir, "WRF", "{}_1km_3338".format(wrf_var))
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)

    wrf_dir = os.path.join(scratch_dir, "WRF", "{}_1km_3338".format(wrf_var))

    with rio.open(template_fp) as tmp:
        temp_arr = np.empty_like(tmp.read())
        temp_meta = tmp.meta
    # setup multiband dirs
    out_multi_dir = out_dir + "_multiband"
    if not os.path.exists(out_multi_dir):
        _ = os.makedirs(out_multi_dir)

    for directory in os.listdir(wrf_dir):
        wrf_fps = sorted(glob.glob(os.path.join(wrf_dir, directory, "*.tif")))
        set_dir = os.path.join(out_dir, directory)
        if not os.path.exists(set_dir):
            _ = os.makedirs(set_dir)

        for fp in wrf_fps:
            fn = os.path.basename(fp)
            print("working on", fn)
            out_fp = os.path.join(set_dir, fn)
            warp_wrf(temp_arr, temp_meta, fp, out_fp)

        # assemble files into geotiffs and netcdf
        print("Making geotiffs and netcdfs")
        bands_fps = sorted(glob.glob(os.path.join(out_dir, directory, "*.tif")))
        bands_arr = np.array([open_raster(fn) for fn in bands_fps])
        with rio.open(bands_fps[0]) as tmp:
            new_meta = tmp.meta.copy()
            new_meta.update(count=bands_arr.shape[0])
            # for netcdf
            idx = np.arange(tmp.width)
            idy = np.arange(tmp.height)
            xc = tmp.xy(np.repeat(0, idx.shape), idx)[0]
            yc = tmp.xy(idy, np.repeat(0, idy.shape))[1]

        # geotiff
        out_tif_fp = os.path.join(
            out_multi_dir, "{}_{}_multiband.tif".format(wrf_var, directory),
        )
        with rio.open(out_tif_fp, "w", **new_meta,) as out:
            out.write(bands_arr)

        tmp = None
        out = None
        del tmp, out

        # netcdf
        dates = get_dates(bands_fps)
        ds = xr.Dataset(
            {wrf_var: (["date", "yc", "xc"], bands_arr)},
            coords={"xc": xc, "yc": yc, "date": dates,},
        )
        attrs = {
            "title": "1km MODIS-ailgned WRF TSK",
            "creation_time": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "epsg": 3338,
            "nc_source": "WRFDS_2000-02-21_serdp.nc",
            "SNAP_version": "0.1.0",
        }
        ds.attrs = attrs
        out_nc_fp = out_tif_fp.replace(".tif", ".nc")
        ds.to_netcdf(out_nc_fp)
        print(out_tif_fp)
        print(out_nc_fp)

