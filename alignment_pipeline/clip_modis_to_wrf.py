"""
Clip the mosaicked, reprojected and rescaled MODIS files 
to the 1km WRF grid, producing the final MODIS data
"""

import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
import os, glob, subprocess, itertools, datetime
from helpers import check_env


def make_cutline(template_fp, temp_dir):
    shp_fp = os.path.join(temp_dir, "clip_modis.shp")
    _ = subprocess.call(
        ["gdaltindex", shp_fp, template_fp]
    )
    return shp_fp


def clip_modis(shp_fp, fp, out_fp):
    _ = subprocess.call(
        ["gdalwarp", "-cutline", shp_fp, "-crop_to_cutline", "-overwrite", fp, out_fp]
    )
    with rio.open(out_fp, mode="r+") as out:
        arr = out.read(1)
        arr[arr == 0] = -9999
        out.write(arr, 1)


def open_raster(fn, band=1):
    with rio.open(fn) as rst:
        arr = rst.read(band).copy()
    rst = None
    return arr


def get_dates(fps):
    jd_lst = [os.path.basename(fp).split("_")[1][1:] for fp in fps]
    dt = [datetime.datetime.strptime(jd_str, "%Y%j") for jd_str in jd_lst]
    dt_arr = np.array(dt, dtype="datetime64")
    return dt_arr


if __name__ == "__main__":
    # check environment
    wrf_env_var = check_env()
    if not wrf_env_var:
        exit("Environment variables incorrectly setup, check README for requirements")
    elif wrf_env_var == "WRF_20KM_DIR":
        exit("this script does not apply to 20km WRF")

    wrf_var = "tsk"
    mod_var = "lst"
    scratch_dir = os.getenv("SCRATCH_DIR")
    temp_dir = os.path.join(scratch_dir, "temp")
    out_base_dir = os.getenv("OUTPUT_DIR")
    modis_dir = os.path.join(scratch_dir, "MODIS", "rescaled")
    out_dir = os.path.join(out_base_dir, "MODIS", "{}_1km_3338".format(mod_var))
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)
    wrf_dir = os.path.join(
        scratch_dir, "WRF", "{}_1km_3338".format(wrf_var), "era_mean_MOD11A2"
    )
    template_fp = os.path.join(wrf_dir, os.listdir(wrf_dir)[0])
    modis_fps = sorted(
        [
            os.path.join(modis_dir, fn)
            for fn in os.listdir(modis_dir)
            if fn.endswith("01.tif")
        ]
    )

    # make cutline fp
    shp_fp = make_cutline(template_fp, temp_dir)

    for fp in modis_fps:
        fn = os.path.basename(fp)
        print("working on", fn[:16])
        out_fp = os.path.join(out_dir, fn)
        clip_modis(shp_fp, fp, out_fp,)

    _ = os.unlink(shp_fp)

    # assemble files into geotiffs and netcdf
    print("Making geotiffs and netcdfs")
    sensors = ["MOD11A2", "MYD11A2"]
    for sensor in sensors:
        print(os.path.join(out_dir, "*{}*.tif".format(sensor)))
        bands_fps = sorted(glob.glob(os.path.join(out_dir, "*{}*.tif".format(sensor))))
        bands_arr = np.array([open_raster(fn) for fn in bands_fps])
        with rio.open(bands_fps[0]) as tmp:
            new_meta = tmp.meta.copy()
            new_meta.update(count=bands_arr.shape[0])
            # for netcdf
            idx = np.arange(tmp.width)
            idy = np.arange(tmp.height)
            xc = tmp.xy(np.repeat(0, idx.shape), idx)[0]
            yc = tmp.xy(idy, np.repeat(0, idy.shape))[1]

        out_dir_multi = out_dir + "_multiband"
        if not os.path.exists(out_dir_multi):
            _ = os.makedirs(out_dir_multi)

        # geotiff
        out_tif_fp = os.path.join(
            out_dir_multi,
            "{}_{}_InteriorAK_006_3338_multiband.tif".format(sensor, mod_var),
        )
        with rio.open(out_tif_fp, "w", **new_meta,) as out:
            out.write(bands_arr)

        tmp = None
        out = None
        del tmp, out

        # netcdf
        dates = get_dates(bands_fps)
        ds = xr.Dataset(
            {mod_var: (["date", "yc", "xc"], bands_arr)},
            coords={"xc": xc, "yc": yc, "date": dates,},
        )
        out_nc_fp = out_tif_fp.replace(".tif", ".nc")
        ds.to_netcdf(out_nc_fp)
        print(out_tif_fp)
        print(out_nc_fp)

