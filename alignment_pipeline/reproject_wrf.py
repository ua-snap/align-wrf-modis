"""
Reproject the 1km WRF TSK resampled to 8-day means, 
producing the final WRF data
"""

import xarray as xr
import numpy as np
import pandas as pd
import rasterio as rio
import os, glob, subprocess, itertools
from helpers import check_env, get_model_groups
from datetime import datetime


def get_output_filepaths(ds, wrf_fp, out_dir):
    times_fp = wrf_fp.replace(".nc", "_times.csv")
    times = pd.read_csv(times_fp)
    times = pd.to_datetime(times.time)
    jdates = list(times.map(lambda x: x.strftime("%Y%j")))
    out_fps = [
        os.path.join(
            out_dir,
            "{}_8Day_daytime_wrf_{}_{}_{}_{}.tif".format(
                variable, group, metric, sensor, jdate
            ),
        )
        for jdate in jdates
    ]
    return out_fps


def concat_bands(ds, variable):
    # concatenate separate bands from warping 1km WRF
    arr = []
    varnames = list(ds.variables.keys())
    # 3 non-band variabls
    for i in np.arange(len(varnames) - 3):
        arr.append(ds.variables[varnames[i]].values)
    arr = np.flip(np.array(arr), axis=1)
    arr[arr == 0] = -9999
    ds_new = xr.Dataset(
        {variable: (["time", "yc", "xc"], arr)},
        coords={
            "xc": ("xc", ds["lon"].values),
            "yc": ("yc", np.flip(ds["lat"].values)),
            "time": np.arange(len(varnames) - 3),
        }
    )
    return ds_new


def get_metadata(ds, wrf_fp, wrf_env_var, variable):
    if wrf_env_var == "WRF_20KM_DIR":
        da = ds[variable]
        bands, rows, cols = da.shape
        res = (20000.0, 20000.0)
        transform = rio.transform.from_origin(
            ds.xc.data.min() - (res[0] / 2.0),
            ds.yc.data.max() + (res[0] / 2.0),
            res[0],
            res[1],
        )
        proj4 = "+a=6370000 +b=6370000 +k=1 +lat_0=90 +lat_ts=64 +lon_0=-152 +no_defs +proj=stere +units=m +x_0=0 +y_0=0"
        crs = rio.crs.CRS.from_string(proj4)
    elif wrf_env_var == "WRF_1KM_DIR":
        rows, cols = ds[variable].shape[1:]
        src = rio.open("netcdf:{}:Band1".format(wrf_fp))
        transform = src.transform
        crs = src.crs
    meta = {
        "count": 1,
        "crs": crs,
        "driver": "GTiff",
        "dtype": "float32",
        "height": rows,
        "nodata": -9999,
        "transform": transform,
        "width": cols,
    }
    return meta


def make_gtiff(arr, meta, out_fp):
    shape = arr.shape
    if len(shape) == 2:
        count = 1
        arr = arr[np.newaxis, ...]
    else:
        count = shape[0]
    meta.update(count=count)
    with rio.open(out_fp, "w", **meta) as out:
        out.write(arr)
    return out_fp


def reproject_wrf_to_modis_20km(fn, out_fn):
    print(out_fn)
    _ = subprocess.call(["gdalwarp", "-q", "-r", "lanczos", fn, out_fn])
    with rio.open(out_fn, mode="r+") as out:
        arr = out.read(1)
        arr[arr == 0] = -9999
        out.write(arr, 1)
    return out_fn


def reproject_wrf_to_3338(fp, out_fp):
    print(out_fp)
    _ = subprocess.call(["gdalwarp", "-t_srs", "epsg:3338", "-q", "-overwrite", fp, out_fp])
    with rio.open(out_fp, mode="r+") as out:
        arr = out.read(1)
        arr[arr == 0] = -9999
        out.write(arr, 1)
    return out_fp


def run_reproject(arr, meta, out_fp, template_fp, wrf_env_var):
    import copy

    fp = make_gtiff(arr, meta, out_fp)
    out_fp = copy.copy(fp).replace(".tif", "_3338.tif")

    if wrf_env_var == "WRF_20KM_DIR":
        with rio.open(template_fp) as tmp:
            arr = np.empty_like(tmp.read())
            tmp_meta = tmp.meta
            tmp_meta.update(compress="lzw")
        with rio.open(out_fp, "w", **tmp_meta) as rst:
            rst.write(arr)
        _ = reproject_wrf_to_modis_20km(fn, out_fp)
    elif wrf_env_var == "WRF_1KM_DIR":
        _ = reproject_wrf_to_3338(fp, out_fp)
    os.unlink(fp)
    return out_fp


def wrapper(x):
    return run_reproject(*x)


def open_raster(fn, band=1):
    with rio.open(fn) as rst:
        arr = rst.read(band).copy()
    rst = None
    return arr


if __name__ == "__main__":
    # check environment
    wrf_env_var = check_env()
    if not wrf_env_var:
        exit("Environment variables incorrectly setup, check README for requirements")
    print("env ok, wrf_env_var:", wrf_env_var)

    # import multiprocessing as mp

    # setup dirs
    variable = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    if wrf_env_var == "WRF_20KM_DIR":
        out_dir = os.path.join(
            os.getenv("OUTPUT_DIR"), "WRF", "{}_20km_3338_lanczos".format(variable)
        )
        template_fp = glob.glob(os.path.join(scratch_dir,'MODIS','modis_20km_3338_lanczos', '*.tif'))[0]
    elif wrf_env_var == "WRF_1KM_DIR":
        out_dir = os.path.join(
            os.getenv("OUTPUT_DIR"), "WRF", "{}_1km_3338".format(variable)
        )
        template_fp = ""

    out_dir_multi = out_dir + "_multiband"
    if not os.path.exists(out_dir_multi):
        _ = os.makedirs(out_dir_multi)
    
    wrf_dir = os.path.join(scratch_dir, "WRF", "WRF_day_hours", variable)
    # out_dir = os.path.join(scratch_dir, "WRF", "WRF_day_hours", variable)

    # base_path = '/workspace/Shared/Tech_Projects/SERDP_Fish_Fire/project_data'

    # out_path = os.path.join(out_dir, "WRF" )
    # out_dir =
    # template_fp = glob.glob(os.path.join(base_path,'MODIS_DATA','modis_20km_3338_lanczos', '*.tif'))[0]
    # groups = ['ERA-Interim_historical','GFDL-CM3_historical','GFDL-CM3_rcp85','NCAR-CCSM4_historical','NCAR-CCSM4_rcp85',]
    groups = get_model_groups(wrf_env_var)
    sensors = ["MOD11A2", "MYD11A2"]
    metrics = ["mean", "min", "max"]

    for group, sensor, metric in itertools.product(groups, sensors, metrics):
        print("working on", group, sensor, metric)
        # use subdir for each set of gtiffs
        set_dir = os.path.join(out_dir, "{}_{}_{}".format(group, metric, sensor))
        if not os.path.exists(set_dir):
            _ = os.makedirs(set_dir)    
        # open modisified data
        wrf_fp = glob.glob(
            os.path.join(wrf_dir, "*{}_{}_{}.nc".format(group, metric, sensor))
        )[0]
        ds = xr.open_dataset(wrf_fp)

        if wrf_env_var == "WRF_1KM_DIR":
            # in case of 1km WRF, the bands (dates) are separate subdatasets from gdalwarp
            print("concatenating bands...", end="")
            ds = concat_bands(ds, variable) 
            print("done")            

        meta = get_metadata(ds, wrf_fp, wrf_env_var, variable)
        out_fps = get_output_filepaths(ds, wrf_fp, set_dir, wrf_env_var)
        # make args
        da = ds[variable]
        args = [
            (a, meta, out_fp, template_fp, wrf_env_var)
            for a, out_fp in zip(list(da.values), out_fps)
        ]

        # serial processing
        # out = [wrapper(x) for x in args]

        # # multicore process is giving errors, but this is how it was deployed. somethings up.
        # pool = mp.Pool( 32 )
        # out = pool.map( wrapper, args )
        # pool.close()
        # pool.join()

        # now put these new GTiffs into a (1) Geotiff
        bands_fps = sorted(glob.glob(os.path.join(set_dir, "*.tif")))
        bands_arr = np.array([open_raster(fn) for fn in bands_fps])

        with rio.open(bands_fps[0]) as tmp:
            new_meta = tmp.meta.copy()
            new_meta.update(compress="lzw", count=bands_arr.shape[0])
            # coords for mb ds
            yc = [tmp.xy(row, 0)[1] for row in np.arange(new_meta["height"])]
            xc = [tmp.xy(0, col)[0] for col in np.arange(new_meta["width"])]

        out_mb_fp = os.path.join(
            out_dir_multi,
            "{}_8Day_daytime_wrf_{}_{}_{}_3338_multiband.tif".format(
                variable, group, metric, sensor
            ),
        )
        # with rio.open(
        #     out_mb_fp,
        #     "w",
        #     **new_meta,
        # ) as out:
        #     out.write(bands_arr)

        # and (2) netcdf
        print("writing netcdf: {}".format(out_mb_fp.replace(".tif", ".nc")), end="...")
        dates = np.array([datetime.strptime(fp.split("_")[-2], "%Y%j") for fp in bands_fps])
        mb_ds = xr.Dataset({variable: (["date", "yc", "xc"], bands_arr)},
            coords={"date": dates, "yc": yc, "xc": xc})
        mb_ds.to_netcdf(out_mb_fp.replace(".tif", ".nc"), "w", "NETCDF4")      
        print("done")

        tmp = None
        out = None
        del bands_arr, tmp, out


