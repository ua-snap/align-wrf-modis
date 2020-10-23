"""
Reproject the resampled 1km WRF TSK to epsg:3338
"""

import os, glob, subprocess, argparse, itertools, copy
import xarray as xr
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.windows import Window
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
            "{}_8Day_daytime_wrf_{}_max_{}.tif".format(
                variable, group, jdate
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
        },
    )
    return ds_new


def get_metadata(ds, wrf_fp, variable):
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


def reproject_wrf_to_3338(fp, out_fp):

    print(out_fp)
    _ = subprocess.call(
        ["gdalwarp", "-t_srs", "epsg:3338", "-q", "-overwrite", fp, out_fp]
    )
    # open array clipped to bounds and write
    with rio.open(out_fp, mode="r+") as out:
        arr = out.read(1)
        arr[arr == 0] = -9999
        out.write(arr, 1)
    return out_fp


def run_reproject(arr, meta, out_fp):
    fp = make_gtiff(arr, meta, out_fp)
    out_fp = copy.copy(fp).replace(".tif", "_3338.tif")
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

    # parse args
    parser = argparse.ArgumentParser(
        description="reproject the 8-day WRF data to EPSG:3338"
    )
    parser.add_argument(
        "-r",
        "--year_range",
        action="store",
        dest="year_range",
        help="WRF years to work on ('2000-2018', '2037-2047', or '2067-2077')",
    )
    # unpack the args here
    args = parser.parse_args()
    year_range = args.year_range
    # check period, set model groups
    model_groups = ["gfdl", "ccsm"]
    valid_ranges = ["2000-2018", "2037-2047", "2067-2077"]
    if year_range not in valid_ranges:
        exit("Invalid year range specified")
    elif year_range == "2000-2018":
        model_groups = ["era"] + model_groups

    # setup dirs
    variable = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    out_dir = os.path.join(
        # os.getenv("OUTPUT_DIR"), "WRF", "{}_1km_3338".format(variable)
        scratch_dir,
        "WRF",
        "{}_1km_3338-slim".format(variable),
    )
    # wrf_dir = os.path.join(scratch_dir, "WRF", "WRF_day_hours", variable)
    wrf_dir = os.path.join(scratch_dir, "WRF", "WRF_day_hours-slim", variable)
    # groups = get_model_groups(wrf_env_var)
    # sensors = ["MOD11A2", "MYD11A2"]
    # metrics = ["mean", "min", "max"]

    
    # for group, sensor, metric in itertools.product(groups, sensors, metrics):
    for group in model_groups:
        # print("working on", group, sensor, metric)
        print("Working on", group)

        # use subdir for each set of GeoTIFFs
        # set_dir = os.path.join(out_dir, "{}_{}_{}".format(group, metric, sensor))
        if (group in ["ccsm", "gfdl"]) & (year_range == "2000-2018"):
            begin_year, end_year = "2007", "2017"
        else:
            years = year_range.split("-")
            begin_year, end_year = years[0], years[1]

        set_dir = os.path.join(out_dir, f"{group}_max_{begin_year}-{end_year}")
        if not os.path.exists(set_dir):
            _ = os.makedirs(set_dir)
        # open modisified data
        wrf_fp = glob.glob(
            os.path.join(wrf_dir, f"*{group}_{begin_year}-{end_year}.nc")
        )[0]
        ds = xr.open_dataset(wrf_fp)
        # concatenate bands from multiband GeoTIFF
        ds = concat_bands(ds, variable)
        # get metadata
        meta = get_metadata(ds, wrf_fp, variable)
        # get output file paths
        out_fps = get_output_filepaths(ds, wrf_fp, set_dir)
        # make args
        da = ds[variable]
        args = [(a, meta, out_fp) for a, out_fp in zip(list(da.values), out_fps)]
        # serial processing
        out = [wrapper(x) for x in args]

