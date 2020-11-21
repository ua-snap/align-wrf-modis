"""Reproject the resampled WRF to epsg:3338"""

import argparse, glob, os, subprocess, time
import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
from datetime import datetime
from helpers import check_env
from multiprocessing import Pool
from rasterio.windows import Window


def get_output_filepaths(ds, wrf_fp, out_dir):
    """Make output filepaths for saving warped data"""
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
    """Concatenate separate bands produced from warping"""
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
    """Make metadata using resampled WRF"""
    rows, cols = ds[variable].shape[1:]
    with rio.open("netcdf:{}:Band1".format(wrf_fp)) as src:
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
    """make GeoTIFF band arr to gdalwarp"""
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
    """gdalwarp to epsg:3338"""
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
    """run the reprojection for a band"""
    fp = make_gtiff(arr, meta, out_fp)
    out_fp = fp.replace(".tif", "_3338.tif")
    _ = reproject_wrf_to_3338(fp, out_fp)
    os.unlink(fp)
    return out_fp


def wrap_reproject(args):
    """Wrapper for reprojecting in parallel"""
    return run_reproject(*args)


if __name__ == "__main__":
    _ = check_env()
    parser = argparse.ArgumentParser(
        description="Extract the date ranges of historical 8-day MODIS data"
    )
    parser.add_argument(
        "-n",
        "--ncpus",
        action="store",
        dest="ncpus",
        type=int,
        help="Number of cores to use with multiprocessing",
    )
    args = parser.parse_args()
    ncpus = args.ncpus
    # setup dirs
    variable = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    out_dir = os.path.join(scratch_dir, "WRF", "reprojected", variable)
    wrf_dir = os.path.join(scratch_dir, "WRF", "resampled", variable)
    periods = [(2007, 2017), (2037, 2047), (2067, 2077)]
    groups = ["era", "gfdl", "ccsm"] 
    periods = [periods[i] for i in [0,0,0,1,1,2,2]]
    groups = [groups[i] for i in [0,1,2,1,2,1,2]]
    for period, group in zip(periods, groups):
        if (period[0] == 2007) & (group == "era"):
            period = (2000, 2018)
        print(f"Working on {group}, {period}", end="...")
        tic = time.perf_counter()
        set_dir = os.path.join(out_dir, f"{group}_{period[0]}-{period[1]}")
        if not os.path.exists(set_dir):
            _ = os.makedirs(set_dir)
        # open modisified data
        wrf_fp = glob.glob(
            os.path.join(wrf_dir, f"*{group}_{period[0]}-{period[1]}.nc")
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
        # process in parallel
        args = [(a, meta, out_fp) for a, out_fp in zip(list(da.values), out_fps)]
        pool = Pool(ncpus)
        out = pool.map(wrap_reproject, args)
        pool.close()
        pool.join()
        duration = round(time.perf_counter() - tic, 1)
        print(f"done, duration: {duration}s")
        print(f"Reprojected files saved to {set_dir}")
