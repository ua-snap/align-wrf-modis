"""Re-project WRF (warp and downsample) to match clipped MODIS"""

import argparse, glob, os, subprocess, time
import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
from datetime import datetime
from helpers import check_env
from multiprocessing import Pool


def make_bands_fps(wrf_fp, wrf_var, group, bands_dir):
    """Make output filepaths for saving warped data"""
    times_fp = wrf_fp.replace(".nc", "_times.csv")
    dates = pd.read_csv(times_fp).time.values
    bands_fps = [
        os.path.join(bands_dir, f"{wrf_var}_max_{group}_{date}_resampled.tif")
        for date in dates
    ]
    return bands_fps


def reproject_wrf(band_arr, fp, temp_meta, band_meta, temp_arr):
    """Warp WRF to prepped MODIS grid raster with gdalwarp"""
    # modify meta for band GeoTIFF, write
    band_arr[band_arr == 0] = -9999
    with rio.open(fp, "w", **band_meta) as src:
        src.write(band_arr, 1)
    # write blank template array to GeoTIFF
    out_fp = fp.replace("resampled", "reprojected")
    with rio.open(out_fp, "w", **temp_meta) as src:
        src.write(temp_arr)
    _ = subprocess.call(["gdalwarp", "-s_srs", "epsg:4326", "-q", fp, out_fp])
    return out_fp


def wrap_reproject(args):
    """Wrapper for warping in parallel"""
    return reproject_wrf(*args)


def read_band(fn, band=1):
    """Read a GeoTIFF band's data, for reprojected WRF"""
    with rio.open(fn) as src:
        return src.read(band).copy()


def get_dates(fps):
    """Get datetime array of dates from filenames, for NetCDFs"""
    dates = [os.path.basename(fp).split("_")[-2] for fp in fps]
    dt = [datetime.strptime(date, "%Y-%m-%d") for date in dates]
    dt_arr = np.array(dt, dtype="datetime64")
    return dt_arr


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
    wrf_var = "tsk"
    mod_var = "lst"
    scratch_dir = os.getenv("SCRATCH_DIR")
    modis_dir = os.path.join(scratch_dir, "MODIS", "clipped", mod_var)
    res_dir = os.path.join(scratch_dir, "WRF", "resampled", wrf_var)
    out_dir = os.path.join(scratch_dir, "WRF", "reprojected", wrf_var)
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)
    out_nc_dir = os.path.join(out_dir, "netcdf")
    if not os.path.exists(out_nc_dir):
        _ = os.makedirs(out_nc_dir)
    template_fp = glob.glob(os.path.join(modis_dir, "*/*"))[0]
    # get clipped MODIS metadata
    with rio.open(template_fp) as src:
        temp_meta = src.meta
        temp_arr = np.empty_like(src.read())
    res_fps = sorted(glob.glob(os.path.join(res_dir, "*.nc")))
    # get metadata for the WRF bands
    with rio.open(f"netcdf:{res_fps[0]}:Band1") as src:
        band_meta = src.meta
    band_meta.update({"driver": "GTiff", "nodata": -9999})
    for fp in res_fps:
        # set up file paths
        fn_meta = os.path.basename(fp).split("_")
        period = fn_meta[-1].split(".")[0]
        group = fn_meta[-2]
        print(f"Working on {group}, {period}", end="...")
        tic = time.perf_counter()
        # for writing the intermediate bands in current projection
        bands_dir = os.path.join(res_dir, "bands", f"{group}_{period}")
        if not os.path.exists(bands_dir):
            _ = os.makedirs(bands_dir)
        bands_fps = make_bands_fps(fp, wrf_var, group, bands_dir)
        # for writing reprojected bands
        out_bands_dir = os.path.join(out_dir, "bands", f"{group}_{period}")
        if not os.path.exists(out_bands_dir):
            _ = os.makedirs(out_bands_dir)
        # combine bands with fps and meta
        with xr.open_dataset(fp) as ds:
            args = [
                (ds[band].values, fp, temp_meta, band_meta, temp_arr)
                for band, fp in zip(ds.variables, bands_fps)
            ]
        # run reprojection
        pool = Pool(ncpus)
        out_bands_fps = pool.map(wrap_reproject, args)
        pool.close()
        pool.join()

        # assemble files into NetCDF
        print("Making NetCDF", end="...")
        pool = Pool(ncpus)
        bands = pool.map(read_band, out_bands_fps)
        pool.close()
        pool.join()
        bands_arr = np.array(bands)
        with rio.open(out_bands_fps[0]) as src:
            idx = np.arange(src.width)
            idy = np.arange(src.height)
            xc = src.xy(np.repeat(0, idx.shape), idx)[0]
            yc = src.xy(idy, np.repeat(0, idy.shape))[1]
        dates = get_dates(out_bands_fps)
        ds = xr.Dataset(
            {wrf_var: (["date", "yc", "xc"], bands_arr)},
            coords={"xc": xc, "yc": yc, "date": dates,},
        )
        out_nc_fp = os.path.join(
            out_nc_dir, f"{wrf_var}_max_{group}_{period}_reprojected.nc"
        )
        ds.to_netcdf(out_nc_fp)
        duration = round(time.perf_counter() - tic, 1)
        print(f"done, duration: {duration}s")
        print(f"Data saved to {out_nc_fp}")
