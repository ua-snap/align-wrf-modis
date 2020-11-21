"""Re-grid (downsample) WRF to the clipped MODIS"""


import argparse, datetime, glob, os, subprocess, time
import numpy as np
import pandas as pd
import rasterio as rio
import xarray as xr
from helpers import check_env
from multiprocessing import Pool


def warp_wrf(temp_arr, temp_meta, fp, out_fp):
    """Warp WRF to prepped MODIS grid raster with gdalwarp"""
    with rio.open(out_fp, "w", **temp_meta) as src:
        src.write(temp_arr)
    _ = subprocess.call(["gdalwarp", "-s_srs", "epsg:3338", "-q", fp, out_fp])
    with rio.open(out_fp, mode="r+") as src:
        arr = src.read(1)
        # values close to zero form from reasmpling
        arr[np.isclose(arr, 0)] = -9999
        src.write(arr, 1)
    return out_fp


def wrap_warp(args):
    """Wrapper for warping in parallel"""
    return warp_wrf(*args)


def read_band(fn, band=1):
    """Read a GeoTIFF band's data, for reprojected WRF"""
    with rio.open(fn) as rst:
        return rst.read(band).copy()


def get_dates(fps):
    """Get datetime array of dates from filenames, for NetCDFs"""
    jd_lst = [os.path.basename(fp).split("_")[-2] for fp in fps]
    dt = [datetime.datetime.strptime(jd_str, "%Y%j") for jd_str in jd_lst]
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
    wrf_var = "tsk"
    mod_var = "lst"
    scratch_dir = os.getenv("SCRATCH_DIR")
    out_base_dir = os.getenv("OUTPUT_DIR")
    template_fp = glob.glob(
        os.path.join(scratch_dir, "MODIS", "clipped", mod_var, "*/*")
    )[0]
    out_dir = os.path.join(scratch_dir, "WRF", "regridded", wrf_var)
    out_nc_dir = os.path.join(out_dir, "netcdf")
    if not os.path.exists(out_nc_dir):
        _ = os.makedirs(out_nc_dir)
    wrf_dir = os.path.join(scratch_dir, "WRF", "reprojected", wrf_var)
    # directory paths containing reprojected GeoTIFFs
    wrf_dps = glob.glob(os.path.join(wrf_dir, "*"))
    with rio.open(template_fp) as src:
        temp_arr = np.empty_like(src.read())
        temp_meta = src.meta
    # operate on each directory of GeoTIFFs
    for dp in wrf_dps:
        print(f"Working on {dp}", end="...")
        tic = time.perf_counter()
        wrf_fps = sorted(glob.glob(os.path.join(dp, "*")))
        # dir to store regridded files for particular period/model group
        bands_dn = os.path.basename(dp)
        bands_dir = os.path.join(out_dir, "bands", bands_dn)
        if not os.path.exists(bands_dir):
            _ = os.makedirs(bands_dir)
        warp_args = []
        for fp in wrf_fps:
            fn = os.path.basename(fp).replace("3338.tif", "regridded.tif")
            out_fp = os.path.join(bands_dir, fn)
            warp_args.append((temp_arr, temp_meta, fp, out_fp))
        pool = Pool(ncpus)
        bands_fps = pool.map(wrap_warp, warp_args)
        pool.close()
        pool.join()

        # assemble files into geotiffs and netcdf
        print("Making NetCDF", end="...")
        pool = Pool(ncpus)
        out_arrs = pool.map(read_band, bands_fps)
        pool.close()
        pool.join()
        bands_arr = np.array(out_arrs)
        with rio.open(bands_fps[0]) as tmp:
            new_meta = tmp.meta.copy()
            new_meta.update(count=bands_arr.shape[0])
            # for netcdf
            idx = np.arange(tmp.width)
            idy = np.arange(tmp.height)
            xc = tmp.xy(np.repeat(0, idx.shape), idx)[0]
            yc = tmp.xy(idy, np.repeat(0, idy.shape))[1]
        out_nc_fp = os.path.join(out_nc_dir, f"{wrf_var}_max_{bands_dn}_regridded.nc")
        dates = get_dates(bands_fps)
        ds = xr.Dataset(
            {wrf_var: (["date", "yc", "xc"], bands_arr)},
            coords={"xc": xc, "yc": yc, "date": dates,},
        )
        ds.to_netcdf(out_nc_fp)
        duration = round(time.perf_counter() - tic, 1)
        print(f"done, duration: {duration}s")
        print(f"Data saved to {out_nc_fp}")
