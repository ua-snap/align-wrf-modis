"""
Resample hourly WRF data to 8-day averages consistent with MODIS data
Operates on 1km WRF data developed for SERDP
"""


import os, glob, calendar, shutil, itertools
import numpy as np
import xarray as xr
import pandas as pd
import datetime
import multiprocessing as mp
import subprocess
import time
from helpers import check_env, get_model_groups


def get_files(wrf_dir, group, begin_year, end_year):
    """get filepaths to raw WRF output in wrf_dir based on year range"""
    # get all files from all years - directory structure important, see README!
    files = sorted(
        [
            fp
            for sublist in [
                glob.glob(os.path.join(wrf_dir, group, s, "*"))
                for s in [
                    s
                    for s in os.listdir(os.path.join(wrf_dir, group))
                    if s in [str(yr) for yr in np.arange(begin_year, end_year + 1)]
                ]
            ]
            for fp in sublist
        ]
    )
    # NOTE FOR FUTURE VERSIONS:
    # code below omits duplicate July 2 files (gfdl/ccsm), taking first 
    #   occurrence as desired file. This should be removed when duplicate files
    #   ar removed from raw data directories.
    fns = [os.path.basename(fp) for fp in files]
    unique_fns, idx = np.unique(fns, return_index=True)
    all_idx = np.arange(len(fns))
    rm_idx = np.setdiff1d(all_idx, idx) - 1
    [i for i in all_idx if i not in rm_idx]
    files = list(np.array(files)[idx])
    return files


def get_drop_vars(fp, keep_vars):
    """get list of variables to drop when reading 1km WRF"""
    ds = xr.open_dataset(fp)
    drop_vars = [var for var in list(ds.data_vars) if var not in keep_vars]
    ds.close()
    ds = None
    return drop_vars


def make_wrf_like_modis(
    ds_sel, ds_coords, meta_df, metric, product, variable
):
    """ resample (temporal) WRF hourlies to the MODIS ranges used in 8-day compositing. """
    # sort the values -- probably unnecessary, but is here.
    meta_df = meta_df.sort_values("date")
    # make the slices we need from the metadata file
    slices = meta_df.apply(
        lambda row: slice(
            row["RANGEBEGINNINGDATE"] + "T" + row["RANGEBEGINNINGTIME"],
            row["RANGEENDINGDATE"] + "T" + row["RANGEENDINGTIME"],
        ),
        axis=1,
    ).tolist()

    if metric == "mean":
        out_arr = np.array(
            [ds_sel.sel(time=sl).mean("time")[variable].values for sl in slices]
        )
    elif metric == "min":
        out_arr = np.array(
            [ds_sel.sel(time=sl).min("time")[variable].values for sl in slices]
        )
    elif metric == "max":
        out_arr = np.array(
            [ds_sel.sel(time=sl).max("time")[variable].values for sl in slices]
        )

    # make dataset from aggregate values and original lat/lon arrays
    # flip vertically to orient array intuitively
    new_times = pd.DatetimeIndex(
        meta_df["date"]
        .astype(str)
        .apply(lambda x: datetime.datetime.strptime(str(x), "%Y%j"))
        .tolist()
    )
    new_ds = xr.Dataset(
        {
            "lat": (["south_north", "west_east"], ds_coords.lat.data),
            "lon": (["south_north", "west_east"], ds_coords.lon.data),
            "aggr": (["time", "south_north", "west_east"], out_arr),
        },
        coords={
            "west_east": ds_coords.west_east.values,
            "south_north": ds_coords.south_north.values,
            "time": new_times,
        },
    )
    new_ds["aggr"].encoding = ds_sel[variable].encoding
    return new_ds


def make_warp_files(temp_dir, nt):
    """make VRT files needed for gdalwarp with -geoloc"""
    # make lat/lon VRT files
    coords = ["lat", "lon"]
    with open("WRF_Processing/vrt_template_latlon.txt", "r") as template_file:
        template = template_file.read()
        for coord in coords:
            vrt_file = open(os.path.join(temp_dir, "{}.vrt".format(coord)), "w")
            vrt_file.write(template.format(coord))
            vrt_file.close()
    # make dataset VRT file
    with open("WRF_Processing/vrt_template_meta.txt", "r") as template_file:
        meta_text = template_file.read()
    with open("WRF_Processing/vrt_template_band.txt", "r") as template_file:
        band_text = template_file.read()
    bands_text = "".join([band_text.format(i + 1, i + 1) for i in np.arange(nt)])
    vrt_text = meta_text.format(temp_dir, temp_dir, bands_text)
    vrt_fp = os.path.join(temp_dir, "aggr.vrt")
    f = open(vrt_fp, "w")
    f.write(vrt_text)
    f.close()
    return vrt_fp


def warp_with_geoloc(temp_dir, nt, out_fp):
    """use gdalwarp to rasterize using geolocation arrays present in 1k WRF, warping to 3338"""
    # first create necessary files
    make_warp_files(temp_dir, nt)
    _ = subprocess.call(
        [
            "gdalwarp",
            "-geoloc",
            "-t_srs",
            "-q",
            "EPSG:4326",
            "-overwrite",
            os.path.join(temp_dir, "aggr.vrt"),
            out_fp,
        ]
    )
    # cleanup files
    temp_files = [os.path.join(temp_dir, fn) for fn in os.listdir(temp_dir)]
    # _ = [os.unlink(fp) for fp in temp_files]
    # os.rmdir(temp_dir)


def modisify_wrf(
    files, temp_dir, out_dir, sensor, meta_df, metrics, variable, ncpus=7,
):
    """ make the WRF data temporally look like MODIS LST 8-Day composite Daytime"""

    def get_day_hours(month):
        return (month >= 9) & (month <= 17)

    def get_day_hours_ds(ds):
        ds_sel = ds.sel(Time=get_day_hours(ds["time.hour"]))
        return ds_sel

    # resample the data following begin / end times from the MODIS file metadata for 8-day composites
    print("reading all files as one dataset to aggregate", end="...")
    tic = time.perf_counter()
    # only need to read in certain variables
    drop_vars = get_drop_vars(files[0], [variable])
    # load without lat/lon vars initially
    ds_mf = xr.open_mfdataset(
        files,
        drop_variables=drop_vars,
        concat_dim="Time",
        combine="nested",
        preprocess=get_day_hours_ds,
    ).compute()
    print("done,", round(time.perf_counter() - tic, 1), "s")
    # re-shape dataset for aggregating by time
    arr = ds_mf[variable.upper()].values
    ds_sel = xr.Dataset(
        {variable: (["time", "south_north", "west_east"], arr)},
        coords={
            "west_east": ds_mf.west_east.values,
            "south_north": ds_mf.south_north.values,
            "time": ds_mf.time.values,
        },
    )
    # open file with lat/lon info
    drop_vars_coords = get_drop_vars(files[0], [variable, "lat", "lon"])
    ds_coords = xr.open_dataset(
        files[0], drop_variables=drop_vars_coords
    )
    # get model
    model = files[0].split("_")[0].split("/")[-3]
    # iterate through metrics
    for metric in metrics:
        fn = "{}_8Day_daytime_wrf_{}_{}_{}.nc".format(
            variable, model, metric, sensor
        )
        out_fp = os.path.join(out_dir, fn)

        ds_res = make_wrf_like_modis(
            ds_sel, ds_coords, meta_df, metric, sensor, variable
        )
        # need length of time dimension for warping later
        nt = len(ds_res["time"].values)

        # write to temp folder for gdalwarp, cleanup
        print("writing temp 8-Day netcdf, aggregate:{}".format(metric), end="...")
        tic = time.perf_counter()
        ds_res.to_netcdf(
            os.path.join(temp_dir, "aggr.nc"), mode="w", format="NETCDF4"
        )

        # move data timeseries to separate file
        out_ts_fp = out_fp.replace(".nc", "_times.csv")
        times = ds_res.time.to_index().to_frame()
        times.to_csv(out_ts_fp)

        ds_res.close()
        del ds_res
        print("done,", round(time.perf_counter() - tic, 1), "s")

        # rasterize using geolocation arrays present in WRF files
        print("warp with geoloc", end="...")
        tic = time.perf_counter()
        warp_with_geoloc(temp_dir, nt, out_fp)
        print("done,", round(time.perf_counter() - tic, 1), "s")
        print(out_fp)

return out_fp


def run(x):
    return modisify_wrf(**x)


if __name__ == "__main__":
    # check environment
    wrf_env_var = check_env()
    if not wrf_env_var:
        exit("Environment variables incorrectly setup, check README for requirements")
    print("env ok, wrf_env_var:", wrf_env_var)

    wrf_dir = os.getenv(wrf_env_var)
    variable = "tsk"
    scratch_dir = os.getenv("SCRATCH_DIR")
    temp_dir = os.path.join(scratch_dir, "temp")
    if not os.path.exists(temp_dir):
        _ = os.makedirs(temp_dir)
    out_dir = os.path.join(scratch_dir, "WRF", "WRF_day_hours", variable)
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)

    # read in a metadata file built from metadata entries within each MOD11A2/MYD11A2 HDF file
    #    --> see "MODIS_Processing/get_begin_end_range_8day_composite_MODIS_hdf_metadata.py" to build this file
    meta_fp = os.path.join(
        scratch_dir,
        "ancillary",
        "MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv",
    )
    if os.path.exists(meta_fp):
        print("reading in metadata", end="...")
        meta_df = pd.read_csv(meta_fp)
        print("done")
    else:
        exit("MODIS metadata not present, see README")

    meta_df["year"] = meta_df["date"].apply(lambda x: int(str(x)[:4]))

    sensors = ["MOD11A2", "MYD11A2"]
    metrics = ["mean", "max", "min"]
    # get names of models used in WRF based on the scale (1km v 20km)
    model_groups = get_model_groups(wrf_env_var)
    begin_year, end_year = 2000, 2019

    # implmenting without iterating over metirc
    # for sensor, metric, group in itertools.product(sensors, metrics, model_groups):
    for sensor, group in itertools.product(sensors, model_groups):
        print("Working on:", sensor, group)

        # subset the metadata DataFrame
        meta_df_sub = meta_df[meta_df["product"] == sensor].sort_values("date")

        # convert the RANGE TIMES to Pandas date objects for easier querying.
        meta_times_end = meta_df_sub["RANGEENDINGDATE"].apply(
            lambda x: pd.to_datetime(x + " 23:00:00", format="%Y-%m-%d %H:%M:%S")
        )
        meta_times_begin = meta_df_sub["RANGEBEGINNINGDATE"].apply(
            lambda x: pd.to_datetime(x + " 23:00:00", format="%Y-%m-%d %H:%M:%S")
        )

        # get filepaths
        print("getting filepaths to base data in", wrf_dir, end="...")
        files = get_files(wrf_dir, group, begin_year, end_year)
        print("done. \nBegin:", files[0], "\nEnd:", files[-1])

        # pull bounding times from file list
        drop_vars = get_drop_vars(files[0], [variable, "time"])
        with xr.open_mfdataset(
            [files[0], files[-1]],
            drop_variables=drop_vars,
            concat_dim="Time",
            combine="nested",
        ) as ds:
            times = ds.time.to_index()

        # slice meta_df_sub to the overlapping times
        meta_df_sub = meta_df_sub[
            (meta_times_begin > times[0]) & (meta_times_end < times[-1])
        ].copy(deep=True)
        print("date ranges:", meta_df_sub.iloc[0].date, meta_df_sub.iloc[-1].date)

        # run
        names = [
            "files",
            "temp_dir",
            "out_dir",
            "sensor",
            "metrics",
            "meta_df",
            "variable",
        ]
        args = dict(
            zip(
                names,
                [
                    files,
                    temp_dir,
                    out_dir,
                    sensor,
                    metrics,
                    meta_df_sub,
                    variable,
                ],
            )
        )
        print("running modisification")
        tic1 = time.perf_counter()
        out = run(args)
        print(
            sensor, group, "modisified,", round(time.perf_counter() - tic1, 1), "s",
        )

