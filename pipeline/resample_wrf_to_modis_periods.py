"""Resample hourly WRF data to 8-day MODIS periods"""

import argparse, glob, os, subprocess, time
import numpy as np
import pandas as pd
import xarray as xr
from helpers import check_env
from multiprocessing import Pool


def get_files(wrf_dir, group, begin_year, end_year):
    """get filepaths to raw WRF output in wrf_dir within years"""
    # directory structure important, see README!
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
    return files


def get_drop_vars(fp, keep_vars):
    """get list of variables to drop when reading 1km WRF (optimization)"""
    ds = xr.open_dataset(fp)
    drop_vars = [var for var in list(ds.data_vars) if var not in keep_vars]
    ds.close()
    ds = None
    return drop_vars


def read_wrf(args):
    """read in the relevant data from a raw WRF files"""
    def get_day_hours(month):
        return (month >= 9) & (month <= 17)

    def get_day_hours_ds(ds):
        ds_sel = ds.sel(Time=get_day_hours(ds["time.hour"]))
        return ds_sel

    fp, drop_vars = args[0], args[1]
    with xr.open_dataset(fp, drop_variables=drop_vars) as ds:
        return get_day_hours_ds(ds)


def make_wrf_like_modis(ds_sel, ds_coords, meta_df, variable):
    """resample (temporal) WRF hourly data (maximum) to the 8-day MODIS periods"""
    # sort the values (probably unnecessary)
    meta_df = meta_df.sort_values("RANGEBEGINNINGDATE")
    # make the slices we need from the metadata file
    slices = meta_df.apply(
        lambda row: slice(
            row["RANGEBEGINNINGDATE"] + "T" + row["RANGEBEGINNINGTIME"],
            row["RANGEENDINGDATE"] + "T" + row["RANGEENDINGTIME"],
        ),
        axis=1,
    ).tolist()
    out_arr = np.array(
        [ds_sel.sel(time=sl).max("time")[variable].values for sl in slices]
    )
    new_times = pd.DatetimeIndex(meta_df["RANGEBEGINNINGDATE"].values)
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
    with open(
        "pipeline/vrt_template/vrt_template_latlon.txt", "r"
    ) as template_file:
        template = template_file.read()
        for coord in coords:
            vrt_file = open(os.path.join(temp_dir, "{}.vrt".format(coord)), "w")
            vrt_file.write(template.format(coord))
            vrt_file.close()
    # make dataset VRT file
    with open(
        "pipeline/vrt_template/vrt_template_meta.txt", "r"
    ) as template_file:
        meta_text = template_file.read()
    with open(
        "pipeline/vrt_template/vrt_template_band.txt", "r"
    ) as template_file:
        band_text = template_file.read()
    bands_text = "".join([band_text.format(i + 1, i + 1) for i in np.arange(nt)])
    vrt_text = meta_text.format(temp_dir, temp_dir, bands_text)
    vrt_fp = os.path.join(temp_dir, "aggr.vrt")
    f = open(vrt_fp, "w")
    f.write(vrt_text)
    f.close()
    return vrt_fp


def warp_with_geoloc(temp_dir, nt, out_fp):
    """use gdalwarp to rasterize using geolocation arrays, output to epsg:4326"""
    # first create necessary files
    vrt_fp = make_warp_files(temp_dir, nt)
    # run
    _ = subprocess.call(
        [
            "gdalwarp",
            "-geoloc",
            "-t_srs",
            "EPSG:4326",
            "-q",
            "-overwrite",
            vrt_fp,
            out_fp,
        ]
    )
    return None


def modisify_wrf(files, temp_dir, out_dir, meta_df, variable, ncpus):
    """read WRF data, resample to MODIS periods, rasterize"""
    # resample the data following begin / end times from the extracted ranges
    # get model
    model = files[0].split("_")[0].split("/")[-3]
    # only need to read in certain variables (ignore lat, lon here)
    drop_vars = get_drop_vars(files[0], [variable.upper(), "time"])
    # make args for reading data via Pool
    args = [(file, drop_vars) for file in files]
    # load data, concatenate
    p = Pool(ncpus)
    ds_lst = p.map(read_wrf, args)
    p.close()
    p.join()
    ds_mf = xr.concat(ds_lst, dim="Time")
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
    drop_vars_coords = get_drop_vars(files[0], [variable.upper(), "lat", "lon"])
    ds_coords = xr.open_dataset(files[0], drop_variables=drop_vars_coords)
    # make output filepath
    period = "-".join([meta_df["RANGEBEGINNINGDATE"].values[i][:4] for i in (0, -1)])
    fn = "{}_max_8Day_daytime_wrf_{}_{}.nc".format(variable, model, period)
    out_fp = os.path.join(out_dir, fn)
    print(f"resampling WRF {model}")
    ds_res = make_wrf_like_modis(ds_sel, ds_coords, meta_df, variable)
    # need length of time dimension for warping later
    nt = len(ds_res["time"].values)
    # write to temp folder for gdalwarp, cleanup
    ds_res.to_netcdf(os.path.join(temp_dir, "aggr.nc"))
    ds_res.close()
    # move data timeseries to separate file
    out_ts_fp = out_fp.replace(".nc", "_times.csv")
    times = ds_res.time.to_index().to_frame()
    times.to_csv(out_ts_fp)
    del ds_res

    # rasterize using geolocation arrays present in WRF files
    _ = warp_with_geoloc(temp_dir, nt, out_fp)
    return out_fp


if __name__ == "__main__":
    # check environment
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
    wrf_dir = os.getenv("WRF_DIR")
    scratch_dir = os.getenv("SCRATCH_DIR")
    temp_dir = os.path.join(scratch_dir, "temp")
    if not os.path.exists(temp_dir):
        _ = os.makedirs(temp_dir)
    variable = "tsk"
    out_dir = os.path.join(scratch_dir, "WRF", "resampled", variable)
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)
    periods = [(2000, 2018), (2037, 2047), (2067, 2077)]
    groups = ["era", "gfdl", "ccsm"] 
    periods = [periods[i] for i in [0,0,0,1,1,2,2]]
    groups = [groups[i] for i in [0,1,2,1,2,1,2]]
    meta_fn_mem = None
    for period, group in zip(periods, groups):
        # read metadata
        if period[0] == 2000:
            meta_fn = "historical_modis_range_metadata.csv"
        else:
            meta_fn = f"future_modis_range_metadata_{period[0]}-{period[1]}.csv"
        if meta_fn_mem != meta_fn:
            meta_fp = os.path.join(scratch_dir, "ancillary", meta_fn)
            meta_df = pd.read_csv(meta_fp)
            meta_fn_mem = meta_fn
        print(f"Working on: {group}, {period}", end="...")
        tic = time.perf_counter()
        # convert the RANGE TIMES to Pandas date objects for easier querying.
        meta_times_end = meta_df["RANGEENDINGDATE"].apply(
            lambda x: pd.to_datetime(x + " 23:00:00", format="%Y-%m-%d %H:%M:%S")
        )
        meta_times_begin = meta_df["RANGEBEGINNINGDATE"].apply(
            lambda x: pd.to_datetime(x + " 23:00:00", format="%Y-%m-%d %H:%M:%S")
        )
        # get raw WRF filepaths
        files = get_files(wrf_dir, group, period[0], period[1])
        # pull bounding times from file list
        drop_vars = get_drop_vars(files[0], ["time"])
        with xr.open_mfdataset(
            [files[0], files[-1]],
            drop_variables=drop_vars,
            concat_dim="Time",
            combine="nested",
        ) as ds:
            times = ds.time.to_index()
        # subset meta_df to the times overlapping available WRF
        meta_df_sub = meta_df[
            (meta_times_begin > times[0]) & (meta_times_end < times[-1])
        ].copy(deep=True)
        out_fp = modisify_wrf(files, temp_dir, out_dir, meta_df_sub, variable, ncpus)
        duration = round(time.perf_counter() - tic, 1)
        print(f"modisified, duration: {duration}s")
        print(f"Resampled WRF saved to {out_fp}")
