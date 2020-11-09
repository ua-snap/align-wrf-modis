# Compute deltas for future WRF TSK and save to sctatch dir

import itertools, os
import numpy as np
import pandas as pd
import xarray as xr
from helpers import check_env


def get_modis_period(np_dt):
    """
    determine the MODIS period for particular date,
    in terms of periods since 3-29
    """
    # origin date, year is arbitrary
    origin = pd.to_datetime("2007-03-29")
    pd_dt = pd.to_datetime(np_dt)
    datestrs = [f"2007-{str(dt.month).zfill(2)}-{dt.day}" for dt in pd_dt]
    date = pd.to_datetime(datestrs)
    return (date - origin).days // 8


def compute_clim(ds):
    """
    compute climatology for aligned WRF dataset
    """
    out_ds = ds.assign(date=lambda x: get_modis_period(ds.date.values))
    out_ds = out_ds.rename({"date": "modis_period"})
    out_da = out_ds.tsk.groupby("modis_period").mean()
    return out_da.where(out_da != -9999)


if __name__ == "__main__":
    # setup file pathing
    scratch_dir = os.getenv("SCRATCH_DIR")
    output_dir = os.getenv("OUTPUT_DIR")
    # eventually change to tsk_aligned/
    wrf_dir = os.path.join(scratch_dir, "WRF", "regridded")
    deltas_dir = os.path.join(scratch_dir, "WRF", "deltas")
    clim_dir = os.path.join(scratch_dir, "WRF", "climatologies")
    if not os.path.exists(deltas_dir):
        _ = os.makedirs(deltas_dir)
    if not os.path.exists(clim_dir):
        _ = os.makedirs(clim_dir)
    out_dir = os.path.join(output_dir, "WRF", "tsk_aligned_adjusted")
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)
    wrf_fn = "tsk_{}_max_{}.nc"
    wrf_fp = os.path.join(wrf_dir, wrf_fn)
    deltas_fp = os.path.join(deltas_dir, "tsk_{}_max_{}_deltas.nc")
    clim_fp = os.path.join(clim_dir, "tsk_{}_max_{}_climatology.nc")
    out_fp = os.path.join(out_dir, wrf_fn)
    # read ERA data and compute climatology. Subset to 2008-2017, to match GCM period.
    era_fp = os.path.join(wrf_dir, wrf_fp.format("era", "2000-2018"))
    with xr.open_dataset(era_fp) as ds:
        era_ds = ds.where(
            (ds.date >= np.datetime64("2008-03-30"))
            & (ds.date <= np.datetime64("2017-11-01")),
            drop=True,
        )
    era_clim_da = compute_clim(era_ds)
    era_clim_fp = clim_fp.format("era", "2008-2017")
    era_clim_da.to_netcdf(era_clim_fp)
    print(f"era climatology saved to {era_clim_fp}")

    gcms = ["gfdl", "ccsm"]
    year_ranges = ["2037-2047", "2067-2077"]
    wrf_period = "2008-2017"
    # for gcm, year_range in itertools.product(gcms, year_ranges):
    for gcm in gcms:
        # compute and save climatologies and deltas for GCM
        with xr.open_dataset(wrf_fp.format(gcm, "2007-2017")) as ds:
            bias_ds = ds.where(ds.date >= np.datetime64("2008-03-30"), drop=True)
        bias_clim_da = compute_clim(bias_ds)
        bias_clim_fp = clim_fp.format(gcm, wrf_period)
        bias_clim_da.to_netcdf(bias_clim_fp)
        print(f"Climatology for {gcm}, {wrf_period} saved to {bias_clim_fp}")
        # compute deltas, save
        deltas_da = bias_clim_da - era_clim_da
        gcm_deltas_fp = deltas_fp.format(gcm, wrf_period)
        deltas_da.to_netcdf(gcm_deltas_fp)
        print(f"Deltas for {gcm}, {wrf_period} saved to {gcm_deltas_fp}")
        # tile deltas to match shape (modis periods) of aligned unadjusted data
        deltas_arr = np.tile(deltas_da, (10, 1, 1))
        modis_period = np.tile(deltas_da.modis_period.values, 10)
        deltas_da = xr.DataArray(
            deltas_arr,
            coords=[modis_period, deltas_da.yc.values, deltas_da.xc.values],
            dims=["modis_period", "yc", "xc"],
        )

        for year_range in year_ranges:
            # apply deltas to aligned WRF data
            cut_str = "2038-03-30"
            if year_range == year_ranges[1]:
                cut_str = "2068-03-29"
            with xr.open_dataset(wrf_fp.format(gcm, year_range)) as ds:
                # filter out beginning periods of
                bias_ds = ds.where(ds.date >= np.datetime64(cut_str), drop=True)
            adj_ds = bias_ds.copy()
            adj_ds.tsk.values = bias_ds.tsk.values + deltas_da.values
            adj_fp = out_fp.format(gcm, year_range.replace("7-", "8-"))
            adj_ds.to_netcdf(adj_fp)
            print(
                f"Delta-adjusted aligned {gcm} data for {year_range} saved to {adj_fp}"
            )

