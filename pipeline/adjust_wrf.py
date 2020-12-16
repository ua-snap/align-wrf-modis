"""Compute and apply deltas for future WRF"""
# saves climatologies and deltas to $SCRATCH_DIR
# saves adjusted WRF to $OUTPUT_DIR

import itertools, os
import numpy as np
import pandas as pd
import xarray as xr
from helpers import check_env, add_metadata


def get_modis_period(np_dt):
    """determine the MODIS period for particular date"""
    # The MODIS period is an integer value represneting the
    #   number of 8-day MODIS periods since 3/29
    # origin date, year is arbitrary
    origin = pd.to_datetime("2007-03-29")
    pd_dt = pd.to_datetime(np_dt)
    datestrs = [f"2007-{str(dt.month).zfill(2)}-{dt.day}" for dt in pd_dt]
    date = pd.to_datetime(datestrs)
    return (date - origin).days // 8


def compute_clim(ds):
    """Compute climatology for an aligned WRF dataset"""
    out_ds = ds.assign(date=lambda x: get_modis_period(ds.date.values))
    out_ds = out_ds.rename({"date": "modis_period"})
    return out_ds.tsk.groupby("modis_period").mean()


if __name__ == "__main__":
    _ = check_env()
    # setup file pathing
    scratch_dir = os.getenv("SCRATCH_DIR")
    output_dir = os.getenv("OUTPUT_DIR")
    wrf_var = "tsk"
    wrf_dir = os.path.join(scratch_dir, "WRF", "reprojected", wrf_var, "netcdf")
    deltas_dir = os.path.join(scratch_dir, "WRF", "deltas", wrf_var)
    if not os.path.exists(deltas_dir):
        _ = os.makedirs(deltas_dir)
    clim_dir = os.path.join(scratch_dir, "WRF", "climatologies", wrf_var)
    if not os.path.exists(clim_dir):
        _ = os.makedirs(clim_dir)
    out_dir = os.path.join(output_dir, "aligned-WRF-MODIS", "WRF")
    if not os.path.exists(out_dir):
        _ = os.makedirs(out_dir)
    wrf_fn = "{}_max_{}_{}_{}.nc"
    wrf_fp = os.path.join(wrf_dir, wrf_fn)
    deltas_fp = os.path.join(deltas_dir, wrf_fn)
    clim_fp = os.path.join(clim_dir, wrf_fn)
    out_fp = os.path.join(out_dir, wrf_fn)
    
    # stuff for metadata
    gcm_lu = {"gfdl": "GFDL-CM3", "ccsm": "NCAR-CCSM4"}
    src_str = "{} downscaled via WRF 4.0"
    title = "1km MODIS-aligned WRF TSK"
    long_name = "Surface skin temperature"
    
    # read ERA data and compute climatology.
    #   Subset to 2008-2017, to match GCM period.
    era_fp = os.path.join(
        wrf_dir, wrf_fn.format(wrf_var, "era", "2000-2018", "reprojected")
    )
    # first, add metadata and save
    era_out_fp = os.path.join(
        out_dir, wrf_fn.format(wrf_var, "era", "2000-2018", "aligned")
    )
    with xr.open_dataset(era_fp) as ds:
        # ds for clim, susbet
        era_out_ds = add_metadata(
            ds, wrf_var, long_name, title, src_str.format("ERA-Interim")
        )
    era_out_ds.to_netcdf(era_out_fp)
    print(f"Adjusted ERA data saved to {era_out_fp}")
    
    # then, subset and compute climatology
    era_ds = era_out_ds.where(
        (era_out_ds.date >= np.datetime64("2008-03-29"))
        & (era_out_ds.date <= np.datetime64("2017-11-01")),
        drop=True,
    )
    era_clim_da = compute_clim(era_ds)
    era_clim_fp = clim_fp.format(wrf_var, "era", "2008-2017", "climatology")
    era_clim_da.to_netcdf(era_clim_fp)
    print(f"era climatology saved to {era_clim_fp}")
    
    # iterate delta adjustment over gcms
    gcms = ["gfdl", "ccsm"]
    year_ranges = ["2037-2047", "2067-2077"]
    wrf_period = "2008-2017"
    for gcm in gcms:
        # compute and save climatologies and deltas for GCM
        wrf_hist_fp = wrf_fp.format(wrf_var, gcm, "2007-2017", "reprojected")
        with xr.open_dataset(wrf_hist_fp) as ds:
            bias_ds = ds.where(ds.date >= np.datetime64("2008-03-29"), drop=True)
        bias_clim_da = compute_clim(bias_ds)
        bias_clim_fp = clim_fp.format(wrf_var, gcm, wrf_period, "climatology")
        bias_clim_da.to_netcdf(bias_clim_fp)
        print(f"Climatology for {gcm}, {wrf_period} saved to {bias_clim_fp}")
        
        # compute deltas, save
        deltas_da = bias_clim_da - era_clim_da
        gcm_deltas_fp = deltas_fp.format(wrf_var, gcm, wrf_period, "deltas")
        deltas_da.to_netcdf(gcm_deltas_fp)
        print(f"Deltas for {gcm}, {wrf_period} saved to {gcm_deltas_fp}")\
        
        # tile deltas to match shape (modis periods) of aligned unadjusted data,
        #   make new deltas data array
        deltas_arr = np.tile(deltas_da, (10, 1, 1))
        modis_period = np.tile(deltas_da.modis_period.values, 10)
        deltas_da = xr.DataArray(
            deltas_arr,
            coords=[modis_period, deltas_da.yc.values, deltas_da.xc.values],
            dims=["modis_period", "yc", "xc"],
        )
        
        # apply deltas to historical GCM dataset
        adj_ds = bias_ds.copy()
        adj_ds[wrf_var].values = bias_ds[wrf_var].values - deltas_da.values
        adj_ds = add_metadata(
            adj_ds, wrf_var, long_name, title, src_str.format(gcm_lu[gcm])
        )
        adj_fp = out_fp.format(wrf_var, gcm, wrf_period, "aligned")
        adj_ds.to_netcdf(adj_fp)
        print(f"Adjusted {gcm} data for {wrf_period} saved to {adj_fp}")
        
        # apply deltas to each future GCM dataset,
        #   indexed by year range
        for year_range in year_ranges:
            # for later filtering of unused WRF data
            cut_str = "2038-03-30"
            if year_range == year_ranges[1]:
                cut_str = "2068-03-29"
            with xr.open_dataset(
                wrf_fp.format(wrf_var, gcm, year_range, "reprojected")
            ) as ds:
                # read WRF and filter out beginning period
                bias_ds = ds.where(ds.date >= np.datetime64(cut_str), drop=True)
            adj_ds = bias_ds.copy()
            adj_ds[wrf_var].values = bias_ds[wrf_var].values - deltas_da.values
            adj_ds = add_metadata(
                adj_ds, wrf_var, long_name, title, src_str.format(gcm_lu[gcm])
            )
            adj_fp = out_fp.format(
                wrf_var, gcm, year_range.replace("7-", "8-"), "aligned"
            )
            adj_ds.to_netcdf(adj_fp)
            print(f"Adjusted {gcm} data for {year_range} saved to {adj_fp}")
