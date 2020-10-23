"""
Create the begin/end dates for the Future MODIS-based 8-day aggregate periods
"""


import os, datetime
import pandas as pd
import numpy as np
from helpers import check_env


def make_write_future_ranges(df, out_fp, periods):
    """
    Make and save data frame of future ranges for specified year
    """
    prev_date = df["RANGEBEGINNINGDATE"].values[-1]
    begin_time = df["RANGEBEGINNINGTIME"].values[-1]
    end_time = df["RANGEENDINGTIME"].values[-1]
    begin_range = pd.date_range(prev_date, "2078-01-08", freq="8D")
    begin_dates = pd.to_datetime(begin_range.values)
    end_dates = (begin_dates - pd.Timedelta(seconds=1))[1:]
    begin_dates = begin_dates[:-1]
    times_df = pd.DataFrame(
        {
            "RANGEBEGINNINGDATE": begin_dates.strftime("%Y-%m-%d"),
            "RANGEBEGINNINGTIME": begin_time,
            "RANGEENDINGDATE": end_dates.strftime("%Y-%m-%d"),
            "RANGEENDINGTIME": end_time,
        },
        index=begin_range[:-1],
    )

    # find start and end days-of-year for subsetting to relevant
    #  time of year for models (growing season: April - October)
    start_doy = datetime.date(2020, 3, 25).timetuple().tm_yday
    # finding yday of last day of October should be done with
    #  leap year to ensure inclusion of all possible periods starting
    #  in October
    end_doy = datetime.date(2020, 10, 31).timetuple().tm_yday
    times_df = times_df.loc[
        (times_df.index.dayofyear >= start_doy) & (times_df.index.dayofyear <= end_doy)
    ]

    # slice and save for each future period
    for years in periods:
        print("saving file as ", out_fp.format(years[0], years[1]))
        times_df[years[0] : years[1]].to_csv(out_fp.format(years[0], years[1]), index=False)


if __name__ == "__main__":
    env_ok = check_env(wrf_env_var=False)
    if not env_ok:
        exit("Environment variables incorrectly setup, check README for requirements")
    # read historical meta data to get last date of MODIS to continue series into future
    scratch_dir = os.getenv("SCRATCH_DIR")
    ancillary_dir = os.path.join(scratch_dir, "ancillary")
    meta_fp = os.path.join(
        ancillary_dir, "historical_modis_range_metadata.csv"
    )
    df = pd.read_csv(meta_fp)

    # Future WRF periods are 2037-07-02 - 2047-12-31, and 2067-07-02 - 2077-12-31
    future_periods = [("2037", "2047"), ("2067", "2077")]
    out_fp = os.path.join(ancillary_dir, "future_modis_range_metadata_{}-{}.csv")

    make_write_future_ranges(df, out_fp, future_periods)
