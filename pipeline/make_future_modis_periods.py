"""
Create the begin/end dates for the Future MODIS-based 8-day aggregate periods
"""


import os, datetime
import pandas as pd
import numpy as np
from helpers import check_env


def make_write_future_ranges(temp_fp, period):
    """
    Make and save data frame of future ranges for specified years. 
    Note: this does not need to be based on historical period, because a 
    new 8-day period is started on the first of every year. 
    """
    def make_dates(year):
        """
        Make starting dates for a given year.
        Use end cutoff of 11/1. Using 10/31 as end date results in extra period
        for leap years. This final period will always be only days from November
        excep on leap years.
        """
        begin_dates = pd.date_range(f"{year}-01-01", f"{year}-11-01", freq="8D")
        return begin_dates.values

    years = np.arange(int(period[0]), int(period[1]) + 1)
    begin_dates = np.array([make_dates(year) for year in years]).flatten()
    begin_dates = pd.to_datetime(begin_dates)
    end_dates = begin_dates + pd.Timedelta(days=7)
    times_df = pd.DataFrame(
        {
            "RANGEBEGINNINGDATE": begin_dates.strftime("%Y-%m-%d"),
            "RANGEBEGINNINGTIME": "00:00:00",
            "RANGEENDINGDATE": end_dates.strftime("%Y-%m-%d"),
            "RANGEENDINGTIME": "23:59:59",
        },
        # index=begin_range[:-1],
        index=begin_dates
    )

    # find start day-of-year for subsetting to relevant
    # time of year for models (growing season starts in April)
    start_doy = datetime.date(2020, 3, 25).timetuple().tm_yday
    times_df = times_df.loc[(times_df.index.dayofyear >= start_doy)]
    out_fp = temp_fp.format(period[0], period[1])
    times_df.to_csv(out_fp, index=False)

    return out_fp


if __name__ == "__main__":
    env_ok = check_env()
    # read historical meta data to get last date of MODIS to continue series into future
    scratch_dir = os.getenv("SCRATCH_DIR")
    ancillary_dir = os.path.join(scratch_dir, "ancillary")
    # meta_fp = os.path.join(
    #     ancillary_dir, "historical_modis_range_metadata.csv"
    # )
    # df = pd.read_csv(meta_fp)

    # Future WRF periods are 2037-07-02 - 2047-12-31, and 2067-07-02 - 2077-12-31
    future_periods = [("2037", "2047"), ("2067", "2077")]
    temp_fp = os.path.join(ancillary_dir, "future_modis_range_metadata_{}-{}.csv")
    for period in future_periods:
        out_fp = make_write_future_ranges(temp_fp, period)
        print(f"{period} date ranges saved as {out_fp}")
