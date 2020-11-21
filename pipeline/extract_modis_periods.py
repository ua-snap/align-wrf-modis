"""Extract dates for the MODIS 8-day aggregate periods, hisotrical and future"""


import os, glob, time, argparse, datetime
import gdal
import multiprocessing as mp
import pandas as pd
import numpy as np
from helpers import check_env


def get_ranges(fp):
    """get date ranges for single MODIS file"""
    fn = os.path.basename(fp)
    f = gdal.Open(fp)
    meta = f.GetMetadata_Dict()
    out = {
        "RANGEBEGINNINGDATE": meta["RANGEBEGINNINGDATE"],
        "RANGEBEGINNINGTIME": meta["RANGEBEGINNINGTIME"],
        "RANGEENDINGDATE": meta["RANGEENDINGDATE"],
        "RANGEENDINGTIME": meta["RANGEENDINGTIME"],
        "fn": fn
    }
    f = None
    return out


def get_future_ranges(period):
    """Make data frame of future ranges for specified period""" 
    # Note: this does not need to be based on historical period, because a 
    #   new 8-day period is started on the first of every year. 
    def make_dates(year):
        """Make starting dates for a given year"""
        # Use end cutoff of 11/1. Using 10/31 as end date results in extra period
        #   for leap years. This final period will always be only days from November
        #   except on leap years, where it would contain 10/31.
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
        index=begin_dates
    )
    # find start day-of-year for subsetting to relevant
    #   time of year for models (growing season starts in April)
    start_doy = datetime.date(2020, 3, 25).timetuple().tm_yday
    return times_df.loc[(times_df.index.dayofyear >= start_doy)]


if __name__ == "__main__":
    _ = check_env()
    # 1) extract historical ranges
    print("Extracting historical MODIS date ranges...")
    tic = time.perf_counter()
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
    modis_dir = os.getenv("MODIS_DIR")
    scratch_dir = os.getenv("SCRATCH_DIR")
    files = sorted(
        glob.glob(os.path.join(modis_dir, "MOD11A2", "*h11*.hdf"))
    )
    # Setup output dir
    ancillary_dir = os.path.join(scratch_dir, "ancillary")
    if not os.path.exists(ancillary_dir):
        _ = os.makedirs(ancillary_dir)
    # run it
    pool = mp.Pool(ncpus)
    out = pool.map(get_ranges, files)
    pool.close()
    pool.join()
    # make a table out of it and subset to relevant periods (growing season)
    df = pd.DataFrame(out)
    df = df.set_index(pd.to_datetime(df["RANGEBEGINNINGDATE"]))
    # find start and end days-of-year for subsetting to relevant
    #  time of year for models (growing season: April - October)
    start_doy = datetime.date(2020, 3, 25).timetuple().tm_yday
    # finding yday of last day of October is done using
    #  leap year to ensure inclusion of all possible periods starting
    #  in October
    end_doy = datetime.date(2020, 10, 31).timetuple().tm_yday
    df = df.loc[
        (df.index.dayofyear >= start_doy) & (df.index.dayofyear <= end_doy)
    ]
    out_fp = os.path.join(ancillary_dir, "historical_modis_range_metadata.csv")
    df.to_csv(out_fp, index=False)
    duration = round(time.perf_counter() - tic, 1)
    print(f"Historical ranges saved as {out_fp},\n duration: {duration}s")

    # 2) create future ranges
    # Future WRF periods are 2037-07-02 - 2047-12-31, and 2067-07-02 - 2077-12-31
    print("Creating future ranges...")
    future_periods = [("2037", "2047"), ("2067", "2077")]
    temp_fp = os.path.join(ancillary_dir, "future_modis_range_metadata_{}-{}.csv")
    for period in future_periods:
        out_df = get_future_ranges(period)
        out_fp = temp_fp.format(period[0], period[1])
        out_df.to_csv(out_fp, index=False)
        print(f"Future date ranges for {period[0]}-{period[1]} saved as {out_fp}")
