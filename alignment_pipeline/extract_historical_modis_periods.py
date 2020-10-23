"""
Extract the begin/end dates for the MODIS 8-day aggregate periods
"""

import os, glob, time, argparse, datetime
import gdal
import multiprocessing as mp
import pandas as pd
import numpy as np
from helpers import check_env


def get_ranges(fp):
    # elems = ["product", "date", "location", "version", "production_date"]
    # fn_dict = dict(zip(elems, os.path.basename(fn).split(".hdf")[0].split(".")))
    fn = os.path.basename(fp)
    # fn_dict["date"] = fn_dict["date"].replace("A", "")
    # fn_dict["fn"] = fn
    f = gdal.Open(fp)
    meta = f.GetMetadata_Dict()
    out = {
        "RANGEBEGINNINGDATE": meta["RANGEBEGINNINGDATE"],
        "RANGEBEGINNINGTIME": meta["RANGEBEGINNINGTIME"],
        "RANGEENDINGDATE": meta["RANGEENDINGDATE"],
        "RANGEENDINGTIME": meta["RANGEENDINGTIME"],
        "fn": fn
    }
    # out.update(fn_dict)
    f = None
    return out


if __name__ == "__main__":
    print("Extracting historical MODIS date ranges...")
    tic = time.perf_counter()

    env_ok = check_env(wrf_env_var=False)
    if not env_ok:
        exit("Environment variables incorrectly setup, check README for requirements")
    # args
    # parse args
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

    base_dir = os.getenv("MODIS_DIR")
    scratch_dir = os.getenv("SCRATCH_DIR")
#     files = sorted([
#         os.path.join(r, fn)
#         for r, s, files in os.walk(base_dir)
#         for fn in files
#         # if fn.endswith(".hdf") and "11A2" in fn
#         if fn.endswith(".hdf") and "MOD11A2" in fn
#     ])
    
    # determine ranges from MOD11A2 because we have determined 
    #  that MYD11A2 time ranges are a subset of MOD11A2
    files = sorted(
        glob.glob(os.path.join(base_dir, "MOD11A2", "*h11*.hdf"))
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
    # finding yday of last day of October is done with
    #  leap year to ensure inclusion of all possible periods starting
    #  in October
    end_doy = datetime.date(2020, 10, 31).timetuple().tm_yday
    df = df.loc[
        (df.index.dayofyear >= start_doy) & (df.index.dayofyear <= end_doy)
    ]

    out_fp = os.path.join(ancillary_dir, "historical_modis_range_metadata.csv")
            # ancillary_dir, "MODIS_LST_8dayComposite_begin_end_range_metadata.csv"
    df.to_csv(out_fp, index=False)
    print(f"Ranges saved to CSV, {round(time.perf_counter() - tic, 1)}s")
    print(f"Output path: {out_fp}\n")

    # subset it since we have posited that they use the same dates in each tile agg
    # new_df = df.groupby(["product", "date"]).apply(lambda x: x.iloc[0].drop("location"))
    # new_df.to_csv(
    #     os.path.join(
    #         ancillary_dir, "MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv"
    #     )
    # )

