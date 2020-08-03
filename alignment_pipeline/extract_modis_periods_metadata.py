"""
Extract the begin/end dates for the MODIS 8-day aggregate periods
"""

import os
import gdal
import multiprocessing as mp
import pandas as pd
import numpy as np
from helpers import check_env


def get_ranges(fn):
    elems = ["product", "date", "location", "version", "production_date"]
    fn_dict = dict(zip(elems, os.path.basename(fn).split(".hdf")[0].split(".")))
    fn_dict["date"] = fn_dict["date"].replace("A", "")
    fn_dict["fn"] = fn
    f = gdal.Open(fn)
    meta = f.GetMetadata_Dict()
    out = {
        "RANGEBEGINNINGDATE": meta["RANGEBEGINNINGDATE"],
        "RANGEBEGINNINGTIME": meta["RANGEBEGINNINGTIME"],
        "RANGEENDINGDATE": meta["RANGEENDINGDATE"],
        "RANGEENDINGTIME": meta["RANGEENDINGTIME"],
    }
    out.update(fn_dict)
    f = None
    return out


if __name__ == "__main__":
    env_ok = check_env(wrf_env_var=False)
    if not env_ok:
        exit("Environment variables incorrectly setup, check README for requirements")
    # args
    base_dir = os.getenv("MODIS_DIR")
    scratch_dir = os.getenv("SCRATCH_DIR")
    files = [
        os.path.join(r, fn)
        for r, s, files in os.walk(base_dir)
        for fn in files
        if fn.endswith(".hdf") and "11A2" in fn
    ]

    # Setup output dir
    ancillary_dir = os.path.join(scratch_dir, "ancillary")
    if not os.path.exists(ancillary_dir):
        _ = os.makedirs(ancillary_dir)

    # run it
    pool = mp.Pool(15)
    out = pool.map(get_ranges, files)
    pool.close()
    pool.join()

    # make a table out of it and dump to disk...
    df = pd.DataFrame(out)
    df.to_csv(
        os.path.join(
            ancillary_dir, "MODIS_LST_8dayComposite_begin_end_range_metadata.csv"
        )
    )

    # subset it since we have posited that they use the same dates in each tile agg
    new_df = df.groupby(["product", "date"]).apply(lambda x: x.iloc[0].drop("location"))
    new_df.to_csv(
        os.path.join(
            ancillary_dir, "MODIS_LST_8dayComposite_begin_end_range_metadata_final.csv"
        )
    )

