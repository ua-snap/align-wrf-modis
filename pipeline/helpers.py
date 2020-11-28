"""Helper functions for MODIS LST-WRF TSK alignment"""

import datetime, os
import numpy as np
from pyproj import CRS


def check_env():
    """Check for all env vars, exit with message if not present."""
    key_vars = ["WRF_DIR", "MODIS_DIR", "SCRATCH_DIR", "OUTPUT_DIR"]
    env_vars = os.environ
    absent = [var not in env_vars for var in key_vars]
    if any(absent):
        absent_idx = np.where(absent)
        miss_vars = [key_vars[idx] for idx in absent_idx[0]]
        exit(f"The following env vars are missing: {miss_vars}")
    return True


def convert_date(dt_arr, epoch):
    """Convert datetime64 array to days (integer) since epoch"""
    return (dt_arr - epoch).astype("timedelta64[D]").astype(np.int32)


def add_metadata(ds, var, var_name, title, source, epoch):
    """Add metadata adhering to CF convention standards"""

    def cs_attrs(c):
        """coordinate system attributes dict"""
        return {
            "standard_name": f"projection_{c}_coordinate",
            "long_name": f"{c}-coordinate in projected coordinate system",
            "units": "meters",
        }

    # global attributes
    ds.attrs = {
        "title": title,
        "history": f"{str(datetime.datetime.utcnow())} Python",
        "Conventions": "CF-1.8",
        "institution": "Scenarios Network for Alaska + Arctic Planning",
        "contact": "kmredilla@alaska.edu",
        "source": source,
        "version": "1.0.0",
        "comment": "Intended for use with SERDP research on fish and fire in Alaskan boreal forests",
    }

    # spatial attributes
    crs = CRS.from_epsg(3338)
    ds["crs"] = np.int16()
    ds["crs"].attrs = crs.to_cf()

    # variable attributes
    ds[var].attrs = {
        "long_name": var_name,
        "units": "Kelvin",
        "grid_mapping": "crs",
    }
    ds.xc.attrs = cs_attrs("x")
    ds.yc.attrs = cs_attrs("y")
    ds.date.attrs = {
        "long_name": "modis_period_start_date",
        "units": f"days since {str(epoch)[:10]}"
    }
    return ds