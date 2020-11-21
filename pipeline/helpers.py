"""Helper functions for MODIS LST-WRF TSK alignment"""

import datetime, os
import numpy as np


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
    return (dt_arr - epoch).astype("timedelta64[D]").astype("int")


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
    ds["crs"] = int()
    ds["crs"].attrs = {
        "grid_mapping_name": "albers_conical_equal_area",
        "crs_wkt": 'PROJCS["NAD83 / Alaska Albers",GEOGCS["NAD83",DATUM["North_American_Datum_1983",SPHEROID["GRS 1980",6378137,298.257222101,AUTHORITY["EPSG","7019"]],TOWGS84[0,0,0,0,0,0,0],AUTHORITY["EPSG","6269"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4269"]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["standard_parallel_1",55],PARAMETER["standard_parallel_2",65],PARAMETER["latitude_of_center",50],PARAMETER["longitude_of_center",-154],PARAMETER["false_easting",0],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["X",EAST],AXIS["Y",NORTH],AUTHORITY["EPSG","3338"]]',
    }

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