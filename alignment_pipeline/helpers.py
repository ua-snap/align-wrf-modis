"""
Helper functions for MODIS LST-WRF TSK alignment
"""

import os


def check_env(wrf_env_var=True):
    """
    Will need to check that environment variables are set up properly,
      return either False or name of defined WRF var
    """
    key_vars = ["MODIS_DIR", "SCRATCH_DIR", "OUTPUT_DIR"]
    env_vars = os.environ

    if any(var not in env_vars for var in key_vars):
        return False

    # can't have both env vars defined
    if ("WRF_1KM_DIR" in env_vars) & ("WRF_20KM_DIR" in env_vars):
        return False

    # Need at least one to be defined
    if wrf_env_var:
        if ("WRF_1KM_DIR" not in env_vars) & ("WRF_20KM_DIR" not in env_vars):
            return False
        elif "WRF_1KM_DIR" in env_vars:
            wrf_env_var = "WRF_1KM_DIR"
        else:
            wrf_env_var = "WRF_20KM_DIR"

    # PLACEHOLDER: check that files in directory match expected WRF structure
    # wrf_dir = env_vars[wrf_env_var]

    return wrf_env_var


def get_model_groups(wrf_env_var):
    """get model groups based on chosen WRF type"""
    # model groups
    model_groups = {
        "WRF_1KM_DIR": ["era", "gfdl", "ccsm"],
        "WRF_20KM_DIR": [
            "ERA-Interim_historical",
            "GFDL-CM3_historical",
            "GFDL-CM3_rcp85",
            "NCAR-CCSM4_historical",
            "NCAR-CCSM4_rcp85",
        ],
    }
    return model_groups[wrf_env_var]

