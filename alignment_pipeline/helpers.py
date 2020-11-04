"""
Helper functions for MODIS LST-WRF TSK alignment
"""

import os


def check_env():
    """
    Check that all necessary env vars are present. 
    Exit with error message if not.
    """
    key_vars = ["WRF_DIR", "MODIS_DIR", "SCRATCH_DIR", "OUTPUT_DIR"]
    env_vars = os.environ

    absent = [var not in env_vars for var in key_vars]
    if any(absent):
        absent_idx = np.where(absent)
        miss_vars = [key_vars[idx] for idx in absent_idx[0]]
        exit(f"The following env vars are missing: {miss_vars}")

    return True
