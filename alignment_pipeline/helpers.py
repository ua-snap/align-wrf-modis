"""
Helper functions for MODIS LST-WRF TSK alignment
"""

import os, argparse
import numpy as np


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


def parse_args():
    """
    Validate args for specifying which WRF data to work on, for select scripts
    Args to check are year_range and model_groups.
    """
    parser = argparse.ArgumentParser(
        description="resample the WRF data to the MODIS time scale"
    )
    parser.add_argument(
        "-r",
        "--year_range",
        action="store",
        dest="year_range",
        help="WRF years to work on ('2000-2018', '2037-2047', or '2067-2077')",
    )
    parser.add_argument(
        "-g",
        "--model_groups",
        action="store",
        dest="model_groups",
        help="WRF model groups to work on: ('era', 'ccsm', or 'gfdl'. \
            String multiple together by separation with '-')",
    )

    args = parser.parse_args()
    year_range = args.year_range
    model_groups = args.model_groups

    if model_groups != None:
        model_groups = model_groups.split("-")
    else:
        model_groups = [""]

    valid_ranges = ["2000-2018", "2037-2047", "2067-2077"]
    if year_range not in valid_ranges:
        exit("Invalid year range specified")
    else:
        years = year_range.split("-")

    valid_groups = ["gfdl", "ccsm", "era"]
    group_check = [group in valid_groups for group in model_groups]
    all_groups = None
    if model_groups[0] == "":
        all_groups = True
        model_groups = ["gfdl", "ccsm"]
    elif any(group_check) == False:
        exit("Invalid model group(s) specified")
    elif ("era" in model_groups) & (year_range != valid_ranges[0]):
        exit("'era' model group only allowed with '2000-2018' range")

    if year_range == valid_ranges[0]:
        if all_groups:
            model_groups = ["era"] + model_groups
    
    return year_range, model_groups

