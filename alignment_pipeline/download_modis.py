"""
Download MODIS LST data to $MODIS_DIR, creating subdirs as necessary
"""

import argparse
import datetime
import itertools
import os
import sys
from pymodis import downmodis
from pathlib import Path
from helpers import check_env


def download_modis(user, pw, prod, wdir):
    """
    wrapper for downmodis functions.
    Non-changing vars (e.g., user/pw) left in global scope.
    """

    def list_present_dates(wdir=wdir, tile="1"):
        """
        Compile list of dates already present in MODIS_DIR for particular source/tile.
        Done to avoid redundant downloads, as this script may have to be applied
          more than once to download all files. 
        Takes to minimum of what is available for both tiles
        """
        dates_present = []
        modis_files = list(wdir.glob("*h1{}v02*.hdf".format(tile)))
        # keep = [os.stat(file).st_size > 0 for file in modis_files]
        # files_present = list(itertools.compress(modis_files, keep))
        for file in modis_files:
            """
            files of size 0 are assumed to be unsuccessful downloads, 
              need to be removed for .dayDownload() to re-download
            """
            if os.stat(file).st_size == 0:
                print("Removing empty file: {}".format(str(file).split("/")[-1]))
                os.remove(file)
                continue
            jday = int(str(file).split(".")[1][5:])
            year = int(str(file).split(".")[1][1:5])
            date = datetime.datetime(year, 1, 1) + datetime.timedelta(jday - 1)
            dates_present.append(date.strftime("%Y.%m.%d"))
        return dates_present

    # choose aqua or terra
    product = "M{}D11A2.006".format(prod)
    prod_info = {
        "MOD11A2.006": ["MOLT", "2000-02-24"],
        "MYD11A2.006": ["MOLA", "2002-07-24"],
    }
    path = prod_info[product][0]
    start_date = prod_info[product][1]
    # setup args for downmodis
    url = "https://e4ftl01.cr.usgs.gov"
    tiles = "h11v02,h12v02"
    modis_download = downmodis.downModis(
        destinationFolder=wdir,
        password=pw,
        user=user,
        url=url,
        tiles=tiles,
        product=product,
        path=path,
    )
    modis_download.connect()
    """
    NOTE:
    Using .downloadsAllDay() was failing unpredictably after variable numbers of
    successful downloads. This iterative approach seems promising after testing.
    """

    h11_dates_present = list_present_dates()
    h12_dates_present = list_present_dates(tile="2")
    if len(h11_dates_present) < len(h12_dates_present):
        dates_present = h11_dates_present
    else:
        dates_present = h12_dates_present
    # Download difference of what files have/have not been downloaded or transferred
    avail_dates = modis_download.getAllDays()
    get_dates = [date for date in avail_dates if date not in dates_present]
    if len(get_dates) == 0:
        print("No new dates to get.")
        return None
    n_diff = len(avail_dates) - len(get_dates)
    print("Downloading {} data, skipping {} files".format(product[:7], n_diff))
    for day in get_dates:
        files_list = modis_download.getFilesList(day)
        """
        AttributeError occurs after some unpredictable number of downloads.
          This might just be a system-specific thing, but this exception will 
          retry the download 5 times.
        """
        for attempt in range(5):
            try:
                modis_download.dayDownload(day, files_list)
            except AttributeError:
                print(
                    "Download failed for {} {}. Failure count: {}.".format(
                        product[:7], day, attempt
                    )
                )
                if attempt < 4:
                    print("Trying again...")
            else:
                break
        else:
            sys.exit(
                "Download failed five times for {} {}, download aborted.".format(
                    product[:7], day, attempt
                )
            )
    print("done.")


if __name__ == "__main__":
    env_ok = check_env(wrf_env_var=False)
    if not env_ok:
        exit("Environment variables incorrectly setup, check README for requirements")

    parser = argparse.ArgumentParser(description="Download MODIS LST data")
    parser.add_argument(
        "-u",
        "--user",
        action="store",
        dest="user",
        type=str,
        help="Earthdata login username",
    )
    parser.add_argument(
        "-p",
        "--password",
        action="store",
        dest="pw",
        type=str,
        help="Earthdata login password",
    )
    # unpack args
    args = parser.parse_args()
    user = args.user
    pw = args.pw

    # build target directories
    modis_dir = Path(os.getenv("MODIS_DIR"))
    myd_dir = modis_dir.joinpath("MYD11A2")
    mod_dir = modis_dir.joinpath("MOD11A2")
    # don't overwrite dirs if exist
    if not myd_dir.is_dir():
        myd_dir.mkdir()
    if not mod_dir.is_dir():
        mod_dir.mkdir()

    download_modis(user, pw, "Y", myd_dir)
    download_modis(user, pw, "O", mod_dir)

