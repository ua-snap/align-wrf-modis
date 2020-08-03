"""
Perform the following with MODIS MOD11A2/MYD11A2 data:
  - converts raw MODIS from HDF4 to GeoTIFF for relevant bands (LST Daytime)
  - extracts and mosaics the two tiles overlapping the WRF domain
  - reprojects the mosaicked files to epsg:3338
  - rescales the projected data to correct values in Kelvin
"""


def extract_bands(fn):
    print(fn)
    # translate to GTiff
    return subprocess.call(
        [
            "gdal_translate",
            "-q",
            "-of",
            "GTiff",
            "-sds",
            "-co",
            "COMPRESS=LZW",
            fn,
            os.path.basename(fn).replace(".", "_").replace("_hdf", ".tif"),
        ]
    )


def run_all_extract_bands(files, output_dir, ncpus=16):
    # change the directory
    os.chdir(output_dir)
    # parallel it
    pool = mp.Pool(ncpus)
    out = pool.map(extract_bands, files)
    pool.close()
    pool.join()
    return output_dir


def move_file(x):
    fn, out_fn = x
    return shutil.move(fn, out_fn)


def move_desired_bands(in_dir, out_dir, keep_bands=["01", "03"], ncpus=14):

    # list and filter the files based on the desired bands
    files = glob.glob(os.path.join(in_dir, "*.tif"))
    files = [
        fn
        for fn in files
        if fn.endswith(tuple("_{}.tif".format(band) for band in keep_bands))
    ]
    out_files = [os.path.join(out_dir, os.path.basename(fn)) for fn in files]
    # move the files to a new dir
    pool = mp.Pool(ncpus)
    moved = pool.map(move_file, list(zip(files, out_files)))
    pool.close()
    pool.join()
    return out_dir


def move_undesired_bands(in_dir, out_dir, ncpus=16):
    files = glob.glob(os.path.join(temp_dir, "*.tif"))
    out_files = [os.path.join(out_dir, os.path.basename(fn)) for fn in files]
    pool = mp.Pool(ncpus)
    out = pool.map(move_file, list(zip(files, out_files)))
    pool.close()
    pool.join()
    return out_dir


def make_band_lookup():
    """hardwired name lookup table for M*D11A2 data..."""
    return {
        "01": "LST_Day_1km",
        "02": "QC_Day",
        "03": "Day_view_time",
        "04": "Day_view_angl",
        "05": "LST_Night_1km",
        "06": "QC_Night",
        "07": "Night_view_time",
        "08": "Night_view_angl",
        "09": "Emis_31",
        "10": "Emis_32",
        "11": "Clear_sky_days",
        "12": "Clear_sky_nights",
    }


def make_mosaic_args(files, out_dir):
    """make arguments to pass to mosaic"""

    # make a dataframe that we can group and use in arguments creation for the mosaicking
    colnames = ["product", "date", "tile", "version", "production", "band"]
    df = pd.DataFrame(
        [os.path.basename(fn).split(".")[0].split("_") for fn in files],
        columns=colnames,
    )
    df["fn"] = files

    def make_args(group, out_dir):
        files = sorted(group["fn"].tolist())
        elems = (
            group[["product", "date", "version", "production", "band"]].iloc[0].tolist()
        )
        out_fn = os.path.join(mosaicked_dir, "_".join(elems) + ".tif").replace(
            "_006_", "_InteriorAK_006_"
        )
        return [files] + [out_fn]

    return [
        make_args(group, out_dir)
        for i, group in df.groupby(["product", "date", "band"])
    ]


def mosaic_tiles(files, out_fn):
    command = ["gdal_merge.py", "-n", "0", "-a_nodata", "0", "-o", out_fn,] + files
    _ = subprocess.call(command)
    return out_fn


def wrap_mosaic_tiles(x):
    """a wrapper for the mosaic f(x) for parallelism"""
    files, out_fn = x
    return mosaic_tiles(files, out_fn)


def run_mosaic_tiles(args, ncpus=5):
    # pool = mp.Pool(ncpus)
    # out = pool.map(wrap_mosaic_tiles, args)
    # pool.close()
    # pool.join()
    out = [wrap_mosaic_tiles(arg) for arg in args]


def warp_to_3338(fn, out_fn):
    return subprocess.call(
        [
            "gdalwarp",
            "-q",
            "-overwrite",
            "-t_srs",
            "EPSG:3338",
            "-co",
            "COMPRESS=LZW",
            fn,
            out_fn,
        ]
    )


def wrap_warp_to_3338(x):
    return warp_to_3338(*x)


def run_warp_to_3338(args, ncpus=5):
    pool = mp.Pool(ncpus)
    out = pool.map(wrap_warp_to_3338, args)
    pool.close()
    pool.join()
    return out


def rescale_values(fn):
    with rasterio.open(fn) as rst:
        meta = rst.meta.copy()
        meta.update(compress="lzw", dtype="float32", nodata=0)
        arr = rst.read(1)

    arr_out = np.copy(arr).astype(np.float32)
    ind = np.where(arr != meta["nodata"])

    # scale it:
    if fn.endswith(("_01.tif", "_05.tif")):
        arr_out[ind] = arr[ind] * 0.02

    elif fn.endswith(("_09.tif", "_10.tif")):
        arr_out[ind] = (arr[ind] * 0.0020) + 0.49

    elif fn.endswith("_03.tif"):
        arr_out[ind] = arr[ind] * 0.1

    else:
        raise BaseException("wrong bands")

    # make the output filename and dump to disk
    out_fn = fn.replace("warped", "rescaled")
    with rasterio.open(out_fn, "w", **meta) as out:
        out.write(arr_out.astype(np.float32), 1)
    return out_fn


def run_rescale_values(files, ncpus):
    pool = mp.Pool(ncpus)
    out = pool.map(rescale_values, files)
    pool.close()
    pool.join()
    return out


if __name__ == "__main__":
    import argparse
    import subprocess, glob, os, shutil
    import pandas as pd
    import numpy as np
    import rasterio
    import multiprocessing as mp
    import time

    # parse args
    parser = argparse.ArgumentParser(
        description="run various processing tasks on MODIS data"
    )
    parser.add_argument(
        "-c", "--convert", action="store_true", dest="convert", help="convert to GTiff"
    )
    parser.add_argument(
        "-m", "--mosaic", action="store_true", dest="mosaic", help="mosaic tiles"
    )
    parser.add_argument(
        "-w", "--warp", action="store_true", dest="warp", help="warp to EPSG:3338"
    )
    parser.add_argument(
        "-r", "--rescale", action="store_true", dest="rescale", help="rescale values"
    )
    parser.add_argument(
        "-n",
        "--ncpus",
        action="store",
        dest="ncpus",
        type=int,
        help="number of cpus to use",
    )

    # unpack the args here
    args = parser.parse_args()
    convert = args.convert
    mosaic = args.mosaic
    warp = args.warp
    rescale = args.rescale
    ncpus = args.ncpus

    # LIST THE DATA
    base_dir = os.getenv("MODIS_DIR")
    scratch_dir = os.getenv("SCRATCH_DIR")
    files = [
        os.path.join(r, fn)
        for r, s, files in os.walk(base_dir)
        for fn in files
        if fn.endswith(".hdf")
    ]

    # SETUP OUTPUT DIRECTORIES
    converted_dir = os.path.join(scratch_dir, "MODIS", "converted")
    if not os.path.exists(converted_dir):
        _ = os.makedirs(converted_dir)

    # make a tempdir to dump MOD11A2 bands
    temp_dir = os.path.join(scratch_dir, "MODIS", "TEMP")
    if not os.path.exists(temp_dir):
        _ = os.makedirs(temp_dir)

    # make a directory to store the unwanted extra bands from the HDF files. (might be needed later)
    extra_bands_dir = os.path.join(scratch_dir, "MODIS", "raw_extra_bands_extracted")
    if not os.path.exists(extra_bands_dir):
        _ = os.makedirs(extra_bands_dir)

    # make a directory for the mosaicked outputs
    mosaicked_dir = os.path.join(scratch_dir, "MODIS", "mosaicked")
    if not os.path.exists(mosaicked_dir):
        _ = os.makedirs(mosaicked_dir)

    # make a directory to store the warped to 3338 outputs
    warped_dir = os.path.join(scratch_dir, "MODIS", "warped")
    if not os.path.exists(warped_dir):
        _ = os.makedirs(warped_dir)

    # make a directory to store the final rescaled outputs
    rescaled_dir = os.path.join(scratch_dir, "MODIS", "rescaled")
    if not os.path.exists(rescaled_dir):
        _ = os.makedirs(rescaled_dir)

    # -----
    # PROCESS BAND EXTRACTION AND CONVERSION TO GTiff
    if convert:
        print("Converting files to GTiff", end="... ")
        tic = time.perf_counter()
        out = run_all_extract_bands(files, temp_dir, ncpus=ncpus)
        print("done,", round(time.perf_counter() - tic, 1), "s")

        # move the desired bands to a new location
        print("Moving desired bands", end="... ")
        tic = time.perf_counter()
        _ = move_desired_bands(
            temp_dir, converted_dir, keep_bands=["01", "03"], ncpus=ncpus
        )
        print("done,", round(time.perf_counter() - tic, 1), "s")

        # move the bands we want to the extra_bands_dir
        print("Moving undesired bands", end="... ")
        tic = time.perf_counter()
        _ = move_undesired_bands(temp_dir, extra_bands_dir, ncpus=ncpus)
        print("done,", round(time.perf_counter() - tic, 1), "s")

    # ---
    # MOSAIC TILES: using gdal_merge.py
    if mosaic:
        print("Mosaicking tiles", end="... ")
        tic = time.perf_counter()
        files = glob.glob(os.path.join(converted_dir, "*.tif"))
        args = make_mosaic_args(files, mosaicked_dir)
        mosaicked = run_mosaic_tiles(args, ncpus=ncpus)
        print("done,", round(time.perf_counter() - tic, 1), "s")

    # ---
    # WARP TO EPSG:3338 AK Albers
    if warp:
        print("warping to 3338", end="... ")
        tic = time.perf_counter()
        files = glob.glob(os.path.join(mosaicked_dir, "*.tif"))
        args = [(fn, os.path.join(warped_dir, os.path.basename(fn))) for fn in files]
        warped = run_warp_to_3338(args, ncpus=ncpus)
        print("done,", round(time.perf_counter() - tic, 1), "s")

    # ---
    # RESCALE TO THE PROPER VALUES -- according to user-guide
    #   https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod11_user_guide.pdf
    if rescale:
        print("Rescaling values according to user guide", end="... ")
        tic = time.perf_counter()
        files = glob.glob(os.path.join(warped_dir, "*.tif"))
        rescaled = run_rescale_values(files, ncpus=ncpus)
        print("done,", round(time.perf_counter() - tic, 1), "s")

