# MODIS LST and WRF TSK alignment

Align WRF-downscaled reanalysis (ERA-Interim) and GCM (GFDL-CM3, NCAR-CCSM4) maximum skin temperature data (TSK) with 8-day MODIS LST data (MOD11A2, MYD11A2) in time and space.

## Background

The problem this codebase solves is that of getting 8-day MODIS LST data and hourly, 1km WRF-downscaled TSK data in the same spatial and temporal scale. To do so, both datasets are projected to the same CRS, and the WRF data is temporally resampled to the matching time periods used by MODIS. 

There are multiple periods' worth of WRF data. There is a historical period, 2008-2017, and two future periods: 2038-2047, and 2068-2077. The historical period is contained within the observed MODIS time period, but the future periods are obviously not. To "align" the future WRF temporally, we estimate the future MODIS periods based on the historical pattern: January 1 of every year marks the start of the first eight-day period for that year, and a new period begins every eight days. 

The processing of MODIS LST and 1km WRF TSK data developed for the SERDP Fish and Fire project is done via the scripts in the `./pipeline` directory. Steps for running the pipeline are given in the next section. This pipeline processes 8-day aggregated MODIS LST data to the extent of the WRF domain in EPSG:3338. The WRF data is then aggregated via maximum to the same 8-day periods as MODIS, and then downsampled to same grid. Finally, future WRF data are adjusted via the delta method. 

## Using this codebase

### 0.) software and Earthdata

This codebase is designed to be run with `pipenv`. To set up the environment, run `pipenv install --dev --sequential` (args needed for proper dependency configuration) from within the project directory.

GDAL >= 2.0 is required. `pipenv` does not necessarily install the correct version of the python GDAL bindings when `pipenv install` is run. If that fails, get the GDAL version installed on the system being used via `gdalinfo --version` and run `pipenv install GDAL==<version>` install the correct version of the GDAL bindings and update the `Pipfile`.

A NASA Earthdata account is required to download the MODIS data. Follow the instructions [here](http://www.pymodis.org/info.html#requirements) for enabling an Earthdata account to do so.

### 1.) environment

The following env vars must be set to run the pipeline:  

`WRF_DIR`: path to directory containing 1km WRF data
`MODIS_DIR`: path to directory where M\*D11A2 data will be downloaded
`SCRATCH_DIR`: path to directory for storing intermediate files
`OUTPUT_DIR`: path to directory for writing aligned MODIS and WRF data

You may need to run `unset GDAL_DATA` if you encounter an error (known to interfere with `rasterio.crs.CRS`)

**Note: currently, the 1km WRF data is not available for download, and so the appropriate files must be retrieved manually and arranged in a particular directory structure - those requirements can be found at the end of this README.**

### 2.) MODIS download

*Skip this step if MODIS data are present.* 

Steps to download the MODIS LST data:
1. register with [NASA EARTHDATA.](https://urs.earthdata.nasa.gov/users/new)
2. run `pipenv run pipeline/download_modis -u <username> -p <password>.py `

### 3.) Running the alignment pipeline

The pipeline is run via `pipenv run python pipeline/<script.py>`, where `<script.py>` is the pipeline script to be executed. Script execution is not automated, and this is not planned. The order for executing the scripts is given next, along with an explanation of the arguments used where applicable.  

In many of the scripts, you may specify the number of cores to use for program execution via `-n` or `--ncpus`. Using more cores can dramatically decrease the time needed to run the pipeline from start to finish.  

**Order for executing scripts in the pipeline:**
1. `process_modis.py -cmwr -n 16`. Process the raw modis with flags for all steps (-c, -m, -w, and -r correspond to convert, mosaic, warp, and reproject). Writes to `$SCRATCH_DIR`
2. `modis_periods.py`. Extract 8-day time periods used in the M\*D11A2 LST data, and create future periods. Writes to `$SCRATCH_DIR`.
3. `resample_wrf.py`. Temporally resample and aggregate the WRF TSK to the MODIS time periods via maximum. Writes to `SCRATCH_DIR`.
4. `wrf_template.py`. Create a WRF template GeoTIFF in epsg:3338 for clipping MODIS. Writes to `$SCRATCH_DIR`.
5. `clip_modis.py`. Clip the MODIS LST data to the template WRF extent in EPSG:3338 and create final output data files. Writes to `$OUTPUT_DIR`.
6. `reproject_wrf.py`. Reproject the WRF TSK data to the MODIS grid in EPSG:3338 (downscales to the slightly lower-res MODIS data). Writes to `$SCRATCH_DIR`.
7. `adjust_wrf.py`. Adjust the future WRF TSK data values via the delta method and create final output data. Writes to `$OUTPUT_DIR`.

Aligned data for MODIS and WRF will be saved in: `$OUTPUT_DIR/MODIS` and `$OUTPUT_DIR/WRF`, repsectively. 

### 4.) Other scripts

The `./scripts` directory contains a script for plotting data (in the `R/` subdirectory) along with some jupyter notebooks with supporting info in the `ancillary` subdirectory. Use of the R script is explained below. Here is a brief description of the scripts in the `ancillary` dir:

`explore_deltas.ipynb`: Assessments of the delta correction method applied for some sample data.
`explore_modis_dates.ipynb`: Demonstrates that the MODIS Terra 8-day dates is a superset of the Aqua dates.

### 5.) Using the R script:

This script relies on having the `$OUTPUT_DIR` environmental variable set for finding the data. This needs to be set to the path containing those output directories, `WRF/` and `MODIS/`. E.g., run `export OUTPUT_DIR=/path/to/folder` before running using this script.

#### `plot_point.R`

For specified coordinates (lat/lon), plot the **historical** weekly MODIS LST and/or WRF TSK for **any** source or combination of sources. This script can be run via the command line with the `Rscript` command, e.g. from the project folder, run `Rscript scripts/R/plot_point.R <args>`. Help on these arguments can be displayed via `Rscript scripts/R/plot_point.R --help`. 

An example usage of this, making use of all arguments, is:  

```
 Rscript scripts/R/plot_point.R  --era --mod -x -147.72 -y 64.84 -o /path/to/output/file.png
```

This will create a plot of the aligned WRF ERA-Interim TSK and MODIS LST values from the MOD11A2 sensor plotted for the coordinates specified and saved to the output path.

## WRF data directory structure

#### `WRF_DIR`

`WRF_DIR` should be the path to the directory containing the 1km WRF TSK data. This directory should be structured in the following way:

```
$WRF_1KM_DIR/
  ccsm/
    2007/
    2008/
    ...
    2077/

  era/
    1979/ # same structure for these dirs
      WRFDS_1979-07-02_serdp.nc
      WRFDS_1979-07-03_serdp.nc
      ...
      WRFDS_1979-12-31_serdp.nc
    1980/
    ...
    2014/

  gfld/
    2007/
    2008/
    ...
    2077/
```

Where each `WRFDS_*-*-*_serdp.nc` file contains data for all variables available, at all 24 integer-valued hours of the date in the filename. 

