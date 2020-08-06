# MODIS LST and WRF TSK alignment

Align WRF-downscaled reanalysis (ERA-Interim) and GCM (GFDL-CM3, NCAR-CCSM4) data with MODIS LST data (MOD11A2) in time and space. Currently, only the WRF TSK variable is supported. 

The processing of MODIS LST and 1km WRF TSK data developed for the SERDP Fish and Fire project is done via the `alignment_pipeline` directory. This pipeline downloads and extracts the relevant 8-day aggregated MODIS LST data to the extent of the WRF domain in EPSG:3338, and downsamples the WRF to this grid. This downsampling occurs because although both products are "1km" resolution, the MODIS rasters in this area are less resolute than WRF when reprojected to EPSG:3338.

Legacy code for exploratory analysis and processing 20km WRF/MODIS alignment is currently retained. 

## Using this codebase

### 0.) software

This codebase is designed to be run with `pipenv`.

GDAL >= 2.0 is required. `pipenv` does not necessarily install the correct version of GDAL when `pipenv install` is run. If that fails, get the GDAL version via `gdalinfo --version` and run `pipenv install GDAL==<version>`.

Follow instructions [here](http://www.pymodis.org/info.html#requirements) for enabling an Earthdata account to download MODIS data.

### 1.) env vars

First, the paths to the necessary data need to be available in the following env vars.  

ONE (and *only* one) of the following two variables must be defined:

`WRF_1KM_DIR`: path to directory containing 1km WRF data  
`WRF_20KM_DIR`: path to directory containing 20km WRF data  

If both are defined, the WRF processing steps will throw an error. 

All of these must be set:  

`MODIS_DIR`: path to directory where M\*D11A2 data will be downloaded
`SCRATCH_DIR`: path to directory for storing intermediate files
`OUTPUT_DIR`: path to directory for writing transformed MODIS and WRF data

*Note: currently, the 1km WRF TSK data is not available for download, and so the appropriate files must be retrieved manually and arranged in a particular directory structure - those requirements can be found at the end of this README.*

### 2.) MODIS download

*Skip this step if MODIS data are present.* 

First, register with [NASA EARTHDATA.](https://urs.earthdata.nasa.gov/users/new)

`pipenv run alignment_pipeline/download_modis -u <username> -p <password>.py `

### 3.) Running the alignment pipeline

1. process the raw modis with flags for all steps: `pipenv run python alignment_pipeline/process_modis.py -cmwr`
2. extract 8-day time periods used in aggregated (averaged) MODIS LST data: `pipenv run python alignment_pipeline/extract_modis_periods_metadata.py`
3. temporally resample and aggregate the WRF TSK to the MODIS time periods: `pipenv run python alignment_pipeline/resample_wrf_to_modis_periods.py`.  
4. reproject the WRF TSK data to epsg:3338: `pipenv run python alignment_pipeline/reporoject_wrf.py`
5. clip the MODIS LST data to the reprojected WRF domain: `pipenv run python alignment_pipeline/clip_modis_to_wrf.py`
6. regrid the WRF to the slightly lower resolution of the clipped MODIS data: `pipenv run python alignment_pipeline/regrid_wrf_to_modis.py`

Aligned data for MODIS and WRF will be saved in: `$OUTPUT_DIR/MODIS` and `$OUTPUT_DIR/WRF`, repsectively. 

### 4.) Using the R scripts

These scripts rely on having the `$OUTPUT_DIR` environmental variable set for finding the data. This needs to be set to the path to the location of those output directories, `WRF/` and `MODIS/`.

#### `plot_era_comparison.R`

For specified coordinates (lat/lon), plot the averaged weekly timeseries of any **single** TSK source with MODIS LST values on the aligned grid for overlapping years, and save as a PNG. This script can be run via the command line with the `Rscript` via `Rscript R_scripts/plot_wrf_comparison.R <args>`. Help on these arguments can be displayed via `Rscript R_scripts/plot_era_comparison.R --help`. 

An example usage of this, making use of all arguments, is:  

```
 Rscript R_scripts/plot_era_comparison.R -t gfdl -a max -s MYD11A2 -x -147.72 -y 64.84 -o /path/to/file.png
```

This will create a plot using the WRF GFDL TSK values, resampled to the MODIS 8-day time periods using **maximum values** as the aggregate function, plotted against the LST values derived from the MYD11A2 sensor, all for the coordinates specified and saved to the output path.

#### `plot_gcm_comparison.R`

For specified coordinates (lat/lon), plot the averaged weekly timeseries of **both** GCM-based TSK data with MODIS LST values on the aligned grid for overlapping years, and save as a PNG. This script can be run via the command line with the `Rscript` via `Rscript R_scripts/plot_gcm_comparison.R <args>`. Help on these arguments can be displayed via `Rscript R_scripts/plot_gcom_comparison.R --help`. The main difference between this script and `plot_wrf_comparison.R` is that it removes the choice of TSK source, but plots both GCM TSK sources instead. As such, it has one less argument. Example:

```
 Rscript R_scripts/plot_weekly_timeseries.R -a max -s MOD11A2 -x -147.72 -y 64.84 -o /path/to/file.png
```

## WRF data directory structure

#### `WRF_1KM_DIR`

`WRF_1KM_DIR` should be the path to the directory containing the 1km WRF TSK data. This directory should be structured in the following way:

```
$WRF_1KM_DIR/
  ccsm/
    2007/
    2010/
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
    2010/
    ...
    2077/
```

Where each `WRFDS_*-*-*_serdp.nc` file contains data for all variables available, at all 24 integer-valued hours of the date in the filename. 

#### `WRF_20KM_DIR`

`WRF_20KM_DIR` should be the path to the directory containing the 20km WRF TSK data. This directory should be structured in the following way:
```
$WRF_20KM_DIR/
  tsk_wrf_hourly_ERA-Interim_2000.nc
  tsk_wrf_hourly_ERA-Interim_2001.nc
  ...
  tsk_wrf_hourly_ERA-Interim_2015.nc
  tsk_hourly_wrf_GFDL-CM3_historical_2000.nc 
  tsk_hourly_wrf_GFDL-CM3_historical_2001.nc
  ...
  tsk_hourly_wrf_GFDL-CM3_historical_2006.nc
  tsk_hourly_wrf_GFDL-CM3_rcp85_2006.nc
  tsk_hourly_wrf_GFDL-CM3_rcp85_2007.nc
  ...
  tsk_hourly_wrf_GFDL-CM3_rcp85_2100.nc
  tsk_hourly_wrf_NCAR-CCSM4_historical_2000.nc 
  tsk_hourly_wrf_NCAR-CCSM4_historical_2001.nc
  ...
  tsk_hourly_wrf_NCAR-CCSM4_historical_2006.nc
  tsk_hourly_wrf_NCAR-CCSM4_rcp85_2006.nc
  tsk_hourly_wrf_NCAR-CCSM4_rcp85_2007.nc
  ...
  tsk_hourly_wrf_NCAR-CCSM4_rcp85_2100.nc

```



