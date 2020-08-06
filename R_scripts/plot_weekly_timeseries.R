# Plot weekly MODIS LST and WRF TSK for any single WRF source (gfdl, ccsm, era)

suppressMessages({
  library(ncdf4)
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(sf)
})


convert_lonlat <- function(nc, ncvar, lonlat) {
  get_idx <- function(cv, x) {
    which(abs(cv - x) == min(abs(cv - x), na.rm = TRUE))
  }
  sfc <- st_sfc(st_point(lonlat), crs = 4326)
  new_sfc <- st_transform(sfc, 3338)
  x <- new_sfc[[1]][1]
  y <- new_sfc[[1]][2]
  xc <- nc$var[[ncvar]]$dim[[1]]$vals # xc
  yc <- nc$var[[ncvar]]$dim[[2]]$vals # yc
  c(get_idx(xc, x), get_idx(yc, y))
}


# function to extract data with date based on row,col
# default cell
extr_cell <- function(fp, ncvar, lonlat) {
  nc <- nc_open(fp)
  rc <- convert_lonlat(nc, ncvar, lonlat)
  varr <- ncvar_get(nc, ncvar, c(rc[1], rc[2], 1), c(1,1,-1))
  # THIS IS WHAT IS SHOULD BE, fix after removing 0 values
  # varr[varr == -9999] <- NA
  varr[varr < 1] <- NA
  # convert to C
  varr <- varr - 273.15
  dates <- ncvar_get(nc, "date")
  date_units <- nc$var[[ncvar]]$dim[[3]]$units
  origin <- strsplit(date_units, " ")[[1]][3]
  dates <- as.Date(dates, origin = origin)
  # use filename to get variable
  nc_close(nc)
  data.frame(
    date = dates, value = varr, source = toupper(ncvar),
    stringsAsFactors = FALSE
  )
}


# function to summarise extracted data as avg, min, max by week of year
aggr_woy <- function(tsk_df, lst_df) {
  # ensure use of only years present for both TSK and LST
  tsk_df <- mutate(tsk_df, year = format(date, "%Y"))
  lst_df <- mutate(lst_df, year = format(date, "%Y"))
  keep_yrs <- intersect(tsk_df$year, lst_df$year)
  df <- rbind(tsk_df, lst_df) %>%
    filter(year %in% keep_yrs) %>%
    mutate(
      week = as.factor(format(date, "%V")),
      source = factor(source, levels = c("LST", "TSK"))
    ) %>%
    group_by(week, source) %>%
    summarise(
      mean = mean(value, na.rm = TRUE),
      min = min(value, na.rm = TRUE), 
      max = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    select(week, source, mean, min, max) %>%
    list(c(min(keep_yrs), max(keep_yrs)))
}


# plot lst and tsk aggregate timeseries
plot_aggr_woy <- function(tsk_fp, 
                          lst_fp, 
                          lonlat, 
                          tsk_source, 
                          aggr_type, 
                          sensor) {
  # prep data for plot
  tsk_df <- extr_cell(tsk_fp, "tsk", lonlat)
  lst_df <- extr_cell(lst_fp, "lst", lonlat)
  aggr <- aggr_woy(tsk_df, lst_df)
  df <- aggr[[1]]
  yrs <- aggr[[2]]
  # make plot
  p <- ggplot(df, aes(x = week, y = mean, group = source)) + 
    geom_line(aes(color = source), size = 1) + 
    geom_ribbon(aes(ymin=min, ymax=max, fill = source), alpha=0.2) + 
    scale_color_manual(values = c("#00AFBB", "#E7B800")) +
    scale_fill_manual(values = c("#00AFBB", "#E7B800")) + 
    scale_x_discrete(breaks = c("01", "10", "20", "30", "40", "50")) + 
    ylab("Temperature (°C)") + xlab("Week of Year") + 
    ggtitle(
      paste0(
        "MODIS LST vs WRF TSK weekly summary for ", lat, "N, ", lon, "E"
      ),
      subtitle = paste0(
        sensor, "/", tsk_source, ", ", aggr_type, ", ", yrs[1], "-", yrs[2]
      )
    ) +
    theme_classic() + 
    theme(
      panel.background = element_rect(colour = "black", size=1),
      legend.title = element_blank())
  return(p)
}


#-- Main ----------------------------------------------------------------------
# setup
out_dir = Sys.getenv("OUTPUT_DIR")
wrf_dir <- file.path(out_dir, "WRF")
modis_dir <- file.path(out_dir, "MODIS")

tsk_dir <- file.path(wrf_dir, "tsk_1km_3338_multiband")
lst_dir <- file.path(modis_dir, "lst_1km_3338_multiband")

tsk_source <- NULL
aggr_type <- NULL
sensor <- NULL
out_fp <- NULL

# setup command line arguments
option_list = list(
  make_option(
    c("-t", "--tsk-source"), type = "character", default = "era",
    help = "Source model for TSK (accepts one of: era, gfdl, or ccsm)", 
    metavar = "character"
  ),
  make_option(
    c("-a", "--aggr-type"), type = "character", default = "max",
    help = "Aggregate method (accepts one of: min, mean, or max)", 
    metavar = "character"
  ),
  make_option(
    c("-s", "--sensor"), type = "character", default = "MOD11A2",
    help = "Sensor type (accepts one of: mod, myd)", 
    metavar = "character"
  ),
  make_option(
    c("-x", "--lon"), type = "double", default = -147.72, 
    help = "longitude (default: -147.72", metavar = "character"
  ),
  make_option(
    c("-y", "--lat"), type = "double", default = 64.84, 
    help = "latitude (default: 64.84 N)", metavar = "character"
  ),
  make_option(
    c("-o", "--out-file"), type = "character", default = "", 
    help ="output filepath (default: $OUTPUT_DIR/plots/woy_modis_lst_wrf_tsk_<tsk-source>_<aggr-type>_<lat>N<lon>W.png)", 
    metavar = "character"
  )
)

# parse args and check
opt = parse_args(OptionParser(option_list=option_list))
if (!(opt$`tsk-source` %in% c("era", "gfdl", "ccsm"))) {
  stop("Invalid argument for --tsk-source")
} else {tsk_source <- opt$`tsk-source`}
if (!(opt$`aggr-type` %in% c("min", "mean", "max"))) {
  stop("Invalid argument for --aggr-type")
} else {aggr_type <- opt$`aggr-type`}
if (!(opt$sensor %in% c("MOD11A2", "MYD11A2"))) {
  stop("Invalid argument for --sensor")
} else {sensor <- opt$sensor}

lon <- opt$lon
lat <- opt$lat
# stop if neither out file or $OUTPUT_DIR are set, otherwise set it
if (opt$`out-file` == "") {
  if (out_dir == "") stop("$OUTPUT_DIR must be set if --out-file not used")
  out_fp <- file.path(
    out_dir, "plots", paste0(
      "woy_modis_lst_wrf_tsk_", 
      tsk_source, "_", aggr_type, "_", sensor, "_", 
      lat, "N", abs(lon), "W.png"
    )
  )
} else {out_fp <- opt$`out-file`}

# source WRF file
tsk_fp <- file.path(
  tsk_dir, paste0(
    "tsk_", tsk_source, "_", 
    aggr_type, "_MOD11A2_multiband.nc"
  )
)
lst_fp <- file.path(
  lst_dir, paste0("lst_", sensor, "_InteriorAK_3338_multiband.nc")
)

# create plot
p <- plot_aggr_woy(tsk_fp, lst_fp, c(lon, lat), tsk_source, aggr_type, sensor)

# write
dir.create(dirname(out_fp), showWarnings = FALSE)
ggsave(out_fp, p, "png", height=5)
cat(out_fp, "\n")

#------------------------------------------------------------------------------

