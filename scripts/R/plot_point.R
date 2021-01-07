# Plot weekly MODIS LST and/or WRF TSK for any source

suppressMessages({
  library(ncdf4)
  library(dplyr)
  library(ggplot2)
  library(optparse)
  library(sf)
})


make_args <- function(fn, nc_dir, ncvar, source) {
  args <- list(
    fp = file.path(nc_dir, fn),
    ncvar = ncvar,
    source = source
  )
}


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
extr_cell <- function(args, lonlat) {
  nc <- nc_open(args$fp)
  ncvar <- args$ncvar
  rc <- convert_lonlat(nc, ncvar, lonlat)
  # helper function to query data from NetCDF
  varr <- ncvar_get(nc, ncvar, c(rc[1], rc[2], 1), c(1,1,-1))
  if (all(is.nan(varr))) stop("ERROR: No data present for one or more sources")
  # convert to C
  varr <- varr - 273.15
  dates <- ncvar_get(nc, "date")
  # accessing the NetCDF without helper functions
  date_units <- nc$var[[ncvar]]$dim[[3]]$units
  origin <- strsplit(date_units, " ")[[1]][3]
  dates <- as.Date(dates, origin = origin)
  # use filename to get variable
  nc_close(nc)
  data.frame(
    date = dates, 
    value = varr, 
    source = args$source,
    stringsAsFactors = FALSE
  ) %>%
    mutate(year = format(date, "%Y"))
}


# function to summarise extracted data as avg, min, max by week of year
# aggr_woy <- function(tsk1_df, tsk2_df, lst_df) {
aggr_woy <- function(df_list) {
  # ensure use of only years present for both all data
  yr_list <- lapply(df_list, function(df) df$year)
  sources <- unlist(lapply(df_list, function(df) unique(df$source)))                
  keep_yrs <- Reduce(intersect, yr_list)                         
  df <- do.call("rbind", df_list) %>%
    filter(year %in% keep_yrs) %>%
    mutate(
      week = as.factor(format(date, "%V")),
      source = factor(source, levels = sources)
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
plot_aggr_woy <- function(extr_args, lonlat) {
  aggr <- lapply(extr_args, extr_cell, lonlat) %>%
    aggr_woy()
  df <- aggr[[1]]
  yrs <- aggr[[2]]
  # make plot
  cmap <- c("#FF0000", "#00A08A", "#F2AD00", "#5BBCD6", "#F98400")
  p <- ggplot(df, aes(x = week, y = mean, group = source)) + 
    geom_line(aes(color = source), size = 1) + 
    geom_ribbon(aes(ymin=min, ymax=max, fill = source), alpha=0.2) + 
    scale_color_manual(values = cmap) +
    scale_fill_manual(values = cmap) + 
    scale_x_discrete(breaks = c("01", "10", "20", "30", "40", "50")) + 
    ylab("Temperature (°C)") + xlab("Week of Year") + 
    ggtitle(
      paste0(
        "Aggregated temperature at ", lat, "N, ", lon, "E"
      ),
      subtitle = paste0(
        "Aligned WRF(max) / MODIS(mean) data, ", yrs[1], "-", yrs[2]
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
out_dir = file.path(Sys.getenv("OUTPUT_DIR"))
# OUTPUT_DIR must be set
if (out_dir == "") stop("$OUTPUT_DIR must be set")
tsk_dir <- file.path(out_dir, "aligned-WRF-MODIS", "WRF")
lst_dir <- file.path(out_dir, "aligned-WRF-MODIS", "MODIS")

# setup command line arguments
option_list = list(
  make_option(
    c("-e", "--era"), action = "store_true", default = FALSE,
    help = "Include ERA-Interim in plot"
  ),
  make_option(
    c("-c", "--ccsm"), action = "store_true", default = FALSE,
    help = "Include CCSM4 in plot"
  ),
  make_option(
    c("-g", "--gfdl"), action = "store_true", default = FALSE,
    help = "Include GFDL CM3 in plot"
  ),
  make_option(
    c("-t", "--mod"), action = "store_true", default = FALSE,
    help = "Include MOD11A2 in plot"
  ),
  make_option(
    c("-a", "--myd"), action = "store_true", default = FALSE,
    help = "Include MYD11A2 in plot"
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
    help ="output filepath (default: $OUTPUT_DIR/plots/aligned_data_comparison_<lat>N<lon>W.png)", 
    metavar = "character"
  )
)

# parse args and check
opt = parse_args(OptionParser(option_list=option_list))   
lon <- opt$lon
lat <- opt$lat

# source WRF files with args for later 
nc_args = list()
if (opt$gfdl) {
  nc_args[["gfdl"]] <- make_args(
      "tsk_max_gfdl_2008-2017_aligned.nc",
      tsk_dir, "tsk", "GFDL-TSK"
  )
}
if (opt$ccsm) {
  nc_args[["ccsm"]] <- make_args(
      "tsk_max_ccsm_2008-2017_aligned.nc",
      tsk_dir, "tsk", "CCSM4-TSK"
  )
}
if (opt$era) {
  nc_args[["era"]] <- make_args(
      "tsk_max_era_2000-2018_aligned.nc",
      tsk_dir, "tsk", "ERA-TSK"
  )
}
if (opt$mod) {
  nc_args[["mod"]] <- make_args(
      "lst_MOD11A2_aligned.nc",
      lst_dir, "lst", "MOD-LST"
  )
}
if (opt$myd) {
  nc_args[["myd"]] <- make_args(
      "lst_MYD11A2_aligned.nc",
      lst_dir, "lst", "MYD-LST"
  )
}
                           
if (opt$`out-file` == "") {
  fp_mod_names <- paste(names(nc_args), collapse = "_")
  out_fp <- file.path(
    out_dir, "plots", 
    paste0(
      "aligned_data_comparison_",
      fp_mod_names, "_", lat, "N", abs(lon), "W.png"
    )
  )
} else {out_fp <- opt$`out-file`}
    
# create plot
p <- plot_aggr_woy(nc_args, c(lon, lat))

# write
dir.create(dirname(out_fp), showWarnings = FALSE)
ggsave(out_fp, p, "png", height=5)
cat(out_fp, "\n")

#------------------------------------------------------------------------------

