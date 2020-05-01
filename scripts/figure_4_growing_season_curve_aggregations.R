# Change in variation of NDVI curves with aggregation
# Jakob Assmann j.assmann@ed.ac.uk 8 October 2018

### Preparations ----
# Dependencies
library(dplyr)
library(tidyr)
library(ggplot2)
library(raster)
library(rasterVis)
library(viridisLite)
library(MCMCglmm)
library(cowplot)

# Set global parameters / load site boundaries and meta data
figure_out_path <- "figures/"
log_path <- "log/"
data_out_path <- "data/fig_4_curve_aggregations/"
site_boundaries <- read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
load("data/meta_data.Rda")

# Function to create an extent object form the boundaries data.frame
get_sent_extent <- function(site_veg_id, sent_boundaries) {
  extent_object <- extent(
    c(sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmax,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymax))
  return(extent_object)
}

# Create extent objects
PS1_HER_extent <- get_sent_extent("PS1_HER", site_boundaries)
PS1_KOM_extent <- get_sent_extent("PS1_KOM", site_boundaries)
PS2_HER_extent <- get_sent_extent("PS2_HER", site_boundaries)
PS2_KOM_extent <- get_sent_extent("PS2_KOM", site_boundaries)
PS3_HER_extent <- get_sent_extent("PS3_HER", site_boundaries)
PS3_KOM_extent <- get_sent_extent("PS3_KOM", site_boundaries)
PS4_HER_extent <- get_sent_extent("PS4_HER", site_boundaries)
PS4_KOM_extent <- get_sent_extent("PS4_KOM", site_boundaries)

# Set time-series combos of interest (in 2016 we can only use PS1 and PS2)

ts_combos <- data.frame(
  site_veg = c("PS1_HER",
               "PS1_KOM",
               "PS2_HER",
               "PS2_KOM",
               "PS3_HER",
               "PS3_KOM",
               "PS4_HER",
               "PS4_KOM"),
  year = format(as.Date(c("2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017"),
                        format = "%Y"), "%Y"),
  stringsAsFactors = F)

# Set aggregation levels for drone data
agg_levels <- c(0.5, 1,5,10,20,33.3)
# fill in ts_combos data frame!
ts_combos <- data.frame(
  ts_combos,
  agg_levels = sort(rep(agg_levels, nrow(ts_combos))),
  stringsAsFactors = F
)

### 1) Perform aggregation ----

# Note:
# as the drone rasters all have odd resultions that are just below or around
# 0.05 m we have to resample rather than just aggregate
# by default 'resample' first performs an aggregation if the scale factor is > 2
# and then shifts the position of the cells accurately using the
# 'bilinear' algorithm.
# See https://github.com/cran/raster/blob/master/R/resample.R

# define function to prepare aggregation templates
create_agg_raster <- function(site_veg, agg_level){
  # Build a raster
  # Empty object
  agg_raster <- raster()
  # Set projection
  projection(agg_raster) <- "+proj=utm +zone=7 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  # Set extent
  extent(agg_raster) <- get(paste0(site_veg, "_extent"))
  # Set dimension
  ncol(agg_raster) <- 100 / agg_level
  nrow(agg_raster) <- 100 / agg_level
  # return raster
  return(agg_raster)
}

# Create aggregation templates for all site_veg combos for 0.05 m
list2env(mapply(function(x, y){create_agg_raster(x, y)},
                setNames(unique(ts_combos$site_veg), 
                         make.names(paste0(unique(ts_combos$site_veg), "_agg_raster_","0.05m"))),
                0.05), 
         envir = .GlobalEnv)
# and 33.333... m
list2env(mapply(function(x, y){create_agg_raster(x, y)},
                setNames(unique(ts_combos$site_veg), 
                         make.names(paste0(unique(ts_combos$site_veg), "_agg_raster_","33.3m"))),
                30), 
         envir = .GlobalEnv)

# Define function to calculate NDVI
NDVI <- function(red_band, nir_band) {(nir_band - red_band) / (nir_band + red_band) }

# Define function to reasmple drone rasters and produce summary statistics
resample_agg_rasters <- function(site_veg_es, year_es, agg_level) {
  cat("starting: ", site_veg_es, "_", year_es, "_", agg_level, "m... ", sep = "")
  # subset meta data
  ts_objects <- meta_data %>% 
    filter(site_veg == site_veg_es, 
           format(date, "%Y") == year_es, 
           band != "NDVI")
  # Remove PS4 HER 2017-07-17 (outlier)
  if  (site_veg_es == "PS4_HER" & year_es == format(as.Date("2017", "%Y"),"%Y")) {
    ts_objects <- ts_objects %>% filter(date != as.Date("2017-07-17"))
  }
  # Load raster files
  list2env(
    lapply(
      setNames(ts_objects$file_path, 
               make.names(ts_objects$object_name)),
      raster), 
    envir = .GlobalEnv)
  
  # Crop rasters to site extent
  list2env(
    lapply(
      setNames(ts_objects$object_name, 
               make.names(paste0(ts_objects$object_name, "_cropped"))),
      function(x) { 
        crop(get(x), get(paste0(site_veg_es, "_extent")))
      }), 
    envir = .GlobalEnv)
  ts_objects$object_cropped <- paste0(ts_objects$object_name, "_cropped")
  
  # Re-sample then aggregate rasters
  
  # If not already done resample to 0.05 m
  agg_raster <- paste0(site_veg_es, "_agg_raster_0.05m")
  if(!file.exists(paste0(data_out_path, year_es, "/",
                  ts_objects$object_cropped[1], "_resamp_0.05m.tif"))){
    list2env(
      lapply(
        setNames(ts_objects$object_cropped, 
                 make.names(paste0(ts_objects$object_name, "_resamp_0.05m"))),
        function(x) {resample(get(x), get(agg_raster))}), 
      envir = .GlobalEnv)
    ts_objects$object_resamp_0.05m <- paste0(ts_objects$object_name, 
                                             "_resamp_0.05m")
    lapply(ts_objects$object_resamp_0.05m, 
           function(x){
             writeRaster(get(x), filename = paste0(data_out_path, year_es, "/",
                                   x, ".tif"),
                         overwrite = T)
           })
  }
  
  # Next aggreagate
  list2env(
    lapply(
      setNames(ts_objects$object_resamp_0.05m, 
               make.names(paste0(ts_objects$object_name, "_agg"))),
      
      function(x) { 
        agg_raster <- aggregate(get(x), floor(agg_level / 0.05), expand = F)
        # if agg_level = 33.3 the division is not even as the pixel size 
        # is 33.25, we re sample to 33.3 exactly using nearest neighbour
        # This "nudges" the pixels into place
        if(agg_level == 33.3) {
          agg_raster <- resample(agg_raster,
                                 get(paste0(site_veg_es, "_agg_raster_33.3m")),
                                 method = "ngb")
        }
        return(agg_raster)
      }), 
    envir = .GlobalEnv)
  ts_objects$object_agg <- paste0(ts_objects$object_name, "_agg")
  
  # Calculate NDVI
  red_bands <- (ts_objects %>% filter(band == "RED"))$object_agg
  nir_bands <- (ts_objects %>% filter(band == "NIR"))$object_agg
  ndvi_rasters <- paste0(unique(ts_objects$flight_id), "_", agg_level, "m_ndvi")
  list2env(
    mapply(
      function(x, y) { NDVI(get(x), get(y))}, 
      setNames(red_bands, 
               make.names(ndvi_rasters)),
      nir_bands),
    envir = .GlobalEnv)
  
  # Calculate summary statistics
  ts_stats <- data.frame(
    flight_id = ts_objects[ts_objects$band == "RED",]$flight_id,
    date = ts_objects[ts_objects$band == "RED",]$date, 
    site_veg = site_veg_es,
    site_name = ts_objects[ts_objects$band == "RED",]$site_name,
    veg_type = ts_objects[ts_objects$band == "RED",]$veg_type,
    agg_level = agg_level,
    ndvi_mean = sapply(mget(ndvi_rasters, 
                            envir = .GlobalEnv), 
                       function(x) { cellStats(x, "mean")}), 
    ndvi_sd = sapply(mget(ndvi_rasters, 
                          envir = .GlobalEnv), 
                     function(x) { cellStats(x, "sd")}), 
    ndvi_min = sapply(mget(ndvi_rasters, 
                           envir = .GlobalEnv), 
                      function(x) { cellStats(x, "min")}), 
    ndvi_max = sapply(mget(ndvi_rasters, 
                           envir = .GlobalEnv), 
                      function(x) { cellStats(x, "max")})) 
  ts_stats$ndvi_mean_plus_sd = ts_stats$ndvi_mean + ts_stats$ndvi_sd
  ts_stats$ndvi_mean_minus_sd = ts_stats$ndvi_mean - ts_stats$ndvi_sd
  
  # Plot the data, as extent and resolution match we can stack this time
  # and use the rasterVis violin plots:
  ts_ndvi_stack <- stack(mget(ndvi_rasters, envir = .GlobalEnv))
  names(ts_ndvi_stack) <- format(ts_stats$date, "%b_%d")
  png(filename = paste0(data_out_path, 
                        year_es, "/", site_veg_es, "_", 
                        year_es, "_", agg_level, "m.png"),
      height = 6, width = 9, unit = "in", res = 300)
  print(bwplot(ts_ndvi_stack,
               main = list(label = paste0(site_veg_es, " ", 
                                          year_es, " aggregated to ", 
                                          agg_level, " m"), cex = 2), 
               scales = list(axis.x.text = list(rot = 0)),
               ylab = list(label = "NDVI", cex = 1.5)))
  dev.off()
  
  # Save stack for later use
  writeRaster(ts_ndvi_stack, filename = paste0(
    data_out_path, year_es, "/", 
    site_veg_es, "_NDVI_stack_", agg_level, ".tif"),
    overwrite=TRUE)
  # clean up 
  rm(list = c(ts_objects$object_name,
              ts_objects$object_cropped,
              ts_objects$object_agg,
              ndvi_rasters), envir = .GlobalEnv)
  gc()
  
  # Status output
  cat("Finished processing: ", 
      site_veg_es, "_", year_es, "_", agg_level, "m.\n", sep = "")
  return(ts_stats)
}

# Extract NDVI for all combos and create day difference column in dataframe
stats_list <- mapply(resample_agg_rasters,
                     ts_combos$site_veg, 
                     ts_combos$year, 
                     ts_combos$agg_levels, 
                     SIMPLIFY = F)

# collapse list into dataframe 
stats_df <- bind_rows(stats_list)
# tidy up
rm(stats_list)

# create year and doy collumns
stats_df$doy <- as.integer(format(stats_df$date, "%j"))
stats_df$year <- factor(format(stats_df$date, "%Y"))

# Save stats data frame for later
save(stats_df , file = paste0(data_out_path, "stats_df.Rda"))
#load(paste0(data_out_path, "stats_df.Rda"))

# Define function to fit a simple binominal model
phen_model <- function(site_veg, year, agg_level){
  cat("starting: ", site_veg, "_", year, "_", agg_level, "m... ", sep = "")
  # Collect Meta Data
  doy <- stats_df[stats_df$site_veg == site_veg & 
                    stats_df$year == year & 
                    stats_df$agg_level == agg_level,]$doy
  # load timeseries brick
  phen_brick <- brick(paste0(data_out_path, year, "/" , site_veg, 
                             "_NDVI_stack_", agg_level, ".tif"))
  # assign doy as layer names
  names(phen_brick) <- doy
  # fir quadratic curve
  fit_quadratic <- function(x) {
    lm(x ~ poly(doy,2, raw = T))$coefficients
  }
  phen_curves <- calc(phen_brick, fun = fit_quadratic)
  
  # quality control with random subsample of 9 points 
  cell_id <- sample(phen_brick@ncols * phen_brick@nrows, 9)
  cell_values <- raster::extract(phen_brick, cell_id)
  cell_values <- data.frame(t(cell_values))
  colnames(cell_values) <- cell_id
  cell_values$doy <- doy
  cell_values <- gather(cell_values, "cell_id", "NDVI", -doy)
  
  # calculate model predictions:
  model_coeffs <- data.frame(raster::extract(phen_curves, cell_id))
  model_coeffs$cell_id <- as.character(cell_id)
  model_coeffs  
  preds_df <- data.frame(
    doy = rep(seq(min(doy), max(doy)), length(cell_id)),
    cell_id = sort(rep(as.character(cell_id), length(seq(min(doy), max(doy))))))
  preds_df$pred <- mapply(function(x, cell_id){
    y <- model_coeffs[model_coeffs$cell_id == cell_id, 1] + 
      model_coeffs[model_coeffs$cell_id == cell_id, 2]*x +
      model_coeffs[model_coeffs$cell_id == cell_id, 3]*x^2
  }, preds_df$doy, preds_df$cell_id)
  
  # plot curves and values for quality control
  curve_plot <- ggplot(
    cell_values, 
    aes(x = doy, y = NDVI, group = cell_id, colour = cell_id)) +
    geom_point() +
    geom_smooth(data = preds_df, 
                mapping = aes(x= doy, y= pred, 
                              group = cell_id))
  save_plot(curve_plot, 
            filename = paste0(data_out_path, "QC/curve_fits/",
                              "/", site_veg, "_", agg_level, ".png"))
  
  # NB Notation from now on: a quadratic term, b linear term, c intercept.
  # visualise and plot rasters
  a_coef <- levelplot(phen_curves[[3]])
  b_coef <- levelplot(phen_curves[[2]])
  c_coef <- levelplot(phen_curves[[1]])
  png(filename = paste0(data_out_path, "QC/curve_fits/", 
                        "/", site_veg, "_", agg_level, "_a_coef.png"))
  print(a_coef)
  dev.off()
  png(filename = paste0(data_out_path, "QC/curve_fits/",
                        "/", site_veg, "_", agg_level, "_b_coef.png"))
  print(b_coef)
  dev.off()
  png(filename = paste0(data_out_path, "QC/curve_fits/",
                        "/", site_veg, "_", agg_level, "_c_coef.png"))
  print(c_coef)
  dev.off()
  
  # write raster
  writeRaster(
    phen_curves, 
    filename =  paste0(data_out_path, 
                       "/curve_fits/", year, "/", site_veg, "_", 
                       agg_level, "m_coefs.tif"), 
    overwrite = T)
  
  # create data_frame to return values
  coefs_df <- data.frame(site_veg = site_veg,
                         veg = substr(site_veg, 5, 7),
                         year = year,
                         agg_level = agg_level,
                         cell_id = seq(1, phen_curves@ncols * phen_curves@nrows),
                         a = getValues(phen_curves[[3]]),
                         b = getValues(phen_curves[[2]]),
                         c = getValues(phen_curves[[1]]))
  
  cat("Finished processing: ", site_veg, "_", year, "_", agg_level, "m.\n", sep = "")
  
  # return df
  return(coefs_df)
}

# execute for all sites and agg_levels
coefs_list <- mapply(
  phen_model, 
  ts_combos$site_veg, 
  ts_combos$year, 
  ts_combos$agg_levels, 
  SIMPLIFY = F)
# collapse list into dataframe 
coefs_df <- bind_rows(coefs_list)
# tidy up
rm(coefs_list)
save(coefs_df , file = paste0(data_out_path, "coefs_df.Rda"))

# quick scatter plots to check out patterns
ordination_both <- ggplot(coefs_df, aes(x = a, y = b, group = veg, colour = veg)) + geom_point() + 
  scale_colour_manual(values = c("#440154FF", "#21908CFF")) + facet_wrap(vars(agg_level))
ordination_HER <- ggplot(coefs_df, aes(x = a, y = b, group = veg, colour = veg)) + geom_point() + 
  scale_colour_manual(values = c("#44015400", "#21908CFF")) + facet_wrap(vars(agg_level))
ordination_KOM <- ggplot(coefs_df, aes(x = a, y = b, group = veg, colour = veg)) + geom_point() + 
  scale_colour_manual(values = c("#440154FF", "#21908C00")) + facet_wrap(vars(agg_level))

ordination_plots <- plot_grid(ordination_both, ordination_HER, ordination_KOM)
save_plot(paste0(data_out_path, "ordination_plots_agg_levels.png"), ordination_plots, base_aspect_ratio = 1.5, base_height = 10)

