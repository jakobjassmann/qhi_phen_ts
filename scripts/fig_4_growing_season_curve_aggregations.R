# Change in variation of NDVI curves with aggregation instead of resampling
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
library(pbmcapply)
library(gridExtra)
library(colorspace)

# Set global parameters / load site boundaries and meta data
figure_out_path <- "figures/fig_4_curve_fits_aggregations/"
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
# In this script I first re-sample the rasters to a common resolution of 5 cm
# and then aggregate cleanly using the mean. The approach taken in the figure
# for the main manuscript takes the opposite approach: it aggregates first and
# then re-samples.

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
                         make.names(paste0(unique(ts_combos$site_veg), 
                                           "_agg_raster_","0.05m"))),
                0.05), 
         envir = .GlobalEnv)
# and 33.333... m
list2env(mapply(function(x, y){create_agg_raster(x, y)},
                setNames(unique(ts_combos$site_veg),
                         make.names(paste0(unique(ts_combos$site_veg),
                                           "_agg_raster_","33.3m"))),
                33.3),
         envir = .GlobalEnv)

# Define function to calculate NDVI
NDVI <- function(red_band, nir_band) {
  (nir_band - red_band) / (nir_band + red_band) 
  }

# Define function to reasmple drone rasters and produce summary statistics
resample_agg_rasters <- function(site_veg_es, year_es, agg_level) {
  cat("starting: ", site_veg_es, "_", year_es, "_", 
      agg_level, "m... ", sep = "")
  # subset meta data
  ts_objects <- meta_data %>% 
    filter(site_veg == site_veg_es, 
           format(date, "%Y") == year_es, 
           band != "NDVI")
  # Remove PS4 HER 2017-07-17 (outlier)
  if  (site_veg_es == "PS4_HER" & 
       year_es == format(as.Date("2017", "%Y"),"%Y")) {
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

## 2) Fit phenology models ----
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
  coefs_df <- data.frame(
    site_veg = site_veg,
    veg = substr(site_veg, 5, 7),
    year = year,
    agg_level = agg_level,
    cell_id = seq(1, phen_curves@ncols * phen_curves@nrows),
    a = getValues(phen_curves[[3]]),
    b = getValues(phen_curves[[2]]),
    c = getValues(phen_curves[[1]]))
  
  cat("Finished processing: ", site_veg, "_", year, "_",
      agg_level, "m.\n", sep = "")
  
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

#load(paste0(data_out_path, "coefs_df.Rda"))

coefs_df <- mutate(coefs_df, agg_level_pretty = paste0(agg_level, " m"))
scatter_both <- ggplot(
  coefs_df, aes(x = a, y = b)) + 
  geom_point() + 
  xlab("coef. 'a'") +
  ylab("coef. 'b'") +
  facet_wrap(vars(agg_level_pretty)) +
  theme_cowplot(20) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

save_plot(paste0(figure_out_path, "../fig_s5_coefs_agg_level_aggregation.png"), 
          scatter_both, base_aspect_ratio = 1.6, base_height = 10)

cor(coefs_df$a, coefs_df$b, method = "spearman")
cor.test(coefs_df$a, coefs_df$b, method = "spearman")

## 3) Plot results -----
# Set colours
her_col <- "#1e9148FF"
kom_col <- "#1e5c91FF"

# Create plots showing the reduction in sd with increasing grain size
coefs_df_sd_summary <- coefs_df %>%
  group_by(site_veg, agg_level) %>% 
  summarise(mean_a = mean(a),
            mean_b = mean(b),
            mean_c = mean(c),
            sd_a = sd(a),
            sd_b = sd(b),
            sd_c = sd(c),
            min_a = min(a),
            min_b = min(b),
            min_c = min(c),
            max_a = max(a),
            max_b = max(b),
            max_c = max(c)) %>%
  mutate(site_name = substring(site_veg,1,3),
         veg = substring(site_veg,5,7),
         sd_a_e5 = sd_a *10^5,
         mean_a_e5 = mean_a *10^5,
         agg_level_pretty =  paste0(agg_level, " m"))

coes_df_sd_summary_mean <- coefs_df_sd_summary %>%
  group_by(agg_level) %>%
  summarise(mean_sd_a_e5 = mean(sd_a_e5))

coefs_df_sd_summary$agg_level_pretty <- ordered(
  coefs_df_sd_summary$agg_level_pretty, 
  levels = c("0.5 m",
             "1 m",
             "5 m",
             "10 m",
             "20 m",
             "33.3 m"))

a_sd_plot <- ggplot(coefs_df_sd_summary,
                    aes(x = agg_level,
                        y = sd_a_e5,
                        colour = veg)) +
  geom_smooth(data = coes_df_sd_summary_mean,
              mapping = aes(x = agg_level, y = mean_sd_a_e5),
              method = "lm",
              formula = y ~ log(x),
              inherit.aes = F,
              colour = "darkgrey",
              se = F,
              size = 1.5,
              alpha = 0.5) +
  geom_point(size = 2.5) +
  xlab("Aggregation level (m)") +
  ylab(expression(Variation~ "in"~a~(sigma~x~10^{-5}))) +
  scale_x_continuous(limits = c(0,35),
                     breaks = c(1,5,10,20,33.3),
                     minor_breaks = c(0.5)) +
  scale_y_continuous(limits = c(0,6),
                     breaks = seq(0,6,1)) +
  scale_colour_manual(values = c(her_col,
                                 kom_col)) + 
  annotate("text", x = 33.3, y = 5.75,
           label = "Tussock sedge tundra",
           colour = her_col,
           size = 6, hjust = 1) +
  annotate("text", x = 33.3, y = 5.25,
           label = "Dryas-vetch tundra",
           colour = kom_col,
           size = 6, hjust = 1) +
  theme_cowplot(20) +
  theme(legend.position = "none") 

save_plot(paste0(figure_out_path, "fig_4_panel_a_aggregation.png"),
          a_sd_plot, 
          base_aspect_ratio = 1.3)

a_mean_plot <- ggplot(coefs_df_sd_summary,
                      aes(x = agg_level,
                          y = mean_a_e5,
                          colour = veg)) +
  geom_point(size = 2.5) +
  xlab("Aggregation level (m)") +
  ylab(expression(Mean~ "of"~a~(mu~x~10^{-5}))) +
  scale_x_continuous(limits = c(0,35),
                     breaks = c(1,5,10,20,33.3),
                     minor_breaks = c(0.5)) +
  scale_y_continuous(limits = c(-25,0),
                     breaks = seq(-25,0,5)) +
  scale_colour_manual(values = c(her_col,
                                 kom_col)) + 
  annotate("text", x = 33.3, y = 0,
           label = "Tussock sedge tundra",
           colour = her_col,
           size = 6, hjust = 1) +
  annotate("text", x = 33.3, y = -2,
           label = "Dryas-vetch tundra",
           colour = kom_col,
           size = 6, hjust = 1) +
  theme_cowplot(20) +
  theme(legend.position = "none")

save_plot(paste0(figure_out_path, "../fig_s4_a_mean_agg_aggregation.png"), a_mean_plot, 
          base_aspect_ratio = 1.6)



## Plot PS2 Komakuk coef a aggreation as an example (panel c)
site_name <- "PS2"
veg_type <- "KOM"
rasters_to_plot <- paste0(site_name, "_",
                          veg_type, "_",
                          agg_levels, "_coefs.tif")

# Set scale breaks
scale_min <- coefs_df_sd_summary %>% 
  filter(site_veg == paste0(site_name, "_",
                            veg_type)) %>%
  pull(min_a) %>% 
  min() %>%
  round(4)
scale_max <- coefs_df_sd_summary %>% 
  filter(site_veg == paste0(site_name, "_",
                            veg_type)) %>%
  pull(max_a) %>% 
  max() %>%
  round(4)
scale_min <- scale_min * 10^5
scale_max <- scale_max * 10^5
scale_breaks <- seq(scale_min, scale_max, by = (scale_max - scale_min)/100)

# Set magma theme
magma_no_borders <- rasterTheme(
  region = magma(100, begin = 0, end = 1), # virdis scale with steps
  axis.line = list(lwd = 2),
  par.main.text = list(font = 1, cex = 2)) # normal face title

# Plot function
agg_levels_pretty <- c("0.5 m",
                       "1 m",
                       "5 m",
                       "10 m",
                       "20 m",
                       "33.3 m")
plot_list <- lapply(agg_levels, function(agg_level){
  plot_drone_native <- levelplot(
    (raster(paste0(data_out_path,
                   "curve_fits/2017/",
                   site_name, "_", veg_type, "_", agg_level, "m_coefs.tif"),
            band = 3) * 10^5),
    main = list(
      label = paste0(agg_levels_pretty[agg_levels == agg_level]), 
      cex = 2), 
    margin = F, 
    maxpixels = 6e5,
    colorkey = F,
    par.settings = magma_no_borders, 
    scales = list(draw = F),
    at = scale_breaks)
})

grid_matrix <- rbind(c( 1, 2, 3, 4, 5, 6))
png(paste0(figure_out_path, "fig_4_panel_c_aggregation.png"), 
    width = 12,
    height = 3, 
    units = "in", 
    res = 300)
grid.arrange(grobs = plot_list,
             ncol = 6, 
             nrow = 1,
             layout_matrix = grid_matrix)
dev.off()

# Scale bar plot
plot_drone_native <- levelplot(
  (raster(paste0(data_out_path,
                 "curve_fits/2017/",
                 site_name, "_", veg_type, "_", agg_levels[1], "m_coefs.tif"),
          band = 3) * 10^5),
  margin = F, # no margins
  maxpixels = 6e5,
  colorkey = list(draw = T, axis.line = list(lwd = 2), 
                  axis.text = list(cex = 1),
                  space='bottom'), 
  par.settings = magma_no_borders, 
  scales = list(draw = F),
  at = scale_breaks)

grid_matrix <- rbind(c( 1, 1, 1, 1, 1, 1))
png(paste0(figure_out_path, "fig_4_scale_bar_aggregation.png"), 
    width = 12,
    height = 3, 
    units = "in", 
    res = 300)
grid.arrange(plot_drone_native,
             ncol = 6, 
             nrow = 1,
             layout_matrix = grid_matrix)
dev.off()


## Panel b) curve samples 
# Example PS2 KOM

ndvi_stack_0.5m <- brick(paste0("data/fig_4_curve_resampled/2017/PS2_KOM_NDVI_stack_0.5.tif"))
ndvi_stack_10m <- brick(paste0("data/fig_4_curve_resampled/2017/PS2_KOM_NDVI_stack_10.tif"))
ndvi_stack_30m <- brick(paste0("data/fig_4_curve_resampled/2017/PS2_KOM_NDVI_stack_30.tif"))
doys <- meta_data %>%
  filter(site_veg == "PS2_KOM" & format(date, "%Y") == 2017 & band == "NDVI") %>%
  mutate(doy = format(date, "%j")) %>% 
  pull(doy)
names(ndvi_stack_0.5m) <- paste0("d", doys)
names(ndvi_stack_10m) <- paste0("d", doys)
names(ndvi_stack_30m) <- paste0("d", doys)
ndvi_values <- bind_rows(
  data.frame(getValues(ndvi_stack_0.5m), 
             agg_level = 0.5, 
             cell_id = 1:nrow(getValues(ndvi_stack_0.5m))),
  data.frame(getValues(ndvi_stack_10m), 
             agg_level = 10,
             cell_id = 1:nrow(getValues(ndvi_stack_10m))),
  data.frame(getValues(ndvi_stack_30m), 
             agg_level = 33.3,
             cell_id = 1:nrow(getValues(ndvi_stack_30m)))) 
ndvi_values <- ndvi_values %>%
  pivot_longer(cols = 1:5, 
               names_to = "d_doy",
               values_to = "NDVI") %>%
  mutate(doy = as.numeric(substr(d_doy,2,4))) %>% 
  dplyr::select(cell_id, agg_level, doy, NDVI)

# Calculate predictions:
coefs_df_sub <- coefs_df %>%
  filter(site_veg == "PS2_KOM", agg_level %in% c(0.5,10,33.3)) %>% 
  arrange(agg_level, cell_id)
# Create data frame for predictions
preds_df <- data.frame(
  doy = c(rep(seq(min(ndvi_values$doy), 
                  max(ndvi_values$doy)), 
              length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 0.5])),
          rep(seq(min(ndvi_values$doy), 
                  max(ndvi_values$doy)), 
              length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 10])),
          rep(seq(min(ndvi_values$doy), 
                  max(ndvi_values$doy)), 
              length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 33.3]))),
  cell_id = c(sort(rep(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 0.5], 
                       length(
                         seq(min(ndvi_values$doy),
                             max(ndvi_values$doy))))),
              sort(rep(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 10], 
                       length(
                         seq(min(ndvi_values$doy),
                             max(ndvi_values$doy))))),
              sort(rep(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 33.3], 
                       length(
                         seq(min(ndvi_values$doy),
                             max(ndvi_values$doy)))))),
  agg_level = c(rep(rep(0.5, 
                        length(
                          seq(min(ndvi_values$doy),
                              max(ndvi_values$doy)))), 
                    length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 0.5])),
                rep(rep(10, 
                        length(
                          seq(min(ndvi_values$doy),
                              max(ndvi_values$doy)))), 
                    length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 10])),
                rep(rep(33.3, 
                        length(
                          seq(min(ndvi_values$doy),
                              max(ndvi_values$doy)))), 
                    length(coefs_df_sub$cell_id[coefs_df_sub$agg_level == 33.3]))),
  cell_id_coefs_df = NA,
  a = NA,
  b = NA,
  c = NA,
  stringsAsFactors = F)
preds_df <- preds_df %>% arrange(agg_level, doy, cell_id)
preds_df[preds_df$agg_level == 0.5, 4:7] <-
  coefs_df_sub[coefs_df_sub$agg_level == 0.5, 5:8]
preds_df[preds_df$agg_level == 10, 4:7] <-
  coefs_df_sub[coefs_df_sub$agg_level == 10, 5:8]
preds_df[preds_df$agg_level == 33.3, 4:7] <-
  coefs_df_sub[coefs_df_sub$agg_level == 33.3, 5:8]
# Verify that matching worked
sum(preds_df$cell_id != preds_df$cell_id_coefs_df)


# # Optional subsample for quick plotting to adjust plot layout
# set.seed(5)
# sample_0.5m <- sample(unique(preds_df$cell_id[preds_df$agg_level == 0.5]),
#                       100)
# preds_df <- bind_rows(
#   preds_df %>% filter(agg_level == 0.5 & cell_id %in% sample_0.5m),
#   preds_df %>% filter(agg_level %in% c(10,33.3))
# )
# 
# # Reduce size of NDVI values df
# ndvi_values <- bind_rows(
#   ndvi_values %>% filter(agg_level == 0.5 & cell_id %in% sample_0.5m),
#   ndvi_values %>% filter(agg_level %in% c(10,33.3))
# )

# Calculate predictions
preds_df$preds <- preds_df$a*((preds_df$doy)^2) + 
  preds_df$b*(preds_df$doy) +
  preds_df$c

# Set grouping factors
ndvi_values$agg_level <- ordered(ndvi_values$agg_level,
                                 levels = c(0.5,10,33.3))
preds_df$agg_level <- ordered(preds_df$agg_level,
                              levels = c(0.5, 10, 33.3))

# Set colour ramp 
col_ramp <- sequential_hcl(5, palette = "Blues3")[1:3]
col_ramp <- sequential_hcl(3, palette = "Inferno")
# col_ramp[3] <- "#00dbfe"  
# Plot predictions
curve_plots <- ggplot(ndvi_values,
                      mapping = aes(x= doy, y= NDVI, 
                                    group = cell_id,
                                    colour = agg_level,
                                    fill = agg_level)) +
  geom_line(data = preds_df[preds_df$agg_level == 0.5,], 
            mapping = aes(x= doy, y= preds, 
                          group = cell_id),
            alpha = 0.01,
            size = 1,
            colour = col_ramp[1]) +
  geom_line(data = preds_df[preds_df$agg_level == 10,], 
            mapping = aes(x= doy, y= preds, 
                          group = cell_id),
            alpha = 0.5,
            size = 1,
            colour = col_ramp[2]) +
  geom_line(data = preds_df[preds_df$agg_level == 33.3,], 
            mapping = aes(x= doy, y= preds, 
                          group = cell_id),
            alpha = 0.5,
            size = 1,
            colour = col_ramp[3]) +
  geom_point(data = ndvi_values[ndvi_values$agg_level == 0.5,], 
             mapping = aes(x= doy, y= NDVI, 
                           group = cell_id),
             alpha = 0.05,
             size = 1,
             shape = 16,
             colour = col_ramp[1]) +
  geom_point(data = ndvi_values[ndvi_values$agg_level == 10,], 
             mapping = aes(x= doy, y= NDVI, 
                           group = cell_id),
             alpha = 0.5,
             size = 1,
             colour = col_ramp[2]) +
  geom_point(data = ndvi_values[ndvi_values$agg_level == 33.3,],
             mapping = aes(x= doy, y= NDVI,
                           group = cell_id),
             alpha = 0.7,
             size = 1,
             shape = 21,
             colour = "black",
             fill = col_ramp[3]) +
  scale_y_continuous(limits = c(0.2, 0.9), breaks = seq(0.2, 0.9, 0.1)) +
  xlab("Day of Year") +
  ylab("NDVI") +
  annotate("text", x = 200, y = 0.225,
           color = "black",
           label = parse(text = "'y = a x' ^ 2 * ' + b x + c'"),
           size = 5) +
  # annotate("rect", xmin = 175, xmax = 186, ymin = 0.73, ymax = 0.90,
  #          color = "black",
  #          fill = "black",
  #          alpha = 0.25) +
  annotate("point", x = 176.5, y = 0.865,
           size = 2,
           color = col_ramp[1]) +
  annotate("point", x = 176.5, y = 0.815,
           size = 2,
           color = col_ramp[2]) +
  annotate("point", x = 176.5, y = 0.765,
           size = 2, shape = 21, 
           fill = col_ramp[3]) +
  annotate("text", x = 176.5, y = 0.865,
           size = 5,
           color = col_ramp[1],
           hjust = 0,
           label = "  0.5 m") +
  annotate("text", x = 176.5, y = 0.815,
           size = 5,
           color = col_ramp[2],
           hjust = 0,
           label = "  10 m") +
  annotate("text", x = 176.5, y = 0.765,
           size = 5,
           color = "black",
           hjust = 0,
           label = "  33.3 m") +
  theme_cowplot(15) +
  theme(legend.position = "none")
# curve_plots
system.time(save_plot(paste0(figure_out_path, "fig_4_panel_b_aggregation.png"),
                      curve_plots,
                      base_aspect_ratio = 1.3))
