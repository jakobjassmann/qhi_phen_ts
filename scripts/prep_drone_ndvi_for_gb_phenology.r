# This script determines and extracts drone NDVI data
# for the ground based phenology monitoring plots.
# Needs to be run in perparation for the figure 5 script.
# Jakob Assmann Feburary 2019

## load dependencies
library(raster)
library(dplyr)
library(rgdal)
library(rasterVis)
library(gridExtra)
library(grid)
library(ggplot2)
library(viridis)

## Prepare global variables 
data_path <- "data/ground_based_phenology/"

## Load required datasets
# Ground based phenology plot coordinates
plot_coords <- read.csv("data/ground_based_phenology/gb_plot_coordinates.csv")
# Meta data
load("data/meta_data.Rda")

# grab utm zone 7 prj4 data
epsg <- make_EPSG()
utm_z7N <- as.character(epsg %>% filter(code == 32607) %>% select(prj4))
rm(epsg)

# Define a function to create a polygon from for coordinate pairs
plot_bound_poly <- function(site_veg){
  plot_bound_coords <- as.matrix(
    plot_coords[plot_coords$site_veg == site_veg,4:5])
  plot_bound_poly <- spPolygons(plot_bound_coords, crs=utm_z7N)
  return(plot_bound_poly)
}

# Apply function to all plots
plot_bound_poly_list <- lapply(unique(plot_coords$site_veg), plot_bound_poly)
names(plot_bound_poly_list) <- paste0(unique(plot_coords$site_veg), 
                                      "_gbplot_poly")
list2env(plot_bound_poly_list, envir = .GlobalEnv)
rm(plot_bound_poly_list)

# Plot for quality control
lapply(paste0(unique(plot_coords$site_veg), "_gbplot_poly"),
       function(x){
         plot(get(x), main = x)
       })
# Looking goood!

# Specifiy combos
ts_combos <- data.frame(
  site_veg = c("PS1_HER",
               "PS1_KOM",
               "PS2_HER",
               "PS2_KOM",
               "PS1_HER",
               "PS1_KOM",
               "PS2_HER",
               "PS2_KOM",
               "PS3_HER",
               "PS3_KOM",
               "PS4_HER",
               "PS4_KOM"),
  year = format(as.Date(c("2016",
                          "2016",
                          "2016",
                          "2016",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017",
                          "2017"),
                        format = "%Y"), "%Y"),
  stringsAsFactors = F)

### 1) Extract gb plot NDVI -----

## Function to extract plot NDVI average
# site_veg_es <- "PS1_KOM"
# year_es <- 2017
extract_gb_plot_NDVI <- function(site_veg_es, year_es) {  
  cat("starting: ", site_veg_es, "_", year_es, "\n", sep = "")
  # subset meta data
  ts_objects <- meta_data %>% 
    filter(site_veg == site_veg_es, 
           format(date, "%Y") == year_es, 
           band == "NDVI")
  # filter multiples for the 3 August for PS1 in 2016 
  # and remove PS4 HER 2017-07-17 (outlier)
  if(site_veg_es == "PS1_HER" & 
     year_es == format(as.Date("2016", "%Y"),"%Y")) {
    ts_objects <- ts_objects[c(1:6),]
  } else if (site_veg_es == "PS1_KOM" & 
             year_es == format(as.Date("2016", "%Y"),"%Y")) {
    ts_objects <- ts_objects[c(1:5),]
  } else if  (site_veg_es == "PS4_HER" & 
              year_es == format(as.Date("2017", "%Y"),"%Y")) {
    ts_objects <- ts_objects %>% filter(date != as.Date("2017-07-17"))
  } else if (site_veg_es == "PS1_KOM") {
    return(NULL)
  }
  # Load files
  list2env(
    lapply(
      setNames(ts_objects$file_path, 
               make.names(ts_objects$object_name)),
      raster), 
    envir = .GlobalEnv)
  
  # Crop and Mask Raster
  list2env(
    lapply(
      setNames(ts_objects$object_name, 
               make.names(paste0(ts_objects$object_name, "_gbplot"))),
      function(x) { 
        cropped_raster<- crop(get(x), get(paste0(site_veg_es, "_gbplot_poly")))
        masked_raster <- mask(cropped_raster, 
                              get(paste0(site_veg_es, "_gbplot_poly")))
        return(masked_raster)
      }), 
    envir = .GlobalEnv)
  ts_objects$object_gbplot <- paste0(ts_objects$object_name, "_gbplot")
  
  # Calculate summary stats
  ts_gbplots_stats <- data.frame(
    flight_id = ts_objects$flight_id,
    date = ts_objects$date, 
    doy = format(ts_objects$date, "%j"),
    site_veg = site_veg_es,
    site_name = ts_objects$site_name,
    veg_type = ts_objects$veg_type,
    ndvi_mean = sapply(mget(ts_objects$object_gbplot, 
                            envir = .GlobalEnv), 
                       function(x) { cellStats(x, "mean")}), 
    ndvi_sd = sapply(mget(ts_objects$object_gbplot, 
                          envir = .GlobalEnv), 
                     function(x) { cellStats(x, "sd")}), 
    ndvi_min = sapply(mget(ts_objects$object_gbplot, 
                           envir = .GlobalEnv), 
                      function(x) { cellStats(x, "min")}), 
    ndvi_max = sapply(mget(ts_objects$object_gbplot, 
                           envir = .GlobalEnv), 
                      function(x) { cellStats(x, "max")})) 
  ts_gbplots_stats$ndvi_mean_plus_sd = ts_gbplots_stats$ndvi_mean + 
    ts_gbplots_stats$ndvi_sd
  ts_gbplots_stats$ndvi_mean_minus_sd = ts_gbplots_stats$ndvi_mean - 
    ts_gbplots_stats$ndvi_sd
  
  ggplot(ts_gbplots_stats,aes(x = doy, y = ndvi_mean)) + geom_point()
  
    # Calculate Time-Series grand min and max
  grand_min <- min(ts_gbplots_stats$ndvi_min)
  grand_max <- max(ts_gbplots_stats$ndvi_max)
  # Round down / up to one decimals
  grand_min <- floor(grand_min * 10) / 10 
  grand_max <- ceiling(grand_max * 10) / 10
  
  # Plot
  gbplots <- lapply(ts_objects$object_gbplot, function(x){
    gbplot <- levelplot(
      get(x), 
      margin = FALSE, 
      colorkey= list(
        space='bottom',
        labels=list(at=seq(grand_min, grand_max, 0.1), font=2)
      ),
      par.settings = list(
        axis.line=list(col='transparent'),
        par.main.text = list(font = 2, # make it bold
                             just = "left", 
                             x = grid::unit(5, "mm"))
      ),
      main = paste0("DOY: ", 
                    ts_gbplots_stats$doy[ts_objects$object_gbplot == x]),
      scales=list(draw=FALSE),
      col.regions=viridis(100),                   
      at=seq(grand_min, grand_max, (grand_max - grand_min) / 100)
      
    )
    return(gbplot)
  })

  # Sort list by ascending doy
  gbplots <- lapply(order(ts_gbplots_stats$doy), 
         function(x){
           gbplot <- gbplots[[x]]
           return(gbplot)
         })
  
  # export to png file
  png(paste0(data_path, site_veg_es, "_", year_es, "_gbplot_ts.png"), 
      width = 4 * length(gbplots),
      height = 4, 
      units = "in", 
      res = 200)
  print(marrangeGrob(grobs = gbplots,
                      ncol = length(gbplots), 
                      nrow = 1,
                      top = textGrob(paste0(site_veg_es, "_", year_es), 
                                     gp=gpar(fontsize=14))))
  dev.off()
  
  # return stats dataframe
  return(ts_gbplots_stats)
}

# Execute for all sites and year combinations
gbstats_list <- mapply(
  extract_gb_plot_NDVI,
  ts_combos$site_veg,
  ts_combos$year,
  SIMPLIFY = F)
# collapse list into dataframe 
gbstats_df <- bind_rows(gbstats_list)
# tidy up
rm(gbstats_list)

gbstats_df$year <- format(gbstats_df$date, "%Y")

# Visual quality control, blurry images for:
# PS1_HER_2017 DOY 198
gbstats_df <- gbstats_df[c(-which(gbstats_df$site_veg == "PS1_HER" & 
                                  gbstats_df$year == 2017 &
                                  gbstats_df$doy == 198)),]
# PS1_KOM_2017 all values were discarded as flight path did not extent
# to gb phenology plot. see beginning of script
# Ps2_KOM_2016 DOY 212 & 218, data have v. small "holes" looks good otherwise,
# so retained in analysis
# PS2_KOM_2017 DOY 192
gbstats_df <- gbstats_df[c(-which(gbstats_df$site_veg == "PS2_KOM" & 
                                  gbstats_df$year == 2017 &
                                  gbstats_df$doy == 192)),]
# PS3_HER_2017 DOY 199
gbstats_df <- gbstats_df[c(-which(gbstats_df$site_veg == "PS3_HER" & 
                                  gbstats_df$year == 2017 &
                                  gbstats_df$doy == 199)),]
# PS3_KOM_2017 DOY 199
gbstats_df <- gbstats_df[c(-which(gbstats_df$site_veg == "PS3_KOM" & 
                                  gbstats_df$year == 2017 &
                                  gbstats_df$doy == 199)),]
# PS4_KOM_2017 DOY 191
gbstats_df <- gbstats_df[c(-which(gbstats_df$site_veg == "PS4_KOM" & 
                                  gbstats_df$year == 2017 &
                                  gbstats_df$doy == 191)),]

# Save stats df to file
save(gbstats_df, file = paste0(data_path, "gbstats_drone.Rda"))

     