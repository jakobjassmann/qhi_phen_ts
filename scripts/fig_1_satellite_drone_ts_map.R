# Phenology Time-Series Quick Script to create a map of all sites in the
# time-series.
# Jakob Assmann j.assmann@ed.ac.uk 11 October 2018

# Dependencies
library(dplyr)
library(tidyverse)
library(raster)
library(rgdal)
library(rasterVis)
library(viridisLite)

library(maps)
library(mapdata)
library(sf)

library(cowplot)
library(ggmap)
library(magick)
library(gridExtra)

library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)

library(oce)

### Global prepartions ----

# Set global variables
site_boundaries <- 
  read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
ts_out_path <- "figures/fig_1_ts_plots/"
sentinel_path <- 
  "/Volumes/BowheadRdge/phen_time_series/sentinel_data/qhi_cld_free/"
drone_path <- "/Volumes/BowheadRdge/phen_time_series/final_outputs"
data_out_path <- "data/fig_1_satellite_drone_ts_map/"

# Prepare meta data
meta_data_global <- data.frame(flight_id = NA,
                               site_veg = NA,
                               site_name = NA,
                               veg_type = NA,
                               file_path= NA,
                               object_name = NA,
                               sensor = NA,
                               date = as.Date("2017-05-01"),
                               mean_NDVI= NA)

### 1 Function definitions ----

## Helper function to create an extent object from the sentinel boundaries file.
get_sent_extent <- function(site_veg_id, sent_boundaries) {
  extent_object <- extent(
    c(sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmax,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymax))
  return(extent_object)
}


## Function to create satellite + drone time-series plots for 2016
plot_ts_2016 <- function(site_name, veg_type){
  
  # Prep MODIS data
  site_coords_key <- read.csv("data/modis/Modis_plot_key.csv")
  MODIS6_data <- read.csv("data/modis/PS_centre_plots_v6_2000-2017.csv")
  
  MODIS6_qikiqtaruk <- merge(site_coords_key, MODIS6_data, 
                             by.x = c("ModLon", "ModLat"), 
                             by.y = c("longitude", "latitude"))
  MODIS6_qikiqtaruk$NDVI_percent <- MODIS6_qikiqtaruk$NDVI/10000
  
  MODIS6_site <- MODIS6_qikiqtaruk %>% 
    filter(CentrePoint == site_name & 
             VegType == veg_type & 
             year == 2016, 
           SummaryQA < 3 & SummaryQA > -1) %>%
    mutate(date = format(as.Date(paste0(DOY, "/2016"), format = "%j/%Y"), "%d/%m/%Y"))
  rm(MODIS6_data, MODIS6_qikiqtaruk)
  # prep meta_data frame for MODIS
  MODIS_meta <- data.frame(
    flight_id = paste0(MODIS6_site$date, "_id"),
    site_veg = rep(paste0(site_name, "_", veg_type), 
                   length(MODIS6_site$date)),
    site_name = rep(site_name, length(MODIS6_site$date)),
    veg_type = rep(veg_type, length(MODIS6_site$date)),
    file_path = rep(NA, length(MODIS6_site$date)),
    object_name = rep(NA, length(MODIS6_site$date)),
    mean_NDVI = MODIS6_site$NDVI_percent,
    date = as.Date(MODIS6_site$date, "%d/%m/%Y"),
    sensor = rep("MODIS", length(MODIS6_site$date))
  )
  
  ## Prepare loading of drone files
  
  # generate site-specific path to drone files
  folder_path <- paste0(drone_path, 
                        "/2016/", site_name, "_", veg_type, "/output/ndvi_maps")
  
  # optain list of files
  file_names <- list.files(folder_path, pattern = "*.tif")
  
  # create meta-data df
  drone_meta <- data.frame(
    flight_id = substr(file_names, 1, 13),
    site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
    site_name = rep(site_name, length(file_names)),
    veg_type = rep(veg_type, length(file_names)),
    file_path = paste0(folder_path, "/", file_names),
    object_name = gsub(".tif", "", file_names))
  drone_meta$file_path <- as.character(drone_meta$file_path)
  drone_meta$object_name <- as.character(drone_meta$object_name)
  
  # kick out all multiples for 3 August 2016 for PS1 HER and KOM
  if(site_name == "PS1" & veg_type == "HER"){
    drone_meta <- drone_meta[1:6,]
  } else if(site_name == "PS1" & veg_type == "KOM"){
    drone_meta <- drone_meta[1:5,]
  }
  
  ## Prepare loading of sentinel files
  
  ## obtain list of sentinel files
  sentinel_path <- paste0(sentinel_path, "/2016/")
  sentinel_files <- list.files(sentinel_path, pattern = "ndvi.tif")
  
  # add meta_data
  sentinel_meta <- data.frame(
    flight_id = gsub("*.*(2016[0-9][0-9][0-9][0-9]).*", 
                     "\\1", sentinel_files),
    site_veg = rep(paste0(site_name, "_", veg_type), length(sentinel_files)),
    site_name = rep(site_name, length(sentinel_files)),
    veg_type = rep(veg_type, length(sentinel_files)),
    file_path = paste0(sentinel_path, "/", sentinel_files),
    object_name = gsub(".tif", "", sentinel_files))
  sentinel_meta$file_path <- as.character(sentinel_meta$file_path)
  sentinel_meta$object_name <- as.character(sentinel_meta$object_name)
  
  # merge drone and sentinel meta data data frames
  sensor_id <- c(rep("Sentinel 2A", 
                     nrow(sentinel_meta)), 
                 rep("drone", nrow(drone_meta)))
  meta_data <- rbind(sentinel_meta, drone_meta)
  meta_data$sensor <- sensor_id
  
  # tidy up
  rm(list = c("sentinel_meta", "drone_meta"))
  
  # add date collumn
  meta_data$date <- as.Date(gsub("*.*(2016[0-9][0-9][0-9][0-9]).*", 
                                 "\\1", meta_data$flight_id), format = "%Y%m%d")
  
  ## load drone and sentinel files
  list2env(
    lapply(
      setNames(meta_data$file_path, 
               make.names(gsub(".*/", "", 
                               gsub(".tif$", "", meta_data$file_path)))),
      raster), 
    envir = .GlobalEnv)
  
  
  # create site extent object
  site_bounds <- get_sent_extent(
    paste0(site_name, "_", veg_type), 
    site_boundaries)
  
  ## crop files
  list2env(
    lapply(
      setNames(meta_data$object_name, 
               make.names(paste0(meta_data$object_name, "_cropped"))), 
      function(x){crop(get(x), site_bounds)}), 
    envir = .GlobalEnv)
  
  # calcualte means
  meta_data$mean_NDVI <- sapply(paste0(meta_data$object_name, "_cropped"), 
                                function(x) cellStats(get(x), mean, na.rm = T))
  
  # For PS1 HER JA20160730_02, PS1 KOM JA20160730_03 and PS2 KOM JA20160630_02  
  # have no calibraiton imagery and hence NDVI is inflated. 
  # Let's indicate that with a different colour in the graphs
  if (site_name == "PS1") {
    meta_data[meta_data$date == as.Date("2016-07-30"),]$sensor <- 
      "drone_nocalib"
    no_calib_true <- NA
  } else if (site_name == "PS2" & veg_type == "KOM") {
    meta_data[meta_data$date == as.Date("2016-06-30"),]$sensor <- 
      "drone_nocalib"
    no_calib_true <- NA
  } else {
    no_calib_true <- 0
  }
  
  # Merge with MODIS data
  meta_data <- rbind(
    MODIS_meta[,match(names(meta_data), names(MODIS_meta))],
    meta_data)
  
  ## Load Landsat 8 data
  load("data/landsat8/meta_data_ls8_with_mean.Rda")
  meta_data_ls8_with_mean <- meta_data_ls8_with_mean[
    meta_data_ls8_with_mean$site_name == site_name &
      meta_data_ls8_with_mean$veg_type == veg_type &
      format(as.Date(meta_data_ls8_with_mean$date), "%Y") == "2016",] 
  
  # Merge with sentinel, drone and MODIS data
  meta_data <- rbind(meta_data, meta_data_ls8_with_mean)
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2016", "%d/%m/%Y") & 
             date <= as.Date("30/09/2016", "%d/%m/%Y"))
  
  ## Prepare plotting
  # Set colour scale and sort dataframe
  meta_data$sensor <- ordered(meta_data$sensor,
                              levels = c(
                                         "Landsat8",
                                         "Sentinel 2A",
                                         #"Sentinel 2B",
                                         "MODIS",
                                         "drone",
                                         "drone_nocalib"))
  meta_data <- meta_data[order(meta_data$sensor),]
  plot_scale = viridis(6)[c(1,2,6,5)]
  
  # Set y-axis label
  y_label <- "NDVI"
  
  # Adjust visibility of axis labels according to plot positions
  x.axis.text.colour <- "black"
  y.axis.text.colour <- "black"
  if(veg_type == "HER") {
    x.axis.text.colour <- "white"
  }
  
  ## Generate plot
  ts_plot <- ggplot(
    data = filter(meta_data, sensor != "drone_nocalib"), 
    mapping = aes(x = date, 
                  y = mean_NDVI, 
                  color = sensor, 
                  shape = sensor, 
                  size = sensor, 
                  fill = sensor), 
    aes(x = date, y = mean_NDVI), 
    inherit.aes = F) +
    
    geom_point() +
    
    geom_smooth(
      mapping = aes(x = date, y= mean_NDVI), 
      data = filter(meta_data, sensor != "drone_nocalib"),
      se= F, 
      method = "glm",
      formula= y ~ poly(x,2), 
      colour = "gray50", 
      linetype = "dashed", 
      size = 1.5,
      span = 1, 
      fullrange = T,
      inherit.aes = F) +
    
    geom_point() +
    
    scale_size_manual(
      values = c(rep(5, 5))) +
    scale_colour_manual(
      values = rep("black",5)) +
    scale_shape_manual(
      values = c(21, 21, 21, 21, 21)) +
    scale_fill_manual(
      values = plot_scale) +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(-0.15, 1.05), 
      breaks = seq(0,1,0.25), 
      labels = paste0(format(seq(0,1,0.25), digits = 2), " ") ) +
    scale_x_date(
      expand = c(0,0), 
      limits = c(as.Date("2016-05-01"), as.Date("2016-09-30")),
      breaks = as.Date(
        c("2016-05-17",
          "2016-06-16",
          "2016-07-16",
          "2016-08-16",
          "2016-09-16")),
      
      labels = c("M",
                 "J",
                 "J",
                 "A",
                 "S")) +
    annotate(
      geom="segment", 
      y = -0.15, 
      yend = -0.12, 
      size = 1.2,
      x = as.Date(
        c("2016-05-01",
          "2016-06-01",
          "2016-07-01",
          "2016-08-01",
          "2016-09-01",
          "2016-09-30")), 
      xend = as.Date(
        c("2016-05-01",
          "2016-06-01",
          "2016-07-01",
          "2016-08-01",
          "2016-09-01",
          "2016-09-30"))) +
    annotate(
      geom="text", 
      x = as.Date("2016-07-15"), 
      y = 0.98, 
      label = "2016", 
      size =  10, 
      fontface = "bold") +
    ylab(y_label) +
    xlab("") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 1.2),
      axis.title = element_text(size = 24, face = "bold"),
      axis.text.x = element_text(hjust = 0.5, 
                                 size = 20, 
                                 colour = x.axis.text.colour),
      axis.text.y = element_text(size = 20, 
                                 colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.position = "none") 
  
  # Tidy up removing all the rasters loaded
  rm(list = c(
    meta_data[meta_data$sensor != "MODIS" & 
                meta_data$sensor != "Landsat8" & 
                !is.na(meta_data$object_name),]$object_name, 
    paste0(meta_data[meta_data$sensor != "MODIS" & 
                       meta_data$sensor != "Landsat8" &
                       !is.na(meta_data$object_name),]$object_name, 
           "_cropped")), 
    envir = .GlobalEnv)
  
  # Export meta data to global meta_data df
  meta_data_global <<- rbind(meta_data_global, meta_data)
  
  ## Return plot object
  return(ts_plot)
}

## Function to create satellite + drone time-series plots for 2016
plot_ts_2017 <- function(site_name, veg_type){
  
  # Prepare MODIS data
  site_coords_key <- read.csv("data/modis/Modis_plot_key.csv")
  MODIS6_data <- read.csv("data/modis/PS_centre_plots_v6_2000-2017.csv")
  
  MODIS6_qikiqtaruk <- merge(site_coords_key, MODIS6_data, 
                             by.x = c("ModLon", "ModLat"), 
                             by.y = c("longitude", "latitude"))
  MODIS6_qikiqtaruk$NDVI_percent <- MODIS6_qikiqtaruk$NDVI/10000
  
  MODIS6_site <- MODIS6_qikiqtaruk %>% 
    filter(CentrePoint == site_name & 
             VegType == veg_type & 
             year == 2017, 
           SummaryQA < 3 & SummaryQA > -1) %>%
    mutate(date = format(as.Date(paste0(DOY, "/2017"), 
                                 format = "%j/%Y"), "%d/%m/%Y"))
  rm(MODIS6_data, MODIS6_qikiqtaruk)
  
  # prep meta_data frame for MODIS
  MODIS_meta <- data.frame(flight_id = paste0(MODIS6_site$date, "_id"),
                           site_veg = rep(paste0(site_name, "_", veg_type), 
                                          length(MODIS6_site$date)),
                           site_name = rep(site_name, length(MODIS6_site$date)),
                           veg_type = rep(veg_type, length(MODIS6_site$date)),
                           file_path = rep(NA, length(MODIS6_site$date)),
                           object_name = rep(NA, length(MODIS6_site$date)),
                           mean_NDVI = MODIS6_site$NDVI_percent,
                           date = as.Date(MODIS6_site$date, "%d/%m/%Y"),
                           sensor = rep("MODIS", length(MODIS6_site$date))
  )
  
  ##  Prepare drone file meta data 
  
  # Set folder path
  folder_path <- paste0(drone_path,
                        "/2017/", site_name, "_", veg_type, "/output/ndvi_maps")
  # optain list of files
  file_names <- list.files(folder_path, pattern = "*.tif")
  
  # create meta-data df
  drone_meta <- data.frame(
    flight_id = gsub("-", "", substr(file_names, 1, 18)),
    site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
    site_name = rep(site_name, length(file_names)),
    veg_type = rep(veg_type, length(file_names)),
    file_path = paste0(folder_path, "/", file_names),
    object_name = gsub("-", "", gsub(".tif", "", file_names)))
  drone_meta$file_path <- as.character(drone_meta$file_path)
  drone_meta$object_name <- as.character(drone_meta$object_name)
  
  ## Prepare sentinel meta data
  
  # get list of sentinel files
  sentinel_path <- paste0(sentinel_path, "/2017/")
  sentinel_files <- list.files(sentinel_path, pattern = "ndvi.tif")
  
  # add meta_data
  sentinel_meta <- data.frame(
    flight_id = gsub("*.*(2017[0-9][0-9][0-9][0-9]).*", 
                     "\\1", sentinel_files),
    site_veg = rep(paste0(site_name, "_", veg_type), length(sentinel_files)),
    site_name = rep(site_name, length(sentinel_files)),
    veg_type = rep(veg_type, length(sentinel_files)),
    file_path = paste0(sentinel_path, "/", sentinel_files),
    object_name = gsub(".tif", "", sentinel_files))
  sentinel_meta$file_path <- as.character(sentinel_meta$file_path)
  sentinel_meta$object_name <- as.character(sentinel_meta$object_name)
  
  # set sensor info for sentinel images
  sentinel_id <- c(
    "Sentinel 2A",
    "Sentinel 2A",
    "Sentinel 2A",
    'Sentinel 2B',
    "Sentinel 2A",
    'Sentinel 2B',
    "Sentinel 2A",
    'Sentinel 2A',
    'Sentinel 2B',
    'Sentinel 2A',
    'Sentinel 2B',
    'Sentinel 2A',
    'Sentinel 2B',
    "Sentinel 2B",
    'Sentinel 2A')
  
  # merge two meta data dfs
  sensor_id <- c(sentinel_id, rep("drone", nrow(drone_meta)))
  meta_data <- rbind(sentinel_meta, drone_meta)
  meta_data$sensor <- sensor_id
  # add date collumn
  meta_data$date <- as.Date(gsub("*.*(2017[0-9][0-9][0-9][0-9]).*", 
                                 "\\1", meta_data$flight_id), format = "%Y%m%d")
  
  # tidy up
  rm(list = c("sentinel_meta", "drone_meta"))
  
  ## load files
  list2env(
    lapply(
      setNames(meta_data$file_path, 
               make.names(meta_data$object_name)),
      raster), 
    envir = .GlobalEnv)
  
  # generate site boundary extent object
  site_bounds <- get_sent_extent(
    paste0(site_name, "_", veg_type), site_boundaries)
  
  ## crop files
  list2env(
    lapply(
      setNames(meta_data$object_name, 
               make.names(paste0(meta_data$object_name, "_cropped"))), 
      function(x){crop(get(x), site_bounds)}), 
    envir = .GlobalEnv)
  
  # calcualte means
  meta_data$mean_NDVI <- sapply(paste0(meta_data$object_name, "_cropped"), 
                                function(x) cellStats(get(x), mean, na.rm = T))
  
  # Merge sentinel and drone with MODIS data
  meta_data <- rbind(
    MODIS_meta[,match(names(meta_data), names(MODIS_meta))], 
    meta_data)
  
  ## Load Landsat 8 data
  load("data/landsat8/meta_data_ls8_with_mean.Rda")
  meta_data_ls8_with_mean <- meta_data_ls8_with_mean[
    meta_data_ls8_with_mean$site_name == site_name &
    meta_data_ls8_with_mean$veg_type == veg_type &
      format(as.Date(meta_data_ls8_with_mean$date), "%Y") == "2017",] 
  
  # Merge with sentinel, drone and MODIS data
  meta_data <- rbind(meta_data, meta_data_ls8_with_mean)
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2017", "%d/%m/%Y") & 
             date <= as.Date("30/09/2017", "%d/%m/%Y"))
  
  ## Prepare plotting
  
  # Define sorting and scale for plot
  meta_data$sensor <- ordered(meta_data$sensor,
                              levels = c(#"drone_nocalib",
                                         "Landsat8",
                                         "Sentinel 2A",
                                         "Sentinel 2B",
                                         "MODIS",
                                         "drone"))
  meta_data <- meta_data[order(meta_data$sensor),]
  plot_scale = viridis(6)[c(1,2,3,6,5)]
  
  # Set y axis label and colour colour
  if(site_name == "PS1" | site_name == "PS2"){
    y_label <- "NDVI"
    y_label_colour <- "white"
  } else {
    y_label = "NDVI"
    y_label_colour <- "black"
  }
  
  # Adjust visibility of axis labels according to plot positions
  x.axis.text.colour <- "black"
  y.axis.text.colour <- "black"
  if(site_name == "PS1" | site_name == "PS2") {
    y.axis.text.colour <- "white"
  }
  if(veg_type == "HER") {
    x.axis.text.colour <- "white"
  }
  
  ## Plot full time series for site
  ts_plot <- ggplot(
    data = filter(meta_data, !(date == "2017-07-17" & site_veg == "PS4_HER")), 
    mapping = aes(x = date, 
                  y = mean_NDVI, 
                  color = sensor, 
                  shape = sensor, 
                  size = sensor, 
                  fill = sensor), 
    aes(x = date, y = mean_NDVI), 
    inherit.aes = F) +
    
    geom_point() +
    
    geom_smooth(
      mapping = aes(x = date, y= mean_NDVI), 
      data = filter(meta_data, sensor != "drone_nocalib", 
                    !(date == "2017-07-17" & site_veg == "PS4_HER")),
      se= F, 
      method = "glm",
      formula= y ~ poly(x,2), 
      colour = "gray50", 
      linetype = "dashed", 
      span = 1, 
      size = 1.5,
      fullrange = T,
      inherit.aes = F) +
    
    geom_point() +
    
    scale_size_manual(values = c(rep(5, 5))) +
    scale_colour_manual(values = rep("black",5)) +
    scale_shape_manual(values = c(21, 21, 21, 21, 21)) +
    scale_fill_manual(values = plot_scale) +
    scale_y_continuous(
      expand = c(0,0), 
      limits = c(-0.15, 1.05), 
      breaks = seq(0,1,0.25), 
      labels = paste0(format(seq(0,1,0.25), digits = 2), " ")) +
    scale_x_date(
      expand = c(0,0), 
      limits = c(as.Date("2017-05-01"), as.Date("2017-09-30")),
      breaks = as.Date(
        c("2017-05-17",
          "2017-06-16",
          "2017-07-16",
          "2017-08-16",
          "2017-09-16")),
      
      labels = c("M",
                 "J",
                 "J",
                 "A",
                 "S")) +
    annotate(
      geom="segment", 
      y = -0.15, 
      yend = -0.12, 
      size = 1.2,
      x = as.Date(
        c("2017-05-01",
          "2017-06-01",
          "2017-07-01",
          "2017-08-01",
          "2017-09-01",
          "2017-09-30")), 
      xend = as.Date(
        c("2017-05-01",
          "2017-06-01",
          "2017-07-01",
          "2017-08-01",
          "2017-09-01",
          "2017-09-30"))) +
    annotate(
      geom="text", 
      x = as.Date("2017-07-15"), 
      y = 0.98, 
      label = "2017", 
      size =  10, 
      fontface = "bold") +
    ylab(y_label) +
    xlab("") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black", size = 1.2),
      axis.title = element_text(size = 24, 
                                face = "bold", 
                                colour = y_label_colour),
      axis.text.x = element_text(hjust = 0.5, 
                                 size = 20, 
                                 colour = x.axis.text.colour),
      axis.text.y = element_text(size = 20, colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.position = "none") 
  
  # Tidy up removing all the rasters loaded
  rm(list = c(
    meta_data[meta_data$sensor != "MODIS" & 
                meta_data$sensor != "Landsat8" & 
                !is.na(meta_data$object_name),]$object_name, 
    paste0(meta_data[meta_data$sensor != "MODIS" & 
                       meta_data$sensor != "Landsat8" &
                       !is.na(meta_data$object_name),]$object_name, 
           "_cropped")), 
    envir = .GlobalEnv)
  
  # Export meta data to global meta_data df
  meta_data_global <<- rbind(meta_data_global, meta_data)
  
  # Return plot object
  return(ts_plot)
}

### 2 Create time-series plots ----

# Execute function calls to produce plots
PS1_HER_tsplot_2016 <- plot_ts_2016("PS1", "HER")
PS1_KOM_tsplot_2016 <- plot_ts_2016("PS1", "KOM")
PS2_HER_tsplot_2016 <- plot_ts_2016("PS2", "HER")
PS2_KOM_tsplot_2016 <- plot_ts_2016("PS2", "KOM")

PS1_HER_tsplot_2017 <- plot_ts_2017("PS1", "HER")
PS1_KOM_tsplot_2017 <- plot_ts_2017("PS1", "KOM")
PS2_HER_tsplot_2017 <- plot_ts_2017("PS2", "HER")
PS2_KOM_tsplot_2017 <- plot_ts_2017("PS2", "KOM")
PS3_HER_tsplot_2017 <- plot_ts_2017("PS3", "HER")
PS3_KOM_tsplot_2017 <- plot_ts_2017("PS3", "KOM")
PS4_HER_tsplot_2017 <- plot_ts_2017("PS4", "HER")
PS4_KOM_tsplot_2017 <- plot_ts_2017("PS4", "KOM")

# Save plots
save_plot(paste0(ts_out_path,"PS1_HER_2016_tsplot.png"), 
          PS1_HER_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS1_KOM_2016_tsplot.png"), 
          PS1_KOM_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS2_HER_2016_tsplot.png"), 
          PS2_HER_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS2_KOM_2016_tsplot.png"), 
          PS2_KOM_tsplot_2016, base_aspect_ratio = 1)

save_plot(paste0(ts_out_path,"PS1_HER_2017_tsplot.png"), 
          PS1_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS1_KOM_2017_tsplot.png"), 
          PS1_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS2_HER_2017_tsplot.png"), 
          PS2_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS2_KOM_2017_tsplot.png"), 
          PS2_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS3_HER_2017_tsplot.png"), 
          PS3_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS3_KOM_2017_tsplot.png"), 
          PS3_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS4_HER_2017_tsplot.png"), 
          PS4_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(ts_out_path,"PS4_KOM_2017_tsplot.png"), 
          PS4_KOM_tsplot_2017, base_aspect_ratio = 1)

## 3 Create legend plot ----
# drone Landsat8 MODIS Sentinel 2A Sentinel 2B
plot_scale <- viridis(6)
label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 3)) +
  scale_x_continuous(expand = c(0,0), 
               limits = c(0,100)) +
  annotate("text", x = 8, y = 2, 
           label = "MODIS", colour = "black" , fontface = "bold", 
           size = 15, hjust = 0) +
  annotate("point", x = 5, y = 2, 
           shape = 21, colour = "black", 
           fill = plot_scale[6], size = 15) +
  
  annotate("text", x = 33, y = 2, 
           label = "Sentinel 2A", colour = plot_scale[3] , 
           fontface = "bold",size = 15,  hjust = 0) +
  annotate("point", x = 30, y = 2,
           shape = 21, colour = "black", 
           fill = plot_scale[2], size = 15) +
  
  annotate("text", x = 68, y = 2, 
           label = "Sentinel 2B", colour = plot_scale[4] , 
           fontface = "bold",size = 15,  hjust = 0) +
  annotate("point", x = 65, y = 2, 
           shape = 21, colour = "black", 
           fill = plot_scale[3], size = 15) +
  
  annotate("text", x = 28, y = 1.3, 
           label = "Landsat 8", colour = plot_scale[1] , 
           fontface = "bold", size = 15, hjust = 0) +
  annotate("point", x = 25, y = 1.3, 
           shape = 21, colour = "black", 
           fill = plot_scale[1], size = 15) +
  
  annotate("text", x = 60, y = 1.3, 
           label = "Drone", colour = plot_scale[5] , 
           fontface = "bold", size = 15, hjust = 0) +
  annotate("point", x = 57, y = 1.3, 
           shape = 21, colour = "black", 
           fill = plot_scale[5], size = 15) +
  
  # annotate("text", x = 39, y = 1.3, 
  #          label = "Drone, not calibrated", 
  #          colour = plot_scale[2] , fontface = "bold", 
  #          size = 15, hjust = 0) +
  # annotate("point", x = 36, y = 1.3, 
  #          shape = 21, colour = "black", 
  #          fill = plot_scale[2], size = 15) +
  
  xlab("") +
  ylab("") +
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  coord_cartesian()

save_plot(paste0(ts_out_path,"label_plot.png"), label_plot,
          base_aspect_ratio = 3)

### 4 Create Canada map ----
world <- ne_countries(scale = "large", returnclass = "sf")
qhi_location <- data.frame(y = c(69.58), x = c(-139.05), label = "Qikiatruk") %>%
  st_as_sf(coords = c("x", "y"),crs = "+proj=longlat +datum=WGS84 +no_defs")
canada_map <- ggplot() + geom_sf(data = world, fill = "#ffffffFF", size = 0.5) + 
  geom_sf(data = qhi_location, colour = "white",
             fill = "#1e5c91FF", shape = 21, size = 8) +  
  coord_sf(xlim = c(-170, -100), ylim = c(55.7, 81.3), expand = F) +
  theme(legend.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, "mm"),
        panel.grid.minor = element_line(colour = "white"),
        panel.grid.major = element_line(colour = "white"),
        plot.margin=grid::unit(c(0,0,0,0), "mm"))
# theme(panel.background = element_rect(fill = "lightskyblue"))
save_plot(canada_map, filename = "figures/fig_1_ts_plots/canada_map.png",
          base_aspect_ratio = 1, base_height = 5)  

### 5 Prepare study site boundaries

# Load ps site boundaries
load("data/site_boundaries/ps_sent_site_bounds.Rda")

# Define function to create a four corner definition of the site boundaries
expand_coords <- function(site_veg) {
  x1 <- new_bounds[new_bounds$site_veg == site_veg,]$xmin
  y1 <- new_bounds[new_bounds$site_veg == site_veg,]$ymin
  
  x2 <- new_bounds[new_bounds$site_veg == site_veg,]$xmax
  y2 <- new_bounds[new_bounds$site_veg == site_veg,]$ymin
  
  x3 <- new_bounds[new_bounds$site_veg == site_veg,]$xmax
  y3 <- new_bounds[new_bounds$site_veg == site_veg,]$ymax
  
  x4 <- new_bounds[new_bounds$site_veg == site_veg,]$xmin
  y4 <- new_bounds[new_bounds$site_veg == site_veg,]$ymax
  
  return__df <- data.frame(site_veg = site_veg,
                           point_id = seq(1,4),
                           x = c(x1, x2, x3, x4),
                           y = c(y1, y2, y3, y4))
  
}

# Apply function to loaded site coordiantes
site_coords <- bind_rows(lapply(new_bounds$site_veg,expand_coords))
# Add a veg_type coloumn to the just created data frame
site_coords$veg_type <- substr(site_coords$site_veg,5,7)

# Convert data frame into sf polygon object
sites_coords_sf <- site_coords %>% 
  dplyr::select(site_veg, x, y) %>%
  mutate(x_1 = x, y_1 = y) %>%
  group_by(site_veg) %>%
  st_as_sf( 
    coords = c("x", "y"), 
    crs = "+proj=utm +zone=7 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
sites_coords_sf <- sites_coords_sf %>% summarise() %>% st_convex_hull() 

# Write out shapefile for map
st_write(sites_coords_sf, "data/site_boundaries/ps_site_bounds.shp",
         delete_dsn =  T)

### The rest of the map figure is assembled in QGIS - 
### See QGIS project file and layout in:
### figures/fig_1_ts_plots/fig_1_drone_satellite_ts_map.qgz

### 6) Calculate correlations and offsets between sensors: ----


# sensor_peak_season_mean <- meta_data_global %>% 
#   mutate(year = as.numeric(format.Date(date, "%Y"))) 
# sensor_peak_season_mean$sensor[sensor_peak_season_mean$sensor == "Sentinel 2A"] <- "sentinel"
# sensor_peak_season_mean$sensor[sensor_peak_season_mean$sensor == "Sentinel 2B"] <- "sentinel"
# 
# sensor_peak_season_mean <- sensor_peak_season_mean %>%
#   filter(as.numeric(format.Date(date, "%j")) >= 196 & 
#            as.numeric(format.Date(date, "%j")) <= 217 # &
#            # as.numeric(format.Date(date, "%Y")) == 2017
#          ) %>%
#   group_by(site_name, veg_type, sensor, year) %>%
#   summarise(mean_NDVI = mean(mean_NDVI, na.rm = T)) %>%
#   filter(sensor != "drone_nocalib")
# 
# 
# # combine_sentinel <- function(x) {
# #   if(x == "sentinel") return("sentinel")
# #   else if(x == "Sentinel 2A") return("sentinel")
# #   else if(x == "Sentinel 2B") return("sentinel")
# #   else return(x)
# # }
# # sensor_peak_season_mean$sensor <- 
#   # modify(sensor_peak_season_mean$sensor, combine_sentinel) 
#   # sensor_peak_season_mean$sensor <- 
#   #   modify(sensor_peak_season_mean$sensor, combine_sentinel) 
#   
#   
# sensor_peak_season_mean_wide <- pivot_wider(sensor_peak_season_mean,
#             names_from = "sensor",
#             values_from = "mean_NDVI") %>% na.omit()
# peak_season_cor <- sensor_peak_season_mean_wide %>%
#   ungroup() %>%
#   dplyr::select(drone, sentinel, MODIS, Landsat8) %>% as.data.frame() %>%
#   cor() %>% round(2)
# 
# sensor_peak_season_diff <- sensor_peak_season_mean %>%
#   group_by(site_name, veg_type) %>%
#   group_map(function(subset, groupings){
#     return(data.frame(
#       site_name = groupings[1],
#       veg_type = groupings[2],
#       modis_drone = abs(subset[subset$sensor == "MODIS",]$mean_NDVI) -
#         abs(subset[subset$sensor == "drone",]$mean_NDVI),
#       modis_sentinel = abs(subset[subset$sensor == "MODIS",]$mean_NDVI -
#         subset[subset$sensor == "sentinel",]$mean_NDVI),
#       sentinel_drone = abs(subset[subset$sensor == "sentinel",]$mean_NDVI -
#         subset[subset$sensor == "drone",]$mean_NDVI)
#     ))
#   }) %>% bind_rows()
# sensor_peak_season_diff_mean <- sensor_peak_season_diff %>% 
#   summarise('difference MODIS:drone' = round(mean(modis_drone),3),
#             'difference MODIS:sentinel' = round(mean(modis_sentinel),3),
#             'difference Sentinel:drone' = round(mean(sentinel_drone),3))
# # Save / export data
# write.csv(peak_season_cor,
#           file = paste0(data_out_path, "sensor_peak_season_cor.csv"))
# write.csv(sensor_peak_season_diff_mean,
#           file = paste0(data_out_path, "sensor_peak_season_mean_diff.csv"))
# 
# save(file = paste0(data_out_path, "meta_data_with_means.Rda"),
#      meta_data_global)
# # load(paste0(data_out_path, "meta_data_with_means.Rda"))

# Divide season into 7 day block startin on the Set doy range
first_march <- c(122, 121) # 2016 was a leap year, 2017
last_sept <- c(274, 273) # 2016, 2017
week_blocks <- data.frame(year = c(rep(2016, length(seq(first_march[1], last_sept[1], 7))),
                    rep(2017, length(seq(first_march[2], last_sept[2], 7)))),
  block_no = c(paste0("2016_", 1:length(seq(first_march[1], last_sept[1], 7))),
               paste0("2017_", 1:length(seq(first_march[2], last_sept[2], 7)))),
  block_strat_doy = c(seq(first_march[1], last_sept[1], 7),
                      seq(first_march[2], last_sept[2], 7)),
  block_end_doy = c(seq(first_march[1], last_sept[1], 7) + 6,
                    seq(first_march[2], last_sept[2], 7) + 6),
  stringsAsFactors = F)

# Create doy columns for easier handling
meta_data_global <- meta_data_global %>% mutate(doy = as.numeric(format.Date(date, "%j")),
                            year = as.numeric(format.Date(date, "%Y")))

# Calculate weekly averges
weekly_means <- lapply(week_blocks$block_no, function(block_no){
  start_doy <- week_blocks$block_strat_doy[week_blocks$block_no == block_no]
  end_doy <- week_blocks$block_end_doy[week_blocks$block_no == block_no]
  year <- week_blocks$year[week_blocks$block_no == block_no]
  
  block_mean <- meta_data_global %>% 
    filter(doy >= start_doy, doy <= end_doy, year == year) %>%
    filter(sensor != "drone_nocalib") %>%
    group_by(sensor) %>% 
    summarise(mean_NDVI = mean(mean_NDVI))
  
  if(nrow(block_mean) != 0) {
    block_mean$year <- year
    block_mean$block_no <- block_no
  } else {
    block_mean <- data.frame(
      sensor = NA,
      mean_NDVI = NA,
      year = year,
      block_no = block_no
    )
  }
  return(block_mean)
}) %>% bind_rows %>%
  na.omit

# Determine sensor combinations
sensor_combinations <- expand.grid(unique(weekly_means$sensor),
                                   unique(weekly_means$sensor))

# Calculate correlations
correlations <- bind_rows(apply(
  sensor_combinations, 
  1, 
  function(sensors){
    sensor_1 <- sensors[1]
    sensor_2 <- sensors[2]
    sensor_1_means <- filter(weekly_means, sensor == sensor_1) %>%
      mutate(sensor_1 = sensor,
             mean_NDVI_1 = mean_NDVI) %>% 
      dplyr::select(-sensor, -mean_NDVI)
    sensor_2_means <- filter(weekly_means, sensor == sensor_2) %>%
      mutate(sensor_2 = sensor,
             mean_NDVI_2 = mean_NDVI) %>% 
      dplyr::select(-sensor, -mean_NDVI)
    combo_means <- full_join(sensor_1_means, sensor_2_means) %>% 
      na.omit()
    sensor_cor <- combo_means %>% 
      group_by(sensor_1, sensor_2) %>% 
      summarise(cor = cor(mean_NDVI_1, mean_NDVI_2, 
                          method = "spearman"),
                n = n())
    return(sensor_cor)
}))

# Convert into cor matrix
cor_matrix <- correlations %>% 
  arrange(sensor_1, sensor_2) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = sensor_2, values_from = cor) 
n_pairs <-  correlations %>% 
  arrange(sensor_1, sensor_2) %>%
  dplyr::select(-cor) %>%
  pivot_wider(names_from = sensor_2, values_from = n)

# Prepare for output
cor_matrix_for_output <- cor_matrix
for(i in 2:ncol(cor_matrix_for_output)){
  for(j in 1:nrow(cor_matrix_for_output)){
    cor_matrix_for_output[j,i] <- paste0(formatC(round(as.numeric(cor_matrix[j,i]),2),
                                                 digits = 2,
                                                 format = "f"), 
                                      " (", 
                                      formatC(as.numeric(n_pairs[j,i]),
                                              width = 2), ")")
    }
}

names(cor_matrix_for_output) <- c("Sensor",
                                  "Drone",
                                  "Landsat 8",
                                  "MODIS",
                                  "Sentinel 2A",
                                  "Sentinel 2B")
cor_matrix_for_output$Sensor <- c(
  "Drone",
  "Landsat 8",
  "MODIS",
  "Sentinel 2A",
  "Sentinel 2B")

# Export
write.csv(cor_matrix_for_output, 
          "data/fig_1_satellite_drone_ts_map/sensor_cor.csv",
          row.names = F)

# Calculate mean differences per weekly block
mean_diffs <- bind_rows(apply(
  sensor_combinations, 
  1, 
  function(sensors){
    sensor_1 <- sensors[1]
    sensor_2 <- sensors[2]
    sensor_1_means <- filter(weekly_means, sensor == sensor_1) %>%
      mutate(sensor_1 = sensor,
             mean_NDVI_1 = mean_NDVI) %>% 
      dplyr::select(-sensor, -mean_NDVI)
    sensor_2_means <- filter(weekly_means, sensor == sensor_2) %>%
      mutate(sensor_2 = sensor,
             mean_NDVI_2 = mean_NDVI) %>% 
      dplyr::select(-sensor, -mean_NDVI)
    combo_means <- full_join(sensor_1_means, sensor_2_means) %>% 
      na.omit()
    sensor_cor <- combo_means %>% 
      group_by(sensor_1, sensor_2) %>% 
      summarise(mean_diff = mean(mean_NDVI_1 - mean_NDVI_2),
                n = n())
    return(sensor_cor)
  }))

# Convert into diff matrix
diff_matrix <- mean_diffs %>% 
  arrange(sensor_1, sensor_2) %>%
  dplyr::select(-n) %>%
  pivot_wider(names_from = sensor_2, values_from = mean_diff) 
n_pairs <-  mean_diffs %>% 
  arrange(sensor_1, sensor_2) %>%
  dplyr::select(-mean_diff) %>%
  pivot_wider(names_from = sensor_2, values_from = n)

# Prepare for output
diff_matrix_for_output <- as.data.frame(diff_matrix)
for(i in 2:ncol(diff_matrix_for_output)){
  for(j in 1:nrow(diff_matrix_for_output)){
    diff_matrix_for_output[j,i] <- paste0(formatC(as.numeric(round(diff_matrix[j,i],3)),
                                                 digits = 3,
                                                 flag = " ",
                                                 format = "f"), 
                                         " (", 
                                         formatC(as.numeric(n_pairs[j,i]),
                                                 width = 2), ")")
  }
}


names(diff_matrix_for_output) <- c("Sensor",
                                  "Drone",
                                  "Landsat 8",
                                  "MODIS",
                                  "Sentinel 2A",
                                  "Sentinel 2B")
diff_matrix_for_output$Sensor <- c(
  "Drone",
  "Landsat 8",
  "MODIS",
  "Sentinel 2A",
  "Sentinel 2B")

# Export
write.csv(diff_matrix_for_output, 
          "data/fig_1_satellite_drone_ts_map/sensor_diff.csv",
          row.names = F)


### 7) Data overview tables ----

## Drone flights
drone_flights <- meta_data_global %>% 
  filter(sensor == "drone") %>%
  group_by(site_name, veg_type, flight_id) %>%
  distinct(date) %>%
  ungroup() %>%
  dplyr::select(site_name, veg_type, date, flight_id) %>%
  arrange(site_name, veg_type, date, flight_id)

# Add meta-data from flight logs
flight_log_meta <- read.csv("data/auxillary/flight_log_meta_data.csv", 
                            stringsAsFactors = F)
flight_log_meta$date <- as.Date(flight_log_meta$date, "%d/%m/%Y")
drone_flights <- merge(drone_flights, flight_log_meta,
                       by.x = "flight_id",
                       by.y = "flight_id") 
drone_flights$date.x == drone_flights$date.y
drone_flights$date <- drone_flights$date.x
drone_flights <- dplyr::select(drone_flights, -date.x, -date.y)
# Calculate solar elevation
drone_flights$posix <- as.POSIXct(paste(drone_flights$date, drone_flights$time),
                                  format = "%Y-%m-%d %H:%M:%OS",
                                  tz = "MST7MDT")
drone_flights$posix_utc <- as.POSIXct(format(drone_flights$posix, tz = "UTC"),
                                      tz = "UTC")

# Calculate solar elevation at sea-level
drone_flights$solar_elevation <- sunAngle(drone_flights$posix_utc,
                                          lat = 69.57, lon = -138.91)$altitude
drone_flights$solar_azimuth <- sunAngle(drone_flights$posix_utc,
                                  lat = 69.57, lon = -138.91)$azimuth

# Load time of solar noon and calculate difference in time to solar noon
solar_noon <- read.csv("data/auxillary/solar_noon.csv")
solar_noon$solar_noon_posix <- as.POSIXct(
  paste(solar_noon$date, solar_noon$solar_noon.UTC.7.),
  format = "%Y-%m-%d %H:%M:%OS",
  tz ="America/Vancouver")
solar_noon$solar_noon_posix_utc <- as.POSIXct(
  format(solar_noon$solar_noon_posix, tz = "UTC"),
  tz = "UTC")
solar_noon$date <- as.Date(solar_noon$date)

drone_flights <- merge(drone_flights, solar_noon, by.x = "date", by.y = "date")
drone_flights$solar_noon_diff <- difftime(drone_flights$posix_utc,
                                          drone_flights$solar_noon_posix_utc,
                                          units = "hours")

# Tidy up
drone_flights <- dplyr::select(drone_flights,
                        site_name,
                        veg_type,
                        date,
                        time,
                        solar_noon_diff,
                        solar_elevation,
                        solar_azimuth,
                        Aircraft_ID,
                        Sensor_ID,
                        Skye_Code,
                        Notes)

drone_flights <- arrange(drone_flights,
                         site_name, veg_type, date)

drone_flights$solar_noon_diff <- round(drone_flights$solar_noon_diff, 2)
drone_flights$solar_elevation <- round(drone_flights$solar_elevation)
drone_flights$solar_azimuth <- round(drone_flights$solar_azimuth)

# Calculate mean difference to solar noon
mean(abs(drone_flights$solar_noon_diff), na.rm = T)

# Filter time 
# Helper function to produce pretty names
pretty_name <- function(value){
  site_names_pretty <- data.frame(
    site_name = as.character(unique(drone_flights$site_name)),
    site_pretty = c(
      "Area 1 - Collinson Head",
      "Area 2 - Bowhead Ridge",
      "Area 3 - Hawk Valley",
      "Area 4 - Hawk Ridge"
    ), stringsAsFactors = F)
  if(value == site_names_pretty[1,1]) return_value <- site_names_pretty[1,2]
  if(value == site_names_pretty[2,1]) return_value <- site_names_pretty[2,2]
  if(value == site_names_pretty[3,1]) return_value <- site_names_pretty[3,2]
  if(value == site_names_pretty[4,1]) return_value <- site_names_pretty[4,2]
  if(value == "HER") return_value <- "Tussock Sedge Tundra"
  if(value == "KOM") return_value <- "Dryas-Vetch Tundra"
  return(return_value)
}
drone_flights$site_name <- as.character(drone_flights$site_name)
drone_flights$veg_type <- as.character(drone_flights$veg_type)
drone_flights$site_name <- modify(drone_flights$site_name, pretty_name)
drone_flights$veg_type <- modify(drone_flights$veg_type, pretty_name)

names(drone_flights) <- c("Area Name", 
                          "Vegetation Type", 
                          "Date",
                          "Time (UTC-6)",
                          "Diff. to Solar Noon (h)",
                          "Solar Elevation",
                          "Solar Azimuth",
                          "Drone Platform",
                          "Sensor",
                          "Sky Code",
                          "Notes")

write.csv(file = paste0(data_out_path, "drone_flights.csv"), 
          drone_flights, 
          row.names = F)


# Satellite Data Sensing dates
# These are the same for all site veg combinations
meta_data_global$sensor[meta_data_global$sensor == "sentinel"] <- "Sentinel 2A"
sentinel_scenes <- meta_data_global %>% 
  filter(sensor == "Sentinel 2A" | sensor ==  "Sentinel 2B" |
           sensor == "Sentinel") %>%
  group_by(sensor) %>% 
  distinct(date) %>% 
  dplyr::select(sensor, date) %>%
  arrange(sensor, date)  
write.csv(file = paste0(data_out_path, "sentinel_scenes.csv"), 
          sentinel_scenes, 
          row.names = F)


# Modis Pixels
# Due to quality control the dates are not the same for all sites
meta_data_global %>%
  filter(sensor == "MODIS") %>%
  group_by(site_veg) %>% 
  distinct(date) %>% summarise(n = n())

modis_pixels <- meta_data_global %>% 
  filter(sensor == "MODIS") %>%
  group_by(site_name, veg_type) %>%
  distinct(date) %>%
  ungroup() %>%
  dplyr::select(site_name, veg_type, date) %>%
  arrange(site_name, veg_type, date)

# Helper function to produce pretty names
pretty_name <- function(value){
  site_names_pretty <- data.frame(
    site_name = as.character(unique(modis_pixels$site_name)),
    site_pretty = c(
      "Area 1 - Collinson Head",
      "Area 2 - Bowhead Ridge",
      "Area 3 - Hawk Valley",
      "Area 4 - Hawk Ridge"
    ), stringsAsFactors = F)
  if(value == site_names_pretty[1,1]) return_value <- site_names_pretty[1,2]
  if(value == site_names_pretty[2,1]) return_value <- site_names_pretty[2,2]
  if(value == site_names_pretty[3,1]) return_value <- site_names_pretty[3,2]
  if(value == site_names_pretty[4,1]) return_value <- site_names_pretty[4,2]
  if(value == "HER") return_value <- "Tussock Sedge Tundra"
  if(value == "KOM") return_value <- "Dryas-Vetch Tundra"
  return(return_value)
}

modis_pixels$site_name <- as.character(modis_pixels$site_name)
modis_pixels$veg_type <- as.character(modis_pixels$veg_type)
modis_pixels$site_name <- modify(modis_pixels$site_name, pretty_name)
modis_pixels$veg_type <- modify(modis_pixels$veg_type, pretty_name)

names(modis_pixels) <- c("Area Name", "Vegetation Type", "Date")

write.csv(file = paste0(data_out_path, "modis_pixels.csv"), 
          modis_pixels,
          row.names = F)
