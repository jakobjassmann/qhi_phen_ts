# Phenology Time-Series Quick Script to create a map of all sites in the
# time-series.
# Jakob Assmann j.assmann@ed.ac.uk 11 October 2018

# Dependencies
library(dplyr)
library(raster)
library(rgdal)
library(rasterVis)
library(viridisLite)

library(maps)
library(mapdata)

library(cowplot)
library(ggmap)
library(magick)
library(gridExtra)

library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)


### Global prepartions ----

# Set global variables
site_boundaries <- 
  read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
ts_out_path <- "figures/fig_1_ts_plots/"
sentinel_path <- 
  "/Volumes/BowheadRdge/phen_time_series/sentinel_data/qhi_cld_free/"
drone_path <- "/Volumes/BowheadRdge/phen_time_series/final_outputs"

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
  
  # !!!!!!!!!!!!!!!!!!! The following section needs to be revised !!!!!!!!!!!!!
  # Prep MODIS data
  site_coords_key <- read.csv("data/modis/Modis_plot_key.csv")
  MODIS6_data <- read.csv("data/modis/PS_centre_plots_v6_2000-2017.csv")
  
  MODIS6_qikiqtaruk <- merge(site_coords_key, MODIS6_data, 
                             by.x = c("ModLon", "ModLat"), 
                             by.y = c("longitude", "latitude"))
  MODIS6_qikiqtaruk$NDVI_percent <- MODIS6_qikiqtaruk$NDVI/10000
  
  MODIS6_site <- MODIS6_qikiqtaruk %>% 
    filter(CentrePoint == site_name & VegType == veg_type & year == 2016, SummaryQA < 3 & SummaryQA > -1) %>%
    mutate(date = format(as.Date(paste0(DOY, "/2016"), format = "%j/%Y"), "%d/%m/%Y"))
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
  
  ### !!!!!!!!!!!!!!!! End of revsions needed !!!!!!!!!!!!!!!!!!!!!!!!!
  
  ## Prepare loading of drone files
  
  # generate site-specific path to drone files
  folder_path <- paste0(drone_path, 
                        "/2016/", site_name, "_", veg_type, "/output/ndvi_maps")
  
  # optain list of files
  file_names <- list.files(folder_path, pattern = "*.tif")
  
  # create meta-data df
  drone_meta <- data.frame(flight_id = substr(file_names, 1, 13),
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
    drone_meta <- drone_meta[1:4,]
  }
  
  ## Prepare loading of sentinel files
  
  ## obtain list of sentinel files
  sentinel_path <- paste0(sentinel_path, "/2016/")
  sentinel_files <- list.files(sentinel_path, pattern = "ndvi.tif")
  
  # add meta_data
  sentinel_meta <- data.frame(flight_id = gsub("*.*(2016[0-9][0-9][0-9][0-9]).*", 
                                               "\\1", sentinel_files),
                              site_veg = rep(paste0(site_name, "_", veg_type), length(sentinel_files)),
                              site_name = rep(site_name, length(sentinel_files)),
                              veg_type = rep(veg_type, length(sentinel_files)),
                              file_path = paste0(sentinel_path, "/", sentinel_files),
                              object_name = gsub(".tif", "", sentinel_files))
  sentinel_meta$file_path <- as.character(sentinel_meta$file_path)
  sentinel_meta$object_name <- as.character(sentinel_meta$object_name)
  
  # merge drone and sentinel meta data data frames
  sensor_id <- c(rep("sentinel", nrow(sentinel_meta)), rep("drone", nrow(drone_meta)))
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
    meta_data[meta_data$date == as.Date("2016-07-30"),]$sensor <- "drone_nocalib"
    no_calib_true <- NA
  } else if (site_name == "PS2" & veg_type == "KOM") {
    meta_data[meta_data$date == as.Date("2016-06-30"),]$sensor <- "drone_nocalib"
    no_calib_true <- NA
  } else {
    no_calib_true <- 0
  }
  
  # Merge with MODIS data
  meta_data <- rbind(
    MODIS_meta[,match(names(meta_data), names(MODIS_meta))],
    meta_data)
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2016", "%d/%m/%Y") & 
             date <= as.Date("30/09/2016", "%d/%m/%Y"))
  
  ## Prepare plotting
  # Set colour scale
  plot_scale = viridis(5)[c(5,2,3,4)]
  
  # There is no calibration data for Ps2 HER for some dates, 
  # adjust colour scale and KOM is messed up correct!
  if((site_name == "PS2" & veg_type == "HER")){
    plot_scale = viridis(5)[c(5,2,3)]
  } else if(site_name == "PS2" & veg_type == "KOM"){
    plot_scale = viridis(5)[c(5,2,4,3)]
  }
  
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
    data = meta_data, 
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
      colour = "grey28", 
      linetype = "dashed", 
      size = 2,
      span = 1, 
      fullrange = T,
      inherit.aes = F) +
    
    geom_point() +
    
    scale_size_manual(
      values = c(rep(5, 4))) +
    scale_colour_manual(
      values = rep("black",4)) +
    scale_shape_manual(
      values = c(21, 21, 21, 21)) +
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
      axis.text.x = element_text(hjust = 0.5, size = 20, colour = x.axis.text.colour),
      axis.text.y = element_text(size = 20, colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.position = "none") 
  
  # Tidy up removing all the rasters loaded
  rm(list = c(meta_data[meta_data$sensor != "MODIS",]$object_name, 
              paste0(meta_data[meta_data$sensor != "MODIS",]$object_name, 
                     "_cropped")), 
     envir = .GlobalEnv)
  
  # Export meta data to global meta_data df
  meta_data_global <<- rbind(meta_data_global, meta_data)
  
  ## Return plot object
  return(ts_plot)
}

## Function to create satellite + drone time-series plots for 2016
plot_ts_2017 <- function(site_name, veg_type){
  
  ##### !!!! The following section needs to be revised !!!!
  
  # Prepare MODIS data
  site_coords_key <- read.csv("data/modis/Modis_plot_key.csv")
  MODIS6_data <- read.csv("data/modis/PS_centre_plots_v6_2000-2017.csv")
  
  MODIS6_qikiqtaruk <- merge(site_coords_key, MODIS6_data, 
                             by.x = c("ModLon", "ModLat"), 
                             by.y = c("longitude", "latitude"))
  MODIS6_qikiqtaruk$NDVI_percent <- MODIS6_qikiqtaruk$NDVI/10000
  
  MODIS6_site <- MODIS6_qikiqtaruk %>% 
    filter(CentrePoint == site_name & VegType == veg_type & year == 2017, SummaryQA < 3 & SummaryQA > -1) %>%
    mutate(date = format(as.Date(paste0(DOY, "/2017"), format = "%j/%Y"), "%d/%m/%Y"))
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
  
  
  ### !!!!!!!!!! End of revision section !!!!!
  
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
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2017", "%d/%m/%Y") & 
             date <= as.Date("30/09/2017", "%d/%m/%Y"))
  
  ## Prepare plotting
  
  # Define scale
  plot_scale = viridis(5)[c(5,2,1,3)]
  
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
    data = meta_data, 
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
      colour = "grey28", 
      linetype = "dashed", 
      span = 1, 
      size = 2,
      fullrange = T,
      inherit.aes = F) +
    
    geom_point() +
    
    scale_size_manual(values = c(rep(5, 4))) +
    scale_colour_manual(values = rep("black",4)) +
    scale_shape_manual(values = c(21, 21, 21, 21)) +
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
      axis.title = element_text(size = 24, face = "bold", colour = y_label_colour),
      axis.text.x = element_text(hjust = 0.5, size = 20, colour = x.axis.text.colour),
      axis.text.y = element_text(size = 20, colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.position = "none") 
  
  # Tidy up removing all the rasters loaded
  rm(list = c(
    meta_data[meta_data$sensor != "MODIS",]$object_name, 
    paste0(meta_data[meta_data$sensor != "MODIS",]$object_name, 
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
plot_scale = viridis(5)[c(3,4,2,1,5)]
label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 3)) +
  scale_x_continuous(expand = c(0,0), 
               limits = c(0,100)) +
  annotate("text", x = 8, y = 2, label = "MODIS", colour = "black" , fontface = "bold", size = 15, hjust = 0) +
  annotate("point", x = 5, y = 2, shape = 21, colour = "black", fill = plot_scale[5], size = 15) +
  
  annotate("text", x = 33, y = 2, label = "Sentinel 2A", colour = plot_scale[3] , fontface = "bold",size = 15,  hjust = 0) +
  annotate("point", x = 30, y = 2, shape = 21, colour = "black", fill = plot_scale[3], size = 15) +
  
  annotate("text", x = 68, y = 2, label = "Sentinel 2B", colour = plot_scale[4] , fontface = "bold",size = 15,  hjust = 0) +
  annotate("point", x = 65, y = 2, shape = 21, colour = "black", fill = plot_scale[4], size = 15) +
  
  annotate("text", x = 15, y = 1.3, label = "Drone", colour = plot_scale[1] , fontface = "bold", size = 15, hjust = 0) +
  annotate("point", x = 12, y = 1.3, shape = 21, colour = "black", fill = plot_scale[1], size = 15) +
  
  annotate("text", x = 39, y = 1.3, label = "Drone, not calibrated", colour = plot_scale[2] , fontface = "bold", size = 15, hjust = 0) +
  annotate("point", x = 36, y = 1.3, shape = 21, colour = "black", fill = plot_scale[2], size = 15) +
  
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

### 4 Make Canada map ----
world <- ne_countries(scale = "large", returnclass = "sf")
qhi_location <- data.frame(y = c(69.58), x = c(-139.05), label = "Qikiatruk") %>%
  st_as_sf(coords = c("x", "y"),crs = "+proj=longlat +datum=WGS84 +no_defs")
canada_map <- ggplot() + geom_sf(data = world, fill = "#ffffffFF", size = 0.5) + 
  geom_sf(data = qhi_location, colour = "white",
             fill = "#1e5c91FF", shape = 21, size = 8) +  
  coord_sf(xlim = c(-180, -100), ylim = c(50, 82.3), expand = F) +
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

### 4 Prepare final map figure ----

# Load qhi / yukon boundaries
qhi_boundaries_shp <- readOGR("data/auxillary/yukon_bounds_utmz7.shp")
crs(qhi_boundaries_shp) <- 
  crs("+proj=utm +zone=7 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# convert into a spatial points object to allow for seamless cropping
lines_df <- as(qhi_boundaries_shp, "SpatialLinesDataFrame")
points_df <- as(lines_df, "SpatialPointsDataFrame")
points_df <- data.frame(group = points_df@data$ID,
                        x = points_df@coords[,1],
                        y = points_df@coords[,2]) 

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

sites_coords_sf <- site_coords %>% 
  dplyr::select(site_veg, x, y) %>%
  mutate(x_1 = x, y_1 = y) %>%
  group_by(site_veg) %>%
  st_as_sf( 
    coords = c("x", "y"), 
    crs = "+proj=utm +zone=7 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
sites_coords_sf <- sites_coords_sf %>% summarise() %>% st_convex_hull() 

st_write(sites_coords_sf, "data/site_boundaries/ps_site_bounds.shp",
         delete_dsn =  T)
plot(sites_coords_sf)
# Define global properties for the map figure
plot_width <- 5000
plot_height <- 7500
cropping_margin <- 2500 # this defines how much exess to leave for the polygons

center_x <- 582000
center_y <- 7720750

x_min <- center_x - 0.5 * plot_width - cropping_margin
x_max <- center_x + 0.5 * plot_width + cropping_margin
y_min <- center_y - 0.5 * plot_height - cropping_margin
y_max <- center_y + 0.5 * plot_height + cropping_margin

# Scale bar
scale_bar_start_x <- 583400
scale_bar_start_y <- 7717100
scale_bar_length <- 1000

# Scale breaks for axes
scale_breaks <- data.frame(
  x = c(seq(580000, 584000, 2000)),
  y = c(seq(7717500,7724500, 3000))
)
scale_breaks <- SpatialPoints(
  scale_breaks, 
  proj4string = 
    CRS("+proj=utm +zone=7 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 "))
scale_breaks_latlong <- spTransform(
  scale_breaks, 
  CRS("+proj=longlat +datum=WGS84"))

# Global colours
veg_colours <- c("#1E9148FF", "#1E5C91FF")

## Create plot of site locations for main panel of map
sites_plot <- ggplot() +
  # Add QHI boundary
  geom_polygon(data=points_df, 
               aes(x=x, y=y, group=group),
               lwd = 1, 
               colour = "black", 
               fill = 'white') +
  # Add Site squares
  geom_polygon(data=site_coords, 
               aes(x=x, y=y, 
                   group=site_veg, 
                   fill = veg_type),
               lwd = 0,
               colour = "black",
               inherit.aes = F) +
  # Define scales
  scale_fill_manual(values = veg_colours) +
  scale_colour_manual(values = veg_colours) +
  scale_x_continuous(
    limits = c(x_min, x_max), 
    breaks = scale_breaks$x[1:5],
    labels= paste0(
      format((round(scale_breaks_latlong$x[1:5],2) * -1), 
             digits = 6), 
      "°W"),
    expand = expand_scale(add = c(-cropping_margin, -cropping_margin))) +
  scale_y_continuous(
    limits = c(y_min, y_max), 
    breaks = scale_breaks$y,
    labels= paste0(round(scale_breaks_latlong$y,2), "°N"),
    expand = expand_scale(add = c(-cropping_margin, -cropping_margin))) +
  # Plot scale bar
  annotate("segment", 
           x = scale_bar_start_x, 
           xend = scale_bar_start_x + scale_bar_length,
           y = scale_bar_start_y,
           yend = scale_bar_start_y,
           colour = "black", 
           size = 1.25)+
  annotate("text", 
           label = "1 km",
           x = scale_bar_start_x + 0.5 * scale_bar_length,
           y = scale_bar_start_y + 125,
           size = 3,
           hjust = 0.5) +
  # Add a blank rectangle for some reason
  annotate("rect", xmin = 579600 , xmax = 580300, ymin = 7719590 , ymax = 7719730, colour = NA, fill = "white", alpha = 1) +
  # Add site labels
  annotate("text", label = "Collinson Head", x = 583050, y = 7719650, size = 3, vjust = 0.5, hjust = 1 ) +
  annotate("text", label = "Bowhead Ridge", x = 582250, y = 7721100, size = 3, vjust = 0.5, hjust = 0 ) +
  annotate("text", label = "Hawk Valley", x = 580300, y = 7719650, size = 3, vjust = 0.5, hjust = 1) +
  annotate("text", label = "Hawk Ridge", x = 580300, y = 7721550, size = 3, vjust = 0.5, hjust = 1 ) +
  annotate("text", label = "Veg. Types:", x = 580300, y = 7717100, size = 3, colour = 'black', vjust = 0, hjust = 1) +
  annotate("text", label = "Dryas-Vetch Tundra", x = 580400, y = 7717100, size = 3, colour = veg_colours[2], vjust = 0, hjust = 0 ) +
  annotate("text", label = "Tussock Sedge Tundra", x = 581000, y = 7717100, size = 3, colour = veg_colours[1], vjust = 0, hjust = 0) +
  annotate("text", label = "Qikiqtaruk - \nHerschel Island", x = center_x, y = center_y - 600, size = 4, fontface = "italic", colour = 'black', vjust = 0, hjust = 0.5 ) +
  # Add rectangles
  annotate("rect", xmin = 582700 , xmax = 583100, ymin = 7720350 , ymax = 7720925, colour = "black", fill = "#00000000" ) +
  annotate("rect", xmin = 582900 , xmax = 583300, ymin = 7719800 , ymax = 7720200, colour = "black", fill = "#00000000" ) +
  annotate("rect", xmin = 580825 , xmax = 581250, ymin = 7720550 , ymax = 7721040, colour = "black", fill = "#00000000" ) +
  annotate("rect", xmin = 580500 , xmax = 580950, ymin = 7721060 , ymax = 7721550, colour = "black", fill = "#00000000" ) +
  # Beatuify plot
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 8, colour = "black"),
    axis.ticks = element_line(colour = "black", size = 1),
    legend.position = "none") +
  coord_fixed()

## Genereate insert plot of QHI
plot_width_qhi <- 21000
plot_height_qhi <- 17000
cropping_margin_qhi <- 2500
center_x_qhi <- 575000
center_y_qhi <- 7720000
x_min_qhi <- center_x_qhi - 0.5 * plot_width_qhi - cropping_margin_qhi
x_max_qhi <- center_x_qhi + 0.5 * plot_width_qhi + cropping_margin_qhi
y_min_qhi <- center_y_qhi - 0.5 * plot_height_qhi - cropping_margin_qhi
y_max_qhi <- center_y_qhi + 0.5 * plot_height_qhi + cropping_margin_qhi

qhi <- ggplot() +
  geom_polygon(data=filter(points_df, group == 58 | group == 19), #, group == 897 | group == 898), 
               aes(x=x, y=y),
               lwd = 1, 
               colour = "black", 
               fill = 'white') +
  scale_fill_manual(values = veg_colours) +
  scale_colour_manual(values = veg_colours) +
  scale_x_continuous(limits = c(x_min_qhi, x_max_qhi),
                     expand = expand_scale(add = c(-cropping_margin_qhi, -cropping_margin_qhi))) +
  scale_y_continuous(limits = c(y_min_qhi, y_max_qhi),
                     expand = expand_scale(add = c(-cropping_margin_qhi, -cropping_margin_qhi))) +
  annotate("rect", 
           xmin = x_min + cropping_margin,
           xmax = x_max - cropping_margin, 
           ymin = y_min + cropping_margin, 
           ymax = y_max - cropping_margin, 
           colour = "black",
           fill = "#00000030" ) +
  annotate("text", label = "Qikiqtaruk", x = center_x_qhi, y = center_y_qhi, size = 2.5, colour = 'black', vjust = -0.5, hjust = 0.6) +
  xlab("") +
  ylab("") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0,0,0,0), "cm")) +
  coord_fixed()

## Prepare final map plot
# Set farme boxes for insert plots
ps1_box_xmin <- 0.595
ps1_box_width <- 0.36
ps1_box_ymin <- 0.116
ps1_box_height <- 0.26

ps2_box_xmin <- ps1_box_xmin
ps2_box_width <- ps1_box_width
ps2_box_ymin <- 0.601
ps2_box_height <- ps1_box_height

ps3_box_xmin <- 0.1675
ps3_box_width <- ps1_box_width / 2 
ps3_box_ymin <- 0.116
ps3_box_height <- ps1_box_height

ps4_box_xmin <- 0.1675
ps4_box_width <- ps1_box_width / 2 
ps4_box_ymin <- 0.655
ps4_box_height <- ps1_box_height

box_buffer <- 0.0025

# Load Sentinel brick
sentinel_brick <- brick("figures/fig_1_ts_plots/L2A_QHI_20170728_cldsmskd_10m_brick.tif")
# 5 Generate final map figure ----
map_plot <- ggdraw() +
  #draw_plot(sites_plot + theme(legend.justification = "bottom"), 0, 0, 1, 1) +
  #draw_plot(qhi, x = 0.745, y = 0.82, height = 0.17, width = 0.21) +

map_plot <- ggdraw() +
  ## PS1
  # HER
  draw_image(paste0(ts_out_path, "/PS1_HER_2017_tsplot.png"),
             x = 0, y = 0) +
  draw_image(paste0(ts_out_path, "/PS1_HER_2016_tsplot.png"),
             x = 0.5, y = 0) +
  # KOM
  draw_image(paste0(ts_out_path, "/PS1_KOM_2017_tsplot.png"),
             x = 0, y = 0.5,
             scale = 0.2) +
  draw_image(paste0(ts_out_path, "/PS1_KOM_2016_tsplot.png"),
             x = 0.5, y = 05)

  annotate("rect", 
           xmin = ps1_box_xmin, xmax = ps1_box_xmin + ps1_box_width,
           ymin = ps1_box_ymin - 6 * box_buffer, ymax = ps1_box_ymin + (ps1_box_height / 2) - 2 * box_buffer, 
           colour = veg_colours[1], fill = NA) +
  annotate("rect", 
           xmin = ps1_box_xmin, xmax = ps1_box_xmin + ps1_box_width,
           ymin = ps1_box_ymin + (ps1_box_height / 2), ymax = ps1_box_ymin + ps1_box_height - 2 * box_buffer,
           colour = veg_colours[2], fill = NA) +

  ## PS2
  annotate("rect", 
         xmin = ps2_box_xmin, xmax = ps2_box_xmin + ps2_box_width,
         ymin = ps2_box_ymin , ymax = ps2_box_ymin + ps2_box_height, 
         colour = NA, fill = "white") +
  # HER
  draw_image(paste0(ts_out_path, "/PS2_HER_2017_tsplot.png"),
             x = 0.35, y = 0.28,
             scale = 0.2) +
  draw_image(paste0(ts_out_path, "/PS2_HER_2016_tsplot.png"),
             x = 0.20, y = 0.28,
             scale = 0.2) +
  # KOM
  draw_image(paste0(ts_out_path, "/PS2_KOM_2017_tsplot.png"),
             x = 0.35, y = 0.15,
             scale = 0.2) +
  draw_image(paste0(ts_out_path, "/PS2_KOM_2016_tsplot.png"),
             x = 0.20, y = 0.15,
             scale = 0.2)+
  annotate("rect", 
           xmin = ps2_box_xmin, xmax = ps2_box_xmin + ps2_box_width,
           ymin = ps2_box_ymin - 6 * box_buffer, ymax = ps2_box_ymin + (ps2_box_height / 2) - 2 * box_buffer, 
           colour = veg_colours[1], fill = NA) +
  annotate("rect", 
           xmin = ps2_box_xmin, xmax = ps2_box_xmin + ps2_box_width,
           ymin = ps2_box_ymin + (ps2_box_height / 2), ymax = ps2_box_ymin + ps2_box_height - 2 * box_buffer,
           colour = veg_colours[2], fill = NA) +
  
  # PS3
  annotate("rect", 
           xmin = ps3_box_xmin, xmax = ps3_box_xmin + ps3_box_width + 12* box_buffer,
           ymin = ps3_box_ymin , ymax = ps3_box_ymin + ps3_box_height, 
           colour = NA, fill = "white") +
  # HER
  draw_image(paste0(ts_out_path, "/PS3_HER_2017_tsplot.png"),
             x = -0.2275, y = -0.205,
             scale = 0.2) +
  # KOM
  draw_image(paste0(ts_out_path, "/PS3_KOM_2017_tsplot.png"),
             x = -0.2275, y = -0.335,
             scale = 0.2) +
  
  annotate("rect", 
           xmin = ps3_box_xmin, xmax = ps3_box_xmin + ps3_box_width + 12* box_buffer,
           ymin = ps3_box_ymin - 6 * box_buffer, ymax = ps3_box_ymin + (ps3_box_height / 2) - 2 * box_buffer, 
           colour = veg_colours[1], fill = NA) +
  annotate("rect", 
           xmin = ps3_box_xmin, xmax = ps3_box_xmin + ps3_box_width + 12* box_buffer,
           ymin = ps3_box_ymin + (ps3_box_height / 2), ymax = ps3_box_ymin + ps3_box_height - 2 * box_buffer,
           colour = veg_colours[2], fill = NA) +
  annotate("segment", x = 0.3, y = 0.375, xend = 0.4, yend = 0.495) +

 # PS 4
 annotate("rect", 
         xmin = ps4_box_xmin, xmax = ps4_box_xmin + ps4_box_width + 12* box_buffer,
         ymin = ps4_box_ymin , ymax = ps4_box_ymin + ps4_box_height, 
         colour = NA, fill = "white") +
  # HER
  draw_image(paste0(ts_out_path, "/PS4_HER_2017_tsplot.png"),
             x = -0.2275, y = 0.335,
             scale = 0.2) +
  # KOM
  draw_image(paste0(ts_out_path, "/PS4_KOM_2017_tsplot.png"),
             x = -0.2275, y = 0.205,
             scale = 0.2) +
  
  annotate("rect", 
           xmin = ps4_box_xmin, xmax = ps4_box_xmin + ps4_box_width + 12* box_buffer,
           ymin = ps4_box_ymin - 6 * box_buffer, ymax = ps4_box_ymin + (ps4_box_height / 2) - 2 * box_buffer, 
           colour = veg_colours[1], fill = NA) +
  annotate("rect", 
           xmin = ps4_box_xmin, xmax = ps4_box_xmin + ps4_box_width + 12* box_buffer,
           ymin = ps4_box_ymin + (ps4_box_height / 2), ymax = ps4_box_ymin + ps4_box_height - 2 * box_buffer,
           colour = veg_colours[2], fill = NA) +
  
  draw_image(paste0(ts_out_path, "/label_plot.png"),
             x = 0.16, y = 0.405,
             scale = 0.15) 

# Save plot file  
ggsave("figures/figure_1_satellite_drone_ts_map.png", 
       map_plot, 
       width = 15, 
       height = 20, 
       unit = "cm")
 

