# Phenology time-series landscape / plot scale estimates of NDVI
# Drone data, sentinel and MODIS means
# Jakob Assmann, j.assmann(at)ed.ac.uk

# load dependencies
library(raster)
library(rgdal)
library(rasterVis)
library(dplyr)
library(ggplot2)
library(viridisLite)
library(gridExtra)
library(cowplot)

# Globally important variables
site_boundaries <- read.csv("private/jassmann/phenology_time_series/ps_sent_site_bounds.csv")
script_path <- "private/jassmann/phenology_time_series/plot_scale_ts/"
sentinel_path <- "/Volumes/BowheadRdge/phen_time_series/sentinel_data/qhi_cld_free/"

# Globally important functions

# Funciton to create an site extent object form the sentine site boundary data frame
get_sent_extent <- function(site_veg_id, sent_boundaries) {
  extent_object <- extent(c(sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmin,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmax,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymin,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymax))
  return(extent_object)
}

# Prep globa meta-Data Frame
meta_data_global <- data.frame(flight_id = NA,
                               site_veg = NA,
                               site_name = NA,
                               veg_type = NA,
                               file_path= NA,
                               object_name = NA,
                               sensor = NA,
                               date = as.Date("2017-05-01"),
                               mean_NDVI= NA)


# Function to create time-series for 2016 plots
plot_ts_2016 <- function(site_name, veg_type){
  # Prep MODIS data
  site_coords_key <- read.csv("private/imyerssmith/Modis_plot_key.csv")
  MODIS6_data <- read.csv("private/imyerssmith/PS_centre_plots_v6_2000-2017.csv")
  
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
  
  ##  drone files
  # Set folder path
  folder_path <- paste0("/Volumes/BowheadRdge/phen_time_series/final_outputs/2016/", 
                       site_name, "_", veg_type, "/output/ndvi_maps")
  # optain list of files
  file_names <- list.files(folder_path, pattern = "*.tif")
  
  # create meta-data df
  meta_data <- data.frame(flight_id = substr(file_names, 1, 13),
                          site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
                          site_name = rep(site_name, length(file_names)),
                          veg_type = rep(veg_type, length(file_names)),
                          file_path = paste0(folder_path, "/", file_names),
                          object_name = gsub(".tif", "", file_names))
  meta_data$file_path <- as.character(meta_data$file_path)
  meta_data$object_name <- as.character(meta_data$object_name)
  
  # kick out all multiples for 3 August 2016 for PS1 HER and KOM
  if(site_name == "PS1" & veg_type == "HER"){
    meta_data <- meta_data[1:6,]
  } else if(site_name == "PS1" & veg_type == "KOM"){
    meta_data <- meta_data[1:4,]
  }
  
  ## get list of sentinel files
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
  
  # merge two meta data dfs
  sensor_id <- c(rep("sentinel", nrow(sentinel_meta)), rep("drone", nrow(meta_data)))
  meta_data <- rbind(sentinel_meta, meta_data)
  meta_data$sensor <- sensor_id
  # add date collumn
  meta_data$date <- as.Date(gsub("*.*(2016[0-9][0-9][0-9][0-9]).*", 
                                 "\\1", meta_data$flight_id), format = "%Y%m%d")
  
  ## load files
  list2env(
    lapply(
      setNames(meta_data$file_path, 
             make.names(gsub(".*/", "", 
                             gsub(".tif$", "", meta_data$file_path)))),
      raster), 
    envir = .GlobalEnv)
  

  # get site boundaries (as extent object)
  site_bounds <- get_sent_extent(paste0(site_name, "_", veg_type), site_boundaries)
  
  ## crop files in slick one liner
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
  meta_data <- rbind(MODIS_meta[,match(names(meta_data), names(MODIS_meta))], meta_data)
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2016", "%d/%m/%Y") & 
             date <= as.Date("30/09/2016", "%d/%m/%Y"))
  # plot curve of NDVI means big font version
  plot_scale = viridis(5)[c(5,2,3,4)]
  # No nocalib data for Ps2 HER, adjust colour scale and KOM is messed up correct!
  if((site_name == "PS2" & veg_type == "HER")){
    plot_scale = viridis(5)[c(5,2,3)]
  } else if(site_name == "PS2" & veg_type == "KOM"){
    plot_scale = viridis(5)[c(5,2,4,3)]
  }
  # Set y axis label
  #if(site_name == "PS1"){
    y_label <- "NDVI"
  #} else {
  #  y_label = " \n"
  #}
  
  
  # Adjust visibility of axis labels according to plot positions
  x.axis.text.colour <- "black"
  y.axis.text.colour <- "black"
  if(veg_type == "HER") {
     x.axis.text.colour <- "white"
  }
  
  ts_plot <- ggplot(data = meta_data, 
                    mapping = aes(x = date, 
                                  y = mean_NDVI, 
                                  color = sensor, 
                                  shape = sensor, 
                                  size = sensor, 
                                  fill = sensor), 
                    aes(x = date, y = mean_NDVI), 
                    inherit.aes = F) +
                    geom_point() +
                    geom_smooth(mapping = aes(x = date, y= mean_NDVI), 
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
    
    scale_size_manual(values = c(rep(5, 4))) +
    scale_colour_manual(values = rep("black",4)) +
    scale_shape_manual(values = c(21, 21, 21, 21)) +
    scale_fill_manual(values = plot_scale) +
    scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1.05), 
                       breaks = seq(0,1,0.25), 
                       labels = paste0(format(seq(0,1,0.25), digits = 2), " ") ) +
    scale_x_date(expand = c(0,0), 
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
                            "S")
                 ) +
    annotate(geom="segment", y = -0.15, yend = -0.12, size = 1.2,
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
    annotate(geom="text", 
             x = as.Date("2016-07-15"), 
             y = 0.98, 
             label = "2016", 
             size =  10, 
             fontface = "bold"
             ) +
    ylab(y_label) +
    xlab("") +
    theme_bw() +
    theme(panel.border = element_blank(),
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
  rm(list = c(meta_data[meta_data$sensor != "MODIS",]$object_name, paste0(meta_data[meta_data$sensor != "MODIS",]$object_name, "_cropped")), envir = .GlobalEnv)
  # Export meta data to global meta_data df
  meta_data_global <<- rbind(meta_data_global, meta_data)

  
  return(ts_plot)
}

## Function to plot 2017 plots
#####
plot_ts_2017 <- function(site_name, veg_type){
  # Prep MODIS data
  site_coords_key <- read.csv("private/imyerssmith/Modis_plot_key.csv")
  MODIS6_data <- read.csv("private/imyerssmith/PS_centre_plots_v6_2000-2017.csv")
  
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
  
  ##  drone files
  # Set folder path
  folder_path <- paste0("/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/", 
                        site_name, "_", veg_type, "/output/ndvi_maps")
  # optain list of files
  file_names <- list.files(folder_path, pattern = "*.tif")
  
  # create meta-data df
  meta_data <- data.frame(flight_id = gsub("-", "", substr(file_names, 1, 18)),
                          site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
                          site_name = rep(site_name, length(file_names)),
                          veg_type = rep(veg_type, length(file_names)),
                          file_path = paste0(folder_path, "/", file_names),
                          object_name = gsub("-", "", gsub(".tif", "", file_names)))
  meta_data$file_path <- as.character(meta_data$file_path)
  meta_data$object_name <- as.character(meta_data$object_name)
  
  
  ## get list of sentinel files
  sentinel_path <- paste0(sentinel_path, "/2017/")
  sentinel_files <- list.files(sentinel_path, pattern = "ndvi.tif")
  
  # add meta_data
  sentinel_meta <- data.frame(flight_id = gsub("*.*(2017[0-9][0-9][0-9][0-9]).*", 
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
  sensor_id <- c(sentinel_id, rep("drone", nrow(meta_data)))
  meta_data <- rbind(sentinel_meta, meta_data)
  meta_data$sensor <- sensor_id
  # add date collumn
  meta_data$date <- as.Date(gsub("*.*(2017[0-9][0-9][0-9][0-9]).*", 
                                 "\\1", meta_data$flight_id), format = "%Y%m%d")
  
  ## load files
  list2env(
    lapply(
      setNames(meta_data$file_path, 
               make.names(meta_data$object_name)),
      raster), 
    envir = .GlobalEnv)
  
  
  # get site boundaries (as extent object)
  site_bounds <- get_sent_extent(paste0(site_name, "_", veg_type), site_boundaries)
  
  ## crop files in slick one liner
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
  # if (site_name == "PS1") {
  #   meta_data[meta_data$date == as.Date("2016-07-30"),]$sensor <- "drone_nocalib"
  #   no_calib_true <- NA
  # } else if (site_name == "PS2" & veg_type == "KOM") {
  #   meta_data[meta_data$date == as.Date("2016-06-30"),]$sensor <- "drone_nocalib"
  #   no_calib_true <- NA
  # } else {
  #   no_calib_true <- 0
  # }
  
  # Merge with MODIS data
  meta_data <- rbind(MODIS_meta[,match(names(meta_data), names(MODIS_meta))], meta_data)
  
  # Filter from May to September
  meta_data <- meta_data %>%  
    filter(date >= as.Date("01/05/2017", "%d/%m/%Y") & 
             date <= as.Date("30/09/2017", "%d/%m/%Y"))
  # plot curve of NDVI means big font version
  plot_scale = viridis(5)[c(5,2,1,3)]
  
  # Set y axis label and colour colour
  if(site_name == "PS1" | site_name == "PS2"){
    y_label <- "NDVI"
    y_label_colour <- "white"
  } else {
    y_label = "NDVI"
    y_label_colour <- "black"
  }
  

  # # Adjust visibility of axis labels according to plot positions
    x.axis.text.colour <- "black"
    y.axis.text.colour <- "black"
  if(site_name == "PS1" | site_name == "PS2") {
    y.axis.text.colour <- "white"
  }
  if(veg_type == "HER") {
    x.axis.text.colour <- "white"
  }

  
  ts_plot <- ggplot(data = meta_data, 
                    mapping = aes(x = date, 
                                  y = mean_NDVI, 
                                  color = sensor, 
                                  shape = sensor, 
                                  size = sensor, 
                                  fill = sensor), 
                    aes(x = date, y = mean_NDVI), 
                    inherit.aes = F) +
    geom_point() +
    geom_smooth(mapping = aes(x = date, y= mean_NDVI), 
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
    scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1.05), 
                       breaks = seq(0,1,0.25), 
                       labels = paste0(format(seq(0,1,0.25), digits = 2), " ") ) +
    scale_x_date(expand = c(0,0), 
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
                            "S")
    ) +
    annotate(geom="segment", y = -0.15, yend = -0.12, size = 1.2,
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
    annotate(geom="text", 
             x = as.Date("2017-07-15"), 
             y = 0.98, 
             label = "2017", 
             size =  10, 
             fontface = "bold"
    ) +
    ylab(y_label) +
    xlab("") +
    theme_bw() +
    theme(panel.border = element_blank(),
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
  ts_plot
  rm(list = c(meta_data[meta_data$sensor != "MODIS",]$object_name, paste0(meta_data[meta_data$sensor != "MODIS",]$object_name, "_cropped")), envir = .GlobalEnv)
    meta_data_global <<- rbind(meta_data_global, meta_data)
  return(ts_plot)
}

#####
# Create Plots 
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
#####
# Save plots
save_plot(paste0(script_path,"PS1_HER_2016_tsplot.png"), 
          PS1_HER_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS1_KOM_2016_tsplot.png"), 
          PS1_KOM_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS2_HER_2016_tsplot.png"), 
          PS2_HER_tsplot_2016, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS2_KOM_2016_tsplot.png"), 
          PS2_KOM_tsplot_2016, base_aspect_ratio = 1)

save_plot(paste0(script_path,"PS1_HER_2017_tsplot.png"), 
          PS1_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS1_KOM_2017_tsplot.png"), 
          PS1_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS2_HER_2017_tsplot.png"), 
          PS2_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS2_KOM_2017_tsplot.png"), 
          PS2_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS3_HER_2017_tsplot.png"), 
          PS3_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS3_KOM_2017_tsplot.png"), 
          PS3_KOM_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS4_HER_2017_tsplot.png"), 
          PS4_HER_tsplot_2017, base_aspect_ratio = 1)
save_plot(paste0(script_path,"PS4_KOM_2017_tsplot.png"), 
          PS4_KOM_tsplot_2017, base_aspect_ratio = 1)

# Make Label Plot
plot_scale = viridis(5)[c(3,4,2,1,5)]
label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(0, 5)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-06-20"), as.Date("2017-10-30"))) +
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 4.5, label = "MODIS", colour = "black" , fontface = "bold", size = 13, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 4.5, shape = 21, colour = "black", fill = plot_scale[5], size = 11) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 3.5, label = "Sentinel 2A", colour = plot_scale[3] , fontface = "bold",size = 13,  hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 3.5, shape = 21, colour = "black", fill = plot_scale[3], size = 11) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 2.5, label = "Sentinel 2B", colour = plot_scale[4] , fontface = "bold",size = 13,  hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 2.5, shape = 21, colour = "black", fill = plot_scale[4], size = 11) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 1.5, label = "Drone", colour = plot_scale[1] , fontface = "bold", size = 13, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 1.5, shape = 21, colour = "black", fill = plot_scale[1], size = 11) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.5, label = "Drone, not calibrated", colour = plot_scale[2] , fontface = "bold", size = 13, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.5, shape = 21, colour = "black", fill = plot_scale[2], size = 11) +
  
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

save_plot(paste0(script_path,"label_plot.png"), label_plot,base_aspect_ratio = 1.6)
# Blank dummy plot
#####
blank_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#####

# Make label plots
##### 
# PS1 label
PS1_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("05/08/2017", "%d/%m/%Y"), y = -0.1, label = "Collinson Head (Site 1)", colour = "black" , fontface = "bold", size = 12, hjust = 0.5) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# PS2 label
PS2_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("05/08/2017", "%d/%m/%Y"), y = -0.1, label = "Bowhead Ridge (Site 2)", colour = "black" , fontface = "bold", size = 12, hjust = 0.5) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# PS3 label
plot_scale = viridis(5)[c(3,4,2,1,5)]
PS3_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("05/08/2017", "%d/%m/%Y"), y = -0.1, label = "Hawk Valley (Site 3)", colour = "black" , fontface = "bold", size = 12, hjust = 0.5) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.61, label = "MODIS", colour = "black" , fontface = "bold", size = 10, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.605, shape = 21, colour = "black", fill = plot_scale[5], size = 8) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.53, label = "Sentinel 2A", colour = plot_scale[3] , fontface = "bold",size = 10,  hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.525, shape = 21, colour = "black", fill = plot_scale[3], size = 8) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.45, label = "Sentinel 2B", colour = plot_scale[4] , fontface = "bold",size = 10,  hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.445, shape = 21, colour = "black", fill = plot_scale[4], size = 8) +

  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.37, label = "Drone", colour = plot_scale[1] , fontface = "bold", size = 10, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.365, shape = 21, colour = "black", fill = plot_scale[1], size = 8) +
  
  annotate("text", x = as.Date("05/07/2017", "%d/%m/%Y"), y = 0.29, label = "Drone, not calibrated", colour = plot_scale[2] , fontface = "bold", size = 10, hjust = 0) +
  annotate("point", x = as.Date("26/06/2017", "%d/%m/%Y"), y = 0.285, shape = 21, colour = "black", fill = plot_scale[2], size = 8) +
  
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# PS4 label
PS4_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("05/08/2017", "%d/%m/%Y"), y = -0.1, label = "Hawk Ridge (Site 4)", colour = "black" , fontface = "bold", size = 12, hjust = 0.5) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Herschel label
HER_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("30/07/2017", "%d/%m/%Y"), y = 0.425, label = "Herschel\nVegetation Type", colour = "#1E9148FF" , fontface = "bold", size = 12, hjust = 0.5, vjust = 0.5) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
# Komakuk label
KOM_label_plot <- ggplot() + 
  scale_y_continuous(expand = c(0,0), limits = c(-0.15, 1)) +
  scale_x_date(expand = c(0,0), 
               limits = c(as.Date("2017-05-01"), as.Date("2017-09-30"))) +
  annotate("text", x = as.Date("30/07/2017", "%d/%m/%Y"), y = 0.425, label = "Komakuk\nVegetation Type", colour = "#1E5C91FF" , fontface = "bold", size = 12, hjust = 0.5, vjust = 0.5) +
  xlab("") +
  ylab(" \n") + 
  theme_bw() + 
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
#####

# Create Grid
grid_matrix <- rbind(c( 1,  6, 11, 16, 21), 
                     c( 2,  7, 12, 17, 22),
                     c( 3,  8, 13, 18, 23),
                     c( 4,  9, 14, 19, 24),
                     c( 5, 10, 15, 20, 25))

## Make actual plot
png(paste0(script_path, "/phenology_time_series.png"), width = 40, height = 40, units = "in", res = 300)
marrangeGrob(
  # NB rows and collums swapped for layout no clue why this works that way?
  grobs = list(
         blank_plot,      HER_label_plot,      KOM_label_plot,      HER_label_plot,      KOM_label_plot,
     PS1_label_plot, PS1_HER_tsplot_2016, PS1_KOM_tsplot_2016, PS1_HER_tsplot_2017, PS1_KOM_tsplot_2017,
     PS2_label_plot, PS2_HER_tsplot_2016, PS2_KOM_tsplot_2016, PS2_HER_tsplot_2017, PS2_KOM_tsplot_2017,
         blank_plot,          blank_plot,      PS3_label_plot, PS3_HER_tsplot_2017, PS3_KOM_tsplot_2017,
         blank_plot,          blank_plot,      PS4_label_plot, PS4_HER_tsplot_2017, PS4_KOM_tsplot_2017),
  ncol = 5, 
  nrow = 5,
  top = "",
  layou_matrix = grid_matrix)
dev.off()

# Calculate mean sensor offset between drone and satellites in July.
meta_data_global$month <- format.Date(meta_data_global$date, "%b")
meta_data_global$year <- format.Date(meta_data_global$date, "%Y")
meta_data_global$sensor[meta_data_global$sensor == "Sentinel 2A"] <- "Sentinel"
meta_data_global$sensor[meta_data_global$sensor == "Sentinel 2B"] <- "Sentinel"
meta_data_global$sensor[meta_data_global$sensor == "sentinel"]<- "Sentinel"
july_meta <- meta_data_global %>% 
  filter(!(site_veg == "PS4_HER" & date == as.Date("2017-07-17")) ) %>%
  group_by(site_veg, year, sensor) %>% 
  filter(month == "Jul" & sensor != "drone_nocalib") %>% 
  summarise(mean_NDVI = mean(mean_NDVI))
library(tidyr)
july_meta <- spread(july_meta, sensor, mean_NDVI)
july_meta$MODIS_drone_diff <- july_meta$MODIS - july_meta$drone
july_meta$Sentinel_drone_diff <- july_meta$Sentinel - july_meta$drone
july_meta$MODIS_Sentinel_diff <- july_meta$MODIS - july_meta$Sentinel
july_meta_stats <- july_meta %>% ungroup() %>% summarise(mean_MODIS_drone_diff = mean(MODIS_drone_diff),
                                                         mean_Sentinel_drone_diff = mean(Sentinel_drone_diff),
                                                         mean_MODIS_Sentinel_diff = mean(MODIS_Sentinel_diff),
                                                         sd_MODIS_drone_diff = sd(MODIS_drone_diff),
                                                         sd_Sentinel_drone_diff = sd(Sentinel_drone_diff),
                                                         sd_MODIS_Sentinel_diff = sd(MODIS_Sentinel_diff)) %>% 
  as.data.frame() %>% round(3)
