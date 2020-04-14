# Phenology Time-Series Data Base
# This script is for creating a dataset with an overview of the sentinel and drone data this is available
# Jakob Assmann j.assmann@ed.ac.uk 3 October 2018

# load dependencies
library(dplyr)

#### 1) Drone meta data  ----
# Function to load site_veg meta data
gather_meta_drone <- function(site_name, veg_type){
  # First of 2016
  # Set folder paths
  folder_path_ndvi_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2016/", 
    site_name, "_", 
    veg_type,
    "/output/ndvi_maps")
  folder_path_red_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2016/", 
    site_name, "_", 
    veg_type, 
    "/output/reflec_maps_red")
  folder_path_nir_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2016/", 
    site_name, 
    "_", 
    veg_type, 
    "/output/reflec_maps_nir")
  # optain list of files
  file_names_ndvi <- list.files(folder_path_ndvi_maps, pattern = "*.tif")
  file_names_red <- list.files(folder_path_red_maps, pattern = "*.tif")
  file_names_nir <- list.files(folder_path_nir_maps, pattern = "*.tif")
  # connetate file names vector
  file_names <- c(file_names_ndvi, file_names_red, file_names_nir)
  
  # create meta-data df
  meta_data_site_2016 <- data.frame(
    flight_id = substr(file_names, 1, 13),
    site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
    site_name = rep(site_name, length(file_names)),
    veg_type = rep(veg_type, length(file_names)),
    band = c(rep("NDVI", length(file_names_ndvi)), 
             rep("RED", length(file_names_red)), 
             rep("NIR", length(file_names_nir))),
    file_path = paste0(c(rep(folder_path_ndvi_maps, length(file_names_ndvi)), 
                         rep(folder_path_red_maps, length(file_names_red)), 
                         rep(folder_path_nir_maps, length(file_names_nir))), 
                       "/", file_names),
    object_name = gsub(".tif", "", file_names))
  # Now 2017
  # Set folder path
  folder_path_ndvi_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/", 
    site_name, 
    "_", 
    veg_type, 
    "/output/ndvi_maps")
  folder_path_red_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/", 
    site_name, 
    "_",
    veg_type,
    "/output/reflec_maps_red")
  folder_path_nir_maps <- paste0(
    "/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/", 
    site_name, 
    "_", 
    veg_type, 
    "/output/reflec_maps_nir")
  # optain list of files
  file_names_ndvi <- list.files(folder_path_ndvi_maps, pattern = "*.tif")
  file_names_red <- list.files(folder_path_red_maps, pattern = "*.tif")
  file_names_nir <- list.files(folder_path_nir_maps, pattern = "*.tif")
  # connetate file names vector
  file_names <- c(file_names_ndvi, file_names_red, file_names_nir)
  
  # create meta-data df
  meta_data_site_2017 <- data.frame(
    flight_id = gsub("-", "", substr(file_names, 1, 18)),
    site_veg = rep(paste0(site_name, "_", veg_type), length(file_names)),
    site_name = rep(site_name, length(file_names)),
    veg_type = rep(veg_type, length(file_names)),
    band = c(rep("NDVI", length(file_names_ndvi)), 
             rep("RED", length(file_names_red)), 
             rep("NIR", length(file_names_nir))),
    file_path = paste0(c(rep(folder_path_ndvi_maps, length(file_names_ndvi)), 
                         rep(folder_path_red_maps, length(file_names_red)), 
                         rep(folder_path_nir_maps, length(file_names_nir))),  
                       "/", file_names),
    object_name = gsub("-", "", gsub(".tif", "", file_names)))
  
  # Merge both years
  meta_data_site <- rbind(meta_data_site_2016, meta_data_site_2017)
  
  # Clear up factors into character vectors
  meta_data_site$file_path <- as.character(meta_data_site$file_path)
  meta_data_site$object_name <- as.character(meta_data_site$object_name)
  
  # check whether global meta_data exists?
  if(exists("meta_data_drone")) {
    # Yes? Append new site meta_data
    meta_data_drone <<- rbind(meta_data_drone, meta_data_site)
  } else {
    # No? Create global data frame
    meta_data_drone <<- meta_data_site
  }
}

gather_meta_drone("PS1", "HER")
gather_meta_drone("PS1", "KOM")
gather_meta_drone("PS2", "HER")
gather_meta_drone("PS2", "KOM")
gather_meta_drone("PS3", "HER")
gather_meta_drone("PS3", "KOM")
gather_meta_drone("PS4", "HER")
gather_meta_drone("PS4", "KOM")


#### 2) Sentinel meta data -----

# Start with 2016
sentinel_path <- 
  "/Volumes/BowheadRdge/phen_time_series/sentinel_data/qhi_cld_free/2016"
## get list of sentinel files
sentinel_files <- list.files(sentinel_path, pattern = "brick.tif")

# create meta_data for 2016
sentinel_meta_2016 <- data.frame(
  flight_id = gsub("*.*(2016[0-9][0-9][0-9][0-9]).*", 
                   "\\1", sentinel_files),
  site_veg = rep(NA, length(sentinel_files)),
  site_name = rep(NA, length(sentinel_files)),
  veg_type = rep(NA, length(sentinel_files)),
  band = "BRICK",
  file_path = paste0(sentinel_path, "/", sentinel_files),
  object_name = gsub(".tif", "", sentinel_files))

## Next 2017
sentinel_path <- 
  "/Volumes/BowheadRdge/phen_time_series/sentinel_data/qhi_cld_free/2017"
## get list of sentinel files
sentinel_files <- list.files(sentinel_path, pattern = "brick.tif")

# create meta_data for 2017
sentinel_meta_2017 <- data.frame(
  flight_id = gsub("*.*(2017[0-9][0-9][0-9][0-9]).*", 
                   "\\1", sentinel_files),
  site_veg = rep(NA, length(sentinel_files)),
  site_name = rep(NA, length(sentinel_files)),
  veg_type = rep(NA, length(sentinel_files)),
  band = "BRICK",
  file_path = paste0(sentinel_path, "/", sentinel_files),
  object_name = gsub(".tif", "", sentinel_files))

# Merge data for 2016 and 2017
sentinel_meta <- rbind(sentinel_meta_2016, sentinel_meta_2017)
sentinel_meta$file_path <- as.character(sentinel_meta$file_path)
sentinel_meta$object_name <- as.character(sentinel_meta$object_name)

# Set sensor info
sensor_id <- c(rep("Drone", nrow(meta_data_drone)),
  rep("Sentinel 2A", nrow(sentinel_meta_2016)),
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

##### 3) Combine meta data ----
meta_data <- rbind(meta_data_drone, sentinel_meta)
# add date collumn
meta_data$date <- as.Date(gsub("*.*(201[6-7][0-9][0-9][0-9][0-9]).*", 
                               "\\1", meta_data$flight_id), format = "%Y%m%d")
# some flight ids from 2017 do not follow the naming convention, 
# making the date conversiosn not work. Luckily there is a quick fix:
meta_data[is.na(meta_data$date),]$date <- 
  as.Date(
    gsub("*.*(201[6-7]_[0-9][0-9]_[0-9][0-9]).*",
         "\\1", 
         meta_data[is.na(meta_data$date),]$flight_id), format = "%Y_%m_%d")
# add sensor_id collumn
meta_data$sensor_id <- sensor_id

# Change order of collumns
meta_data <- meta_data[,c(1,8,9,2,3,4,5,6,7)]

# Save dataframe
save(meta_data, file = "data/meta_data.Rda")
