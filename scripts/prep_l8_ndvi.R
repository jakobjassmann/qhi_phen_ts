# QHI Pehnology time-series Landsat 8 prep

# Dependencies
library(dplyr)
library(raster)
library(sf)

# Load scene lists
ls8_dirs <- c(list.dirs("/Volumes/BowheadRdge/phen_time_series/L8/2016"),
         list.dirs("/Volumes/BowheadRdge/phen_time_series/L8/2017"))
ls8_dirs <- ls8_dirs[!grepl(".*201[67]$", ls8_dirs)]
ls8 <- data.frame(scene_id = gsub("^.*/201[76]/(.*)", "\\1", ls8_dirs),
                  date = as.Date(gsub(".*(201[67])([0-9]{2})([0-9]{2}).*", 
                                      "\\1-\\2-\\3", 
                                      ls8_dirs)),
                  cloud_free = NA,
                  folder_path = ls8_dirs,
                  stringsAsFactors = F)
write.csv(ls8, 
          "data/auxillary/ls8_cloud_data.csv",
          row.names = F)

# Set scene boundaries for QA viewing
qhi_site_bounds <- as_Spatial(read_sf("data/site_boundaries/ps_site_bounds.shp"))

qhi_coords <- c(x = 575600, y = 7721888)
qhi_extent <- extent(c(qhi_coords[1] - 10000,
                       qhi_coords[1] + 10000,
                       qhi_coords[2] - 10000,
                       qhi_coords[2] + 10000))
qhi_extent <- as(qhi_extent, "SpatialPolygons")
crs(qhi_extent) <- crs(qhi_site_bounds)

# Set raster options
rasterOptions(progress = "text")


# Check for cloud contamination
lapply(ls8$scene_id, function(scene_id){
  # Load recent database file
  ls8_recent <- read.csv("data/auxillary/ls8_cloud_data.csv",
                         stringsAsFactors = F)
  # Load filenames
  file_list <- list.files(ls8$folder_path[ls8$scene_id == scene_id], 
                          full.names = T)
  
  # Status update
  cat(paste0("Scene date: ", ls8$date[ls8$scene_id == scene_id], "\n"))
  
  # Stack RGB bands
  ls8_rgb <- stack(file_list[grepl("band4", file_list)],
                   file_list[grepl("band3", file_list)],
                   file_list[grepl("band2", file_list)])
  
  # Convert geometries into LS8 projection
  qhi_extent <- spTransform(qhi_extent, crs(ls8_rgb))
  qhi_site_bounds <- spTransform(qhi_site_bounds, crs(ls8_rgb))
  
  # Crop and stretch
  cat("Cropping tile...\n")
  ls8_rgb <- crop(ls8_rgb, qhi_extent)
  cat("Stretching tile...\n")
  ls8_rgb <- stretch(ls8_rgb, minv = 0, maxv= 255,
                     minq = 0, maxq = 0.999)
  
  # Load qa band and create cloud masks (cloud masks are not reliable)
  # qa_band <- raster(file_list[grepl("pixel_qa", file_list)])
  # qa_band <- crop(qa_band, qhi_extent)
  # cloud_mask <- calc(qa_band, fun = function(x) bitwShiftR(bitwAnd(x, bitwShiftL(1,3)),3))
  # cloud_mask <- reclassify(cloud_mask, c(-Inf, 0, 1, 0.1, +Inf, NA))
  # cloud_shadows <- calc(qa_band, fun = function(x) bitwShiftR(bitwAnd(x, bitwShiftL(1,4)),4))
  # cloud_shadows <- reclassify(cloud_shadows, c(-Inf, 0, 1, 1, +Inf, NA))
  # 
  # # Apply masks
  # ls8_rgb <- mask(ls8_rgb, cloud_mask)
  # ls8_rgb <- mask(ls8_rgb, cloud_shadows)
  
  # Plot data
  par(mar = c(0, 0, 0, 0))
  plot(qhi_extent)
  plotRGB(ls8_rgb, add = T)
  plot(qhi_extent, add = T)
  plot(qhi_site_bounds, col = "red", add = T)
  
  # Prompt
  answer <- readline(prompt = "AoI cloud free and inlcuded in scene? [y/n]:")
  while(!sum(answer %in% c("y", "n"))) answer <- 
    readline(prompt = "Invalid answer! - AoI cloud free and inlcuded in scene? [y/n]:")
  
  # Assign answer and write out to file
  ls8_recent$cloud_free[ls8_recent$scene_id == scene_id] <- answer
  write.csv(ls8_recent, "data/auxillary/ls8_cloud_data.csv",
            row.names = F)
  gc()
  return(answer)
  })

# Load results from manual assessment
ls8 <- read.csv("data/auxillary/ls8_cloud_data.csv",
                stringsAsFactors = F)

# produce NDVI rasters
lapply(ls8$scene_id[ls8$cloud_free == "y"], function(scene_id){
  # Load filenames
  file_list <- list.files(ls8$folder_path[ls8$scene_id == scene_id],
                          full.names = T)
  
  # Status update
  cat(paste0("Processing scene date: ", ls8$date[ls8$scene_id == scene_id], "\n",
             scene_id, "\n"))
  
  # Load bands
  ls8_red <- raster(file_list[grepl("band4", file_list)])
  ls8_nir <- raster(file_list[grepl("band5", file_list)])

  # Crop to AoI
  qhi_extent <- spTransform(qhi_extent, crs(ls8_red))
  ls8_red <- crop(ls8_red, qhi_extent)
  ls8_nir <- crop(ls8_nir, qhi_extent)
  
  # Calculate NDVI
  ndvi <- (ls8_nir - ls8_red) / (ls8_nir + ls8_red)

  # Export raster
  writeRaster(ndvi, 
              paste0(ls8$folder_path[ls8$scene_id == scene_id],
                gsub(".*/(LC08_.*)_sr_band4.tif","/\\1_", 
                     file_list[grepl("band4", file_list)]),
                     "ndvi.tif"))
  # return northing
  return(NULL)
})

# Quick QC of NDVI rasters
ndvi_rasters <- list.files("/Volumes/BowheadRdge/phen_time_series/L8/", 
                           pattern = "ndvi.tif",
                           recursive = T, 
                           full.names = T)
lapply(ndvi_rasters,
       function(x){
         cat(paste0(x, "\n"))
         ndvi_raster <- raster(x)
         ndvi_raster <- reclassify(ndvi_raster, c(-Inf, -1, NA, 1, + Inf, NA))
         par(mar = c(0, 0, 0, 0))
         plot(ndvi_raster)
         readline(prompt = "Hit a key for next raster")
       })  

# Get summary stats per season
ls8 %>% filter(cloud_free == "y") %>% 
  mutate(year = format.Date(date, "%Y")) %>%
  group_by(year) %>%
  summarise(n = n())

## Calculate mean NDVI per scene and site ----

# Load site coords
site_boundaries <- read_sf("data/site_boundaries/ps_site_bounds.shp")

meta_data_ls8_with_mean <- bind_rows(lapply(
  ls8$scene_id[ls8$cloud_free == "y"],
  function(scene_id){
    meta_data_ls8 <- data.frame(
      flight_id = scene_id,
      site_veg = site_boundaries$site_veg,
      site_name = substr(site_boundaries$site_veg, 1,3),
      veg_type = substr(site_boundaries$site_veg, 5,7),
      file_path = list.files(ls8$folder_path[ls8$scene_id == scene_id], 
                             pattern = "_ndvi.tif",
                             full.names = T),
      object_name = scene_id,
      sensor = "Landsat8",
      date = ls8$date[ls8$scene_id == scene_id],
      mean_NDVI = NA
    )  
    ndvi_raster <- raster(list.files(ls8$folder_path[ls8$scene_id == scene_id], 
                                     pattern = "_ndvi.tif",
                                     full.names = T))
    meta_data_ls8$mean_NDVI <- raster::extract(ndvi_raster, 
                                               as_Spatial(site_boundaries),
                                       fun = mean,
                                       weighted = T)
    return(meta_data_ls8)
  }))

# Save mean
save(meta_data_ls8_with_mean, 
     file = "data/landsat8/meta_data_ls8_with_mean.Rda")
write.csv(meta_data_ls8_with_mean, 
          "data/landsat8/meta_data_ls8_with_mean.csv")
