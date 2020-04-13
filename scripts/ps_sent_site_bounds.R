#### Determine Site Boundaries based on sentinel grid.


## load dependencies
library(raster)
library(dplyr)
library(rgdal)

# Define output paths
site_boundaries <- read.csv("private/jassmann/phenology_time_series/ps_site_boundaries.csv")
script_path <- "private/jassmann/phenology_time_series/"

# load drone outputs
PS1_HER <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS1_HER/output/ndvi_maps/JA20160729_01_ndvi.tif")  
PS1_KOM <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS1_KOM/output/ndvi_maps/cropped/JA20160729_02_ndvi_cropped.tif")
PS2_HER <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS2_HER/output/ndvi_maps/JA20160730_05_ndvi.tif")
PS2_KOM <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS2_KOM/output/ndvi_maps/cropped/JA20160730_06_ndvi_cropped.tif")
PS4_HER <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS3_HER/output/ndvi_maps/JA20160727_01_nocalib_ndvi.tif")
PS4_KOM <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS3_KOM/output/ndvi_maps/JA20160727_02_nocalib_ndvi.tif")
PS4_HER <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS4_HER/output/ndvi_maps/JA20160724_06_ndvi.tif")
PS4_KOM <- raster("/Volumes/BowheadRdge/phen_time_series/2016/PS4_KOM/output/ndvi_maps/JA20160724_05_ndvi.tif")

# load sentinel ndvi map
sentinel <- raster("/Volumes/BowheadRdge/phen_time_series/sentinel_data/final_outputs/qhi_ndvi/S2A_USER_MSI_L2A_QHI_20160714_cldsmskd_10m_ndvi.tif")

# grab utm zone 7 prj4 data
epsg <- make_EPSG()
utm_z7N <- as.character(epsg %>% filter(code == 32607) %>% select(prj4))
rm(epsg)

# prepare a function that creates a polygon from the four corner coordinates
bound_polygon <- function(s_name, v_type, site_bnds = site_boundaries){
  # collect corner coordinates
  corners <- site_bnds %>% filter(site_name == s_name & veg_type == v_type) %>% select(long, lat) # NB LONG LAT ORDER IS REALLY IMPORTANT
  crdref <- CRS('+proj=longlat +datum=WGS84')
  corn_polygon <- spPolygons(as.matrix(corners), crs=crdref) # NB as.matrix needed as spPolygons can't handle dfs
  # now transform into UTM Zone7
  corn_polygon <- spTransform(corn_polygon, utm_z7N)
  return(corn_polygon)
}

# create spatial polygons for all site veg type combinations
PS1_HER_bounds <- bound_polygon("PS1", "HER")
PS1_KOM_bounds <- bound_polygon("PS1", "KOM")
PS2_HER_bounds <- bound_polygon("PS2", "HER")
PS2_KOM_bounds <- bound_polygon("PS2", "KOM")
PS4_HER_bounds <- bound_polygon("PS3", "HER")
PS4_KOM_bounds <- bound_polygon("PS3", "KOM")
PS4_HER_bounds <- bound_polygon("PS4", "HER")
PS4_KOM_bounds <- bound_polygon("PS4", "KOM")

# crop rasters
PS1_HER <- crop(PS1_HER, PS1_HER_bounds)
PS1_KOM <- crop(PS1_KOM, PS1_KOM_bounds)
PS2_HER <- crop(PS2_HER, PS2_HER_bounds)
PS2_KOM <- crop(PS2_KOM, PS2_KOM_bounds)
PS4_HER <- crop(PS4_HER, PS3_HER_bounds)
PS4_KOM <- crop(PS4_KOM, PS3_KOM_bounds)
PS4_HER <- crop(PS4_HER, PS4_HER_bounds)
PS4_KOM <- crop(PS4_KOM, PS4_KOM_bounds)

PS1_HER_sent <- crop(sentinel, PS1_HER_bounds)
PS1_KOM_sent <- crop(sentinel, PS1_KOM_bounds)
PS2_HER_sent <- crop(sentinel, PS2_HER_bounds)
PS2_KOM_sent <- crop(sentinel, PS2_KOM_bounds)
PS4_HER_sent <- crop(sentinel, PS3_HER_bounds)
PS4_KOM_sent <- crop(sentinel, PS3_KOM_bounds)
PS4_HER_sent <- crop(sentinel, PS4_HER_bounds)
PS4_KOM_sent <- crop(sentinel, PS4_KOM_bounds)

# prep dataframe for new extends based on sentinel grid
new_bounds <- data.frame(site_veg = c("PS1_HER", 
                                      "PS1_KOM",
                                      "PS2_HER",
                                      "PS2_KOM",
                                      "PS3_HER",
                                      "PS3_KOM",
                                      "PS4_HER",
                                      "PS4_KOM"),
                         site = c("PS1",
                                  "PS1",
                                  "PS2",
                                  "PS2",
                                  "PS3",
                                  "PS3",
                                  "PS4",
                                  "PS4"),
                         veg_type = c("HER",
                                      "KOM",
                                      "HER",
                                      "KOM",
                                      "HER",
                                      "KOM", 
                                      "HER", 
                                      "KOM"),
                         xmin = NA,
                         xmax = NA,
                         ymin = NA,
                         ymax = NA)

#####
# plot rasters, extends and shap of plots
# then define bounds

# PS1 HER first
plot(PS1_HER_sent)
plot(extent(PS1_HER), col = "blue", add = T)
plot(PS1_HER_bounds, col = "red", add = T)
# 1 cell to big in each dimension, cut off 1 row of cells to north and east.
new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$xmin <- PS1_HER_sent@extent@xmin + 0
new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$xmax <- PS1_HER_sent@extent@xmax - 10
new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$ymin <- PS1_HER_sent@extent@ymin + 0
new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$ymax <- PS1_HER_sent@extent@ymax - 10
# check output
PS1_HER_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$xmin,
                             new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$xmax,
                             new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$ymin,
                             new_bounds[which(new_bounds$site_veg == "PS1_HER"),]$ymax))
plot(PS1_HER_sent)
plot(PS1_HER_bounds_new, col = "red", add = T)
# nice one!

# PS1 KOM 
plot(PS1_KOM_sent)
plot(extent(PS1_KOM), col = "blue", add = T)
plot(PS1_KOM_bounds, col = "red", add = T)
# 1 cell to big in each dimension, cut off 1 row of cells to north and east.
new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$xmin <- PS1_KOM_sent@extent@xmin + 10
new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$xmax <- PS1_KOM_sent@extent@xmax + 0
new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$ymin <- PS1_KOM_sent@extent@ymin + 10
new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$ymax <- PS1_KOM_sent@extent@ymax - 10
# check output
PS1_KOM_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS1_KOM"),]$ymax))
plot(PS1_KOM_sent)
plot(PS1_KOM_bounds_new, col = "red", add = T)
# nice one!

# PS2 HER 
plot(PS2_HER_sent)
plot(extent(PS2_HER), col = "blue", add = T)
plot(PS2_HER_bounds, col = "red", add = T)
# 1 cell to big in each dimension, cut off 1 row of cells to north and east.
new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$xmin <- PS2_HER_sent@extent@xmin + 20
new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$xmax <- PS2_HER_sent@extent@xmax - 10
new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$ymin <- PS2_HER_sent@extent@ymin + 0
new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$ymax <- PS2_HER_sent@extent@ymax - 10
# check output
PS2_HER_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS2_HER"),]$ymax))
plot(PS2_HER_sent)
plot(PS2_HER_bounds_new, col = "red", add = T)
# nice one!

# PS2 KOM 
plot(PS2_KOM_sent)
plot(extent(PS2_KOM), col = "blue", add = T)
plot(PS2_KOM_bounds, col = "red", add = T)
# remove 1 row in east and one in south
new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$xmin <- PS2_KOM_sent@extent@xmin + 0
new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$xmax <- PS2_KOM_sent@extent@xmax - 10
new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$ymin <- PS2_KOM_sent@extent@ymin + 10
new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$ymax <- PS2_KOM_sent@extent@ymax + 0
# check output
PS2_KOM_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS2_KOM"),]$ymax))
plot(PS2_KOM_sent)
plot(PS2_KOM_bounds_new, col = "red", add = T)
# nice one!

# PS3 HER 
plot(PS3_HER_sent)
plot(extent(PS3_HER), col = "blue", add = T)
plot(PS3_HER_bounds, col = "red", add = T)
# one cell less in the west and one cell less in the south
new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$xmin <- PS3_HER_sent@extent@xmin + 0
new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$xmax <- PS3_HER_sent@extent@xmax - 10
new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$ymin <- PS3_HER_sent@extent@ymin + 10
new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$ymax <- PS3_HER_sent@extent@ymax + 0
# check output
PS3_HER_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS3_HER"),]$ymax))
plot(PS3_HER_sent)
plot(PS3_HER_bounds_new, col = "red", add = T)
# nice one!

# PS3 KOM 
plot(PS3_KOM_sent)
plot(extent(PS3_KOM), col = "blue", add = T)
plot(PS3_KOM_bounds, col = "red", add = T)
# one cell less in the south
new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$xmin <- PS3_KOM_sent@extent@xmin + 0
new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$xmax <- PS3_KOM_sent@extent@xmax - 0
new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$ymin <- PS3_KOM_sent@extent@ymin + 10
new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$ymax <- PS3_KOM_sent@extent@ymax + 0
# check output
PS3_KOM_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS3_KOM"),]$ymax))
plot(PS3_KOM_sent)
plot(PS3_KOM_bounds_new, col = "red", add = T)
# nice one!

# PS4 HER 
plot(PS4_HER_sent)
plot(extent(PS4_HER), col = "blue", add = T)
plot(PS4_HER_bounds, col = "red", add = T)
# 1 cell less in the north
new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$xmin <- PS4_HER_sent@extent@xmin + 0
new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$xmax <- PS4_HER_sent@extent@xmax + 0
new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$ymin <- PS4_HER_sent@extent@ymin + 0
new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$ymax <- PS4_HER_sent@extent@ymax - 10
# check output
PS4_HER_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS4_HER"),]$ymax))
plot(PS4_HER_sent)
plot(PS4_HER_bounds_new, col = "red", add = T)
# nice one!

# PS4 KOM 
plot(PS4_KOM_sent)
plot(extent(PS4_KOM), col = "blue", add = T)
plot(PS4_KOM_bounds,  add = T)
# remove 1 row in the west
new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$xmin <- PS4_KOM_sent@extent@xmin + 10
new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$xmax <- PS4_KOM_sent@extent@xmax + 0
new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$ymin <- PS4_KOM_sent@extent@ymin + 0
new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$ymax <- PS4_KOM_sent@extent@ymax + 0
# check output
PS4_KOM_bounds_new <- extent(c(new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$xmin,
                               new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$xmax,
                               new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$ymin,
                               new_bounds[which(new_bounds$site_veg == "PS4_KOM"),]$ymax))
plot(PS4_KOM_sent)
plot(PS4_KOM_bounds_new, col = "red", add = T)
# nice one!

# write a quick funciton for future use that creates an extent object form the data.frame
get_sent_extent <- function(site_veg_id, sent_boundaries) {
  extent_object <- extent(c(sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmin,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmax,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymin,
                            sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymax))
  return(extent_object)
  }
# test
new_extent_object <- get_sent_extent("PS1_KOM", new_bounds)
PS1_KOM_bounds_new == new_extent_object
# nice one!

# safe file
write.csv(new_bounds, file = paste0(script_path, "ps_sent_site_bounds.csv"))
save(new_bounds, file = paste0(script_path, "ps_sent_site_bounds.Rda"))
