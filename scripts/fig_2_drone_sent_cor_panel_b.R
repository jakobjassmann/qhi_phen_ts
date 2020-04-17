# Phenology Time Series Snapshop PS2 KOM 2017
# comparison between multispectral measurments from drone and sentinel
# Jakob Assmann j.assmann@ed.ac.uk 3 October 2018
# Updated 17 April 2020 

### Preparations ----

# load dependencies
library(raster)
library(rgdal)
library(rasterVis)
library(dplyr)
library(ggplot2)
library(viridisLite)
library(gridExtra)
library(grid)
library(cowplot)
library(png)

# Set global parameters / load site boundaries and meta data
figure_out_path <- "figures/"
site_boundaries <- read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
load("data/meta_data.Rda")

### 1) Plot RGB rasters ----
# Load RGB files
# This step was done on the UoE Workstation. The original RGB tif has since
# disappeared however, the output still exists and can be found in:
# /figures/fig_2_panel_b/PS2_KOM_20170717_RGB_cropped.png

# PS2_KOM_20170717_RGB <- 
#   brick("/Volumes/BowheadRdge/phen_time_series/2017/RGB/PS2_KOM_2017-07-17_RGB.tif")
# PS2_KOM_20170717_RGB_cropped <-  crop(PS2_KOM_20170717_RGB, PS2_KOM_extent)
# 
# png(paste0(figure_out_path, "/figure_2_panel_bPS2_KOM_20170717_RGB_cropped.png"), 
#     width = 2,
#     height = 2, 
#     units = "in", 
#     res = 300)
# plotRGB(PS2_KOM_20170717_RGB_cropped, stretch = "hist")
# dev.off()

### 2) Gather NDVI rasters ----

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
PS2_KOM_extent <- get_sent_extent("PS2_KOM", site_boundaries)

# Load drone and sentinel rasters for PS2 on 17 July 2017
# drone as raster
PS2_drone_rasters <- meta_data %>% 
  filter(date == as.Date("2017-07-17"), 
         site_veg == "PS2_KOM", 
         band != "NDVI")
list2env(
  lapply(
    setNames(PS2_drone_rasters$file_path, 
             make.names(PS2_drone_rasters$object_name)),
    raster), 
  envir = .GlobalEnv)

# sentinel as brick
PS2_sentinel_rasters <- meta_data %>% 
  filter(date == as.Date("2017-07-17"), is.na(site_veg))
list2env(
  lapply(
    setNames(PS2_sentinel_rasters$file_path, 
             make.names(PS2_sentinel_rasters$object_name)),
    brick), 
  envir = .GlobalEnv)



# Crop rasters / brick to site extent
PS2_KOM_20170717_50m_red_cropped <- crop(PS2_KOM_20170717_50m_red, PS2_KOM_extent)
PS2_KOM_20170717_50m_nir_cropped <- crop(PS2_KOM_20170717_50m_nir, PS2_KOM_extent)
PS2_KOM_sentinel <- crop(L2A_QHI_20170717_cldsmskd_10m_brick, PS2_KOM_extent)


# Aggregate drone rasters to sentinel grid and resolution then resample to 
# using nearest neighbour to match senntinel grid
PS2_KOM_20170717_50m_red_cropped_agg <- aggregate(
  PS2_KOM_20170717_50m_red_cropped, 
  floor(10/res(PS2_KOM_20170717_50m_nir_cropped)[1]), 
  expand = F,
  fun = mean)
PS2_KOM_20170717_50m_red_cropped_resamp <- resample(
  PS2_KOM_20170717_50m_red_cropped_agg,
  PS2_KOM_sentinel, method = 'ngb')
PS2_KOM_20170717_50m_nir_cropped_agg<- aggregate(
  PS2_KOM_20170717_50m_nir_cropped, 
  floor(10/res(PS2_KOM_20170717_50m_nir_cropped)[1]), 
  expand = F,
  fun = mean)
PS2_KOM_20170717_50m_nir_cropped_resamp <- resample(
  PS2_KOM_20170717_50m_nir_cropped_agg,
  PS2_KOM_sentinel, method = 'ngb')

# Calcuate NDVI for drone and sentinel data
# NDVI =  (NIR - RED) / (NIR + RED)
NDVI <- function(red_band, nir_band) {(nir_band - red_band) / (nir_band + red_band) }

PS2_KOM_20170717_ndvi <- NDVI(PS2_KOM_20170717_50m_red_cropped, PS2_KOM_20170717_50m_nir_cropped)
PS2_KOM_20170717_ndvi_resamp <- NDVI(PS2_KOM_20170717_50m_red_cropped_resamp, PS2_KOM_20170717_50m_nir_cropped_resamp)
PS2_KOM_20170717_sentinel_ndvi <- NDVI(PS2_KOM_sentinel[[3]], PS2_KOM_sentinel[[4]])
PS2_KOM_20170717_diff <- PS2_KOM_20170717_sentinel_ndvi - PS2_KOM_20170717_ndvi_resamp

### 3) Plot NDVI Rasters ----

# Set global variables
site_name <- "PS2"
veg_type <- "KOM"
observation_date <- "20170717"
# Define Colour Ramp Themes for Levelplot with no borders on plot
viridis_no_borders <- rasterTheme(
  region = viridis(100, begin = 0, end = 1), # virdis scale with steps
  axis.line = list(lwd = 2),
  par.main.text = list(font = 1, cex = 2)) # normal face title
magma_no_borders <- rasterTheme(
  region = magma(100, begin = 0, end = 1), # virdis scale with steps
  axis.line = list(lwd = 2),
  par.main.text = list(font = 1, cex = 2)) # normal face title
purple_no_borders <- rasterTheme(
  region = colorRampPalette(c(viridis(3)[1], "white", viridis(3)[2]))(100), # virdis scale with steps
  axis.line = list(lwd = 2),
  par.main.text = list(font = 1, cex = 2)) # normal face title

# Determine adn set Scale Limits
scale_limits <-  sapply(
  c(paste0(site_name, "_", veg_type, "_", observation_date, "_ndvi"),
    paste0(site_name, "_", veg_type, "_", observation_date, "_ndvi_resamp"),
    paste0(site_name, "_", veg_type, "_", observation_date,"_sentinel_ndvi"),
    paste0(site_name, "_", veg_type, "_", observation_date,"_diff")),
  function(x) {
    cat(paste(x, 
              "min:", round(minValue(get(x)),3), 
              "max:", round(maxValue(get(x)),3)), "\n") 
    return(c(floor(minValue(get(x))*10)/10,
             ceiling(maxValue(get(x))*10)/10))}
)
  
native_low <- 0
native_up <- 1
native_breaks <- (native_up - native_low) / 100
resamp_low <- min(scale_limits[1,2:3]) 
resamp_up <- max(scale_limits[2,2:3])
resamp_breaks <- (resamp_up - resamp_low) / 100
diff_low <- -1* max(c(-1*scale_limits[1,4], scale_limits[2,4]))
diff_up <- max(c(-1*scale_limits[1,4], scale_limits[2,4]))
diff_breaks <- (diff_up - diff_low) / 100
  
scale_breaks_native <- seq(native_low, native_up, by = native_breaks) 
scale_breaks_resampled <- seq(resamp_low, resamp_up, by = resamp_breaks)
scale_breaks_diff <- seq(diff_low, diff_up, by = diff_breaks)
  
  
# plots 
plot_drone_native <- levelplot(
  get(paste0(site_name, "_", veg_type, "_", observation_date, "_ndvi")),
  main = list(label = "Drone 0.05 m    ", cex = 2), 
  margin = F, # no margins
  maxpixels = 6e5,
  xlab = list(label = "NDVI", cex = 1.5),
  colorkey = list(draw = T, axis.line = list(lwd = 2), axis.text = list(cex = 1.5)), 
  par.settings = viridis_no_borders, 
  scales = list(draw = F),
  at = scale_breaks_native)

plot_drone_resamp_magma <- levelplot(
  get(paste0(site_name, "_", veg_type, "_", observation_date, "_ndvi_resamp")),
  main = list(label = "Drone 10 m       ", cex = 2), 
  margin = F, # no margins
  maxpixels = 6e5,
  xlab = list(label = "NDVI", cex = 1.5),
  colorkey = list(draw = T, axis.line = list(lwd = 2), axis.text = list(cex = 1.5)), 
  par.settings = magma_no_borders, 
  scales = list(draw = F),
  at = scale_breaks_resampled)

plot_sentinel_magma <- levelplot(
  get(paste0(site_name, "_", veg_type, "_", observation_date,"_sentinel_ndvi")),
  main = list(label = "Sentinel 10 m      ", cex = 2), 
  margin = F, # no margins
  maxpixels = 6e5,
  xlab = list(label = "NDVI", cex = 1.5),
  colorkey = list(draw = T, axis.line = list(lwd = 2), axis.text = list(cex = 1.5)), 
  par.settings = magma_no_borders, 
  scales = list(draw = F),
  at = scale_breaks_resampled)

plot_diff_purple <- levelplot(
  get(paste0(site_name, "_", veg_type, "_", observation_date, "_diff")),
  main = list(label = "Difference        ", cex = 2), 
  margin = F, # no margins
  maxpixels = 6e5,
  colorkey = list(draw = T, axis.line = list(lwd = 2), axis.text = list(cex = 1.5)), 
  par.settings = purple_no_borders, 
  scales = list(draw = F),
  at = scale_breaks_diff)

### 4) Export plots ----
png(paste0(figure_out_path, "fig_2_panel_b/plot_drone_rgb.png"), 
    width = 4,
    height = 4,
    units = "in",
    res = 300)
par(yaxt = "n",
    xaxt = "n",
    mar = c(3,3,3,3),
    cex.main = 2,
    font.main = 1,
    family = "Helvetica") 
plotRGB(brick("figures/fig_2_panel_b/PS2_KOM_20170717_RGB_cropped.png"),
        axes = T,
        main = "Drone 0.013 m",
        xlab = "RGB")
box(lwd = 2)
dev.off()

png(paste0(figure_out_path, "fig_2_panel_b/plot_drone_native.png"), 
        width = 4,
        height = 4,
        units = "in",
        res = 300)
plot_drone_native
dev.off()

png(paste0(figure_out_path, "fig_2_panel_b/plot_drone_resamp_magma.png"), 
    width = 4,
    height = 4,
    units = "in",
    res = 300)
plot_drone_resamp_magma
dev.off()

png(paste0(figure_out_path, "fig_2_panel_b/plot_sentinel_magma.png"), 
    width = 4,
    height = 4,
    units = "in",
    res = 300)
plot_sentinel_magma
dev.off()

png(paste0(figure_out_path, "fig_2_panel_b/plot_diff_purple.png"), 
    width = 4,
    height = 4,
    units = "in",
    res = 300)
plot_diff_purple
dev.off()
