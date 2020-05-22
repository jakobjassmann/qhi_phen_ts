# Ground-based time-series validation for Figure 5
# Jakob Assmann jakobjassmann@gmail.com Feb 2019

library(dplyr)
library(cowplot)
library(tidyverse)
library(raster)
library(rgdal)
library(rasterVis)
library(gridExtra)
library(grid)
library(viridis)

#### 1) Preparations ----
data_out_path <- "data/fig_5_ground_based_phenology/"
figure_out_path <- "figures/fig_5_ground_based_phenology/"
# Load phenology and NDVI observations
load(paste0(data_out_path, "gb_phenmeans_abs.Rda"))
load(paste0(data_out_path, "gbstats_drone.Rda"))
# Load Meta-Data
load("data/meta_data.Rda")
#### Merge two data frames
# Determine distinct doys each data frame and export for manual matching
doy_distinct_drone <- gbstats_df %>% dplyr::select(site_veg, year, doy) %>% distinct()
doy_distinct_phen <- phen_means %>% ungroup() %>% 
  dplyr::select(site_veg, year, doy) %>% distinct()
#write.csv(doy_distinct_drone, 
#          file = paste0(script_path, "doy_distinct_drone.csv"))
#write.csv(doy_distinct_phen, 
#          file = paste0(script_path, "doy_distinct_phen.csv"))
# Manual matching of date pairs carried out outside r

# Load manually made key
doy_drone_phen_key <- read.csv(file = paste0(data_out_path,
                                             "doy_drone_phen_key.csv"), 
                               stringsAsFactors = F)
# Make clean version of key for output as table
drone_phen_key_pretty <- doy_drone_phen_key %>% 
  mutate(site = substr(site_veg,1,3),
         veg = substr(site_veg,5,7)) %>%
  group_by(site, veg) %>% 
  na.omit() %>% 
  dplyr::select(site,veg, year, drone_doy, phen_doy)
write.csv(drone_phen_key_pretty, 
          paste0(data_out_path, "doy_drone_phen_key_pretty.csv"),
          row.names = F)

# Calculate max and mean time difference between drone and phen obs
max(abs(doy_drone_phen_key$drone_doy - doy_drone_phen_key$key_phen), na.rm =T)
mean(abs(doy_drone_phen_key$drone_doy - doy_drone_phen_key$key_phen), na.rm =T)
# Three days, manual matching worked well!

# Retain only values that are not na in the key
doy_drone_phen_key <- doy_drone_phen_key[!is.na(doy_drone_phen_key$key_phen),]

# Prep matching identifiers for data frames
gbstats_df$site_veg_year_doy <- paste0(gbstats_df$site_veg,
                                       gbstats_df$year,
                                       gbstats_df$doy)
phen_means$site_veg_year_doy <- paste0(phen_means$site_veg,
                                       phen_means$year,
                                       phen_means$doy)
doy_drone_phen_key$key_phen_doy <- paste0(doy_drone_phen_key$site_veg,
                                          doy_drone_phen_key$year,
                                          doy_drone_phen_key$key_phen)
doy_drone_phen_key$key_drone_doy <- paste0(doy_drone_phen_key$site_veg,
                                           doy_drone_phen_key$year,
                                           doy_drone_phen_key$key_drone)
# retain only matching site_veg, year, doy combos in the drone and phen dfs
phen_means_matching <- phen_means[phen_means$site_veg_year_doy %in% 
                                    doy_drone_phen_key$key_phen_doy,]
gbstats_df <- gbstats_df[gbstats_df$site_veg_year_doy %in% 
                           doy_drone_phen_key$key_drone_doy,]
# check that worked
phen_means_matching$site_veg_year_doy %in% doy_drone_phen_key$key_phen_doy
gbstats_df$site_veg_year_doy %in% doy_drone_phen_key$key_drone_doy
# Nice!

# subset by veg type as veg types have different numbers of species.
gbstats_df_HER <- gbstats_df[gbstats_df$veg_type == "HER",]
gbstats_df_KOM <- gbstats_df[gbstats_df$veg_type == "KOM",]
# Expand df by 4 for HER
gbstats_df_HER <- rbind(gbstats_df_HER,gbstats_df_HER,
                        gbstats_df_HER,gbstats_df_HER)
# and by 3 for KOM
gbstats_df_KOM <- rbind(gbstats_df_KOM,gbstats_df_KOM,
                        gbstats_df_KOM)

# remerge dfs into one big df
gbstats_df <- rbind(gbstats_df_HER, gbstats_df_KOM)
# Sort by site_veg_year_doy_colum
gbstats_df <- gbstats_df[order(gbstats_df$site_veg_year_doy),]

# and sort the phen means also
phen_means_matching <- phen_means_matching[order(
  phen_means_matching$site_veg_year_doy),]

# Finally create merged df of observaitons
gb_phen_ndvi <- data.frame(site = substr(gbstats_df$site_veg, 1,3),
                           site_veg = gbstats_df$site_veg,
                           year = gbstats_df$year,
                           site_veg_year = paste0(gbstats_df$site_veg, "_",
                                                  gbstats_df$year),
                           doy = as.integer(gbstats_df$doy),
                           species = phen_means_matching$species,
                           mean_leaf = phen_means_matching$mean_leaf,
                           sd_leaf = phen_means_matching$sd_leaf,
                           mean_ndvi = gbstats_df$ndvi_mean,
                           sd_ndvi = gbstats_df$ndvi_sd)

save(gb_phen_ndvi, file = paste0(data_out_path, "gb_phen_ndvi.Rda"))

### 2) Leav Length vs. NDVI correlation ----

# Calculate spearmans correlation for each time-series
cor_spear <- gb_phen_ndvi %>% 
  group_by(species, site_veg_year) %>%
  group_map(function(x, y) {
    data.frame(
      species = y[1],
      site_veg_year = y[2],
      cor_coef = cor(x$mean_leaf, x$mean_ndvi, method = "spearman"),
      p_value = cor.test(x$mean_leaf, x$mean_ndvi, method = "spearman")$p.value)
  }) %>%
  bind_rows() 

# Derive mean for all non-evergreen species
round(mean(cor_spear$cor_coef[cor_spear$species != "DRYINT"]), 2)

# Calculate mean corellation coef for a species
cor_spear_sum <- cor_spear %>%
  group_by(species) %>% 
  summarise(mean_cor_coef = round(mean(cor_coef),2))

# Add pretty names to data franes
name_lookup <- data.frame(species = c("ARCLAT",
                       "DRYINT",
                       "ERIVAG",
                       "SALARC",
                       "SALPUL"),
           name_pretty = c("Arctagrostis latifolia",
                           "Dryas integrifolia",
                           "Eriophorum vaginatum",
                           "Salix arctica",
                           "Salix pulchra"))
cor_spear <- cor_spear %>%
  merge(name_lookup, by.x = "species", by.y = "species") %>%
  dplyr::select(name_pretty, site_veg_year, cor_coef, p_value)
names(cor_spear) <- c("Species",
                      "Time-Series",
                      "r",
                      "p-value")
cor_spear$r <- round(cor_spear$r, 2)
cor_spear$'p-value' <- round(cor_spear$'p-value', 3)

cor_spear_sum <- cor_spear_sum %>%
  merge(name_lookup, by.x = "species", by.y = "species") %>%
  dplyr::select(name_pretty, mean_cor_coef)
names(cor_spear_sum) <- c("Species",
                          "mean r")

# Export tables
write.csv(cor_spear, 
          paste0(data_out_path, "leaf_length_cor_by_ts.csv"),
          row.names = F)
write.csv(cor_spear_sum, 
          paste0(data_out_path, "leaf_length_cor_by_spp.csv"),
          row.names = F)

### 2) Leaf Length vs. NDVI Plot ----
pretty_leaf_vs_ndvi_plot <- function(species_to_plot) {
  gb_phen_ndvi_sub <- gb_phen_ndvi %>% filter(species == species_to_plot)
  colourkey <- data.frame(species = c("ARCLAT",
                                      "DRYINT",
                                      "ERIVAG",
                                      "SALARC",
                                      "SALPUL"),
                          name_pretty = c("Arctagrostis latifolia",
                                          "Dryas integrifolia",
                                          "Eriophorum vaginatum",
                                          "Salix arctica",
                                          "Salix pulchra"),
                          colour = c("#112F41FF",
                                     "#0894A1FF",
                                     "#47AB6CFF",
                                     "#F2B134FF",
                                     "#ED553BFF"),
                          x_axis_fac = c(50,
                                         5,
                                         20,
                                         10,
                                         5),
                          stringsAsFactors = F)
  species_name <- colourkey[colourkey$species == species_to_plot,]$name_pretty
  species_colour <- colourkey[colourkey$species == species_to_plot,]$colour
  x_axis_fac <- colourkey[colourkey$species == species_to_plot,]$x_axis_fac
  min_x <- floor(min(gb_phen_ndvi_sub$mean_leaf) / x_axis_fac) * x_axis_fac
  max_x <- ceiling( max(gb_phen_ndvi_sub$mean_leaf) / x_axis_fac) * x_axis_fac
  min_y <- floor(min(gb_phen_ndvi_sub$mean_ndvi) / 0.1) * 0.1
  max_y <- ceiling(max(gb_phen_ndvi_sub$mean_ndvi) / 0.1) * 0.1
  print(paste0(min_x, " ", max_x))
  print(paste0(min_y, " ", max_y))
  mean_cor_coef <- cor_spear_sum[cor_spear_sum$Species == species_name,2]
  leaf_vs_ndvi_by_spp <- 
    ggplot(gb_phen_ndvi_sub, aes(x = mean_leaf, y = mean_ndvi, 
                                 group = site_veg_year,
                                 linetype = year)) +
    geom_point(colour = species_colour) +
    geom_smooth(method = 'lm', colour = species_colour, se = F,) +
    ylab("Mean NDVI") +
    xlab("Mean length of longest leaf (cm)") +
    scale_x_continuous(limits = c(min_x, max_x),
                       breaks = seq(min_x, max_x, x_axis_fac)) +
    
    scale_y_continuous(limits = c(min_y, max_y), 
                       breaks = seq(min_y, max_y, 0.1)) +
    scale_linetype_manual(values = c(2,1)) +
    annotate("text", x = ((max_x - min_x) / 2 + min_x), y = max_y, 
             label = species_name, 
             size = 6, 
             colour = species_colour,
             fontface = "italic",
             hjust = 0.5) +
    annotate("text", x = max_x, y = min_y, 
             label = paste0("mean r = ", mean_cor_coef), 
             size = 5.5, 
             colour = "black",
             hjust = 1,
             vjust = 0) +
    theme(legend.position = "none")
  return(leaf_vs_ndvi_by_spp)
}
pretty_leaf_vs_ndvi_plots <- lapply(sort(unique(gb_phen_ndvi$species)), 
                                    pretty_leaf_vs_ndvi_plot)
pretty_leaf_vs_ndvi_plots_grid <- 
  plot_grid(plotlist = pretty_leaf_vs_ndvi_plots, ncol = 5)
save_plot(pretty_leaf_vs_ndvi_plots_grid, 
          filename = paste0(figure_out_path, 
                            "leaf_length_vs_ndvi_pretty_grid_abs.png"),
          base_aspect_ratio = 5,
          base_height = 4)

### 3) Ground-Based Time-Series for PS2  ----
# Now create PS2 HER and KOM time series plots
PS2_gb_phen_ndvi <- gb_phen_ndvi %>% filter(site == "PS2")
PS2_HER_plot <- function(species_to_plot){
  PS2_gb_phen_ndvi_sub <- PS2_gb_phen_ndvi %>% 
    filter(site_veg_year == "PS2_HER_2017", species == species_to_plot)
  colourkey <- data.frame(species = c("ARCLAT",
                                      "DRYINT",
                                      "ERIVAG",
                                      "SALARC",
                                      "SALPUL"),
                          name_pretty = c("Arctagrostis latifolia",
                                          "Dryas integrifolia",
                                          "Eriophorum vaginatum",
                                          "Salix arctica",
                                          "Salix pulchra"),
                          colour = c("#112F41FF",
                                     "#0894A1FF",
                                     "#47AB6CFF",
                                     "#F2B134FF",
                                     "#ED553BFF"),
                          y_axis_fac = c(10,
                                         1,
                                         20,
                                         10,
                                         5),
                          stringsAsFactors = F)
  species_name <- colourkey[colourkey$species == species_to_plot,]$name
  species_colour <- colourkey[colourkey$species == species_to_plot,]$colour
  x_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$x_axis_labels
  y_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$y_axis_labels
  y_step <- colourkey[colourkey$species == species_to_plot,]$y_axis_fac
  y_min <- floor(min(PS2_gb_phen_ndvi_sub$mean_leaf) / y_step) * y_step
  y_max <- ceiling( max(PS2_gb_phen_ndvi_sub$mean_leaf) / y_step) * y_step
  print(paste0(y_min, " ", y_max, " ", y_step))
  ggplot(PS2_gb_phen_ndvi_sub, aes(x = doy, y = mean_leaf)) +
    geom_point(colour = species_colour) +
    scale_x_continuous(limits = c(170, 230), 
                       breaks = seq(170,230, 10)) +
    scale_y_continuous(limits = c(y_min, y_max), 
                       breaks = seq(y_min,y_max, y_step)) +
    geom_smooth(method ="lm", se = F, colour = species_colour) +
    ylab("Mean Leaf-Length (mm)") +
    xlab("Day of Year") +
    annotate("text", x = 200, y = y_max, 
             label = species_name, 
             size = 6, 
             colour = species_colour,
             fontface = "italic") +
    theme(axis.title.x = element_text(colour = x_axis_label),
          axis.title.y = element_text(colour = y_axis_label))
}
PS2_HER_plot_list <- lapply(c("ARCLAT", "DRYINT", "ERIVAG", "SALPUL"),
                            PS2_HER_plot)
PS2_gb_phen_ndvi_sub <- PS2_gb_phen_ndvi %>% 
  filter(site_veg_year == "PS2_HER_2017") %>% distinct(doy, mean_ndvi)
PS2_HER_plot_NDVI <- ggplot(PS2_gb_phen_ndvi_sub, aes(x = doy, y = mean_ndvi)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F, colour = "black") +
  scale_x_continuous(limits = c(170, 230), 
                     breaks = seq(170,230, 10)) +
  scale_y_continuous(limits = c(0.4, 0.8), 
                     breaks = seq(0.4, 0.8, 0.1)) +
  ylab("Mean NDVI") +
  xlab("Day of Year") +
  annotate("text", x = 200, y = 0.8, 
           label = "NDVI", 
           size = 6)
PS2_HER_plot_list$ndvi <- PS2_HER_plot_NDVI
PS2_HER_plot_list_grid <- plot_grid(plotlist = PS2_HER_plot_list, nrow = 1)
save_plot(PS2_HER_plot_list_grid, 
          filename = paste0(figure_out_path, "PS2_HER_2017_plot_list_grid.png"),
          base_aspect_ratio = 5,
          base_height = 4)


# PS2 KOM
# Now create PS2 time series plots
PS2_gb_phen_ndvi <- gb_phen_ndvi %>% filter(site == "PS2")
PS2_KOM_plot <- function(species_to_plot){
  PS2_gb_phen_ndvi_sub <- PS2_gb_phen_ndvi %>% 
    filter(site_veg_year == "PS2_KOM_2017", species == species_to_plot)
  colourkey <- data.frame(species = c("ARCLAT",
                                      "DRYINT",
                                      "ERIVAG",
                                      "SALARC",
                                      "SALPUL"),
                          name_pretty = c("Arctagrostis latifolia",
                                          "Dryas integrifolia",
                                          "Eriophorum vaginatum",
                                          "Salix arctica",
                                          "Salix pulchra"),
                          colour = c("#112F41FF",
                                     "#0894A1FF",
                                     "#47AB6CFF",
                                     "#F2B134FF",
                                     "#ED553BFF"),
                          y_axis_fac = c(10,
                                         1,
                                         20,
                                         2,
                                         1),
                          stringsAsFactors = F)
  species_name <- colourkey[colourkey$species == species_to_plot,]$name
  species_colour <- colourkey[colourkey$species == species_to_plot,]$colour
  x_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$x_axis_labels
  y_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$y_axis_labels
  y_step <- colourkey[colourkey$species == species_to_plot,]$y_axis_fac
  y_min <- floor(min(PS2_gb_phen_ndvi_sub$mean_leaf) / y_step) * y_step
  y_max <- ceiling( max(PS2_gb_phen_ndvi_sub$mean_leaf) / y_step) * y_step
  print(paste0(y_min, " ", y_max, " ", y_step))
  ggplot(PS2_gb_phen_ndvi_sub, aes(x = doy, y = mean_leaf)) +
    geom_point(colour = species_colour) +
    scale_x_continuous(limits = c(170, 230), 
                       breaks = seq(170,230, 10)) +
    scale_y_continuous(limits = c(y_min, y_max), 
                       breaks = seq(y_min,y_max, y_step)) +
    geom_smooth(method ="lm", se = F, colour = species_colour) +
    ylab("Mean Leaf-Length (mm)") +
    xlab("Day of Year") +
    annotate("text", x = 200, y = y_max, 
             label = species_name, 
             size = 6, 
             colour = species_colour,
             fontface = "italic") +
    theme(axis.title.x = element_text(colour = x_axis_label),
          axis.title.y = element_text(colour = y_axis_label))
}
PS2_KOM_plot_list <- lapply(c("ARCLAT", "DRYINT", "SALARC"), PS2_KOM_plot)
PS2_gb_phen_ndvi_sub <- PS2_gb_phen_ndvi %>% 
  filter(site_veg_year == "PS2_KOM_2017") %>% distinct(doy, mean_ndvi)
PS2_KOM_plot_NDVI <- ggplot(PS2_gb_phen_ndvi_sub, 
                            aes(x = doy, y = mean_ndvi)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F, colour = "black") +
  scale_x_continuous(limits = c(170, 230), 
                     breaks = seq(170,230, 10)) +
  scale_y_continuous(limits = c(0.4, 0.8), 
                     breaks = seq(0.4, 0.8, 0.1)) +
  ylab("Mean NDVI") +
  xlab("Day of Year") +
  annotate("text", x = 200, y = 0.8, 
           label = "NDVI", 
           size = 6)
PS2_KOM_plot_list$ndvi <- PS2_KOM_plot_NDVI
PS2_KOM_plot_list_grid <- plot_grid(plotlist = PS2_KOM_plot_list, nrow = 1)
save_plot(PS2_KOM_plot_list_grid, 
          filename = paste0(figure_out_path, "PS2_KOM_2017_plot_list_grid.png"),
          base_aspect_ratio = 4,
          base_height = 4)

### 4) Drone-Based Time-Series for PS2 ----

# Load grond plot coordinates
plot_coords <- read.csv("data/fig_5_ground_based_phenology/gb_plot_coordinates.csv")

# grab utm zone 7 prj4 data
epsg <- make_EPSG()
utm_z7N <- as.character(epsg %>% filter(code == 32607) %>% dplyr::select(prj4))
rm(epsg)

# Define a function to create a polygon from for coordinate pairs
plot_bound_poly <- function(site_veg){
  plot_bound_coords <- 
    as.matrix(plot_coords[plot_coords$site_veg == site_veg,4:5])
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

# Set specifiy combos
ts_combos <- data.frame(
  site_veg = c("PS2_HER",
               "PS2_KOM"),
  year = format(as.Date(c("2017",
                          "2017"),
                        format = "%Y"), "%Y"),
  stringsAsFactors = F)

# Define function for creating plot time-series
pretty_plots <- function(site_veg_es, year_es, agg_level) {  
  cat("starting: ", site_veg_es, "_", year_es, sep = "")
  # subset meta data
  ts_objects <- meta_data %>% 
    filter(site_veg == site_veg_es, 
           format(date, "%Y") == year_es, 
           band == "NDVI")
  # filter multiples for the 3 August for PS1 in 2016 and remove PS4 HER 2017-07-17 (outlier)
  if(site_veg_es == "PS2_KOM" & year_es == format(as.Date("2011", "%Y"),"%Y")) {
    ts_objects <- ts_objects[c(2,3,4,5),]
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
  
  # create doy column
  ts_objects$doy <- format(ts_objects$date, "%j")
  ts_objects <- ts_objects[order(ts_objects$doy),]
  # filter out doys not in ground-based validation time-series
  ts_objects <- ts_objects %>% 
    filter(doy %in% 
             PS2_gb_phen_ndvi[PS2_gb_phen_ndvi$site_veg == site_veg_es,]$doy)
  
  # plot_first in time series
  first_raster <- get(ts_objects$object_gbplot[1])
  # aggregate too account for geo-locaiton error in the time series
  if(agg_level != 1) {first_raster <- aggregate(first_raster, agg_level)}
  ndvi_min <- floor(10 * first_raster@data@min) / 10
  ndvi_max <- ceiling(10 * first_raster@data@max) / 10
  
  first_doy_plot <- levelplot(
    first_raster, 
    margin = FALSE, 
    colorkey= list(
      space='bottom',
      labels=list(at=seq(ndvi_min, ndvi_max, 0.1),
                  font=2, cex = 1.3)
    ),
    par.settings = list(
      axis.line=list(col='transparent'),
      par.main.text = list(font = 2, # make it bold
                           just = "left",
                           cex = 1.3,
                           x = grid::unit(15, "mm"))
    ),
    main = paste0("DOY: ", ts_objects$doy[1]),
    scales=list(draw=FALSE),
    col.regions=viridis(100),                  
    at=seq(ndvi_min, ndvi_max, (ndvi_max - ndvi_min) / 100),
    legend=list(left=list(fun=grid::textGrob("NDVI", y=0.05, x = 1.4, 
                                             gp=gpar(cex=1.3, 
                                                     fontface = "bold"))))
  )
  
  # Calculate raster differences per doy. 
  diff_rasters <- lapply(
    ts_objects$object_gbplot[2:length(ts_objects$object_gbplot)],
    function(x){
      reprojected_raster <- resample(get(x), first_raster)
      diff_raster <- reprojected_raster - first_raster
      names(diff_raster) <- x
      return(diff_raster)
    })
  min_diff <- min(unlist(lapply(diff_rasters, function(x){x@data@min})))
  max_diff <- max(unlist(lapply(diff_rasters, function(x){x@data@max})))
  min_diff <- floor(min_diff * 10) / 10
  max_diff <- ceiling(max_diff * 10) / 10
  diff_plots <- mapply(function(x, doy){
    levelplot(
      x, 
      margin = FALSE,  
      colorkey= list(
        space='bottom',
        labels=list(at=seq(min_diff, max_diff, 0.1),
                    font=2, cex = 1.3)
      ),
      par.settings = list(
        axis.line=list(col='transparent'),
        par.main.text = list(font = 2, # make it bold
                             just = "left",
                             x = grid::unit(26, "mm"),
                             cex = 1.3
        )
      ),
      main = paste0("DOY: ", doy),
      scales=list(draw=FALSE),
      col.regions=magma(100),                  
      
      at=seq(min_diff, max_diff, (max_diff - min_diff) / 100),
      legend=list(left=list(fun=grid::textGrob("Diff. NDVI", y=0.05, x = 1.425, 
                                               gp=gpar(cex=1.3, 
                                                       fontface = "bold")))))
  }, diff_rasters, ts_objects$doy[2:length(ts_objects$doy)], SIMPLIFY = F)
  gbplots <- list(first_doy_plot)
  for(i in 1:length(diff_plots)){gbplots[[1+i]] <- diff_plots[[i]]}
  png(paste0(figure_out_path, 
             site_veg_es, "_", year_es, 
             "_gb_drone_ts_diff_", agg_level * 5, "cm.png"), 
      width = 4 * length(gbplots),
      height = 4, 
      units = "in", 
      res = 200)
  print(marrangeGrob(grobs = gbplots,
                     ncol = length(gbplots), 
                     nrow = 1,
                     top = textGrob("")
                     #top = textGrob(paste0(site_veg_es, "_", year_es), 
                     #                gp=gpar(fontsize=14))
  )
  )
  dev.off()
  
  # return stats dataframe
  return(NULL)
}

# Execute function for time-series
# with native resolution (agg = 1)
pretty_plots("PS2_HER", "2017", 1)
pretty_plots("PS2_KOM", "2017", 1)

# Also available with other aggregation steps
# pretty_plots("PS2_HER", "2017", 3)
# pretty_plots("PS2_KOM", "2017", 3)
# 
# pretty_plots("PS2_HER", "2017", 6)
# pretty_plots("PS2_KOM", "2017", 6)