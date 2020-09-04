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
library(sf)

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
                           mean_leaf_stand = phen_means_matching$mean_leaf_stand,
                           sd_leaf = phen_means_matching$sd_leaf,
                           mean_ndvi = gbstats_df$ndvi_mean,
                           sd_ndvi = gbstats_df$ndvi_sd)

gb_phen_ndvi <- mutate(gb_phen_ndvi, veg_type = substr(site_veg, 5,7))

save(gb_phen_ndvi, file = paste0(data_out_path, "gb_phen_ndvi.Rda"))

### 2) Leaf Length vs. NDVI correlation ----

# Community level mean_ndvi
gb_phen_ndvi_com <- gb_phen_ndvi %>% 
  group_by(site_veg_year, site, veg_type, year, doy) %>% 
  summarise(comm_mean_leaf_stand = mean(mean_leaf_stand),
            mean_ndvi = mean(mean_ndvi)
            )

# Calculate spearmans correlation for each time-series
cor_spear_com <- gb_phen_ndvi_com %>% 
  group_by(site_veg_year) %>%
  group_map(function(x, y) {
    z <- data.frame(
      site_veg_year = y[1],
      cor_coef = cor(x$comm_mean_leaf_stand, x$mean_ndvi, method = "spearman"),
      p_value = cor.test(x$comm_mean_leaf_stand, x$mean_ndvi, method = "spearman")$p.value,
      n = x %>% summarise(n = n()) %>% pull(n))
  }) %>%
  bind_rows() 

# Calculate mean community correlation across all time-series
cor_spear_com_mean <- round(mean(cor_spear_com$cor_coef), 2)

# Export tables
write.csv(cor_spear_com %>%
            mutate(cor_coef = round(cor_coef, 2),
                   p_value = formatC(round(p_value, 3), digits = 3)), 
          paste0(data_out_path, "standard_leaf_length_cor_by_ts.csv"),
          row.names = F)

# Plot time-series for all sites and years
colour_scale_sites <- c("#4A44F2FF",
                        "#F20505FF", 
                        "#F2BE22FF",
                        "#9C9DA6FF")
com_mean_leaf_vs_ndvi_plot <- ggplot(gb_phen_ndvi_com, 
       aes(x = comm_mean_leaf_stand, 
           y = mean_ndvi, 
           colour = site, 
           group = site_veg_year,
           linetype = year,
           shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Mean standardised\nlongest leaf length",
       y = "\nMean NDVI",
       shape = "Vegetation Type",
       linetype = "Year",
       colour = "Site") +
  scale_y_continuous(limits = c(0.4, 0.8), breaks = seq(0.4,0.8,0.1)) +
  scale_x_continuous(limits = c(-1.25,1.25), breaks = seq(-2, 2, 0.5)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Site 1", "Site 2", "Site 3", "Site 4")) +
  scale_linetype_manual(values = c(2,1)) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  annotate("text", x = 1.25, y = 0.41, hjust = 1, size = 5.5,
           label = paste0("mean ρ = ", cor_spear_com_mean)) +
  theme_cowplot(18) +
  theme(legend.position = "none")

save_plot("figures/fig_5_ground_based_phenology/com_mean_leaf_vs_ndvi_plot.png",
          com_mean_leaf_vs_ndvi_plot,
          base_height = 5,
          base_aspect_ratio = 1.3 - 0.455)


# Doy vs. com mean leaf
com_doy_vs_mean_leaf_plot <- ggplot(gb_phen_ndvi_com, 
                                     aes(x = doy, 
                                         y = comm_mean_leaf_stand, 
                                         colour = site, 
                                         group = site_veg_year,
                                         linetype = year,
                                         shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Day of year\n",
       y =  "Mean standardised\nlongest leaf length",
       shape = "Vegetation Type",
       linetype = "Year",
       colour = "Site") +
  scale_y_continuous(limits = c(-1.25,1.25), breaks = seq(-2, 2, 0.5)) +
  scale_x_continuous(limits = c(170,230), breaks = seq(170, 230, 10)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Site 1", "Site 2", "Site 3", "Site 4")) +
  scale_linetype_manual(values = c(2,1)) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  theme_cowplot(18) +
  theme(legend.position = "none")

save_plot("figures/fig_5_ground_based_phenology/com_doy_vs_mean_leaf_plot.png",
          com_doy_vs_mean_leaf_plot,
          base_height = 5,
          base_aspect_ratio = 1.35 - 0.455)

# Doy vs. mean NDVI
com_doy_vs_ndvi_plot <- ggplot(gb_phen_ndvi_com, 
                               aes(x = doy, 
                                   y = mean_ndvi, 
                                   colour = site, 
                                   group = site_veg_year,
                                   linetype = year,
                                   shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Day of year\n",
       y =  "\nMean NDVI",
       shape = "Vegetation Type",
       linetype = "Year",
       colour = "Site") +
  scale_y_continuous(limits = c(0.4,0.8), breaks = seq(0.4, 0.8, 0.1)) +
  scale_x_continuous(limits = c(170,230), breaks = seq(170, 230, 10)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Area 1", "Area 2", "Area 3", "Area 4")) +
  scale_linetype_manual(values = c(2,1)) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  theme_cowplot(18) +
  theme(legend.key.width=unit(0.455,"inch"))

save_plot("figures/fig_5_ground_based_phenology/com_doy_vs_ndvi_plot.png",
          com_doy_vs_ndvi_plot,
          base_height = 5,
          base_aspect_ratio = 1.35)

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
                                         40,
                                         10,
                                         5),
                          stringsAsFactors = F)
  y_label <- ""
  if(species_to_plot == "ARCLAT") y_label <- "Mean NDVI" 
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
    ylab(y_label) +
    xlab("Mean length of longest leaf (mm)") +
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
    theme_cowplot(20) +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 15))
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
    ylab("Mean length of longest leaf (mm)") +
    xlab("Day of Year") +
    annotate("text", x = 200, y = y_max, 
             label = species_name, 
             size = 6, 
             colour = species_colour,
             fontface = "italic") +
    theme(axis.title.x = element_text(colour = x_axis_label),
          axis.title.y = element_text(colour = y_axis_label),
          legend.position = "none")
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
    ylab("Mean length of longest leaf (mm)") +
    xlab("Day of Year") +
    annotate("text", x = 200, y = y_max, 
             label = species_name, 
             size = 6, 
             colour = species_colour,
             fontface = "italic") +
    theme(axis.title.x = element_text(colour = x_axis_label),
          axis.title.y = element_text(colour = y_axis_label),
          legend.position = "none")
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

## For all sites comprehensively all in one common panel

plot_species <- function(species_to_plot){
  gb_phen_sub <- gb_phen_ndvi %>% filter(species == species_to_plot)
  colourkey <- data.frame(species = c("ARCLAT",
                                      "DRYINT",
                                      "ERIVAG",
                                      "SALARC",
                                      "SALPUL"),
                          name_pretty = c("A. latifolia",
                                          "D. integrifolia",
                                          "E. vaginatum",
                                          "S. arctica",
                                          "S. pulchra"),
                          colour = c("#112F41FF",
                                     "#0894A1FF",
                                     "#47AB6CFF",
                                     "#F2B134FF",
                                     "#ED553BFF"),
                          y_axis_fac = c(50,
                                         5,
                                         20,
                                         5,
                                         5),
                          stringsAsFactors = F)
  y_label <- ""
  if(species_to_plot == "ARCLAT") y_label <- "Mean length of longest leaf (mm)"
  species_name <- colourkey[colourkey$species == species_to_plot,]$name
  species_colour <- colourkey[colourkey$species == species_to_plot,]$colour
  x_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$x_axis_labels
  y_axis_label <- 
    colourkey[colourkey$species == species_to_plot,]$y_axis_labels
  y_step <- colourkey[colourkey$species == species_to_plot,]$y_axis_fac
  y_min <- floor(min(gb_phen_sub$mean_leaf) / y_step) * y_step
  y_max <- ceiling( max(gb_phen_sub$mean_leaf) / y_step) * y_step
  print(paste0(y_min, " ", y_max, " ", y_step))
  ggplot(gb_phen_sub, aes(x = doy, y = mean_leaf, group = site_veg_year,
                      linetype = year)) +
    geom_point(colour = species_colour) +
    scale_x_continuous(limits = c(170, 230), 
                       breaks = seq(170,230, 20)) +
    scale_y_continuous(limits = c(y_min, y_max), 
                       breaks = seq(y_min,y_max, y_step)) +
    scale_linetype_manual(values = c(2,1)) +
    geom_smooth(method ="lm", se = F, colour = species_colour) +
    ylab(y_label) +
    xlab("Day of Year") +
    annotate("text", x = 200, y = y_max, 
             label = species_name, 
             size = 7, 
             colour = species_colour,
             fontface = "italic") +
    theme_cowplot(24) + 
    theme(axis.title.x = element_text(colour = x_axis_label),
          axis.title.y = element_text(colour = y_axis_label, size = 18),
          legend.position = "none")
}

gb_phen_ts_all_list <- lapply(unique(gb_phen_ndvi$species), plot_species)
gb_phen_ndvi_sub <- gb_phen_ndvi %>% 
  group_by(site_veg_year, year) %>% 
  distinct(doy, mean_ndvi)
NDVI_ts_all <- ggplot(gb_phen_ndvi, 
                            aes(x = doy, y = mean_ndvi, group = site_veg_year,
                                linetype = year)) +
  geom_point() +
  geom_smooth(method = 'lm', se = F, colour = "black") +
  scale_x_continuous(limits = c(170, 230), 
                     breaks = seq(170,230, 20)) +
  scale_y_continuous(limits = c(0.4, 0.8), 
                     breaks = seq(0.4, 0.8, 0.1)) +
  ylab("Mean NDVI") +
  xlab("Day of Year") +
  annotate("text", x = 200, y = 0.8, 
           label = "NDVI", 
           size = 7) +
  theme_cowplot(24) +
  theme(legend.position = "none")
gb_phen_ts_all_list$ndvi <- NDVI_ts_all
ts_all_plot_list_grid <- plot_grid(plotlist = gb_phen_ts_all_list, nrow = 1)
save_plot(ts_all_plot_list_grid, 
          filename = paste0(figure_out_path, "ts_gb_phen_all.png"),
          base_aspect_ratio = 5,
          base_height = 5)

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

# Create meta_data rgb rasters
rgb_meta <- data.frame(
  site_veg = c("PS2_HER",
               "PS2_KOM"),
  file_path = c("/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/PS2_HER/output/rgb/PS2_HER_2017-07-17_RGB.tif",
                "/Volumes/BowheadRdge/phen_time_series/final_outputs/2017/PS2_KOM/output/rgb/PS2_KOM_2017-07-17_RGB.tif"),
  date = as.Date(c("2017-07-17", "2017-07-17")),
  stringsAsFactors = F
)
# Define function for creating plot time-series
pretty_plots <- function(site_veg_es, year_es, agg_level, short = F) {  
  cat("starting: ", site_veg_es, "_", year_es, sep = "")
  # subset meta data
  ts_objects <- meta_data %>% 
    filter(site_veg == site_veg_es, 
           format(date, "%Y") == year_es, 
           band == "NDVI")
  # filter multiples for the 3 August for PS1 in 2016 and remove PS4 HER 2017-07-17 (outlier)
  if(site_veg_es == "PS2_KOM" & year_es == format(as.Date("2017", "%Y"),"%Y")) {
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
  # perhaps aggregate to account for geo-locaiton error in the time series
  if(agg_level != 1) {first_raster <- aggregate(first_raster, agg_level)}
  ndvi_min <- floor(10 * first_raster@data@min) / 10
  ndvi_max <- ceiling(10 * first_raster@data@max) / 10
  
  first_doy_plot <- levelplot(
    first_raster, 
    margin = FALSE, 
    colorkey= list(
      space='bottom',
      labels=list(at=seq(ndvi_min, ndvi_max, 0.1),
                  font=1, cex = 1.5)
    ),
    par.settings = list(
      axis.line=list(col='transparent'),
      par.main.text = list(font = 1,
                           just = "left",
                           cex = 1.5,
                           x = grid::unit(15, "mm"))
    ),
    main = paste0("DOY: ", ts_objects$doy[1]),
    scales=list(draw=FALSE),
    col.regions=viridis(100),                  
    at=seq(ndvi_min, ndvi_max, (ndvi_max - ndvi_min) / 100),
    legend=list(left=list(fun=grid::textGrob("NDVI", y=0.05, x = 1.4, 
                                             gp=gpar(cex=1.5))))
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
  
  # For final figure 5 plot only 3 of the ndvi rasters. (option short = T)
  if(site_veg_es == "PS2_HER" & year_es == format(as.Date("2017", "%Y"),"%Y") & 
     short == T) {
    ts_objects <- ts_objects[c(1,3,5),]
    diff_rasters <- list(diff_rasters[[2]], diff_rasters[[4]])
  }
  
  diff_plots <- mapply(function(x, doy){
    levelplot(
      x, 
      margin = FALSE,  
      colorkey= list(
        space='bottom',
        labels=list(at=seq(min_diff, max_diff,0.2),
                    font=1, cex = 1.5)
      ),
      par.settings = list(
        axis.line=list(col='transparent'),
        par.main.text = list(font = 1, 
                             just = "left",
                             x = grid::unit(26, "mm"),
                             cex = 1.5
        )
      ),
      main = paste0("DOY: ", doy),
      scales=list(draw=FALSE),
      col.regions=magma(100),                  
      
      at=seq(min_diff, max_diff, (max_diff - min_diff) / 100),
      legend=list(left=list(fun=grid::textGrob("Δ NDVI", y=0.05, x = 1.425, 
                                               gp=gpar(cex=1.5)))))
  }, diff_rasters, ts_objects$doy[2:length(ts_objects$doy)], SIMPLIFY = F)
  gbplots <- list(first_doy_plot)
  for(i in 1:length(diff_plots)){gbplots[[1+i]] <- diff_plots[[i]]}
  
  ## Add RGB plot
  # Load brick
  if(site_veg_es == "PS2_HER" ) rgb_brick <- brick(rgb_meta$file_path[rgb_meta$site_veg == site_veg_es])  # throw out alpha band
  if(site_veg_es == "PS2_KOM" ) rgb_brick <- brick(rgb_meta$file_path[rgb_meta$site_veg == site_veg_es])  # throw out alpha band
  # crop brick
  rgb_brick <- crop(rgb_brick, get(paste0(site_veg_es, "_gbplot_poly")))
  rgb_brick <- mask(rgb_brick, get(paste0(site_veg_es, "_gbplot_poly")))
  # discard alpha band
  rgb_brick <- rgb_brick[[1:3]]
  rgb_values <- getValues(rgb_brick)
  # Prep colour layer values with NAs 
  cols <- rgb_values[,1]
  # Create HEX GB color for cell values and update colour layer values
  cols[!is.na(cols)] <- rgb(rgb_values[!is.na(cols),], maxColorValue=255)
  cols <- factor(cols)
  # Creat single band temp raster
  rgb_raster <- raster(rgb_brick)
  # re-assign cell values
  rgb_raster[] <- cols
  
  # Plot
  rgb_plot <- levelplot(
    rgb_raster, 
    margin=FALSE,  # don't plot margins
    scales=list(draw=FALSE), # suppress axis labels
    col.regions=as.character(levels(cols)),
    colorkey=FALSE,
    par.settings = list(
      axis.line=list(col='transparent'),
      par.main.text = list(font = 1,
                           just = "left",
                           x = grid::unit(26, "mm"),
                           cex = 1.5
      )
    ),
    main = paste0("RGB"),
    #xlab.top = list("ffadsf", col = "white", cex = 0.8), # Quick work around to make sure spacing is the same
    xlab = list("df", col = "white", cex = 4.6), # Quick work around to make sure spacing is the same
    ylab = list("sd", col = "white", cex = 8) # This is a quick work around to replace the color key bar...
  ) 
  gbplots[[length(gbplots) + 1]] <- rgb_plot
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
pretty_plots("PS2_HER", "2017", 1, short = T)
pretty_plots("PS2_KOM", "2017", 1)

# Also available with other aggregation steps
# pretty_plots("PS2_HER", "2017", 3)
# pretty_plots("PS2_KOM", "2017", 3)
# 
# pretty_plots("PS2_HER", "2017", 6)
# pretty_plots("PS2_KOM", "2017", 6)


### 6) Sentinel-2 Ground-based phenology analysis ----

# retrieve distinct drone dates
doy_distinct_phen <- phen_means %>% ungroup() %>% 
  dplyr::select(site_veg, year, doy) %>% distinct() %>%
  mutate(doy_plus3 = as.numeric(doy) + 3,
         doy_minus3 = as.numeric(doy) - 3)
# ... adding columns with doy +- 3 days to match with sentinel data
# 3 days is the maximum day difference between drone and ground based phenology
# boservations

# retrieve distinc sentinel-2 dates
doy_distinct_sentinel <- meta_data %>%
  filter(sensor_id == "Sentinel 2A" | sensor_id == "Sentinel 2B") %>%
  mutate(doy = format.Date(date, "%j"),
         year = format.Date(date, "%Y")) %>%
  distinct()

gb_phen_sent_combos <- bind_rows(
  lapply(doy_distinct_sentinel$object_name, 
         function(scene_id){
           cat("Scene id: ", scene_id, "\n")
           scene_doy <- as.numeric(doy_distinct_sentinel$doy[doy_distinct_sentinel$object_name == scene_id])
           scene_year <- as.numeric(doy_distinct_sentinel$year[doy_distinct_sentinel$object_name == scene_id])
           scene_stats <- doy_distinct_sentinel[doy_distinct_sentinel$object_name == scene_id,]
           phen_date <- doy_distinct_phen %>%
             filter(year == scene_year & 
                      scene_doy >= doy_minus3 & 
                      scene_doy <= doy_plus3) 
           
           if(nrow(phen_date) == 0) {
             cat("Matching gb obs in +- 3 days: 0\n")
             phen_date[1,] <- NA
             names(phen_date) <- paste0("gb_", names(phen_date))
             names(scene_stats) <- paste0("sentinel_", names(scene_stats))
             phen_combo <- bind_cols(scene_stats,
                                     phen_date)
           } else{
             cat("Matching gb obs in +- 3 days: ", nrow(phen_date), "\n")
             phen_date$diff <- as.numeric(phen_date$doy) - scene_doy
             # For each site select the matching gb phen entry with the lowest
             # absolute date difference, if there are multiple ones, arrange by absolute difference
             # and take the earlier one.
             phen_date <- bind_rows(
               lapply(
                 unique(phen_date$site_veg),
                 function(site_id){
                   phen_date_sub <- phen_date %>% filter(site_veg == site_id)
                   phen_date_sub <- phen_date_sub[phen_date_sub$diff == min(abs(phen_date_sub$diff)),] %>%
                     arrange(diff)
                   phen_date_sub <- phen_date_sub[1,]
                   return(phen_date_sub)
                 }))
             names(phen_date) <- paste0("gb_", names(phen_date))
             names(scene_stats) <- paste0("sentinel_", names(scene_stats))
             scene_stats <- do.call("rbind", replicate(nrow(phen_date), 
                                                       scene_stats, 
                                                       simplify = FALSE))
             phen_combo <- bind_cols(scene_stats,
                                     phen_date) 
             
           }
           return(phen_combo)
         }))

# Next check how many sentinel scenes don't ahve a matching gb phen record:
gb_phen_sent_combos %>% filter(is.na(gb_site_veg)) %>%
  distinct(sentinel_doy) %>%
  summarise(n = n())
# 18 / 24 scenes don't have gb phen data! i.e. 6 do

# Throw out scenes without gb phen data.
gb_phen_sent_combos <- gb_phen_sent_combos %>% filter(!is.na(gb_site_veg))


## Load gb_plot geometries
# Ground based phenology plot coordinates
plot_coords <- read.csv("data/fig_5_ground_based_phenology//gb_plot_coordinates.csv")
# grab utm zone 7 prj4 data
epsg <- make_EPSG()
utm_z7N <- as.character(epsg %>% filter(code == 32607) %>% dplyr::select(prj4))
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

# And create one sf containing all the geometries
gb_polys <- lapply(plot_bound_poly_list, st_as_sf)
gb_polys <- do.call(rbind, gb_polys)
gb_polys$id <- unique(plot_coords$site_veg)

# NDVI helper function
NDVI <- function(red_band, nir_band) {
  (nir_band - red_band) / (nir_band + red_band) 
}

## Visualise the position of the plots on the sentinel grid
# load asentinel scene
sentinel_brick <- brick(
  unique(gb_phen_sent_combos$sentinel_file_path[1]))
sentinel_brick <- crop(sentinel_brick,
                       as_Spatial(st_buffer(gb_polys, 1000)))
# Calculate NDVI
sentinel_ndvi <- NDVI(sentinel_brick[[3]], sentinel_brick[[4]])

lapply(gb_polys$id,
       function(site){
         site_poly <- gb_polys[gb_polys$id == site,]
         sentinel_crop <- crop(sentinel_ndvi, as_Spatial(st_buffer(site_poly,10)))
         png(paste0("data/fig_5_ground_based_phenology/gb_phen_sentinel/", 
                    site, ".png"), 
             width = 4,
             height = 4,
             units = "in",
             res = 300)
         plot(sentinel_crop, main = site)
         plot(site_poly, add = T)
         dev.off()
       })
# Not all polygons fall into a single pixel -> weighted mean extraction.

# Extract sentinel NDVI gb plot ndvi mean for each scene
gb_phen_sent_combos <- bind_rows(lapply(unique(gb_phen_sent_combos$sentinel_object_name),
       function(scene_id){
         gb_phen_sent_combos_sub <- gb_phen_sent_combos %>% 
           filter(sentinel_object_name == scene_id)
         # load sentinel scene
         sentinel_brick <- brick(
           unique(gb_phen_sent_combos_sub$sentinel_file_path))
         # Crop to something managable
         sentinel_brick <- crop(sentinel_brick,
                                as_Spatial(st_buffer(gb_polys, 1000)))
         # Calculate NDVI
         sentinel_ndvi <- NDVI(sentinel_brick[[3]], sentinel_brick[[4]])
         # Extract mean weighted NDVI
         gb_plot_polys <- gb_polys %>% 
           filter(id %in% gb_phen_sent_combos_sub$gb_site_veg)
         gb_sentinel_ndvi <- raster::extract(sentinel_ndvi, as_Spatial(gb_plot_polys),
                         fun = mean, weights = T, sp = T) %>% st_as_sf %>%
           setNames(c("gb_site_veg", "mean_NDVI", "geometry")) %>%
           st_drop_geometry() %>% 
           dplyr::select(gb_site_veg, mean_NDVI)
         gb_phen_sent_combos_sub <- full_join(gb_phen_sent_combos_sub, gb_sentinel_ndvi)
       }))

# Save for later use
save(gb_phen_sent_combos, file = "data/fig_5_ground_based_phenology/gb_phen_sentinel/gb_phen_sent_combos.Rda")
load("data/fig_5_ground_based_phenology/gb_phen_sentinel/gb_phen_sent_combos.Rda")
ggplot(gb_phen_sent_combos, aes(x = gb_doy, y = mean_NDVI, color = gb_site_veg)) +
  geom_point() +
  facet_wrap(~gb_year)
# Combine with phenology data
names(phen_means)
gb_phen_sent <- phen_means %>% inner_join(gb_phen_sent_combos, by = c("year" = "gb_year",
                                                      "doy" = "gb_doy", "site_veg" = "gb_site_veg"))

# summarise to community mean
gb_phen_sent_com <- gb_phen_sent %>% 
  mutate(site_veg_year = paste0(site_veg, year),
         site = substr(site_veg, 1,3),
         veg_type = substr(site_veg, 5,7)) %>%
  group_by(site_veg_year, site, veg_type, year, doy) %>% 
  summarise(comm_mean_leaf_stand = mean(mean_leaf_stand),
            mean_ndvi = mean(mean_NDVI),
            sd_ndvi = sd(mean_NDVI)) 

# Check data points available
gb_phen_sent_com %>% group_by(site_veg_year) %>%
  summarise(n = n())

# Filter out all site site_veg_year with less than 4 data points
gb_phen_sent_com <- gb_phen_sent_com%>% 
  filter(site_veg_year %in% (gb_phen_sent_com %>% group_by(site_veg_year) %>%
                               summarise(n = n()) %>% filter(n >= 4) %>% pull(site_veg_year)))

# Remove single outlier for site 3 above 0.75
# PS3 KOM 2017 DOY 185 - 0.790 
gb_phen_sent_com <- filter(gb_phen_sent_com, comm_mean_leaf_stand <= 0.75)
gb_phen_sent_com %>% filter(site_veg_year == "PS3_KOM2017") %>% summarise(mean = mean(comm_mean_leaf_stand))
# Calculate spearmans correlation for those time-series with more than 5 data
# points
cor_spear_com_sent <- gb_phen_sent_com %>%
  group_by(site_veg_year) %>%
  group_map(function(x, y) {
    data.frame(
      site_veg_year = y[1],
      cor_coef = cor(x$comm_mean_leaf_stand, x$mean_ndvi, method = "spearman"),
      p_value = cor.test(x$comm_mean_leaf_stand, x$mean_ndvi, method = "spearman")$p.value,
      n = x %>% summarise(n = n()) %>% pull(n))
  }) %>%
  bind_rows()

# Calculate mean community correlation across all time-series
cor_spear_com_mean_sent <- round(mean(cor_spear_com_sent$cor_coef), 2)

# Export tables
write.csv(cor_spear_com_sent %>%
            mutate(cor_coef = round(cor_coef, 2),
                   p_value = formatC(round(p_value, 3), digits = 3)), 
          paste0(data_out_path, "standard_leaf_length_cor_by_ts_sent.csv"),
          row.names = F)

# Plot time-series for all sites and years
colour_scale_sites <- c("#4A44F2FF",
                        "#F20505FF", 
                        "#F2BE22FF",
                        "#9C9DA6FF")
com_mean_leaf_vs_ndvi_plot_sent <- ggplot(gb_phen_sent_com, 
                                     aes(x = comm_mean_leaf_stand, 
                                         y = mean_ndvi, 
                                         colour = site, 
                                         group = site_veg_year,
                                         linetype = veg_type,
                                         shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Mean standardised\nlongest leaf length",
       y = "\nMean NDVI",
       shape = "Vegetation Type",
       linetype = "Vegeation Type",
       colour = "Site") +
  scale_y_continuous(limits = c(0.5, 0.8), breaks = seq(0.4,0.8,0.1)) +
  scale_x_continuous(limits = c(-1.25,1.25), breaks = seq(-2, 2, 0.5)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Site 1", "Site 2", "Site 3", "Site 4")) +
  scale_linetype_manual(values = c(1,2),
                        labels = c("Tussock Sedge", "Dryas-Vetch")) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  annotate("text", x = 1.25, y = 0.51, hjust = 1, size = 5.5,
            label = paste0("mean ρ = ", cor_spear_com_mean_sent)) +
  theme_cowplot(18) +
  theme(legend.position = "none")

save_plot("figures/fig_5_ground_based_phenology/sent_com_mean_leaf_vs_ndvi_plot.png",
          com_mean_leaf_vs_ndvi_plot_sent,
          base_height = 5,
          base_aspect_ratio = 1.35 - 0.455)


# Doy vs. com mean leaf
com_doy_vs_mean_leaf_plots_sent <- ggplot(gb_phen_sent_com, 
                                    aes(x = as.numeric(doy), 
                                        y = comm_mean_leaf_stand, 
                                        colour = site, 
                                        group = site_veg_year,
                                        linetype = veg_type,
                                        shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Day of year\n",
       y =  "Mean standardised\nlongest leaf length",
       shape = "Vegetation Type",
       linetype = "Vegetation Type",
       colour = "Site") +
  scale_y_continuous(limits = c(-1.25,1.25), breaks = seq(-2, 2, 0.5)) +
  scale_x_continuous(limits = c(170,230), breaks = seq(170, 230, 10)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Site 1", "Site 2", "Site 3", "Site 4")) +
  scale_linetype_manual(values = c(1,2),
                        labels = c("Tussock Sedge", "Dryas-Vetch")) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  theme_cowplot(18) +
  theme(legend.position = "none")

save_plot("figures/fig_5_ground_based_phenology/sent_com_doy_vs_mean_leaf_plot.png",
          com_doy_vs_mean_leaf_plots_sent,
          base_height = 5,
          base_aspect_ratio = 1.35 - 0.455)

# Doy vs. mean NDVI
com_doy_vs_ndvi_plot_sent <- ggplot(gb_phen_sent_com, 
                               aes(x = as.numeric(doy), 
                                   y = mean_ndvi, 
                                   colour = site, 
                                   group = site_veg_year,
                                   linetype = veg_type,
                                   shape = veg_type)) + 
  geom_point(size = 3) + 
  geom_smooth(method = "lm", se = F) +
  labs(x = "Day of year\n",
       y =  "\nMean NDVI",
       shape = "Vegetation Type",
       linetype = "Vegetation Type",
       colour = "Site") +
  scale_y_continuous(limits = c(0.4,0.8), breaks = seq(0.4, 0.8, 0.1)) +
  scale_x_continuous(limits = c(170,230), breaks = seq(170, 230, 10)) +
  scale_shape_manual(values = c(16, 17),
                     labels = c("Tussock Sedge", "Dryas-Vetch"))+
  scale_color_manual(values = colour_scale_sites,
                     labels = c("Area 1", "Area 2", "Area 3", "Area 4")) +
  scale_linetype_manual(values = c(1,2),
                        labels = c("Tussock Sedge", "Dryas-Vetch")) +
  guides(color = guide_legend(order = 1,
                              title = NULL),
         linetype = guide_legend(order = 2,
                                 title = NULL,
                                 override.aes = list(color = "black")),
         shape = guide_legend(order = 3, title = NULL)) +
  theme_cowplot(18) +
  theme(legend.key.width=unit(0.455,"inch"))

save_plot("figures/fig_5_ground_based_phenology/sent_com_doy_vs_ndvi_plot.png",
          com_doy_vs_ndvi_plot_sent,
          base_height = 5,
          base_aspect_ratio = 1.35)


