# Phenology Time-Series Drone and Sentinel correlation (Fig 2 Panel a)
# Jakob Assmann j.assmann@ed.ac.uk 4 October 2018
# Updated 15 April 2020

### Preparations ----
# Dependencies
library(dplyr)
library(ggplot2)
library(raster)
library(rasterVis)
library(viridisLite)
library(MCMCglmm)
library(cowplot)

# Set global parameters / load site boundaries and meta data
data_out_path <- "data/fig_2_drone_sent_cor/"
figure_out_path <- "figures/"
site_boundaries <- read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
load("data/meta_data.Rda") # (created by pre_gather_meta_data.R)

### 1) Extract pixel by pixel sentinel and drone data ----

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
PS1_HER_extent <- get_sent_extent("PS1_HER", site_boundaries)
PS1_KOM_extent <- get_sent_extent("PS1_KOM", site_boundaries)
PS2_HER_extent <- get_sent_extent("PS2_HER", site_boundaries)
PS2_KOM_extent <- get_sent_extent("PS2_KOM", site_boundaries)
PS3_HER_extent <- get_sent_extent("PS3_HER", site_boundaries)
PS3_KOM_extent <- get_sent_extent("PS3_KOM", site_boundaries)
PS4_HER_extent <- get_sent_extent("PS4_HER", site_boundaries)
PS4_KOM_extent <- get_sent_extent("PS4_KOM", site_boundaries)

# Set Drone - Sentinel Date combos
dro_sent_combos <- data.frame(
  drone_date = as.Date(c("2017-07-17",
                        "2017-07-28",
                        "2017-07-09",
                        "2017-07-11",
                        "2017-07-18",
                        "2017-06-26",
                        "2017-07-06",
                        "2017-07-26",
                        "2017-08-07")),
  sentinel_date = as.Date(c("2017-07-17",
                            "2017-07-28",
                            "2017-07-08",
                            "2017-07-12",
                            "2017-07-19",
                            "2017-06-28",
                            "2017-07-04",
                            "2017-07-28",
                            "2017-08-05"))
)
dro_sent_combos$diff <- 
  dro_sent_combos$sentinel_date - dro_sent_combos$drone_date

# Set up empty data frame for output
pixel_combos <- data.frame(
  combo_id = as.character(NA), 
  drone_date = as.Date(NA), 
  sentinel_date = as.Date(NA), 
  diff = as.difftime("0", "%d"), 
  sentinel_id = as.character(NA), 
  site_name = as.character(NA), 
  veg_type = as.character(NA), 
  pixel_id = as.character(NA), 
  drone_ndvi = as.double(NA), 
  sentinel_ndvi = as.double(NA))

# function to calculate NDVI
NDVI <- function(red_band, nir_band) {
  (nir_band - red_band) / (nir_band + red_band) 
  }

# define funciton to extract ndvi for a single combo
extract_ndvi <- function(site_name_en, 
                         veg_type_en, 
                         drone_sites, 
                         sentinel_brick){
  # grab dates
  drone_date <- drone_sites$date[1]
  sentinel_date <- sentinel_brick$date[1]
  
  # Load Drone bands
  drone_bands <- drone_sites %>% 
    filter(site_name == site_name_en, veg_type == veg_type_en)
  list2env(
    lapply(
      setNames(drone_bands$file_path, 
               make.names(drone_bands$object_name)),
      raster), 
    envir = .GlobalEnv)
  # Load Sentinel Brick
  list2env(
    lapply(
      setNames(sentinel_brick$file_path, 
               make.names(sentinel_brick$object_name)),
      brick), 
    envir = .GlobalEnv)
  
  # Crop rasters and brick
  list2env(
    lapply(
      setNames(drone_bands$object_name, 
               make.names(paste0(drone_bands$object_name, "_cropped"))),
      function(x) { 
        crop(get(x), get(paste0(site_name_en, "_", veg_type_en, "_extent")))
      }), 
    envir = .GlobalEnv)
  list2env(
    lapply(
      setNames(sentinel_brick$object_name, 
               make.names(paste0(sentinel_brick$object_name, "_cropped"))),
      function(x) { 
        crop(get(x), get(paste0(site_name_en, "_", veg_type_en, "_extent")))
      }), 
    envir = .GlobalEnv)
  
  # aggregate and resample drone rasters to sentinel
  list2env(
    lapply(
      setNames(paste0(drone_bands$object_name, "_cropped"), 
               make.names(paste0(drone_bands$object_name, "_cropped_resamp"))),
      function(x) { 
        resample(get(x),
                 get(paste0(sentinel_brick$object_name, "_cropped")),
                 method = "bilinear")
      }), 
    envir = .GlobalEnv)
  
  # prep name vectors for drone band objects
  red_band <- paste0(drone_bands[drone_bands$band == "RED",]$object_name, 
                     "_cropped_resamp")
  nir_band <- paste0(drone_bands[drone_bands$band == "NIR",]$object_name, 
                     "_cropped_resamp")
  
  # calculated NDVI
  drone_ndvi <- NDVI(get(red_band), get(nir_band))
  sentinel_ndvi <- NDVI(get(paste0(sentinel_brick[1,]$object_name, 
                                   "_cropped"))[[3]],
                        get(paste0(sentinel_brick[1,]$object_name, 
                                   "_cropped"))[[4]])
  #extract NDVI value pairs
  ndvi_value_pairs <- data.frame(
    combo_id = paste0("d", drone_date, '_s', sentinel_date),
    drone_date = drone_date,
    sentinel_date = sentinel_date,
    sentinel_id = sentinel_brick[1,]$sensor_id,
    site_veg = paste0(site_name_en, veg_type_en),
    site_name = site_name_en,
    veg_type = veg_type_en,
    pixel_id = paste0(site_name_en, "_", veg_type_en, "_x", 
                      rep(seq(1,10),10), "y", sort(rep(seq(1,10),10))),
    drone_ndvi = getValues(drone_ndvi),
    sentinel_ndvi = getValues(sentinel_ndvi),
    stringsAsFactors = F
    )
  
  # Tidy up
  rm(list = c(drone_bands$object_name,
              sentinel_brick$object_name,
              paste0(drone_bands$object_name, "_cropped"),
              paste0(sentinel_brick$object_name, "_cropped"),
              paste0(drone_bands$object_name, "_cropped_resamp")
  ), 
  envir = .GlobalEnv)
  return(ndvi_value_pairs)
}


extract_combo <- function(drone_date, sentinel_date){
  # first identify sites that match for combo
  drone_sites <- meta_data %>% 
    filter(sensor_id == "Drone", 
           date == drone_date, 
           (band == "RED" | band == "NIR"))
  
  # and the matching sentinel brick
  sentinel_brick <- meta_data %>% 
    filter(sensor_id == "Sentinel 2A" | sensor_id == "Sentinel 2B", 
           date == sentinel_date)
  
  # Extract pixel by pixel ndvi pairs for all site_name veg_type combos
  combo_ndvi_list <- mapply(function(x,y) {
    extract_ndvi(x, y,  drone_sites, sentinel_brick) 
    },
    drone_sites$site_name, drone_sites$veg_type, SIMPLIFY = F)
 # collapse list into data frame
  combo_ndvi_df <- bind_rows(combo_ndvi_list)
  # return data frame
  return(combo_ndvi_df)
}

# Extract NDVI for all combos and create day difference column in final df
pixel_combos_list <- mapply(function(x,y,z) {
  combo <- extract_combo(x, y)
  combo$diff<- z
  return(combo)
  }, 
  dro_sent_combos$drone_date, 
  dro_sent_combos$sentinel_date, 
  dro_sent_combos$diff, SIMPLIFY = F)
# collapse list into dataframe 
pixel_combos <- bind_rows(pixel_combos_list)
# tidy up
rm(pixel_combos_list)

# Save to Rda file
save(pixel_combos, file = paste0(data_out_path, "pixel_combos.Rda"))

load(paste0(data_out_path, "pixel_combos.Rda"))

#### 2) Garphical exploration of relationship -----

# PS4_HER on the 17 July 2017 is a big outlier let's remove it!
pixel_combos$site_veg_date <- paste0(pixel_combos$site_veg, "_", pixel_combos$drone_date)
pixel_combos <- filter(pixel_combos, site_veg_date != "PS4HER_2017-07-17")

# Plots to get to know the data
ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = veg_type)) + 
  geom_point()

ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = site_veg_date)) + 
  geom_point()

ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, y = drone_ndvi, colour = site_veg)) +
  geom_point()
ggplot(filter(pixel_combos, site_name == "PS1"), 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = site_veg_date)) + 
  geom_point()

ggplot(filter(pixel_combos, site_name == "PS2"), 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = site_veg_date)) + 
  geom_point()

ggplot(filter(pixel_combos, site_name == "PS3"), 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = site_veg_date)) + 
  geom_point()
ggplot(filter(pixel_combos, site_name == "PS4"), 
       mapping = aes(x = sentinel_ndvi,
                     y = drone_ndvi, 
                     colour = site_veg_date)) + geom_point()
ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, colour = combo_id)) +
  geom_point()

ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi,
                     colour = sentinel_id)) +
  geom_point()

ggplot(pixel_combos, 
       mapping = aes(x = sentinel_ndvi, 
                     y = drone_ndvi, 
                     colour = factor(diff))) + 
  geom_point()

# looks nice... almost like
# slopes are consistent within sites, satellites
# one oultier is the PS3 HER early season site, where the slope seem ing
# to be slightly different (maybe not a linear repsonse?)
# difference in days seems to have a clear influence on the intercept

##### 3) Linear Models ----

# Define helper function to process model outputs
mcmc_output_to_table <- function(mcmc_model, file_name, digits = 4){
  solutions <- data.frame(
    effect = row.names(as.data.frame(summary(mcmc_model)$solutions)),
    as.data.frame(round(summary(mcmc_model)$solutions, digits)))
  row.names(solutions) <- NULL
  solutions$eff.samp <- round(solutions$eff.samp)
  colnames(solutions) <- c("Effect Name", "Effect Size", "Lower 95% CI", "Upper 95% CI",
                           "Effective Sample Size", "pMCMC")
  random <- data.frame(summary(mcmc_model)$Rcovariances)
  random <- round(random, digits+1)
  random$name <- rownames(random)
  random$pMCMC <- " "
  random <- random[, c("name", "post.mean", "l.95..CI", "u.95..CI", "eff.samp", "pMCMC")]  
  random$eff.samp <- round(random$eff.samp)
  colnames(random) <- names(solutions)
  rownames(random) <- NULL
  
  model_formula <- paste(as.character(summary(mcmc_model)$fixed.formula)[2], 
                         summary(mcmc_model)$fixed.formula[1], 
                         summary(mcmc_model)$fixed.formula[3])
  final_table <- rbind(
    setNames(rep(" ", length(names(solutions))), names(solutions)),
    setNames(c(model_formula, rep(" ", length(names(solutions))-1)), names(solutions)),
    setNames(rep(" ", length(names(solutions))), names(solutions)),
    setNames(as.character(names(solutions)), as.character(names(solutions))),
    data.frame(apply(solutions, 2, as.character), stringsAsFactors = F),
    setNames(rep(" ", length(names(solutions))), names(solutions)))
  colnames(final_table) <- colnames(solutions)
  final_table = rbind(final_table, random)
  write.csv(file = file_name, final_table, row.names = F, col.names = F)
}


# Take random subsample of 10% for each sentinel / drone data combo
set.seed(5)
pixel_combos_sample <- pixel_combos %>% group_by(combo_id) %>% sample_frac(0.1)

# Model simple relationship
sentinel_drone_model_simple <- MCMCglmm(sentinel_ndvi ~ drone_ndvi, 
                                 data = as.data.frame(pixel_combos_sample),
                                 nitt = 50000,
                                 pr = T) 
summary(sentinel_drone_model_simple)
save(sentinel_drone_model_simple, file = paste0(data_out_path, "model_simple.rda"))
#load(paste0(data_out_path, "model_simple.rda"))
mcmc_output_to_table(sentinel_drone_model_simple, 
                     paste0(data_out_path, "model_simple_table.csv"))

# Claculate R2 - based on https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model
# Variance in fitted values
var_fitted <- var(as.vector(apply(sentinel_drone_model_simple$Sol,2,mean) %*% 
                              t(sentinel_drone_model_simple$X)))
# Caluclate marginal r2 
r2 <- var_fitted/(var_fitted+sum(apply(sentinel_drone_model_simple$VCV,2,mean)))

# Let's test whether there is an effect of doy, sentinel, 
# veg_type + interaction and on the relationship
sentinel_drone_model_full <- MCMCglmm(
  sentinel_ndvi ~ drone_ndvi + 
    diff + veg_type + sentinel_id + 
    diff:drone_ndvi + veg_type:drone_ndvi + sentinel_id:drone_ndvi, 
  data = as.data.frame(pixel_combos_sample),
  nitt = 50000,
  pr = T) 
summary(sentinel_drone_model_full)
save(sentinel_drone_model_full, file = paste0(data_out_path, "model_full.rda"))
#load(paste0(data_out_path, "model_full.rda"))
mcmc_output_to_table(sentinel_drone_model_full, 
                     paste0(data_out_path, "model_full_table.csv"))

# Claculate R2 - based on https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model
# Variance in fitted values
var_fitted <- var(as.vector(apply(sentinel_drone_model_full$Sol,2,mean) %*% 
                              t(sentinel_drone_model_full$X)))
# Caluclate marginal r2 
r2 <- var_fitted/(var_fitted+sum(apply(sentinel_drone_model_full$VCV,2,mean)))


# The date difference intercept is the only tested variable that does not have
# a significant effect on the relationship, however it affects the slope 
# significantly

##### 4) Final plot ---- 

## Pretty plot for publication

# Extract model predictions and confidence intervals for simple model
#load(paste0(data_out_path, "model_simple.rda"))

preds <- predict.MCMCglmm(sentinel_drone_model_simple, 
                          interval = "confidence", 
                          type = "response")
preds <- cbind(preds, 
               data.frame(
                 drone_ndvi = sentinel_drone_model_simple$X[,2],
                 stringsAsFactors = F))

sentinel_drone_plot <- ggplot(data = pixel_combos_sample, 
                              mapping = aes(x = drone_ndvi, 
                                            y = sentinel_ndvi)) +
  geom_abline(intercept = 0, slope = 1, 
              color="black", linetype="dashed", 
              size= 1, alpha = 1 ) +
  geom_point(size = 2, color = "black") +
  geom_line(data = preds,
            mapping = aes(x = drone_ndvi,
                          y = fit
            ),
            colour = "blue",
            inherit.aes = FALSE,
            size = 1) +
  geom_ribbon(data = preds,
              mapping = aes(x = drone_ndvi,
                            ymin = lwr,
                            ymax = upr
              ),
              fill = "blue",
              alpha = 0.5,
              inherit.aes = FALSE) + 
  scale_x_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_y_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_colour_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_fill_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  xlab('Drone NDVI')+
  ylab("Sentinel NDVI") +
  theme_cowplot(18)

save_plot(paste0(figure_out_path, "fig_2_panel_a/drone_sent_cor.png"), 
       plot = sentinel_drone_plot, base_aspect_ratio =  1.3)
       

### EOF 