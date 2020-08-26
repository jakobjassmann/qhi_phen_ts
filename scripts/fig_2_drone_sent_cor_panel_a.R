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
library(pbmcapply)

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

# Identify Drone - Sentinel Date combos
meta_data$doy <- as.numeric(format.Date(meta_data$date, "%j"))
meta_data$year <- as.numeric(format.Date(meta_data$date, "%Y"))
meta_data_drone <- meta_data %>% 
  filter(sensor_id == "Drone", 
         band == "NDVI", 
         flight_id != "PS4_HER_20170717",
         !grepl("nocalib", meta_data$file_path)) 

meta_data_drone$doy_plus2 <- meta_data_drone$doy + 2
meta_data_drone$doy_minus2 <- meta_data_drone$doy - 2

# Grab distinct sentinel meta scene dates
meta_data_sentinel_distinct <- meta_data %>%
  filter(sensor_id == "Sentinel 2A" |
           sensor_id == "Sentinel 2B") %>% 
  distinct(date , .keep_all = T)

sentinel_scene_combos <- bind_rows(
  lapply(meta_data_drone$flight_id,
         function(flight){
           cat(paste0(flight, "\n"))
           doy_minus2 <- meta_data_drone[meta_data_drone$flight_id == flight,
                                         c("doy_minus2")]
           doy_plus2 <- meta_data_drone[meta_data_drone$flight_id == flight,
                                        c("doy_plus2")]
           year <- meta_data_drone[meta_data_drone$flight_id == flight,
                                   c("year")]
           sentinel_scenes <- meta_data_sentinel_distinct[
             meta_data_sentinel_distinct$year == rep(year, nrow(meta_data_sentinel_distinct)) & 
               meta_data_sentinel_distinct$doy >= rep(doy_minus2, nrow(meta_data_sentinel_distinct)) &
               meta_data_sentinel_distinct$doy <= rep(doy_plus2, nrow(meta_data_sentinel_distinct)),]
           
           cat(paste0("Sentinel scenes:", nrow(sentinel_scenes), "\n"))
           
           if(nrow(sentinel_scenes) == 1){
             sentinel_scenes$site_veg <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                         c("site_veg")] 
             sentinel_scenes$site_name <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                          c("site_name")] 
             sentinel_scenes$veg_type <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                         c("veg_type")] 
             scene_combos <- 
               bind_cols(
                 setNames(sentinel_scenes, paste0("sentinel_",
                                                  names(sentinel_scenes))),
                 setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                          paste0("drone_", 
                                 names(meta_data_drone))))
           } else if(nrow(sentinel_scenes) > 1) {
             sentinel_scenes$site_veg <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                         c("site_veg")] 
             sentinel_scenes$site_name <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                          c("site_name")] 
             sentinel_scenes$veg_type <- meta_data_drone[meta_data_drone$flight_id == flight,
                                                         c("veg_type")] 
             scene_combos <- setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                                      paste0("drone_", 
                                             names(meta_data_drone)))
             scene_combos <- do.call("rbind", replicate(nrow(sentinel_scenes), 
                                                        scene_combos, 
                                                        simplify = FALSE))
             scene_combos <- bind_cols(
               setNames(sentinel_scenes, paste0("sentinel_",
                                                names(sentinel_scenes))),
               scene_combos)
           } else if(nrow(sentinel_scenes) == 0) {
             scene_combos <- setNames(meta_data_drone[meta_data_drone$flight_id == flight, 
                                                      1:ncol(sentinel_scenes)],
                                      paste0("sentinel_", 
                                             names(sentinel_scenes)))
             scene_combos[,1:ncol(scene_combos)] <- NA 
             scene_combos <- bind_cols(scene_combos, 
                                       setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                                                paste0("drone_", 
                                                       names(meta_data_drone))))
             
           }
           return(scene_combos)
         }))
sentinel_scene_combos %>% na.omit() %>% 
  distinct(drone_date, sentinel_date)

# Check how many drone flights do not have a ls 8 scene
sentinel_scene_combos %>% filter(is.na(sentinel_flight_id)) %>%
  distinct(drone_flight_id) %>%
  summarise(n = n())
# 47 Drone scenes don't have a matching Sentinel scene

# Throw out nas
sentinel_scene_combos <- sentinel_scene_combos %>% filter(!is.na(sentinel_flight_id))

# Filter drone flights with more than two matching ls8 scenes
sentinel_scene_combos %>% 
  na.omit() %>%
  group_by(drone_flight_id) %>% 
  summarise(n = n()) %>%
  filter(n > 1)
sentinel_scene_combos %>% na.omit() %>% 
  distinct(drone_date, sentinel_date) %>%
  arrange(drone_date)
# Multiple flights have multiple combinatins, select those ones that are
# unique
sentinel_scene_combos <- filter(sentinel_scene_combos,
                             !(drone_date == as.Date("2017-07-06") &
                                 sentinel_date == as.Date("2017-07-08")) &
                             !(drone_date == as.Date("2017-07-10") &
                                 sentinel_date == as.Date("2017-07-08")) &
                               !(drone_date == as.Date("2017-07-10") &
                                   sentinel_date == as.Date("2017-07-12")) &
                             !(drone_date == as.Date("2017-07-17") &
                                 sentinel_date == as.Date("2017-07-19")) &
                              !(drone_date == as.Date("2017-07-18") &
                                  sentinel_date == as.Date("2017-07-17")) )
# check that worked
sentinel_scene_combos %>% 
  na.omit() %>%
  group_by(drone_flight_id) %>% 
  summarise(n = n()) %>%
  filter(n > 1)

# Check how many scene combos exist
nrow(sentinel_scene_combos)
# 35 -> that's not bad.

# For legacy reasons this was previously set manually.
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
# The dro_sent_combos and sentinel_scene_combos distinct dates are identical.
# distinct_combos <- sentinel_scene_combos %>% na.omit() %>% 
#   distinct(drone_date, sentinel_date) %>%
#   arrange(drone_date) %>%
#   as.data.frame() 
# distinct_combos$drone_date == arrange(dro_sent_combos, drone_date)$drone_date
# distinct_combos$sentinel_date == arrange(dro_sent_combos, drone_date)$sentinel_date

# We keep using the dro_sent_combos data frame for future reasons. 
rm(sentinel_scene_combos)

# Add date diff colum to the dro_sent_combos data frame
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

#load(paste0(data_out_path, "pixel_combos.Rda"))

#### 2) Drone-Sentinel Garphical exploration of relationship -----

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

##### 3) Drone-Sentinel - Linear Models ----

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
# Get number of pxiels for stats
nrow(pixel_combos_sample)
# 740 

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

##### 4) Drone-Sentinel - Final plot ---- 

## Pretty plot for publication

# Set colours
her_col <- "#1e9148FF"
kom_col <- "#1e5c91FF"

# Extract model predictions and confidence intervals for simple model
#load(paste0(data_out_path, "model_simple.rda"))

preds <- predict.MCMCglmm(sentinel_drone_model_simple, 
                          interval = "confidence", 
                          type = "response")
preds <- cbind(preds, 
               data.frame(
                 drone_ndvi = sentinel_drone_model_simple$X[,2],
                 stringsAsFactors = F))

# Specify difference in days for plotting (little = big, a lot = small)
pixel_combos_sample$diff_plot <- abs(pixel_combos_sample$diff)
max(pixel_combos_sample$diff_plot)
pixel_combos_sample$diff_plot[pixel_combos_sample$diff_plot == 0] <- 3
pixel_combos_sample$diff_plot[pixel_combos_sample$diff_plot == 2] <- -5
pixel_combos_sample$diff_plot[pixel_combos_sample$diff_plot == 1] <- 2
pixel_combos_sample$diff_plot[pixel_combos_sample$diff_plot == -5] <- 1

sentinel_drone_plot <- ggplot(data = pixel_combos_sample, 
                              mapping = aes(x = drone_ndvi, 
                                            y = sentinel_ndvi,
                                            colour = veg_type,
                                            shape = sentinel_id,
                                            size = diff_plot)) +
  geom_abline(intercept = 0, slope = 1, 
              color="black", linetype="dashed", 
              size= 1, alpha = 1 ) +
  geom_point() +
  scale_shape_manual(values = c(21,22)) +
  geom_ribbon(data = preds,
              mapping = aes(x = drone_ndvi,
                            ymin = lwr,
                            ymax = upr
              ),
              fill = "black",
              alpha = 0.4,
              inherit.aes = FALSE) + 
  geom_line(data = preds,
            mapping = aes(x = drone_ndvi,
                          y = fit
            ),
            colour = "black",
            inherit.aes = FALSE,
            size = 1) +
  scale_size(range = c(2, 4)) +
  scale_x_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_y_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_colour_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_fill_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  annotate("text", x = 0.65, y = 0.3, 
           label = "  Dryas-vetch Tundra", hjust = 0,
           colour = kom_col, size = 4) +
  annotate("text", x = 0.65, y = 0.34, 
           label = "  Tussock Sedge Tundra", hjust = 0,
           colour = her_col, size = 4) +
  annotate("point", x = 0.65, y = 0.3, 
           size = 3, shape = 21,
           colour = kom_col) +
  annotate("point", x = 0.65, y = 0.34, 
           size = 3, shape = 21,
           colour = her_col) +
  annotate("text", x = 0.65, y = 0.40, 
           label = "  Sentinel-2 B", hjust = 0,
           colour = "black", size = 4) +
  annotate("text", x = 0.65, y = 0.44, 
           label = "  Sentinel-2 A", hjust = 0,
           colour = "black", size = 4) +
  annotate("point", x = 0.65, y = 0.40, 
           size = 3, shape = 22,
           colour = "black") +
  annotate("point", x = 0.65, y = 0.44, 
           size = 3, shape = 21,
           colour = "black") +
  annotate("text", x = 0.3, y = 0.9, 
           label = "   0 days difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("text", x = 0.3, y = 0.86, 
           label = "   1 day difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("text", x = 0.3, y = 0.82, 
           label = "   2 days difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("point", x = 0.3, y = 0.9, 
           size = 4, shape = 21,
           colour = "black") +
  annotate("point", x = 0.3, y = 0.86, 
           size = 3, shape = 21,
           colour = "black") +
  annotate("point", x = 0.3, y = 0.82, 
           size = 2, shape = 21,
           colour = "black") +
  xlab('Drone NDVI')+
  ylab("Sentinel NDVI") +
  theme_cowplot(15) +
  theme(legend.position = "none")
sentinel_drone_plot
save_plot(paste0(figure_out_path, "fig_2_panel_a/drone_sent_cor.png"), 
       plot = sentinel_drone_plot, base_aspect_ratio =  1.3)
       

#### 5) Drone-Landsat 8 - Extract pxiel by pixel drone ----

# Load Landsat8 meta data
load("data/landsat8/meta_data_ls8_with_mean.Rda")
meta_data_ls8_distinct <- meta_data_ls8_with_mean %>%
  distinct(flight_id, file_path, date) %>%
  mutate(doy = as.numeric(format.Date(date, "%j")),
         year = as.numeric(format.Date(date, "%Y")))

# Identify drone landsat combinations with two day buffer
meta_data$doy <- as.numeric(format.Date(meta_data$date, "%j"))
meta_data$year <- as.numeric(format.Date(meta_data$date, "%Y"))
meta_data_drone <- meta_data %>% 
  filter(sensor_id == "Drone", 
         band == "NDVI", 
         flight_id != "PS4_HER_20170717",
         !grepl("nocalib", meta_data$file_path)) 

meta_data_drone$doy_plus2 <- meta_data_drone$doy + 2
meta_data_drone$doy_minus2 <- meta_data_drone$doy - 2

# Find matching landsat scenes for each flight
ls8_scene_combos <- bind_rows(lapply(meta_data_drone$flight_id,
       function(flight){
         cat(paste0(flight, "\n"))
         doy_minus2 <- meta_data_drone[meta_data_drone$flight_id == flight,
                                        c("doy_minus2")]
         doy_plus2 <- meta_data_drone[meta_data_drone$flight_id == flight,
                                       c("doy_plus2")]
         year <- meta_data_drone[meta_data_drone$flight_id == flight,
                                 c("year")]
         landsat8_scenes <- meta_data_ls8_distinct[
           meta_data_ls8_distinct$year == rep(year, nrow(meta_data_ls8_distinct)) & 
             meta_data_ls8_distinct$doy >= rep(doy_minus2, nrow(meta_data_ls8_distinct)) &
             meta_data_ls8_distinct$doy <= rep(doy_plus2, nrow(meta_data_ls8_distinct)),]
         cat(paste0("landsat 8 scenes:", nrow(landsat8_scenes), "\n"))
         if(nrow(landsat8_scenes) == 1){
           scene_combos <- 
             bind_cols(
               setNames(landsat8_scenes, paste0("landsat8_",
                                                names(landsat8_scenes))),
               setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                        paste0("drone_", 
                               names(meta_data_drone))))
           } else if(nrow(landsat8_scenes) > 1) {
             scene_combos <- setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                                      paste0("drone_", 
                                             names(meta_data_drone)))
             scene_combos <- do.call("rbind", replicate(nrow(landsat8_scenes), 
                                                        scene_combos, 
                                                        simplify = FALSE))
             scene_combos <- bind_cols(
               setNames(landsat8_scenes, paste0("landsat8_",
                                                names(landsat8_scenes))),
               scene_combos)
           } else if(nrow(landsat8_scenes) == 0) {
             scene_combos <- setNames(meta_data_drone[meta_data_drone$flight_id == flight, 
                                                      1:ncol(landsat8_scenes)],
                                  paste0("landsat8_", 
                                         names(landsat8_scenes)))
             scene_combos[,1:ncol(scene_combos)] <- NA 
             scene_combos <- bind_cols(scene_combos, 
                         setNames(meta_data_drone[meta_data_drone$flight_id == flight,],
                                  paste0("drone_", 
                                         names(meta_data_drone))))
             
           }
         return(scene_combos)
       }))

# Check how many drone flights do not have a ls 8 scene
ls8_scene_combos %>% filter(is.na(landsat8_flight_id)) %>%
  distinct(drone_flight_id) %>%
  summarise(n = n())
# 56 Drone scenes don't have a matching LS8 scene

# Throw out nas
ls8_scene_combos <- na.omit(ls8_scene_combos)

# Filter drone flights with more than two matching ls8 scenes
ls8_scene_combos %>% 
  na.omit() %>%
  group_by(drone_flight_id) %>% 
  summarise(n = n()) %>%
  filter(n > 1)
# Two flights on 24 June 2017 which have a landsat scene on the day as well as
# two days before. Remove the two days before scene
ls8_scene_combos <- filter(ls8_scene_combos,
                           !(drone_date == as.Date("2017-06-24") &
                             landsat8_date == as.Date("2017-06-22")))
# check that worked
ls8_scene_combos %>% 
  na.omit() %>%
  group_by(drone_flight_id) %>% 
  summarise(n = n()) %>%
  filter(n > 1)

# Check how many scene combos exist
nrow(ls8_scene_combos)
# 35 -> that's not bad.

# Create day difference colum
ls8_scene_combos$diff <- ls8_scene_combos$drone_doy -
  ls8_scene_combos$landsat8_doy
ls8_scene_combos$diff
mean(abs(ls8_scene_combos$diff))

# define funciton to extract ndvi for a single combo
ls8_extract_ndvi <- function(landsat8_scene_id, 
                         drone_flight_id){
  cat(drone_flight_id)
  # grab dates
  drone_date <- ls8_scene_combos$drone_date[
    ls8_scene_combos$drone_flight_id == drone_flight_id]
  landsat8_date <- ls8_scene_combos$landsat8_date[
    ls8_scene_combos$landsat8_flight_id == landsat8_scene_id]
  landsat8_date <- landsat8_date[1]
  
  # and site veg
  site_veg <- ls8_scene_combos$drone_site_veg[
    ls8_scene_combos$drone_flight_id == drone_flight_id]
  
  # Load bands
  drone_bands <- meta_data %>% 
    filter(flight_id == drone_flight_id,
           band %in% c("RED", "NIR"))
  list2env(
    lapply(
      setNames(drone_bands$file_path, 
               make.names(drone_bands$object_name)),
      raster), 
    envir = .GlobalEnv)
  # Load Landsat8 NDVI raster
  landsat8_ndvi_raster <- raster(unique(ls8_scene_combos$landsat8_file_path[
    ls8_scene_combos$landsat8_flight_id == landsat8_scene_id]))
  
  # Crop landsat raster to site extent
  site_veg_extent <- as(get(paste0(site_veg, "_extent")), "SpatialPolygons")
  crs(site_veg_extent) <- crs(get(drone_bands$object_name[1]))
  # Transform to Landsat 8 crs if needed
  if(as.character(crs(site_veg_extent)) != as.character(crs(landsat8_ndvi_raster))) {
  site_veg_extent <- spTransform(site_veg_extent, crs(landsat8_ndvi_raster))
  }
  landsat8_ndvi_raster_cropped <- 
        crop(landsat8_ndvi_raster, 
             site_veg_extent,
             snap = "in")
  
  # aggregate and resample drone rasters to landsat8
  list2env(
    lapply(
      setNames(drone_bands$object_name, 
               make.names(paste0(drone_bands$object_name, "_cropped_resamp"))),
      function(x) { 
        if(as.character(as.character(crs(get(x))) != 
                        as.character(crs(landsat8_ndvi_raster_cropped)))){
          output_raster <- projectRaster(get(x), landsat8_ndvi_raster_cropped)
        } else {
          output_raster <- resample(get(x),
                  landsat8_ndvi_raster_cropped,
                  method = "bilinear")
        }
        return(output_raster)
      }), 
    envir = .GlobalEnv)
  
  # prep name vectors for drone band objects
  red_band <- paste0(drone_bands[drone_bands$band == "RED",]$object_name, 
                     "_cropped_resamp")
  nir_band <- paste0(drone_bands[drone_bands$band == "NIR",]$object_name, 
                     "_cropped_resamp")
  
  # calculated NDVI
  drone_ndvi <- NDVI(get(red_band), get(nir_band))

  #extract NDVI value pairs
  ndvi_value_pairs <- data.frame(
    combo_id = paste0("d", drone_date, '_s', landsat8_date),
    drone_date = as.Date(drone_date),
    landsat8_date = as.Date(landsat8_date),
    landsat8_id = landsat8_scene_id,
    site_veg = site_veg,
    site_name = substr(site_veg, 1, 3),
    veg_type = substr(site_veg, 5, 7),
    pixel_id = paste0(site_veg, "_",
                      expand.grid(paste0("x", 1:ncol(landsat8_ndvi_raster_cropped)),
                                  paste0("y", 1:nrow(landsat8_ndvi_raster_cropped)))[,1],
                      expand.grid(paste0("x", 1:ncol(landsat8_ndvi_raster_cropped)),
                                         paste0("y", 1:nrow(landsat8_ndvi_raster_cropped)))[,2]),
    drone_ndvi = getValues(drone_ndvi),
    landsat8_ndvi = getValues(landsat8_ndvi_raster_cropped),
    stringsAsFactors = F
  )
  
  # Tidy up
  rm(list = c(drone_bands$object_name,
              paste0(drone_bands$object_name, "_cropped"),
              paste0(drone_bands$object_name, "_cropped_resamp")
  ), 
  envir = .GlobalEnv)
  return(ndvi_value_pairs)
}

# turn id filed into character
ls8_scene_combos$landsat8_flight_id <- as.character(ls8_scene_combos$landsat8_flight_id)
ls8_scene_combos$drone_flight_id <- as.character(ls8_scene_combos$drone_flight_id)

# Extract NDVI for all combos and create day difference column in final df
ls8_pixel_combos <- pbmcmapply(ls8_extract_ndvi,
                           ls8_scene_combos$landsat8_flight_id,
                           ls8_scene_combos$drone_flight_id,
                           SIMPLIFY = F,
                           mc.cores = 4)
ls8_pixel_combos <- bind_rows(ls8_pixel_combos)
ls8_pixel_combos$diff <- as.numeric(format.Date(ls8_pixel_combos$landsat8_date, "%j")) -
  as.numeric(format.Date(ls8_pixel_combos$drone_date, "%j"))
# Save to Rda file
save(ls8_pixel_combos, file = paste0(data_out_path, "ls8_pixel_combos.Rda"))

# Get number of pixels for stats
nrow(ls8_pixel_combos)
# 198

##### 6) Drone-Landsat 8 - Linear Models ----

# # Take random subsample of 10% for each sentinel / drone data combo
# set.seed(5)
# pixel_combos_sample <- pixel_combos %>% group_by(combo_id) %>% sample_frac(0.1)

# Model simple relationship
landsat8_drone_model_simple <- MCMCglmm(landsat8_ndvi ~ drone_ndvi, 
                                        data = as.data.frame(ls8_pixel_combos),
                                        nitt = 50000,
                                        pr = T) 
summary(landsat8_drone_model_simple)
save(landsat8_drone_model_simple, file = paste0(data_out_path, "ls8_model_simple.rda"))
#load(paste0(data_out_path, "model_simple.rda"))
mcmc_output_to_table(landsat8_drone_model_simple, 
                     paste0(data_out_path, "ls8_model_simple_table.csv"))

# Claculate R2 - based on https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model
# Variance in fitted values
var_fitted <- var(as.vector(apply(landsat8_drone_model_simple$Sol,2,mean) %*% 
                              t(landsat8_drone_model_simple$X)))
# Caluclate marginal r2 
r2 <- var_fitted/(var_fitted+sum(apply(landsat8_drone_model_simple$VCV,2,mean)))


# Let's test whether there is an effect of doy, sentinel, 
# veg_type + interaction and on the relationship
landsat8_drone_model_full <- MCMCglmm(
  landsat8_ndvi ~ drone_ndvi + 
    diff + veg_type  + 
    diff:drone_ndvi + veg_type:drone_ndvi, 
  data = as.data.frame(ls8_pixel_combos),
  nitt = 50000,
  pr = T) 
summary(landsat8_drone_model_full)
save(landsat8_drone_model_full, file = paste0(data_out_path, "ls8_model_full.rda"))
#load(paste0(data_out_path, "model_full.rda"))
mcmc_output_to_table(landsat8_drone_model_full, 
                     paste0(data_out_path, "ls8_model_full_table.csv"))

# Claculate R2 - based on https://www.researchgate.net/post/How_can_I_calculate_R2_for_an_Bayesian_MCMC_multilevel_model
# Variance in fitted values
var_fitted <- var(as.vector(apply(landsat8_drone_model_full$Sol,2,mean) %*% 
                              t(landsat8_drone_model_full$X)))
# Caluclate marginal r2 
r2 <- var_fitted/(var_fitted+sum(apply(landsat8_drone_model_full$VCV,2,mean)))


# The date difference intercept is the only tested variable that does not have
# a significant effect on the relationship, however it affects the slope 
# significantly

##### 4) Drone-Landsat 8 - Final plot ---- 

## Pretty plot for publication

# Set colours
her_col <- "#1e9148FF"
kom_col <- "#1e5c91FF"

# Extract model predictions and confidence intervals for simple model
#load(paste0(data_out_path, "model_simple.rda"))

ls8_preds <- predict.MCMCglmm(landsat8_drone_model_simple, 
                          interval = "confidence", 
                          type = "response")
ls8_preds <- cbind(ls8_preds, 
               data.frame(
                 drone_ndvi = landsat8_drone_model_simple$X[,2],
                 stringsAsFactors = F))

# Specify difference in days for plotting (little = big, a lot = small)
ls8_pixel_combos$diff_plot <- abs(ls8_pixel_combos$diff)
max(ls8_pixel_combos$diff_plot)
ls8_pixel_combos$diff_plot[ls8_pixel_combos$diff_plot == 0] <- 3
ls8_pixel_combos$diff_plot[ls8_pixel_combos$diff_plot == 2] <- -5
ls8_pixel_combos$diff_plot[ls8_pixel_combos$diff_plot == 1] <- 2
ls8_pixel_combos$diff_plot[ls8_pixel_combos$diff_plot == -5] <- 1

landsat8_drone_plot <- ggplot(data = ls8_pixel_combos, 
                              mapping = aes(x = drone_ndvi, 
                                            y = landsat8_ndvi,
                                            colour = veg_type,
                                            size = diff_plot)) +
  geom_abline(intercept = 0, slope = 1, 
              color="black", linetype="dashed", 
              size= 1, alpha = 1 ) +
  geom_point(shape = 21) +
  scale_shape_manual(values = c(21,22)) +
  geom_ribbon(data = ls8_preds,
              mapping = aes(x = drone_ndvi,
                            ymin = lwr,
                            ymax = upr
              ),
              fill = "black",
              alpha = 0.4,
              inherit.aes = FALSE) + 
  geom_line(data = ls8_preds,
            mapping = aes(x = drone_ndvi,
                          y = fit
            ),
            colour = "black",
            inherit.aes = FALSE,
            size = 1) +
  scale_size(range = c(2, 4)) +
  scale_x_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_y_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_colour_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_fill_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  annotate("text", x = 0.65, y = 0.3, 
           label = "  Dryas-vetch Tundra", hjust = 0,
           colour = kom_col, size = 4) +
  annotate("text", x = 0.65, y = 0.34, 
           label = "  Tussock Sedge Tundra", hjust = 0,
           colour = her_col, size = 4) +
  annotate("point", x = 0.65, y = 0.3, 
           size = 3, shape = 21,
           colour = kom_col) +
  annotate("point", x = 0.65, y = 0.34, 
           size = 3, shape = 21,
           colour = her_col) +
  annotate("text", x = 0.3, y = 0.9, 
           label = "   0 days difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("text", x = 0.3, y = 0.86, 
           label = "   1 day difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("text", x = 0.3, y = 0.82, 
           label = "   2 days difference", hjust = 0,
           colour = "black", size = 4) +
  annotate("point", x = 0.3, y = 0.9, 
           size = 4, shape = 21,
           colour = "black") +
  annotate("point", x = 0.3, y = 0.86, 
           size = 3, shape = 21,
           colour = "black") +
  annotate("point", x = 0.3, y = 0.82, 
           size = 2, shape = 21,
           colour = "black") +
  xlab('Drone NDVI')+
  ylab("Landsat 8 NDVI") +
  theme_cowplot(15) +
  theme(legend.position = "none")
landsat8_drone_plot
save_plot(paste0(figure_out_path, "fig_2_panel_a/drone_ls8_cor.png"), 
          plot = landsat8_drone_plot, base_aspect_ratio =  1.3)


### EOF 