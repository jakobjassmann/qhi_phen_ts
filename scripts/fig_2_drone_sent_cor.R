# Phenology Time-Series Drone and Sentinel correlation 
# Jakob Assmann j.assmann@ed.ac.uk 4 October 2018

# Dependencies
library(dplyr)
library(ggplot2)
library(raster)
library(rasterVis)
library(viridisLite)
library(MCMCglmm)

# Set global parameters / load site boundaries and meta data
data_out_path <- "data/fig_2_drone_sent_cor"
figure_out_path <- "figures"
site_boundaries <- read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
load("data/meta_data.Rda") # (created by pre_gather_meta_data.R)

### 1) Extract pixel by pixel sentinel and drone data ----
#####

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
  
  # resample drone rasters to sentinel
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
#####

## Plot relationship
#####

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

ggsave(paste0(figure_out_path, "sentinel_drone_plot_diff.png"), 
       width = 12, height = 8)

# looks nice... almost like
# slopes are consistent within sites, satellites
# one oultier is the PS3 HER early season site, where the slope seem to be slightly different (maybe not a linear repsonse?)
# difference in days has a clear influence on the intercept

#####

## Let's model it!
sentinel_drone_model_simple <- MCMCglmm(drone_ndvi~ sentinel_ndvi, 
                                 data = pixel_combos,
                                 nitt = 50000,
                                 pr = T) 
summary(sentinel_drone_model_simple)
save(sentinel_drone_model_diff, file = paste0(data_out_path, "model_simple.rda"))
load(paste0(data_out_path, "model_simple.rda"))

# Let's test whether there is an effect of doy, sentinel, veg_type + interaction and on the relationship
sentinel_drone_model_full <- MCMCglmm(
  drone_ndvi~ sentinel_ndvi + 
    diff + veg_type + sentinel_id + 
    diff:sentinel_ndvi + veg_type:sentinel_ndvi + sentinel_id:sentinel_ndvi, 
  data = pixel_combos,
  nitt = 50000,
  pr = T) 
summary(sentinel_drone_model_full)
save(sentinel_drone_model_full, file = paste0(data_out_path, "model_full.rda"))
load(paste0(data_out_path, "model_full.rda"))
out <- capture.output(summary(sentinel_drone_model_full))
cat(out, file=paste0(data_out_path, "model_full.txt"), sep = "\n", append=F)

sentinel_drone_model_veg <- MCMCglmm(drone_ndvi~ sentinel_ndvi + 
                                        veg_type + 
                                        veg_type:sentinel_ndvi, 
                                      data = pixel_combos,
                                      nitt = 50000,
                                      pr = T) 
summary(sentinel_drone_model_veg)
save(sentinel_drone_model_veg, file = paste0(data_out_path, "model_veg.rda"))
load(paste0(data_out_path, "model_veg.rda"))
out <- capture.output(summary(sentinel_drone_model_veg))
cat(out, file=paste0(data_out_path, "model_veg.txt"), sep="\n", append=F)


a <- 10000
pa_prior <- list(
  R=list(V=diag(1), nu=0.002),
  G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))
sentinel_drone_model_diff <- MCMCglmm(drone_ndvi~ sentinel_ndvi, 
                                 random = ~  diff, 
                                 data = pixel_combos,
                                 nitt = 50000,
                                 prior = pa_prior,
                                 pr = T) 
summary(sentinel_drone_model_diff)
plot(sentinel_drone_model_diff)
save(sentinel_drone_model_diff, file = paste0(data_out_path, "model_diff.rda"))

sentinel_drone_model_sentinel_id <- MCMCglmm(drone_ndvi~ sentinel_ndvi, 
                                      random = ~  sentinel_id, 
                                      data = pixel_combos,
                                      nitt = 19000,
                                      prior = pa_prior,
                                      pr = T) 
summary(sentinel_drone_model_sentinel_id)
plot(sentinel_drone_model_sentinel_id)
a <- 10000
pa_prior <- list(
  R=list(V=diag(1), nu=0.002),
  G=list(G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a),
         G1=list(V=diag(1), nu=1, alpha.mu=c(0), alpha.V=diag(1)*a)))


sentinel_drone_model_diff_n_sent_id <- sentinel_drone_model
sentinel_drone_model_diff_n_sent_id <- MCMCglmm(drone_ndvi~ sentinel_ndvi, 
                                 random = ~  diff + sentinel_id, 
                                 data = pixel_combos,
                                 nitt = 50000,
                                 prior = pa_prior,
                                 pr = T) 
summary(sentinel_drone_model_diff_n_sent_id)
plot(sentinel_drone_model_diff)
save(sentinel_drone_model_diff_n_sent_id, 
     file = paste0(data_out_path, "model_diff_n_sent_id.rda"))

#####

## Pretty plot for publication

# Extract veg model predictions and confidence intervals
load(paste0(data_out_path, "model_veg.rda"))

# PS4_HER on the 17 July 2017 is a big outlier let's remove it!
pixel_combos$site_veg_date <- paste0(pixel_combos$site_veg, 
                                     "_", pixel_combos$drone_date)
pixel_combos <- filter(pixel_combos, site_veg_date != "PS4HER_2017-07-17")

preds <- predict.MCMCglmm(sentinel_drone_model_veg, 
                          interval = "confidence", 
                          type = "response")
preds <- cbind(preds, 
               data.frame(
                 entinel_ndvi = sentinel_drone_model_veg$X[,2],
                 veg_type = as.character(sentinel_drone_model_veg$X[,3]),
                 interac = as.character(sentinel_drone_model_veg$X[,4]), 
                 stringsAsFactors = F))
preds[preds$veg_type == "0",]$veg_type <- "HER"
preds[preds$veg_type == "1",]$veg_type <- "KOM"

sentinel_drone_plot <- ggplot(data = pixel_combos, 
                              mapping = aes(x = sentinel_ndvi, 
                                            y = drone_ndvi, 
                                            colour = veg_type)) +
  geom_abline(intercept = 0, slope = 1, 
              color="black", linetype="dashed", 
              size= 1, alpha = 1 ) +
  geom_point(alpha = 0.5, size = 2) +
  geom_line(data = preds,
            mapping = aes(x = sentinel_ndvi,
                          y = fit,
                          group = veg_type,
                          linetype = veg_type
                          #colour = veg_type
            ),
            colour = "black",
            inherit.aes = FALSE,
            size = 2) +


  scale_x_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_y_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_colour_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_fill_manual(values = c("#1E9148FF", "#1E5C91FF")) +
  scale_linetype_manual(values = c("solid", "twodash")) +
  xlab('\nSentinel NDVI')+
  ylab("Drone NDVI\n") +
  annotate("text", x = 0.79, y = 0.38, 
           label = 'Herschel', colour = "#1E9148FF", 
           size = 12, hjust = 0) +
  annotate("text", x = 0.79, y = 0.33, 
           label = 'Komakuk', colour = "#1E5C91FF", 
           size = 12, hjust = 0) +
  annotate("line", x = c(0.73, 0.78), y = c(0.38, 0.38),
           size = 2, linetype = "solid") +
  annotate("line", x =  c(0.73, 0.78), y = c(0.33, 0.33),
           size = 2, linetype = "twodash") +
  
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 28, face = "bold"),
        axis.text.x = element_text(hjust = 0.5, size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        legend.position = "none") 
sentinel_drone_plot
save(sentinel_drone_plot, 
     file = paste0(data_out_path,
                   "sentinel_drone_plot_veg.Rda"))
ggsave(paste0(figure_out_path, "sentinel_drone_plot_veg.png"), 
       plot = sentinel_drone_plot, width = 12, height = 8)
       
## Plot with dark colours for models for Isla:

load(paste0(data_out_path, "pixel_combos.Rda"))
load(paste0(data_out_path, "model_veg.rda"))

# PS4_HER on the 17 July 2017 is a big outlier let's remove it!
pixel_combos$site_veg_date <- paste0(pixel_combos$site_veg, 
                                     "_", 
                                     pixel_combos$drone_date)
pixel_combos <- filter(pixel_combos, site_veg_date != "PS4HER_2017-07-17")


preds$veg_type2 <- paste0(as.character(preds$veg_type), "_2")
preds$veg_type2 <- factor(
  preds$veg_type2, 
  levels = c(unique(as.character(pixel_combos$veg_type)), 
             unique(as.character(preds$veg_type2))))
pixel_combos$veg_type2 <- factor(
  as.character(pixel_combos$veg_type), 
  levels = c(unique(as.character(pixel_combos$veg_type)), 
             unique(as.character(preds$veg_type2))))
levels(preds$veg_type2)
levels(pixel_combos$veg_type2)

plot_colours <- c("#1E9148FF", # Herschel Pixel Combos
                  "#1E5C91FF", # Komakuk Pixel Combos
                  "#004421FF",  # Herschel Predictions
                  "#002a6dFF") # Komakuk Predictions
names(plot_colours) <- levels(pixel_combos$veg_type2)
plot_scale <- scale_colour_manual(name = "veg_type2",values = plot_colours)

sentinel_drone_plot <- ggplot(
  data = pixel_combos, 
  mapping = aes(x = sentinel_ndvi, y = drone_ndvi, colour = veg_type2)) +
  geom_abline(intercept = 0, slope = 1, 
              color="black", linetype="dashed", 
              size= 1.5, alpha = 1 ) +
  geom_point(alpha = 0.5, size = 2) +
  geom_ribbon(data = preds,
              mapping = aes(x = sentinel_ndvi,
                            ymin = lwr,
                            ymax = upr, 
                            fill = veg_type2),
              inherit.aes = FALSE) +
  geom_line(data = preds,
            mapping = aes(x = sentinel_ndvi,
                          y = fit,
                          group = veg_type2, 
                          colour = veg_type2
            ),
            inherit.aes = FALSE,
            size = 0.5) +
  
  
  scale_x_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_y_continuous(limits = c(0.28,0.9), breaks = seq(0.3,0.9,0.1)) +
  scale_colour_manual(values = plot_colours, aesthetics = c("colour","fill")) +
  xlab('\nSentinel NDVI')+
  ylab("Drone NDVI\n") +
  annotate("text", x = 0.75, y = 0.38, 
           label = 'Herschel', colour = plot_colours[1], size = 12, hjust = 0) +
  annotate("text", x = 0.75, y = 0.33, 
           label = 'Komakuk', colour = plot_colours[2], size = 12, hjust = 0) +
  
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 28, face = "bold"),
        axis.text.x = element_text(hjust = 0.5, size = 24, colour = "black"),
        axis.text.y = element_text(size = 24, colour = "black"),
        legend.position = "none") 
sentinel_drone_plot
ggsave(paste0(figure_out_path, "sentinel_drone_plot_veg_isla.png"), 
       plot = sentinel_drone_plot, width = 12, height = 8)

### EOF 