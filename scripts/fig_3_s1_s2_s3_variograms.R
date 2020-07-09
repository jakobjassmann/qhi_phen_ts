# Varigrams for the TS in 2017
# Jakob Assmann j.assmann@ed.ac.uk 8 October 2017 updated 20 April 2020

### Preparations ----
# depenencies
library(dplyr)
library(ggplot2)
library(raster)
library(rasterVis)
library(viridisLite)
library(usdm)
library(parallel)
library(gstat)
library(cowplot)

# Set global parameters / load site boundaries and meta data
figure_out_path <- "figures/"
log_path <- "log/"
data_out_path <- "data/fig_3_variograms/"
site_boundaries <- read.csv("data/site_boundaries/ps_sent_site_bounds.csv")
load("data/meta_data.Rda")

# Prepare extent objects
get_sent_extent <- function(site_veg_id, sent_boundaries) {
  extent_object <- extent(
    c(sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$xmax,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymin,
      sent_boundaries[which(sent_boundaries$site_veg == site_veg_id),]$ymax))
  return(extent_object)
}

PS1_HER_extent <- get_sent_extent("PS1_HER", site_boundaries)
PS1_KOM_extent <- get_sent_extent("PS1_KOM", site_boundaries)
PS2_HER_extent <- get_sent_extent("PS2_HER", site_boundaries)
PS2_KOM_extent <- get_sent_extent("PS2_KOM", site_boundaries)
PS3_HER_extent <- get_sent_extent("PS3_HER", site_boundaries)
PS3_KOM_extent <- get_sent_extent("PS3_KOM", site_boundaries)
PS4_HER_extent <- get_sent_extent("PS4_HER", site_boundaries)
PS4_KOM_extent <- get_sent_extent("PS4_KOM", site_boundaries)

# filter meta data
meta_sub <- meta_data %>% 
  filter(site_name == "PS1" | site_name == "PS2",
         as.character(date) == "2017-06-26" | 
           as.character(date) == "2017-07-26" |
           as.character(date) == "2017-08-09", 
         band == "NDVI")
meta_sub <- meta_data %>% 
  filter(site_name == "PS4",
         as.character(date) == "2017-07-26" |
           as.character(date) == "2017-07-28", 
         band == "NDVI") %>%
  bind_rows(meta_sub) 

meta_sub <- meta_data %>% 
  filter(site_name == "PS3",
         as.character(date) == "2017-07-18", 
         band == "NDVI") %>%
  bind_rows(meta_sub) 

# Load raster files
list2env(
  lapply(
    setNames(meta_sub$file_path, 
             make.names(meta_sub$object_name)),
    raster), 
  envir = .GlobalEnv)

# crop rasters
list2env(
  mapply(
    function(x,y){crop(get(x, envir = .GlobalEnv), get(y, envir = .GlobalEnv))},
    setNames(meta_sub$object_name, 
             make.names(paste0(meta_sub$object_name, "_cropped"))), 
    paste0(meta_sub$site_veg, '_extent')
  ),
  envir = .GlobalEnv)
meta_sub$object_cropped <- paste0(meta_sub$object_name, "_cropped")

# 1) Convert to spatial dfs ----

# Convert rasters into spatial pixels data frames
cl <- makeCluster(3)
clusterExport(cl=cl, varlist=meta_sub$object_cropped)
list2env(
  parLapply(cl, 
            setNames(meta_sub$object_cropped, 
                     make.names(paste0(meta_sub$object_cropped, "_spdf"))),
            function(x){ 
              library(gstat)
              library(raster)
              start.time <- Sys.time()
              cat(as.character(start.time),": Starting ", x, "...")
              spdf <- as(get(x), "SpatialPixelsDataFrame" ) 
              end.time <- Sys.time()
              time.elapsed <- end.time - start.time
              cat(" finished after: ", as.character(time.elapsed), " hh:mm:ss\n")
              return(spdf)}
  ),
  envir = .GlobalEnv)
stopCluster(cl)
meta_sub$object_spdf <- paste0(meta_sub$object_cropped, "_spdf")


#check wehether they are square 
sapply(sort(meta_sub$object_spdf), 
       function(x){ get(x)@grid@cellsize[1]== get(x)@grid@cellsize[2]})

# adjust round grids of non-square rasters to nearest 0.001 mm
# due to some weird reason with floating point precision and R this
# will have to be done by hand

PS1_HER_20170626_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_HER_20170626_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_HER_20170626_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04786
PS1_HER_20170626_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04786

PS1_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04947
PS1_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04947

PS1_HER_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_HER_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_HER_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04823
PS1_HER_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04823

PS1_KOM_20170626_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_KOM_20170626_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_KOM_20170626_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04754
PS1_KOM_20170626_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04754

PS1_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04806
PS1_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04806

PS1_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS1_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2]
PS1_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.0471
PS1_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.0471

PS2_HER_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[1] 
PS2_HER_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS2_HER_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.05597
PS2_HER_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.05597

PS2_KOM_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS2_KOM_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS2_KOM_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.05512
PS2_KOM_2017_06_26_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.05512

PS2_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS2_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS2_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.05369
PS2_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.05369

PS2_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS2_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS2_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.05296
PS2_KOM_20170809_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.05296

PS3_HER_20170718_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS3_HER_20170718_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS3_HER_20170718_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04574
PS3_HER_20170718_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04574

PS3_KOM_20170718_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS3_KOM_20170718_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS3_KOM_20170718_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04978
PS3_KOM_20170718_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04978

PS4_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS4_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS4_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04594
PS4_HER_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04594

PS4_KOM_20170728_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS4_KOM_20170728_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS4_KOM_20170728_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.05132
PS4_KOM_20170728_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.05132

# Export to gstat dataframes to scratch folder this will save memory later on.
sapply(meta_sub$object_spdf, function(object_spdf) {
  save(list = object_spdf, file = paste0(data_out_path, object_spdf, ".Rda"))
})

# clear memory
rm(list = ls()[grep("_spdf", ls())])
gc()

### 2) Sample variograms ----
# extract variance data
cl <- makeCluster(3)
clusterExport(cl=cl, "data_out_path")
clusterEvalQ(cl, library(gstat))

# Open connection to log file and export
parallelLog <- paste0(log_path,"parallel_log.txt")
clusterExport(cl=cl,"parallelLog")

# Define function to calculuate a variogram 
sample_variogram <- function(object_spdf, thin = 20 , cut_off, bin_width){ 
  # Time the operation
  start.time <- Sys.time()
  cat(as.character(start.time),": Starting ", object_spdf, "... \n", sep="", 
      file = parallelLog, append=TRUE)
  
  # load spdf object 
  load(paste0(data_out_path, object_spdf, ".Rda"))
  
  # Wrapper to access raster layer in variogram function
  raster_layer_name <- names(get(object_spdf))[1]
  
  # Sample the variogram (this can take ages)
  vario <- variogram(get(raster_layer_name) ~ 1, 
                     get(object_spdf)[sample(nrow(get(object_spdf)) / thin),],
                     cutoff = cut_off,
                     width = bin_width,
                     verbose = T) 
  # Stop timer
  end.time <- Sys.time()
  time.elapsed <- end.time - start.time
  cat(as.character(end.time), ": ... finished ", object_spdf, ".\n", sep="", 
      file = parallelLog, append=TRUE)
  
  # Save variogram to harddrive
  copy_of_vario <- vario
  assign(paste0(object_spdf,"_", cut_off, "m_vario"), copy_of_vario)
  local_env <- environment()
  save(list = paste0(object_spdf,"_", cut_off, "m_vario"), 
       file=paste0(data_out_path, object_spdf, "_", cut_off,"m_vario.Rda"), 
       envir = local_env)
  
  # clean memory
  gc()
  
  # Return variogram
  return(get(paste0(object_spdf,"_", cut_off, "m_vario")))
}

# Export to cluster
clusterExport(cl=cl, "sample_variogram")

# Calculate the variograms (this will likely take some time)
list2env(
  parLapply(cl, 
            setNames(meta_sub$object_spdf, 
                     make.names(paste0(meta_sub$object_spdf, "_15m_vario"))),
            function(x){
              sample_variogram(x, thin = 20, cut_off = 15, bin_width = 0.15)
            }),
  envir = .GlobalEnv)
meta_sub$object_15m_vario <- paste0(meta_sub$object_spdf, "_15m_vario")

list2env(
  parLapply(cl, setNames(meta_sub$object_spdf, 
                         make.names(paste0(meta_sub$object_spdf, 
                                           "_45m_vario"))),
            function(x){
              sample_variogram(x, thin = 20, cut_off = 45, bin_width = 3)
            }),
  envir = .GlobalEnv)
meta_sub$object_45m_vario <- paste0(meta_sub$object_spdf, "_45m_vario")

stopCluster(cl)

### 3) Fit variogram models ----

# Define helper function
fit_vario <- function(object_vario){
  
  # Time the operation
  start.time <- Sys.time()
  cat(as.character(start.time),": Starting ", object_vario, "... \n", sep="", 
      file=parallelLog, append=TRUE)
  
  # Fit variogram model
  vario_fit <- fit.variogram(get(object_vario), vgm(c("Exp", "Mat", "Sph")))
  
  # Stop timer
  end.time <- Sys.time()
  cat(as.character(end.time), ": ... finished ", object_vario, ".\n", sep="", 
      file=parallelLog, append=TRUE)
  
  # Save variogram to harddrive
  copy_of_vario_fit <- vario_fit
  assign(paste0(object_vario,'_fit'), copy_of_vario_fit)
  save(list = paste0(object_vario, "_fit"), 
       file=paste0(data_out_path, object_vario, "_fit.Rda"))
  
  # clean memory
  gc()
  
  # Restore output to console
  
  # Return variogram fit
  return(vario_fit)
}


# Execute extraction in parallel
cl <- makeCluster(3)
clusterExport(cl=cl, "parallelLog")
clusterExport(cl=cl, "data_out_path")
clusterEvalQ(cl, library(gstat))

clusterExport(cl=cl, varlist=meta_sub$object_15m_vario)
list2env(
  parLapply(
    cl, 
    setNames(meta_sub$object_15m_vario, 
             make.names(paste0(meta_sub$object_15m_vario , "_fit"))),
    fit_vario),
  envir = .GlobalEnv)

clusterExport(cl=cl, varlist=meta_sub$object_45m_vario)
list2env(
  parLapply(
    cl, 
    setNames(meta_sub$object_45m_vario, 
             make.names(paste0(meta_sub$object_45m_vario, "_fit"))),
    fit_vario),
  envir = .GlobalEnv)

### 4) Export variance data ----

# Extract data from variograms
varios_15m_df <- bind_rows(
  mapply(
    function(x, q, y, z) {
      var_dist <- get(x)$dist
      var_gamma <- get(x)$gamma
      var_np <- get(x)$np
      return_df <- data.frame(flight_id = q,
                              site_veg = y,
                              date = z,
                              dist = var_dist,
                              gamma = var_gamma,
                              np = var_np,
                              max_dist = "15m",
                              stringsAsFactors = F)
      return(return_df)
    }, 
    paste0(meta_sub$object_15m_vario),
    meta_sub$flight_id,
    meta_sub$site_veg,
    meta_sub$date,
    SIMPLIFY = F
  ))
varios_45m_df <- bind_rows(
  mapply(
    function(x, q, y, z) {
      var_dist <- get(x)$dist
      var_gamma <- get(x)$gamma
      var_np <- get(x)$np
      return_df <- data.frame(flight_id = q,
                              site_veg = y,
                              date = z,
                              dist = var_dist,
                              gamma = var_gamma,
                              np = var_np,
                              max_dist = "45m",
                              stringsAsFactors = F)
      return(return_df)
    }, 
    paste0(meta_sub$object_45m_vario),
    meta_sub$flight_id,
    meta_sub$site_veg,
    meta_sub$date,
    SIMPLIFY = F
  ))

varios_all_df <- rbind(varios_15m_df, varios_45m_df)
varios_all_df$site_date <- 
  paste0(varios_all_df$site_veg, "_" ,varios_all_df$date)
varios_all_df$flight_mdist <- 
  paste0(varios_all_df$flight_id, "_" ,varios_all_df$max_dist)

# Extract perpare variogram fits
vario_fits_15m_df <- bind_rows(mapply(function(x,q,y,z){
  vario_line <- variogramLine(get(paste0(x, "_fit")), 15)
  fit_model <- get(paste0(x, "_fit"))$model[2]
  fit_psill <- get(paste0(x, "_fit"))$psill[2]
  fit_range <- get(paste0(x, "_fit"))$range[2]
  return_df <- data.frame(flight_id = q,
                          site_veg = y,
                          date = z,
                          dist = vario_line$dist,
                          gamma = vario_line$gamma,
                          model = fit_model,
                          sill = fit_psill,
                          range = fit_range,
                          max_dist = "15m",
                          stringsAsFactors = F)
  return(return_df)
}, 
paste0(meta_sub$object_15m_vario),
meta_sub$flight_id,
meta_sub$site_veg,
meta_sub$date,
SIMPLIFY = F
))

vario_fits_45m_df <- bind_rows(mapply(function(x,q,y,z){
  vario_line <- variogramLine(get(paste0(x, "_fit")), 45)
  fit_model <- get(paste0(x, "_fit"))$model[2]
  fit_psill <- get(paste0(x, "_fit"))$psill[2]
  fit_range <- get(paste0(x, "_fit"))$range[2]
  return_df <- data.frame(flight_id = q,
                          site_veg = y,
                          date = z,
                          dist = vario_line$dist,
                          gamma = vario_line$gamma,
                          model = fit_model,
                          sill = fit_psill,
                          range = fit_range,
                          max_dist = "45m",
                          stringsAsFactors = F)
  return(return_df)
}, 
paste0(meta_sub$object_15m_vario),
meta_sub$flight_id,
meta_sub$site_veg,
meta_sub$date,
SIMPLIFY = F
))

varios_fits_all_df <- rbind(vario_fits_15m_df, vario_fits_45m_df)
varios_fits_all_df$site_date <- paste0(varios_fits_all_df$site_veg, "_" ,varios_fits_all_df$date)
varios_fits_all_df$flight_mdist <- paste0(varios_fits_all_df$flight_id, "_" ,varios_fits_all_df$max_dist)

varios_all_df$flight_id <- factor(varios_all_df$flight_id)
varios_fits_all_df$flight_id <- factor(varios_fits_all_df$flight_id, levels = sort(unique(varios_fits_all_df$flight_id)))

levels(varios_all_df$flight_id)
levels(varios_fits_all_df$flight_id)

### 5) Plot data ----

# Prepare plots
# Set colours
her_col <- "#1e9148FF"
kom_col <- "#1e5c91FF"
# Set site names
site_names <- data.frame(
  site_name = c("PS1", "PS2", "PS3", "PS4"),
  site_name_full = c("Collinson Head",
                     "Bowhead Ridge",
                     "Hawk Valley",
                     "Hawk Ridge")
  )

# Set overall peak season range mean
range_mean <- varios_fits_all_df %>% 
  filter(max_dist == "15m") %>%
  filter(date == as.Date("2017-07-26") |
           date == as.Date("2017-07-28")) %>%
  filter(site_veg != "PS3_KOM") %>%
  pull(range) %>%
  unique() %>%
  mean(na.rm = T)

##  _5a) Figure 3a ---- 
## Paired Komakuk and Herschel variograms fo PS2 (Figure 3 panel a)
site_name <- "PS2"
varios <- varios_all_df %>% 
  filter(site_veg == paste0(site_name, "_HER") 
         | site_veg == paste0(site_name, "_KOM")) %>% 
  filter(date == as.Date("2017-07-26") |
           date == as.Date("2017-07-28")) %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels = c("HER", "KOM")))
vario_fits <- varios_fits_all_df %>% 
  filter(site_veg == paste0(site_name, "_HER") 
         | site_veg == paste0(site_name, "_KOM")) %>% 
  filter(date == as.Date("2017-07-26") |
           date == as.Date("2017-07-28")) %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels =c("HER", "KOM")))

her_sill <- unique(round(vario_fits$sill[vario_fits$veg == "HER"], 3))
kom_sill <- unique(round(vario_fits$sill[vario_fits$veg == "KOM"], 3))
if(her_sill > kom_sill){
  her_sill <- her_sill + 0.001
  kom_sill <- kom_sill - 0.001
}
if(her_sill < kom_sill){
  her_sill <- her_sill - 0.001
  kom_sill <- kom_sill + 0.001
}

site_name_full <- site_names$site_name_full[site_names$site_name == site_name] 

max_gamma <- round(max(varios$gamma),3) + 0.001
# Plot 4m plot
vario_plot_4m <- ggplot(filter(varios, max_dist == "15m" & dist < 4),
                        aes(x = dist, y = gamma, 
                            group = veg,
                            colour = veg)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, max_dist == "15m" & dist < 4), 
            mapping = aes(x=dist, y = gamma, 
                          group = veg, 
                          colour = veg,
                          alpha = 0.6),
            size = 1.5,
            inherit.aes = F) +
  geom_vline(xintercept = range_mean, colour = "black", size = 1, alpha = 0.6) +
  scale_colour_manual(values = c(her_col, kom_col)) +
  scale_x_continuous(breaks = seq(0,4,0.5)) +
  scale_y_continuous(breaks = seq(0.000, max_gamma,0.001)) +
  #ggtitle(paste0("Site ", substr(site_name, 3, 3), " - ", site_name_full)) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  annotate("text", x = 4, y = her_sill, label = "Tussock Sedge Tundra", 
           hjust = 1, colour = her_col, size = 6) +    
  annotate("text", x = 4, y = kom_sill, label = "Dryas-Vetch Tundra", 
           hjust = 1, colour = kom_col, size = 6) +
  annotate("text", x = range_mean + 0.1, y = 0, label = paste0(round(range_mean, 2), " m"), 
           hjust = 0, vjust = 0, colour = "black", size = 5) +
  annotate("segment", x = rep(-0.2, 13), xend = rep(-0.23, 13), yend = seq(0, max_gamma, 0.0005), y = seq(0, max_gamma, 0.0005)) +
  annotate("segment", x = seq(0, 4, 0.1), xend = seq(0, 4, 0.1), yend = rep(-0.0003, 41), y = rep(-0.000375, 41)) +
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) +
  coord_cartesian(xlim = c(0, 4),
                  ylim = c(0, max_gamma),
                  clip = 'off') 

save_plot(vario_plot_4m, 
          filename = paste0(figure_out_path, 
                            "/fig_3_variograms/",
                            site_name, "_vario_4m.png"),
          base_aspect_ratio = 1.6)
# 45 m plot 
vario_plot_45m <- ggplot(varios,
                        aes(x = dist, y = gamma, 
                            group = veg,
                            colour = veg)) +
  geom_point(size = 2) + 
  geom_line(data = vario_fits, 
            mapping = aes(x=dist, y = gamma, 
                          group = veg, 
                          colour = veg,
                          alpha = 0.6),
            size = 1.5,
            inherit.aes = F) +
  scale_colour_manual(values = c(her_col, kom_col)) +
  scale_x_continuous(breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,max_gamma), 
                     breaks = seq(0.000, max_gamma,0.001)) +
  #ggtitle(paste0("Site ", substr(site_name, 3, 3), " - ", site_name_full)) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  annotate("text", x = 45, y = her_sill, label = "Tussock Sedge Tundra", 
           hjust = 1, colour = her_col, size = 6) +    
  annotate("text", x = 45, y = kom_sill, label = "Dryas-Vetch Tundra", 
           hjust = 1, colour = kom_col, size = 6) +
  annotate("segment", x = rep(-2.3, 13), xend = rep(-2.6, 13), yend = seq(0, max_gamma, 0.0005), y = seq(0, max_gamma, 0.0005)) +
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) +
  coord_cartesian(xlim = c(0, 45),
                  clip = 'off') 
save_plot(vario_plot_45m, 
          filename = paste0(figure_out_path, 
                            "/fig_3_variograms/",
                            site_name, "_vario_45m.png"),
          base_aspect_ratio = 1.6)

## _5b) Figure S1 ----
## Peak season variograms for all sites (figure_s1 panel a her and b kom)
varios <- varios_all_df %>% 
  filter(date == as.Date("2017-07-26") |
           date == as.Date("2017-07-28")) %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels = c("HER", "KOM")))
vario_fits <- varios_fits_all_df %>% 
  filter(date == as.Date("2017-07-26") |
           date == as.Date("2017-07-28")) %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels =c("HER", "KOM")))

palette <- sequential_hcl(5, palette = "Greens 3")[c(-4,-5)]

vario_plot_her <- ggplot(filter(varios, veg == "HER"),
                         aes(x = dist, y = gamma, 
                             group = site_date,
                             colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "HER"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.006), 
                     breaks = seq(0.000, 0.006,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Peak-Season Variograms", subtitle ="Tussock Sedge Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "Site 1", 
           hjust = 1, colour = palette[1]	, size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "Site 2", 
           hjust = 1, colour = palette[2]	, size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "Site 4", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 

palette <- sequential_hcl(5, palette = "Blues 3")[c(-4,-5)]

vario_plot_kom <- ggplot(filter(varios, veg == "KOM" & site_veg != "PS3_KOM"),
                         aes(x = dist, y = gamma, 
                             group = site_date,
                             colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "KOM" & site_veg != "PS3_KOM"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.006), 
                     breaks = seq(0.000, 0.006,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Peak-Season Variograms", subtitle ="Dryas-Vetch Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "Site 1", 
           hjust = 1, colour = palette[1], size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "Site 2", 
           hjust = 1, colour = palette[2]	, size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "Site 4", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 
save_plot(paste0(figure_out_path, "fig_s1_peak_season_varios.png"),
          plot_grid(vario_plot_her, vario_plot_kom, ncol = 2, labels = "auto",
                    label_size = 26),
          base_aspect_ratio = 2.6)

## _5c) Figure S2 ----
## Cross season progression of variograms for ps 1 and ps2
varios <- varios_all_df %>% 
  mutate(site_name = substr(site_veg, 1,3)) %>%
  filter(site_name == "PS1" | site_name == "PS2") %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels = c("HER", "KOM")))
vario_fits <- varios_fits_all_df %>% 
  mutate(site_name = substr(site_veg, 1,3)) %>%
  filter(site_name == "PS1" | site_name == "PS2") %>%
  mutate(veg = ordered(substr(site_veg, 5, 8), levels =c("HER", "KOM")))

palette <- sequential_hcl(4, palette = "Plasma")[-4]

vario_plot_ps1_her <- ggplot(filter(varios, veg == "HER" & site_name == "PS1"),
                         aes(x = dist, y = gamma, 
                             group = site_date,
                             colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "HER" & site_name == "PS1"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.008), 
                     breaks = seq(0.000, 0.008,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Site 1 - Collinson Head", subtitle ="Tussock Sedge Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "2017-06-26", 
           hjust = 1, colour = palette[1]	, size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "2017-07-26", 
           hjust = 1, colour = palette[2], size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "2017-08-09", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 

vario_plot_ps1_kom <- ggplot(filter(varios, veg == "KOM" & site_name == "PS1"),
                             aes(x = dist, y = gamma, 
                                 group = site_date,
                                 colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "KOM" & site_name == "PS1"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.008), 
                     breaks = seq(0.000, 0.008,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Site 1 - Collinson Head", subtitle ="Dryas-Vetch Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "2017-06-26", 
           hjust = 1, colour = palette[1]	, size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "2017-07-26", 
           hjust = 1, colour = palette[2]	, size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "2017-08-09", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 

vario_plot_ps2_her <- ggplot(filter(varios, veg == "HER" & site_name == "PS1"),
                             aes(x = dist, y = gamma, 
                                 group = site_date,
                                 colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "HER" & site_name == "PS1"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.008), 
                     breaks = seq(0.000, 0.008,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Site 2 - Bowhead Ridge", subtitle ="Tussock Sedge Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "2017-06-26", 
           hjust = 1, colour = palette[1], size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "2017-07-26", 
           hjust = 1, colour = palette[2], size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "2017-08-09", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 

vario_plot_ps2_kom <- ggplot(filter(varios, veg == "KOM" & site_name == "PS1"),
                             aes(x = dist, y = gamma, 
                                 group = site_date,
                                 colour = site_date)) +
  geom_point(size = 2) + 
  geom_line(data = filter(vario_fits, veg == "KOM" & site_name == "PS1"), 
            mapping = aes(x=dist, y = gamma, 
                          group = site_date, 
                          colour = site_date),
            size = 1.5,
            alpha = 0.6,
            inherit.aes = F) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.008), 
                     breaks = seq(0.000, 0.008,0.002)) +
  scale_colour_manual(values = palette) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  ggtitle("Site 2 - Bowhead Ridge", subtitle ="Dryas-Vetch Tundra") +
  annotate("text", x = 45, y = 0.0018, label = "2017-06-26", 
           hjust = 1, colour = palette[1]	, size = 6) +    
  annotate("text", x = 45, y = 0.001, label = "2017-07-26", 
           hjust = 1, colour = palette[2]	, size = 6) +    
  annotate("text", x = 45, y = 0.0002, label = "2017-08-09", 
           hjust = 1, colour = palette[3]	, size = 6) +    
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 

cross_season_plots <- plot_grid(vario_plot_ps1_her,
                                vario_plot_ps1_kom,
                                vario_plot_ps2_her,
                                vario_plot_ps2_kom,
                                nrow = 2, ncol = 2,
                                labels = "auto",
                                label_size = 26)
save_plot(paste0(figure_out_path, "fig_s2_cross_season_varios.png"),
          cross_season_plots,
          base_height = 8,
          base_aspect_ratio = 1.3)

## _5d) Figure S3 ----
### Plot Site 3 variograms 
site_name <- "PS3"
site_name_full <- site_names$site_name_full[site_names$site_name == site_name] 

her_sill <- 0.005
kom_sill <- 0.0025

varios <- varios_all_df %>% 
  filter(site_veg == paste0(site_name, "_HER") 
         | site_veg == paste0(site_name, "_KOM")) %>% 
  mutate(veg = ordered(substr(site_veg, 5, 8), levels = c("HER", "KOM")))
vario_fits <- varios_fits_all_df %>% 
  filter(site_veg == paste0(site_name, "_HER") 
         | site_veg == paste0(site_name, "_KOM")) %>% 
  mutate(veg = ordered(substr(site_veg, 5, 8), levels =c("HER", "KOM")))

vario_plot_45m_PS3 <- ggplot(varios,
                         aes(x = dist, y = gamma, 
                             group = veg,
                             colour = veg)) +
  geom_point(size = 2) + 
  geom_line(data = vario_fits, 
            mapping = aes(x=dist, y = gamma, 
                          group = veg, 
                          colour = veg,
                          alpha = 0.6),
            size = 1.5,
            inherit.aes = F) +
  scale_colour_manual(values = c(her_col, kom_col)) +
  scale_x_continuous(limits = c(0,45), 
                     breaks = seq(0,45,5)) +
  scale_y_continuous(limits = c(0,0.005), 
                     breaks = seq(0.000, 0.005,0.001)) +
  ggtitle(paste0("Site ", substr(site_name, 3, 3), " - ", site_name_full, " 18 July 2020")) +
  ylab("NDVI semivariance") +
  xlab("Distance (m)") +
  annotate("text", x = 45, y = her_sill, label = "Tussock Sedge Tundra", 
           hjust = 1, colour = her_col, size = 6) +    
  annotate("text", x = 45, y = kom_sill, label = "Dryas-Vetch Tundra", 
           hjust = 1, colour = kom_col, size = 6) +
  theme_cowplot(18) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0)) 
save_plot(vario_plot_45m_PS3, 
          filename = paste0(figure_out_path, 
                            "/fig_s3_site3_variograms.png"),
          base_aspect_ratio = 1.6)

# Back up commands to allow loading of already calculated variograms and fits 
# cat(unlist(lapply(meta_sub$object_15m_vario, function(x) paste0("load(\"data/fig_3_variograms/", x, "_fit.Rda\")"))), sep = "\n")
# cat(unlist(lapply(meta_sub$object_45m_vario, function(x) paste0("load(\"data/fig_3_variograms/", x, "_fit.Rda\")"))), sep = "\n")
# cat(unlist(lapply(meta_sub$object_15m_vario, function(x) paste0("load(\"data/fig_3_variograms/", x, ".Rda\")"))), sep = "\n")
# cat(unlist(lapply(meta_sub$object_45m_vario, function(x) paste0("load(\"data/fig_3_variograms/", x, ".Rda\")"))), sep = "\n")
