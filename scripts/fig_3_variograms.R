# Varigrams for the TS in 2017
# Jakob Assmann j.assmann@ed.ac.uk 8 October 2017

# depenencies
library(dplyr)
library(ggplot2)
library(raster)
library(rasterVis)
library(viridisLite)
library(usdm)
library(parallel)
library(gstat)

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
  filter(site_name == "PS3" | site_name == "PS4",
         as.character(date) == "2017-07-26" |
           as.character(date) == "2017-07-28", 
         band == "NDVI") %>%
  bind_rows(meta_sub) 

# Load files
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

# Create variograms
# Calculate variograms
# Convert rasters into spatial pixels data frames
# parallel for this one
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

PS3_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1]
PS3_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] 
PS3_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[1] <- 0.04794
PS3_KOM_20170726_50m_ndvi_cropped_spdf@grid@cellsize[2] <- 0.04794

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

# extract variance data
cl <- makeCluster(3)
clusterExport(cl=cl, "data_out_path")
clusterEvalQ(cl, library(gstat))

# Open connection to log file and export
parallelLog <- paste0(log_path,"parallel_log.txt")
clusterExport(cl=cl,"parallelLog")

# Definte function to calculuate the variogram (ready for parallel processing)
sample_variogram <- function(x, thin = 20 , cut_off, bin_width){ 
  # Time the operation
  start.time <- Sys.time()
  cat(as.character(start.time),": Starting ", x, "... \n", sep="", 
      file = parallelLog, append=TRUE)
  
  # load spdf object 
  load(paste0(scratch_folder, object_spdf, ".Rda"))
  
  # Wrapper to access raster layer in variogram function
  raster_layer_name <- names(get(object_spdf))[1]
  
  # Sample the variogram (this can take ages)
  vario <- variogram(get(raster_layer_name) ~ 1, 
                     get(x)[sample(nrow(get(x)) / thin),],
                     cutoff = cut_off,
                     width = bin_width,
                     verbose = T) 
  # Stop timer
  end.time <- Sys.time()
  time.elapsed <- end.time - start.time
  cat(as.character(end.time), ": ... finished ", x, ".\n", sep="", 
      file = parallelLog, append=TRUE)
  
  # Save variogram to harddrive
  copy_of_vario <- vario
  assign(paste0(x,"_", cut_off, "m_vario"), copy_of_vario)
  local_env <- environment()
  save(list = paste0(x,"_", cut_off, "m_vario"), 
       file=paste0(data_out_path, x, "_", cut_off,"m_vario.Rda"), 
       envir = local_env)
  
  # clean memory
  gc()
  
  # Return variogram
  return(get(paste0(x,"_", cut_off, "m_vario")))
}

# Export to cluster
clusterExport(cl=cl, "sample_variogram")

# go fo the variograms (this will likely take ages)
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


# fit variogram models

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
       file=paste0(scratch_folder, object_vario, "_fit.Rda"))
  
  # clean memory
  gc()
  
  # Restore output to console
  
  # Return variogram fit
  return(vario_fit)
}

# Execute extraction in parallel
clusterExport(cl=cl, varlist=meta_sub$object_15m_vario)
list2env(
  parLapply(cl, 
            setNames(meta_sub$object_15m_vario, 
                     make.names(paste0(meta_sub$object_15m_vario , "_fit"))),
            fit_vario,
            envir = .GlobalEnv)
  clusterExport(cl=cl, varlist=meta_sub$object_45m_vario)
  list2env(
    parLapply(cl, 
              setNames(meta_sub$object_45m_vario, 
                       make.names(paste0(meta_sub$object_45m_vario , "_fit"))),
              fit_vario,
  envir = .GlobalEnv)


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
# Plot all data
plot_colours <- c("#1E914870","#1E9148AA", "#1E9148FF", "#1E5C9170", "#1E5C91AA", "#1E5C91FF")
varios_all_plot <- ggplot(varios_all_df, aes(x = dist, y = gamma, group = flight_mdist , colour = flight_id)) +
  geom_point(size = 2) + 
  geom_line(data = filter(varios_fits_all_df, max_dist == "45m"), 
            mapping = aes(x=dist, y = gamma, group = flight_mdist , colour = flight_id),
            size = 1.5,
            inherit.aes = F) +
  scale_colour_manual(values = plot_colours) +
  scale_x_continuous(limits = c(0,45), breaks =  seq(0,45, 5)) +
  scale_y_continuous(limits = c(0,0.011), breaks = seq(0,0.011, 0.002)) +
  # annotate("text", x = 35,  y = 0.0095, label = paste0("Range: ",round(range,2), " m"), colour = "black", hjust = 0, size = 5) +
  # annotate("text", x = 35,  y = 0.009, label = paste0("Sill: ", round(sill, 4)), colour = "black", hjust = 0, size = 5) +
  annotate("text", x= 17.5, y = 0.0009, label = "Komakuk 26 June 2017", colour = plot_colours[4], hjust = 0, size = 7) +
  annotate("text", x= 17.5, y = 0.0005, label = "Komakuk 27 July 2017", colour = plot_colours[5], hjust = 0, size = 7) +
  annotate("text", x= 17.5, y = 0.0001, label = "Komakuk 09 Aug. 2017", colour = plot_colours[6], hjust = 0, size = 7) +
  annotate("text", x= 32.5, y = 0.0009, label = "Herschel  26 June 2017", colour = plot_colours[1], hjust = 0, size = 7) +
  annotate("text", x= 32.5, y = 0.0005, label = "Herschel  27 July 2017", colour = plot_colours[2], hjust = 0, size = 7) +
  annotate("text", x= 32.5, y = 0.0001, label = "Herschel  09 Aug. 2017", colour = plot_colours[3], hjust = 0, size = 7) +
  
  ylab("NDVI semivariance\n") +
  xlab("\n Distance (m)") +
  theme_bw()  +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(hjust = 0.5, size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        legend.position = "none") 

varios_all_plot


# Close up 
# Calculate range mean
range_mean <- mean(unique(varios_fits_all_df$range))

varios_1m_plot <- ggplot(filter(varios_all_df, max_dist == "15m"), aes(x = dist, y = gamma, group = flight_mdist , colour = flight_id)) +
  geom_vline(xintercept = range_mean, colour = "black", size = 0.5, alpha = 0.6) +
  geom_point(size = 4) + 
  geom_line(data = filter(varios_fits_all_df,max_dist == "15m"), 
            mapping = aes(x=dist, y = gamma, group = flight_mdist , 
                          colour = flight_id),
            size = 1.5,
            inherit.aes = F) +
  scale_colour_manual(values = plot_colours) +
  scale_x_continuous(limits = c(0,1), breaks =  seq(0,1,  0.2)) + # originally breaks at 0.1
  scale_y_continuous(limits = c(0,0.008), breaks = seq(0,0.01, 0.002)) +
  annotate("text", x = 0.55, y = 0.0004, label = paste0("Mean range: ", round(range_mean, 2)),
           hjust = 0, size = 5, colour = "black", alpha = 0.6) + # font size originally 4
  #annotate("text", x = 0.55, y = 0.0003, label = paste0(round(range_mean, 2)),
  #         hjust = 0, size = 5, colour = "black", alpha = 0.6) + # font size originally 4
  ylab("Semivariance") +
  xlab("Distance (m)") +
  theme_bw()  +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 16, face = "bold"), # font size originally 12
        axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.text.x = element_text(hjust = 0.5, size = 16, colour = "black"),# font size originally 10
        axis.text.y = element_text(size = 16, colour = "black"), # font size originally 10
        legend.position = "none",
        plot.margin = margin(0,0.5,0,0, "cm")) 

varios_1m_plot
library(cowplot)
ggdraw() +
  draw_plot(varios_all_plot, 0, 0, 1, 1) +
  # geom_rect(xmin = 0.73, xmax = 1, ymin = 0.73, ymax = 1, colour = "black", fill = "#FFFFFF00") +
  draw_plot(varios_1m_plot, 0.61, 0.64, 0.375, 0.35) +
  draw_plot_label(c("A", "B"), c(0, 0.57), c(1, 1), size = 28)
ggsave(paste0(figure_out_path, "vario_all.png"), width = 12, height = 8)

# finaly calculate min and mean sample size for each distance bin in the variograms
varios_all_df %>% group_by(flight_id, max_dist) %>% summarise(min_np = min(np), mean_np = mean(np))
# 1 PS2_HER_2017_06_26      15m   1575137  34997691
# 2 PS2_HER_2017_06_26      45m 430577883 590919868
# 3   PS2_HER_20170726      15m   1620711  37070002
# 4   PS2_HER_20170726      45m 459589096 626514466
# 5   PS2_HER_20170809      15m   1728699  42141992
# 6   PS2_HER_20170809      45m 520144345 711877532
# 7 PS2_KOM_2017_06_26      15m   1624349  37257191
# 8 PS2_KOM_2017_06_26      45m 461933804 629582394
# 9   PS2_KOM_20170726      15m   1712948  41371323
# 10   PS2_KOM_20170726      45m 510486983 698535398
# 11   PS2_KOM_20170809      15m   2110116  43722034
# 12   PS2_KOM_20170809      45m 539193736 738115478
# Looks fairly similar for all variograms. Let's report overall stats
varios_all_df %>% group_by(max_dist) %>% summarise(min_np = min(np), mean_np = mean(np))
# 1      15m   1575137  39426706
# 2      45m 430577883 665924189
# That seems sensible, a minimum smaple size of 1.6m point pairs and mean of 39 m for the 14 m ones
# and min sample of 430m and mean of 666m for the big ones.
