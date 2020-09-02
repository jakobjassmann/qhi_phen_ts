# Phenology Time-Series Script of Plant Composition
# Isla Myers-Smith 1st September 2020

# Dependencies
library(tidyverse)
library(ggpubr)

data <- read.csv("data/plant_cover/qhi_cover_ITEX_1999_2017.csv")

data2 <- data %>%
  filter(year > 2015, !is.na(gfnarrowarft)) %>%
  mutate(subsite = dplyr::recode(
    sub_name, 
    'QHI:HE' = "Tussock Sedge Tundra", 
    'QHI:KO' = "Dryas-Vetch Tundra")) %>%
  rename(species = name) %>% 
  mutate(functional_group = dplyr::recode(
    gfnarrowarft, 
    'OTHER' = "Other", 
    'ROCK' = "Other", 
    'WATER' = "Other", 
    'SOIL' = "Other", 
    'DUNG' = "Other", 
    'FUNGI' = "Other", 
    'FUNGUS' = "Other",
    'FORBSV' = "Forbs", 
    'ESHRUB' = "Evergreen shrub", 
    'DSHRUB' = "Deciduous shrub", 
    'GRAMINOID' = "Grasses and sedges", 
    'LICHEN' = 'Lichen', 
    'MOSS' = "Moss", 
    'SHRUBU' = "Deciduous shrub",
    'SHRUB' = "Other vegetation", 
    'LIVER' = "Other vegetation", 
    'LITTER' = "Litter"))

data3 <- data2 %>%
  filter(functional_group != "Other vegetation") %>%
  group_by(year, subsite, plot, functional_group) %>%
  summarise(plot_cover = sum(cover)) %>%
  ungroup() %>%
  mutate(plot_cover_max = max(plot_cover), 
         plot_cover_100 = (plot_cover/plot_cover_max)*100)

data4 <- data3 %>% 
  group_by(year, subsite, functional_group) %>% 
  summarise(mean_cover = mean(plot_cover_100), 
            sd_cover = sd(plot_cover_100)) %>% 
  group_by(subsite, functional_group) %>% 
  summarise(mean_cover2 = mean(mean_cover), 
            sd_cover2 = mean(sd_cover)) %>% 
  mutate(functional_group = fct_reorder(functional_group, -mean_cover2))

# figure

# Set y axis label and colour colour
y_label_colour <- "black"

# Adjust visibility of axis labels according to plot positions
x.axis.text.colour <- "black"
y.axis.text.colour <- "black"

# plot a

(composition_plot <- ggplot(data4) +
    geom_col(aes(x = functional_group, 
                 y = mean_cover2, 
                 colour = subsite,
                 fill = subsite,
                 group = subsite),
             width = 0.8,
             position = position_dodge(0.88)) +
    geom_errorbar(aes(x = functional_group,
                      ymin = mean_cover2 - sd_cover2,
                      ymax = mean_cover2 + sd_cover2,
                      colour = subsite,
                      group = subsite),
                  width = 0.4,
                  position = position_dodge(0.88),
                  colour = "black") +
    #annotate("text", x = 2.3, y = 100, label = "A. Percent cover of functional groups", size = 6) +
    scale_colour_manual(values = c("#1e9148FF", "#1e5c91FF")) +
    scale_fill_manual(values = c("#1e9148FF", "#1e5c91FF")) +
    coord_cartesian(ylim = c(0,100)) +
    labs(fill = "", colour = "") +
    ylab("Functional Group % Cover\n") +
    xlab("") +
    guides(fill = guide_legend(title = NULL),
           colour = "none") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 0, 1, "cm"),
      axis.line = element_line(colour = "black", size = 1.2),
      axis.title = element_text(size = 16, 
                                face = "bold", 
                                colour = y_label_colour),
      axis.text.x = element_text(size = 14, 
                                 colour = x.axis.text.colour,
                                 angle = 45,
                                 hjust= 0.95,
                                 vjust= 1),
      axis.text.y = element_text(size = 14, colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 16, colour = x.axis.text.colour),
      legend.text = element_text(size = 14, colour = x.axis.text.colour)))

# plot b

data5 <- data2 %>%
  filter(functional_group != "Other vegetation") %>%
  filter(species == "Arctagrostis latifolia" | 
           species == "Dryas integrifolia" |
           species == "Eriophorum vaginatum" | 
           species == "Salix arctica" | 
           species == "Salix pulchra") %>% 
  group_by(year, subsite, plot, species) %>%
  summarise(plot_cover = mean(cover)) %>%
  ungroup() %>%
  mutate(plot_cover_max = max(plot_cover), 
         plot_cover_100 = (plot_cover/plot_cover_max)*100)

data6 <- data5 %>%
  group_by(year, subsite, species) %>% 
  summarise(mean_cover = mean(plot_cover_100), 
            sd_cover = sd(plot_cover_100)) %>% 
  group_by(subsite, species) %>% 
  summarise(mean_cover2 = mean(mean_cover), 
            sd_cover2 = mean(sd_cover)) %>% 
  ungroup() %>% 
  add_row(subsite = "Tussock Sedge Tundra", 
          species = "Salix arctica",
          mean_cover2 = 0, 
          sd_cover2 = NA) %>% 
  add_row(subsite = "Dryas-Vetch Tundra", 
          species = "Eriophorum vaginatum", 
          mean_cover2 = 0, 
          sd_cover2 = NA) %>% 
  mutate(species = fct_reorder(species, -mean_cover2))

(species_plot <- ggplot(data6) +
    geom_col(aes(x = species, 
                 y = mean_cover2, 
                 colour = subsite,
                 fill = subsite,
                 group = subsite),
             width = 0.8,
             position = position_dodge(0.88)) +
    geom_errorbar(aes(x = species,
                      ymin = mean_cover2 - sd_cover2,
                      ymax = mean_cover2 + sd_cover2,
                      colour = subsite,
                      group = subsite),
                  width = 0.4,
                  position = position_dodge(0.88),
                  colour = "black") +
    #annotate("text", x = 2.3, y = 100, label = "B. Percent cover of focal species", size = 6) +
    scale_colour_manual(values = c("#1e9148FF", "#1e5c91FF")) +
    scale_fill_manual(values = c("#1e9148FF", "#1e5c91FF")) +
    coord_cartesian(ylim = c(0,100)) +
    ylab("Species % Cover\n") +
    xlab("") +
    guides(fill = guide_legend(title = NULL),
           colour = "none") +
    theme_bw() +
    theme(
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(1, 1, 0, 1, "cm"),
      axis.line = element_line(colour = "black", size = 1.2),
      axis.title = element_text(size = 16, 
                                face = "bold", 
                                colour = y_label_colour),
      axis.text.x = element_text(size = 14, 
                                 colour = x.axis.text.colour,
                                 angle = 45,
                                 hjust= 0.95,
                                 vjust= 1,
                                 face = "italic"),
      axis.text.y = element_text(size = 14, colour = y.axis.text.colour),
      axis.ticks = element_line(size = 1.2),
      axis.ticks.length = unit(0.4, "cm"),
      axis.ticks.x = element_blank(),
      legend.title = element_text(size = 16, colour = x.axis.text.colour),
      legend.text = element_text(size = 14, colour = x.axis.text.colour)))

figure <- ggarrange(
  composition_plot, 
  species_plot, 
  ncol = 2, 
  nrow = 1, 
  widths = c(0.58, 0.42), 
  labels=c("a  Percent cover of functional groups",
           "b  Percent cover of focal species"),
  font.label = list(size = 24),
  label.x = 0.01,
  label.y = 0.99,
  vjust = c(1, 1), 
  hjust = c(0, 0),
  common.legend = TRUE,
  legend = "bottom") 

ggsave("Figures/fig_x_plant_cover/fig_x_plant_cover.png", 
       figure, width = 14, height = 8)
