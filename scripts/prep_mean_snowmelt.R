# Quick script to derive mean snowmelt for 2016 and 2017 for the
# long-term phenology plots on Qikiqtaruk 

# Jakob Assmann 3/9/2020

library(dplyr)

# Load phenology data
qhi_phen <- read.csv("data/auxillary/qiki_phen_2019.csv")

# Extract mean snowmelt
snowmelt <- qhi_phen %>% 
  group_by(Year) %>%
  filter(Year %in% c(2016,2017)) %>%
  select(Year, P1) %>%
  summarise(mean_snowmelt = mean(P1)) %>%
  mutate(date = as.Date(paste0(Year, "-", mean_snowmelt), "%Y-%j")) %>%
  setNames(c("year", "doy_snowmelt", "date_snowmelt"))

# Exprot
write.csv(snowmelt, "data/auxillary/snowmelt_qhi.csv",
          row.names = F)
