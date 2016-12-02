rm(list=ls())
library(rgdal)
library(raster)
library(rgeos)
library(SDMTools)
library(maptools)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)

# Script to take in data sources and create county master list
# 0. Read in county name and id 
# 1. Mosquito abundace
# 2. County GDP
# 3. Importation 
# 4. Clean Texas County Temperature 


county_parms <- read.csv("../csvs/county_id.csv")

county_parms  %>% 
  separate(col = Geography, into = c("county", "state"), sep = -15) %>%
  mutate(county = tolower(county)) %>%
  dplyr::select(id, county) -> county_parms

# 1. Mosquito Abundance

# Function for calculating county averages
mean.county <- function(raster, shapefile, name) { 
  raster.polygon <- extract(raster, shapefile)
  raster.mean <- data.frame(unlist(lapply(raster.polygon, FUN = mean, na.rm = TRUE)))
  colnames(raster.mean) <- name
  raster.mean <- data.frame(raster.mean)
  return(raster.mean)
}

# Reading in texas aegypti occurence data and texas shape files 
texas.aegypti <- raster("../data/texas.aegypti.asc")
texas.county <- readOGR(dsn = "../TexasCountyShapeFiles", layer = "texas.county")

# calculating mean mosquito occurrence for each county, the shapefiles are in order of county id 
mosquito.occurrence <- mean.county(raster = texas.aegypti, shapefile = texas.county, name = "occurrence") 

# Converting to proxy for mosquito abundance 
mosquito.abundance <- -log(1-mosquito.occurrence)
mosquito.abundance <- mosquito.abundance
county_parms$mosquito.abundance <- mosquito.abundance$occurrence


# 2. County GDP 
economic_data <- read.csv("../csvs/US_calc_xi_090110.csv")
economic_data %>% 
  filter(state == "Texas") %>% 
  dplyr::select(Area.Name, GDP_PC105) %>% 
  mutate(county = tolower(Area.Name)) %>% 
  mutate(gdp = as.numeric(paste(GDP_PC105)) * 10e5) %>%
  dplyr::select(county, gdp) -> economic_data_tx

county_parms <- left_join(county_parms, economic_data_tx)

# 3. County Importation
importation <- read.csv("../csvs/zika_importation_data.csv")
importation %>%
  dplyr::select(Geography, probability) %>%
  mutate(county = tolower(Geography)) %>%
  mutate(importation.projected = round(round(probability, digits = 3)*.9, digits = 3)) %>% 
  dplyr::select(county, importation.projected) -> importation.prob

county_parms <- left_join(county_parms, importation.prob)
write.csv(county_parms, file = "../csvs/county_r0_parameters.csv", row.names = FALSE)
  
# 4. Texas Weather 
temps <- read_tsv("../data/weather/tx_county_temps.txt")

temps <- temps %>% filter(is.na(Notes)) %>%
  rename(avg_max_temp = `Avg Daily Max Air Temperature (C)`,
         avg_min_temp=`Avg Daily Min Air Temperature (C)`) %>%
  mutate(avg_temp = (avg_max_temp + avg_min_temp) / 2,
         county = tolower(str_replace_all(County, pattern = " County, TX", ""))) %>%
  arrange(`Month Code`) %>%
  mutate(Month = factor(Month, levels = unique(Month))) %>%
  dplyr::select(county, month=Month, avg_temp) %>%
  spread(key = month, value = avg_temp)

write.csv(temps, file = "../csvs/tx_county_temps.csv", row.names = FALSE)
