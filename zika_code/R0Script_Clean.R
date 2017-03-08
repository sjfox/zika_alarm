library(plyr)
library(tidyverse)
library(rgdal)
library(raster)
library(rgeos)
library(SDMTools)
library(sp)
library(maptools)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)
library(scam)
library(scales)
library(stringr)
###################################
# This script is to use after going through calculate_perkinsr0.
# Different starting points based on whether you're normalizing rnots or just plotting
# Script is separated into two parts: (1) Multiple months and (2) Single Month

#######################################################
# Part 1 - MULTIPLE MONTHS 
# Step 2: Script to read csvs of rnots and plot across 
# Texas based on hypothetical ranges 
########################################################
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")
county.parms <- read.csv("../csvs/county_r0_parameters.csv")
fig_path <- "../ExploratoryFigures/"
csv_path <- "../csvs/rnot_temp/"
fileList <- list.files(path=csv_path, pattern="r0.csv")

# Read in csvs 
allData <- lapply(fileList, function(.file){
  dat <-read.csv(paste0(csv_path,.file), header=T)
  dat$month <-as.character(.file)
  dat    # return the dataframe
})
# combine into a single dataframe
rnot.temps <- do.call(rbind, allData)
colnames(rnot.temps)[1:1000] <- seq(1:1000)
colnames(rnot.temps)[1001:1002] <- c("county", "file")

#---------------------Normalizing Across the 6 Months 
# Calculate Median of Each across the Months 
# Avg for each rep across the 6 months

rnot.temps %>% 
  separate(file, c("month", "file"), "_") %>%
  gather(key = rep, value = rnot, 1:1000) %>%
  group_by(county, month) %>%
  summarise(monthly.median = median(rnot)) %>%
  ungroup() %>%
  mutate(normalized.rnot = (monthly.median - min(monthly.median))/(max(monthly.median)-min(monthly.median))) %>%
  mutate(hypo.max.1.5 = normalized.rnot *1.5) -> rnot.norm.hypo

#---------------------Plotting  

rnot.norm.hypo %>% 
  gather(key = metric, value = rnot, -county, -month) %>%
  left_join(y = county.parms, by = "county") -> plot.data
plot.data$month <- factor(plot.data$month, 
                                      levels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct"))

texas.county <- readOGR(dsn = "../TexasCountyShapeFiles", layer = "texas.county")
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, plot.data, by = "id", all.x = TRUE)
rnot.data <- merge.texas.county[order(merge.texas.county$id),]

rnot.data %>% filter(metric == "normalized.rnot") %>%
  ggplot(aes(x=long, y = lat)) +
  geom_polygon(aes(group = group, fill = rnot), color = "grey", size = .1) +
  facet_wrap(facets = "month", nrow = 2) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_gradient(name = "Relative Transmission Risk", low = "white", high = "darkseagreen") +
  theme_cowplot() + 
  theme(strip.background = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(), 
        line = element_blank(),
        strip.text.x = element_text(size = 14),
        legend.position = "bottom", 
        legend.key.size =  unit(0.4, "in")) -> plot.months

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/plot.months.pdf", plot = plot.months,  base_height = 8, base_aspect_ratio = 1.5) 

# ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))
# geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
# geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
# size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +


######################################################
# Section 2: ONE MONTH 
###################################################
rnot.ests$county = county.parms$county

rnot.ests %>% gather(key = rep, value = rnot, -county) %>%
  group_by(county) %>%
  summarise(median = median(rnot),
            lowerbound = quantile(rnot, probs = .025),
            upperbound = quantile(rnot, probs = .975)) -> rnot.ests.summary

# ------------- Normalization ---------------------
# Normalizing just the median

rnot.ests.summary %>%
  mutate(normalized.rnot = (median - min(median))/(max(median)-min(median))) %>%
  mutate(hypo.max.1.5 = normalized.rnot *1.5) -> rnot.norm.1.5

#-------------- Plotting ---------------------
# Choose whether plotting normalized or raw values
left_join(county.parms, rnot.ests.summary, by = "county") %>%
  gather(quantile, value = rnot, 6:8) -> plot.data
texas.county <- readOGR(dsn = "../TexasCountyShapeFiles", layer = "texas.county")
texas.county.f <- fortify(texas.county, region = "ID")


merge.texas.county <- merge(texas.county.f, plot.data , by = "id", all.x = TRUE)
rnot.data <- merge.texas.county[order(merge.texas.county$id),]
colorends <- c("white", "blue", "yellow", "red")
gradientends <- c(0,1,1.01,max(rnot.data$rnot))

plot.range <- ggplot(rnot.data, aes(x=long, y = lat)) +
  geom_polygon(data = rnot.data, aes(group = group, fill = round(rnot, digits = 1)), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  facet_wrap(~quantile) +
  scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) + 
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  theme(legend.position = "bottom", legend.key.size =  unit(0.5, "in")) +
  theme(strip.background = element_blank(),strip.text.x = element_blank())

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/plot.range.pdf", plot = plot.range,base_height = 4, base_aspect_ratio = 2.0) 

# Code for including points on map 
#ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))
# geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
# geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
# size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +



####################### Econ Plots ###################################################
#### Ploting relationship for SCAM functions 
econ <- seq(6, 15, by=0.1)
plot(econ, predict(scam.est.list[[1]], newdata=data.frame(econ=econ)), type="l", ylim=c(-2, 5))
for(ii in 2:1000){
  lines(econ, predict(scam.est.list[[ii]], newdata=data.frame(econ=econ)), type="l")
}


