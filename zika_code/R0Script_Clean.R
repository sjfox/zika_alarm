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

#######################################################
# Step 2: # Script to read csvs of rnots and plot across 
# Texas based on hypothetical ranges 
########################################################
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")

county_parms <- read.csv("../csvs/county_r0_parameters.csv")

fig_path <- "../ExploratoryFigures/"
csv_path <- "../csvs/rnot_temp/"

fileList <- list.files(path=csv_path, pattern="r0.csv")

# Read in csvs - Currently just all available csvs 
allData <- lapply(fileList, function(.file){
  dat <-read.csv(paste0(csv_path,.file), header=T)
  dat$month <-as.character(.file)
  dat    # return the dataframe
})

# combine into a single dataframe
r0.temps <- do.call(rbind, allData)
colnames(r0.temps)[1:1000] <- seq(1:1000)

# For each county calculate the average of each rep across the months
# Avg for each rep across the 6 months
r0.temps %>% 
  gather(key = rep, value = rnot, 1:1000) %>%
  group_by(county, rep) %>%
  summarise(rep.6median.rnot = median(rnot)) -> rnot.rep.sixmonth 
# Avg over the reps 
rnot.rep.sixmonth %>%
  group_by(county) %>%
  summarise(final.rnot = median(rep.6median.rnot)) -> rnot.sixmonth

r0.plot <- left_join(county_parms, rnot.sixmonth)

# Normalized absolute rnot values and calculated hypothetical rnots based on 
# range 
r0.plot %>%
  mutate(normalized.rnot = (final.rnot - min(final.rnot))/(max(final.rnot)-min(final.rnot))) %>%
  mutate(hypo.max.2 = normalized.rnot * 2) %>% 
  mutate(hypo.max.1.5 = normalized.rnot *1.5) -> r0.hypo

#########################
# PLOTTING 
#########################
texas.county <- readOGR(dsn = "../TexasCountyShapeFiles", layer = "texas.county")
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, r0.hypo, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]


colorends <- c("white", "darkseagreen", "yellow", "red")
gradientends <- c(0,1,1.01,max(rnott.data.decrease$hypo.max.1.5))
gradientends <- c(0,.25,.5,.75)

ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))

## 
plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = round(normalized.rnot, digits = 1)), color = "grey", size = .1) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
  #  scale_fill_gradient(name = expression("R"[0]), low = "white", high = "darkseagreen") +
    # geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
    # geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  # size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right", legend.key.size =  unit(0.5, "in")) +
    theme(strip.text.x = element_text(size = 14)) 

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/relativernot.pdf", plot = plot.rnott,  base_height = 8, base_aspect_ratio = 1.1) 

####################### Econ Plots ###################################################
#### Ploting relationship for SCAME functions 
econ <- seq(6, 15, by=0.1)
plot(econ, predict(scam.est.list[[1]], newdata=data.frame(econ=econ)), type="l", ylim=c(-2, 5))
for(ii in 2:1000){
  lines(econ, predict(scam.est.list[[ii]], newdata=data.frame(econ=econ)), type="l")
}

#######################R0 Sensitivity R0 Calculations #################################

# logMF.high <- -.9-(.07/.5)*GDP #High
# logMF.low <- -2.6-(.07/.5)*GDP #Low 
# logMF.high.slope <- -1.35-(.9/.5)*GDP



# # Function to calculate the mean value of all cells 
# # in a shapefile 
# mean.county <- function(raster, shapefile, name) { 
#   raster.polygon <- extract(raster, shapefile)
#   raster.mean <- data.frame(unlist(lapply(raster.polygon, FUN = mean, na.rm = TRUE)))
#   colnames(raster.mean) <- name
#   raster.mean <- data.frame(raster.mean)
#   return(raster.mean)
# }

###########################################################################
# #########################################################
# # Step 1: Getting mosquito abundance data for each county 
# #########################################################
# 
# # Importing the aegypti occurrence projections from Kraemer 
# tx.aegypti <- raster("../data/mosquito/texas.aegypti.asc")  
# 
# texas.county <- readShapeSpatial('../TexasCountyShapeFiles/texas.county.shp',
#                                  proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# ## setting the projection to match Texas shape files  
# # bb <- extent(-180, 180, -90, 90)
# # extent(world.aegypti) <- bb
# # world.aegypti <- setExtent(world.aegypti, bb, keepres=TRUE)
# # projection(world.aegypti) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")
# # 
# # # Cutting the world layer to Texas 
# # world.cut <- crop(world.aegypti, extent(texas.county))
# # texas.cut <- mask(world.cut, texas.county)
# 
# # Calculating Mosquito Occurrence for each county as the 
# # average occurrence for each cell in a county 
# mosquito.occurrence <- mean.county(raster = tx.aegypti, shapefile = texas.county, name = "occurrence") 
# 
# # Converting to proxy 
# mosquito.abundance <- -log(1-mosquito.occurrence)
# mosquito.abundance <- as.numeric(mosquito.abundance)


#######################################################
# Step 3: Calculate R0 
#######################################################

# R0 equation is set to use the August eip from the county master 
# Baseline parameters 

#mos.hum.transmission = .634
#c_over_r = 9 *.77
#alpha = .63
#mos.mortality = 1/14
#eip = county_master$aug.eip


#rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
#                               exp(-mos.mortality*eip))/mos.mortality

#county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected





