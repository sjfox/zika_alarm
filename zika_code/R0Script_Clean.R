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

county_parms <- read.csv("../csvs/county_r0_parameters.csv")


#######################################################
# Step 2: Read in Perkins functions and parameters
# Perkins 2016 Nature Microbiology
########################################################
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")
rnot_calc <- function(mosq_abundance, gdp, temperature, 
                      a, b, c.r, mort.fun, eip.fun, scam.est){
  # Function that returns a single rnot estimate from single values
  require(scam)
  # Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
  g <- 1 / mort.fun(temperature)
  e <- eip.fun(temperature)
  gdp_scaling <- predict(scam.est, newdata=data.frame(econ=log(gdp)))
  as.numeric(gdp_scaling * mosq_abundance * a ^ 2 * b * c.r * exp(-g * e) / g)
}

rnot_calc_dist <- function(mosq_abundance, gdp, temperature,
                      a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function that returns the full distribution of R0s for single set of abundances gdp and temperature
  require(scam)
  # Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
  args <- list(mort.fun = mort.fun.list,
                eip.fun = eip.fun.list,
                scam.est = scam.est.list)
  unlist(purrr::pmap(args, rnot_calc, a=a, b=b, c.r=c.r, 
                     mosq_abundance=mosq_abundance, gdp=gdp, temperature=temperature))
}

rnot_calc_counties <- function(mosq_abundance, gdp, temperature,
                               a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function returns median, hi, low r0 for vector of abundances/gdp/temperatures
  args <- list(mosq_abundance = as.list(mosq_abundance),
               gdp = as.list(gdp),
               temperature = as.list(temperature))
  
  cty_rnot_dists <- purrr::pmap(args, rnot_calc_dist, a=a, b=b, c.r=c.r, 
              mort.fun.list=mort.fun.list, eip.fun.list=eip.fun.list,scam.est.list=scam.est.list)
  # browser()
  quantile_df <- function(x, probs, na.rm =F, names = F){
    z <- quantile(x, probs, na.rm, names)
    df <- data.frame(low_rnot=numeric(1), median_rnot = numeric(1), high_rnot=numeric(1))
    df[1,] <- z
    df
  }
  cty_rnot_dists %>% purrr::map( ~ quantile_df(.x, probs = c(0.025, 0.5, 0.975))) %>%
    bind_rows() %>%
    mutate(low_rnot = ifelse(low_rnot<0, 0, low_rnot))
}


rnot_ests <- rnot_calc_counties(county_parms$mosquito.abundance, 
                   county_parms$GDP, 
                   county_parms$August, 
                   a=a, b=b, c.r=c.r, 
                   mort.fun.list=mort.fun, 
                   eip.fun.list=eip.fun, 
                   scam.est.list=scam.est.list)

rnot_data <- cbind(county_parms, rnot_ests)

#######################################################
# Step 3: Calculate R0 
#######################################################

# R0 equation is set to use the August eip from the county master 
# Baseline parameters 

mos.hum.transmission = .634
c_over_r = 9 *.77
alpha = .63
mos.mortality = 1/14
eip = county_master$aug.eip


rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
                                     exp(-mos.mortality*eip))/mos.mortality

county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected



#####################################
# Plotting 
####################################

fig_path <- "../ExploratoryFigures/"


texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, rnot_data, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]


colorends <- c("white", "darkseagreen", "yellow", "red")
gradientends <- c(0,1,1.01,max(rnott.data$median_rnot))

ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))

## 
plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = low_rnot), color = "grey", size = .1) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
    # geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
    # geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  # size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right", legend.key.size =  unit(0.5, "in")) +
    theme(strip.text.x = element_text(size = 14)) 

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/rnott.pdf", plot = plot.rnott,  base_height = 8, base_aspect_ratio = 1.1) 


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
# 
