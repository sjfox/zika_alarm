library(rgdal)
library(raster)
library(rgeos)
library(SDMTools)
library(sp)
library(maptools)



# Cutting World Layer to Texas
setwd(dir = "Downloads/")
world.aegypti <- raster("aegypti.tif.ovr")
bb <- extent(-180, 180, -90, 90)
extent(world.aegypti) <- bb
world.aegypti <- setExtent(world.aegypti, bb, keepres=TRUE)
projection(world.aegypti) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")


world.cut <- crop(world.aegypti, extent(texas.county))
perkins.tx <- mask(world.cut, texas.county)

castro.tx_raw <- read.asc("aedes_aegypti_avg.asc")
castro.tx_raw <- raster.from.asc(castro.tx_raw)
plot(castro.tx_raw)
texas

plot(log(castro.tx_raw))

castro.tx_raw <- crop(castro.tx_raw, extent(texas.county))
castro.tx_raw <- mask(castro.tx_raw, texas.county)
cellStats(x = castro.tx_raw, stat = "sum")

castro.mean_raw <- mean.county(raster = castro.tx_raw, shapefile = texas.county, name = "habitat.castro")

castro.county <- sum.county(raster = castro.tx, shapefile = texas.county, name = "mosquito.county")



castro.mosquito <- -log(1-castro.mean)




# Calculating Mosquito Abundance
mosquito.occurrence <- mean.county(raster = perkins.tx, shapefile = texas.county, name = "occurrence")
mosquito.abundance <- -log(1-mosquito.occurrence)

#already on the county_ids
county_ids$mosquito.abundance <- mosquito.abundance ; colnames(county_ids$mosquito.abundance) <- ""



### SEPARATING RO INTO METRO AREAS 
mosquito_abundance = ddply(county_ids,.variables = 'Metro', summarise, mean(mosquito.abundance))
mosquito_abundance = mosquito_abundance[-16, ]
colnames(mosquito_abundance) = c("MetroArea", "mosquito.abundance")


GDP_metro= ddply(.data = county_ids, .variables = 'Metro', summarise, mean(GDP_capita))
GDP_metro = GDP_metro[-16, ]


# Calculating Ranges of GDP to set Ranges of Multiplication Factor 
GDP <- log(county_ids$GDP_capita)-10
logMF <- -1.7-(.3/1)*GDP
MF <- exp(logMF)
m.to.h.full <- MF * county_ids$mosquito.abundance
m.to.h.castro <- MF * castro.mosquito

county_ids$m.to.h.full <- m.to.h.full
county_ids$m.to.h.castro <- m.to.h.castro

GDP_metro <- log(GDP_metro$..1)-10

logMF = -1.7 - (.3/1)*(GDP_metro)

metro_MF <- exp(logMF)

m.to.h.metro <- mosquito_abundance$mosquito.abundance*metro_MF



metro <- seq(1:15)
m.t.h.metro <- data.frame(cbind(metro, m.to.h.metro))

##### Function to calculate R0 
calculateR0 <- function(metro, m.to.h) {
  mos.hum.transmission = .4
  c_over_r = 3.5
  alpha = .67
  mos.mortality = 1/25
  ext.incub = 5
  R0 = c_over_r * (m.to.h.castro*mos.hum.transmission*alpha^2*
                   exp(-mos.mortality*ext.incub))/mos.mortality
  metro = cbind(metro, R0)
  return(df)
}



metroR0 <- cbind(metro, R0)

#Merging Metro Area
county_ids$MetroR0 <- NA


for (i in 1:nrow(metroR0)) {
  indices = which(county_ids$Metro == metroR0[i,1])
  county_ids$MetroR0[indices] <- metroR0[i,2]
}


indices = which(is.na(county_ids$Metro)) 
county_ids[indices, 13] <- county_ids[indices, 11]
county_ids$metro_round <- round(county_ids$MetroR0, digits = 1)

#### Plotting 

#Importation
county_ids$castro_R0 <- R0
county_ids$castro_combined <- R0 * county_ids$probability

county_ids <- cbind(county_ids, castro.county)
county_ids$relativeR0 <-county_ids$mosquito.county/county_ids$Population.Proportion
county_ids$Combined.Castro <- county_ids$relativeR0*county_ids$probability


county_ids$raw_combined <- county_ids$probability * county_ids$habitat.castro




texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_ids, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]

hist(county_ids$Combined.Castro, breaks = 100)
breaks <- c(seq(from = 0, to = .1, by =  .01), .2, .3, .4)
plot <- ggplot()+geom_polygon(data = final.plot, aes(x=long, y = lat, group = group, fill = Combined.Castro), color = "black", size = .25) + 
  coord_map() +   scale_fill_gradient(name = "Combined", low = "yellow", high = "red")
  
  #scale_fill_gradient2(name = "R0", low = "yellow", midpoint = 1, mid = "red" ,high = "dark red",  na.value = "white")
plot
 
sort(unique(county_ids$relativeR0))

plot
sort(unique(county_ids$metro_round))
  