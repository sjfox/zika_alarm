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

plot(log(castro.tx_raw))

castro.tx_raw <- crop(castro.tx_raw, extent(texas.county))
castro.tx_raw <- mask(castro.tx_raw, texas.county)
cellStats(x = castro.tx_raw, stat = "sum")

castro.mean_raw <- mean.county(raster = castro.tx_raw, shapefile = texas.county, name = "habitat.castro")

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
  R0 = c_over_r * (m.to.h.metro*mos.hum.transmission*alpha^2*
                   exp(-mos.mortality*ext.incub))/mos.mortality
  metro = cbind(metro, R0)
  return(df)
}



metroR0 <-metro

#Merging Metro Area
county_ids$MetroR0 <- NA


for (i in 1:nrow(metroR0)) {
  indices = which(county_ids$Metro == metroR0[i,1])
  county_ids$MetroR0[indices] <- metroR0[i,2]
}


indices = which(is.na(county_ids$Metro)) 
county_ids[indices, 13] <- county_ids[indices, 11]
county_ids$metro_round <- round(county_ids$MetroR0, digits = 1)
write.csv(county_ids, file = "county_ids.csv")

#### Plotting 

#Importation
county_ids$castro_R0 <- R0
county_ids$castro_combined <- R0 * county_ids$probability

relativeR0 <- castro.mean_raw/as.data.frame(county_ids$Population.Proportion)
county_ids$relativeR0 <- NA

county_ids <- cbind(county_ids, relativeR0)



county_ids$raw_combined <- county_ids$probability * county_ids$habitat.castro


texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_ids, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]


plot <- ggplot()+geom_polygon(data = final.plot, aes(x=long, y = lat, group = group, fill = metro_round), color = "black", size = .25) + 
  coord_map() + 
  scale_fill_gradient2(name = "R0", low = "yellow", midpoint = 1, mid = "red" ,high = "dark red",  na.value = "white") +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  labs(x=NULL, y = NULL)+
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) 
  #geom_polygon(data=lamb_county, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)

plot

#lamb_county <- final.plot[which(final.plot$rownames == 141),]
  
 


sort(unique(county_ids$metro_round))
  