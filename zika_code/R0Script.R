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



# Calculating Ranges of GDP to set Ranges of Multiplication Factor 

# First calculate the R0 for each county
GDP <- log(county_ids$GDP_capita)-10
#logMF <- -1.5-(.5/1)*GDP #High
logMF <- -1.7-(.3/1)*GDP #Expected
#logMF <- -1.9-(.3/1.5)*GDP #Low 

MF <- exp(logMF)
m.to.h.full <- MF * r0_low$mosquito.abundance

county_ids$m.to.h.full <- m.to.h.full
r0_expected$m.to.h.full <- m.to.h.full


### SEPARATING RO INTO METRO AREAS 
#When calculating the metro R0 
mosquito_abundance.metro = ddply(r0_high,.variables = 'Metro', summarise, mean(mosquito.abundance))
mosquito_abundance.metro = mosquito_abundance[-16, ]
colnames(mosquito_abundance.metro) = c("MetroArea", "mosquito.abundance")


GDP_metro = ddply(.data = r0_low, .variables = 'Metro', summarise, mean(GDP_capita))
GDP_metro = GDP_metro[-16, ]


GDP_metro <- log(GDP_metro$..1)-10
logMF.metro = -1.7 - (.3/1)*(GDP_metro)
metro_MF <- exp(logMF.metro)
m.to.h.metro <- mosquito_abundance.metro$mosquito.abundance*metro_MF
metro <- seq(1:15)
m.t.h.metro <- data.frame(cbind(metro, m.to.h.metro))




##### Function to calculate R0 
calculateR0 <- function(metro, m.to.h) {
  mos.hum.transmission = .4
  c_over_r = 3.5
  alpha = .67
  mos.mortality = 1/31 # changing this to 31 
  ext.incub = 5
  R0 = c_over_r * (m.to.h.metro*mos.hum.transmission*alpha^2*
                   exp(-mos.mortality*ext.incub))/mos.mortality
  metro = cbind(metro, R0)
  return(df)
}

r0_expected$R0 <- R0

metroR0 <-metro
r0_expected$MetroR0 <- NA

#Merging Metro Area
county_ids$MetroR0 <- NA


for (i in 1:nrow(metroR0)) {
  indices = which(r0_expected$Metro == metroR0[i,1])
  r0_expected$MetroR0[indices] <- metroR0[i,2]
}


indices = which(is.na(county_ids$Metro)) 
indices = which(is.na(r0_expected$Metro))
county_ids[indices, 13] <- county_ids[indices, 11]
r0_expected[indices, "MetroR0"] <- r0_expected[indices, "R0"]
county_ids$metro_round <- round(county_ids$MetroR0, digits = 1)
r0_expected$MetroR0_round <- round(r0_expected$MetroR0, digits = 1)

write.csv(r0_low, file = "r0_expected.csv")


#### Plotting 
#comparing R0s 
r0_compare <- data.frame(cbind(r0_expected$Old_id, r0_expected$id, r0_expected$rownames, r0_expected$MetroR0_round, 
                    r0_low$MetroR0_round, r0_high$MetroR0_round))
colnames(r0_compare) <- c("old_id", "id", "rownames", "expectedR0", "lowR0", "highR0")

#Importation When comparing to Castro suitability maps
#county_ids$castro_R0 <- R0
#county_ids$castro_combined <- R0 * county_ids$probability

#relativeR0 <- castro.mean_raw/as.data.frame(county_ids$Population.Proportion)
#county_ids$relativeR0 <- NA

#county_ids <- cbind(county_ids, relativeR0)



#county_ids$raw_combined <- county_ids$probability * county_ids$habitat.castro

hist(r0_compare$highR0); max(r0_compare$highR0)
hist(r0_compare$lowR0)
hist(r0_compare$expectedR0)
breaks = seq(0,3,1)
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
merge.texas.county <- merge(texas.county.f, r0_compare, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]

breaks = seq(from = 0, to = 2, by = .2)
plot <- ggplot()+geom_polygon(data = final.plot, aes(x=long, y = lat, group = group, fill = R0_round), color = "black", size = .25) + 
  coord_map() + 
  scale_fill_gradient2(name = expression("R"[0]), low = "yellow", mid = "red",
                      high = "dark red",  midpoint = 1, na.value = "white", breaks = breaks, limits = c(0,2)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12, margin = margin(), debug = FALSE), legend.title = element_text(size = 18)) +
  theme(legend.key.size =  unit(0.5, "in")) 
  #geom_polygon(data=lamb_county, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)

plot


sort(r0_compare$expectedR0)
#lamb_county <- final.plot[which(final.plot$rownames == 141),]
  
 


sort(unique(county_ids$metro_round))
  