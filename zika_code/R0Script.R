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


# Calculating Mosquito Abundance
mosquito.occurrence <- mean.county(raster = perkins.tx, shapefile = texas.county, name = "occurrence")
mosquito.abundance <- -log(1-mosquito.occurrence)

#already on the county_ids
county_ids$mosquito.abundance <- mosquito.abundance ; colnames(county_ids$mosquito.abundance) <- ""

mosquito_range = ddply(county_ids,.variables = 'Metro', summarise, mean(mosquito.abundance))
mosquito_range = mosquito_range[-16, ]
colnames(mosquito_range) = c("MetroArea", "mosquito.abundance")

GDP_metro= ddply(.data = county_ids, .variables = 'Metro', summarise, mean(GDP_capita))
GDP_metro = GDP_metro[-16, ]


# Calculating Ranges of GDP to set Ranges of Multiplication Factor 
GDPrange <- range(log(county_ids$GDP_capita))
GDP <- log(county_ids$GDP_capita)-10

GDP_log <- log(GDP_metro$..1)-10

MF = -(1.7+.5/1.5*GDP_log)
metro_MF <- exp(MF)
m.to.h.metro <- mosquito_range$mosquito.abundance*metro_MF
metro <- seq(1:15)
m.t.h.metro <- data.frame(cbind(metro, m.to.h.metro))

#Assign the Highest MF to the lowest GDP County 
county_ids$MF <- exp(MF)

m.to.h <- county_ids$mosquito.abundance*county_ids$MF 

##### Function to calculate R0 
calculateR0 <- function(df, m.to.h) {
  mos.hum.transmission = .4
  c_over_r = 3.5
  alpha = .67
  mos.mortality = 1/25
  ext.incub = 5
  R0 = c_over_r * (m.to.h.metro*mos.hum.transmission*alpha^2*
                   exp(-mos.mortality*ext.incub))/mos.mortality
  df = cbind(df, R0)
  return(df)
}

R0_metro<- data.frame(cbind(metro, R0))
colnames(R0_metro) <- 

county_ids <- calculateR0(df = county_ids, m.to.h)


colnames(county_ids)[10] <- "R0"
county_ids$R0
hist(county_ids$R0)


#Mering Metro Area
county_ids$MetroR0 <- NA
for (i in 1:nrow(R0_metro)) {
  indices = which(county_ids$Metro == R0_metro$Metro[i])
  county_ids$MetroR0[indices] <- R0_metro$..1[i]
}


texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_ids, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]


plot <- ggplot()+geom_polygon(data = final.plot, aes_string(x="long", y = "lat", group = "group", fill = "MetroR0" ), color = "black", size = .25)+coord_map() +
  scale_fill_gradient2(name = "R0", limits = c(0,2), low = "yellow", midpoint = 1, mid = "red" ,high = "dark red",  na.value = "white")
 plot
  