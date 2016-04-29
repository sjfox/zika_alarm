library(rgdal)
library(raster)
library(rgeos)
library(SDMTools)
library(sp)
library(maptools)


# Cutting World Layer to Texas
county_master <- read.csv("../csvs/county_master.csv")
world.aegypti <- raster("~/Downloads/aegypti.tif.ovr")
bb <- extent(-180, 180, -90, 90)
extent(world.aegypti) <- bb
world.aegypti <- setExtent(world.aegypti, bb, keepres=TRUE)
projection(world.aegypti) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

world.cut <- crop(world.aegypti, extent(texas.county))
perkins.tx <- mask(world.cut, texas.county)


# Calculating Mosquito Abundance
mosquito.occurrence <- mean.county(raster = perkins.tx, shapefile = texas.county, name = "occurrence")
mosquito.abundance <- -log(1-mosquito.occurrence)
mosquito.abundance <- as.numeric(mosquito.abundance)
county_master$mosquito.abundance <- mosquito.abundance


mosquito_abundance.metro = ddply(county_master,.variables = 'Metro', summarise, mean(mosquito.abundance))
mosquito_abundance.metro <- mosquito_abundance.metro[-16,]
GDP_metro = ddply(.data = county_master, .variables = 'Metro', summarise, mean(GDP))
GDP_metro <- GDP_metro[-16,]

# Calculating Ranges of GDP to set Ranges of Multiplication Factor 
GDP <- log(county_master$GDP)-10
logMF.high <- -1-(.3/1)*GDP #High
logMF.expected <- -1.7-(.3/1)*GDP #Expected
logMF.low <- -2.1-(.3/1)*GDP #Low 
logMF.high.slope <- -1-(1/1)*GDP


MF.high <- exp(logMF.high)
MF.expected <- exp(logMF.expected)
MF.low <- exp(logMF.low)
MF.high.slope <- exp(logMF.high.slope)

mosquito.high <- MF.high * mosquito.abundance
mosquito.low <- MF.low * mosquito.abundance
mosquito.expected <- MF.expected * mosquito.abundance
mosquito.high.slope <- MF.high.slope * mosquito.abundance


GDP <- log(GDP_metro$..1)-10
logMF.metro <- -1.7-(.3/1)*GDP
MF.metro <- exp(logMF.metro)
mosquito.metro <- MF.metro * mosquito_abundance.metro$..1

##### Calculate R0 
mos.hum.transmission = .4
c_over_r = 3.5
alpha = .67
mos.mortality = 1/25 
ext.incub = 5


#RNOTT SCENARIOS 
rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
                                     exp(-mos.mortality*ext.incub))/mos.mortality

rnott.high = c_over_r * (mosquito.high*mos.hum.transmission*alpha^2*
                                 exp(-mos.mortality*ext.incub))/mos.mortality

rnott.low = c_over_r * (mosquito.low*mos.hum.transmission*alpha^2*
                                exp(-mos.mortality*ext.incub))/mos.mortality
rnott.high.slope = c_over_r * (mosquito.high.slope*mos.hum.transmission*alpha^2*
                                 exp(-mos.mortality*ext.incub))/mos.mortality

rnott.metro = c_over_r * (mosquito.metro*mos.hum.transmission*alpha^2*
                            exp(-mos.mortality*ext.incub))/mos.mortality
metro <- seq(1:15)
metro.rnott <- cbind(metro, round(rnott.metro, digits = 1)); colnames(metro.rnott) <- c("metro", "rnott.metro")
county_master<- merge(x = county_master, y = metro.rnott[,c("metro", "rnott.metro")], 
                       by.x="Metro", by.y="metro", all.x=TRUE, sort=FALSE)

indices = which(is.na(county_master$Metro)) 
county_master[indices, "rnott.metro"] <- county_master[indices, "rnott.expected.round"]



#### Combining to make dataframe of all R0s
county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected
county_master$rnott.low <- NA; county_master$rnott.low <- rnott.low
county_master$rnott.high <- NA; county_master$rnott.high <- rnott.high
county_master$rnott.high.slope <- NA; county_master$rnott.high.slope <- rnott.high.slope
county_master$rnott.expected.round <- round(county_master$rnott.expected, digits = 1)
unique(sort(county_master$rnott.expected.round))
write.csv(county_master, file= "../csvs/county_master.csv")


#### Plotting 
fig_path <- "../ExploratoryFigures/"

texas.county <- readShapeSpatial('../TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
rnott_plot.m <- melt(data = county_master, id.vars = c("id", "Geography"), measure.vars = c("rnott.expected", "rnott.low", "rnott.high", "rnott.high.slope"))
colnames(rnott_plot.m) <- c("id", "Geography", "rnott.scenario", "rnott")


#breaks = seq(0,3,1)
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]
cols <- rev(heat.colors(length(breaks.rnott)))

breaks.rnott <- seq(from = 0, to = 4, by = 1)
labels <- c(rnott.expected = "Expected", rnott.low = "Low", rnott.high = "High", rnott.high.slope = "Mixed")
sp + facet_grid(. ~ sex, labeller=labeller(sex = labels))

plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color = "black", size = .25) +
    facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    scale_fill_gradient(name = expression("R"[0]), low = "yellow", high = "red") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right") +
    theme(legend.text=element_text(size=12, margin = margin(), debug = FALSE), legend.title = element_text(size = 18)) +
    theme(legend.key.size =  unit(0.5, "in")) 
  
save_plot(filename = paste0(fig_path, "supplement_R0.pdf"), plot = plot.rnott, base_height = 8, base_aspect_ratio = 1.1)  




########### script to test the id of a county
county.test <- read.csv("../csvs/county_test.csv")
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_master, by = "id", all.x = TRUE)
test.plot <- merge.texas.county[order(merge.texas.county$id),]

test.sub.county <- test.plot[test.plot$Geography == "Medina County, Texas", ]


plot <- ggplot(test.plot, aes(x = long, y = lat)) + geom_polygon(data = test.plot, aes(group = group, fill = "red"), color = "black", size = .25)
layer <- plot + geom_polygon(data =test.sub.county, aes(x = long, y = lat, group = group, fill = "grey"))
save_plot(filename = paste0(fig_path, "testing_id.pdf"), plot = layer, base_height = 8, base_aspect_ratio = 1.1)


    