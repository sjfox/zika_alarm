library(rgdal)
library(raster)
library(rgeos)
library(SDMTools)
library(sp)
library(maptools)
library(ggplot2)
library(cowplot)
library(reshape2)

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
county_master <- read.csv("Documents/zika_alarm/csvs/county_master.csv")
#mosquito.occurrence <- mean.county(raster = perkins.tx, shapefile = texas.county, name = "occurrence")
#mosquito.abundance <- -log(1-mosquito.occurrence)
#mosquito.abundance <- as.numeric(mosquito.abundance)
#mosquito.abundance <- county_master$mosquito.abundance
#county_master$mosquito.abundance <- mosquito.abundance
mosquito.abundance <- county_master$mosquito.abundance

#mosquito_abundance.metro = ddply(county_master,.variables = 'Metro', summarise, mean(mosquito.abundance))
#mosquito_abundance.metro <- mosquito_abundance.metro[-16,]
#GDP_metro = ddply(.data = county_master, .variables = 'Metro', summarise, mean(GDP))
#GDP_metro <- GDP_metro[-16,]

# Calculating Ranges of GDP to set Ranges of Multiplication Factor 
GDP <- log(county_master$GDP)-10
logMF.high <- -.9-(.07/.5)*GDP #High
logMF.expected <- -1.79-(.07/.5)*GDP #Expected
logMF.low <- -2.6-(.07/.5)*GDP #Low 
logMF.high.slope <- -1.35-(.9/.5)*GDP


MF.high <- exp(logMF.high)
MF.expected <- exp(logMF.expected)
MF.low <- exp(logMF.low)
MF.high.slope <- exp(logMF.high.slope)

mosquito.high <- MF.high * mosquito.abundance
mosquito.low <- MF.low * mosquito.abundance
mosquito.expected <- MF.expected * mosquito.abundance
mosquito.high.slope <- MF.high.slope * mosquito.abundance


#GDP <- log(GDP_metro$..1)-10
#logMF.metro <- -1.7-(.3/1)*GDP
#MF.metro <- exp(logMF.metro)
#mosquito.metro <- MF.metro * mosquito_abundance.metro$..1

##### Calculate R0 
mos.hum.transmission = .634
###c_over_r = 3.5
c_over_r = 9 *.77
alpha = .63
mos.mortality = 1/14



#RNOTT SCENARIOS 
for (i in 1:nrow(county_master)) {
rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
                                     exp(-mos.mortality*county_master$eip))/mos.mortality
}

for (i in 1:nrow(county_master)) {
  rnott.low<- c_over_r *(mosquito.low * mos.hum.transmission*alpha^2*
                                 exp(-mos.mortality*county_master$eip))/mos.mortality
}

for (i in 1:nrow(county_master)) {
  rnott.high <- c_over_r *(mosquito.high * mos.hum.transmission*alpha^2*
                                 exp(-mos.mortality*county_master$eip))/mos.mortality
}

for (i in 1:nrow(county_master)) {
  rnott.high.slope <- c_over_r *(mosquito.high.slope * mos.hum.transmission*alpha^2*
                             exp(-mos.mortality*county_master$eip))/mos.mortality
}



#### Combining to make dataframe of all R0s
county_master.old <- county_master
county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected
#county_master$rnott.temp <- NA; county_master$rnott.temp <- rnott.expected
county_master$rnott.low <- NA; county_master$rnott.low <- rnott.low
county_master$rnott.high <- NA; county_master$rnott.high <- rnott.high
county_master$rnott.high.slope <- NA; county_master$rnott.high.slope <- rnott.high.slope
county_master$rnott.expected.round <- round(county_master$rnott.expected, digits = 1)
unique(sort(county_master$rnott.expected.round))
write.csv(county_master, file= "Documents/zika_alarm/csvs/county_master.csv")
head(county_master)

#### Plotting 
fig_path <- "../ExploratoryFigures/"

texas.county <- readShapeSpatial('Documents/zika_alarm/TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
rnott_plot.m <- melt(data = county_master, id.vars = c("id", "Geography"), measure.vars = c("rnott.expected", "rnott.low", "rnott.high", "rnott.high.slope"))
colnames(rnott_plot.m) <- c("id", "Geography", "rnott.scenario", "rnott")


#breaks = seq(0,3,1)
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
merge.texas.county <- merge(texas.county.f, county_master, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]
head(rnott.data)
#cols <- rev(heat.colors(length(breaks.rnott)))

breaks.rnott <- seq(from = 0, to = 4, by = 1)
labels <- c(rnott.expected = "Medium", rnott.low = "Strong", rnott.high = "Weak", rnott.high.slope = "Mixed")

plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color = "black", size = .25) +
    facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    scale_fill_gradient(name = expression("R"[0]), low = "yellow", high = "dark red") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right") +
    theme(strip.text.x = element_text(size = 16)) 

plot.rnott
  
plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) + 
  geom_polygon(data = rnott.data, aes(group = group, fill = rnott.expected), color = "black", size = .25) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_gradient(name = expression("R"[0]), low = "yellow", high = "red") +
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=12, margin = margin(), debug = FALSE), legend.title = element_text(size = 18)) +
  theme(legend.key.size =  unit(0.5, "in")) 

save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/R0sensitive_temp.pdf", plot = plot.rnott, base_height = 8, base_aspect_ratio = 1.1)  


#### Histograms 
breaks <- c(0,1,6)
mycolors <- c("yellow", "red")

rnott_plot.m.copy <- rnott_plot.m
rnott_plot.m.copy$b <- cut(rnott_plot.m.copy$rnott, breaks)
range(rnott_plot.m.copy$rnott)
head(rnott_plot.m.copy$b)

hist.rnott <- ggplot(rnott_plot.m.copy, aes(rnott, fill = rnott_plot.m.copy$b)) +
  geom_histogram(binwidth = .1) +
  facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
  geom_vline(xintercept = 1, colour="red", linetype = "longdash") +
  scale_x_continuous(name = expression("R"[0]), breaks = seq(0,5,.5)) + 
  scale_y_continuous(name = "Count", breaks = seq(0,100, 10)) +
  scale_fill_manual(breaks = levels(rnott_plot.m.copy$b), values = mycolors,
                                    name = "Risk Level", drop = FALSE) 
+ 
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + 
  theme(strip.text.x = element_text(size = 16)) 

hist.rnott
  #theme(legend.position = "right") +
  theme(legend.text=element_text(size=16, margin = margin(), debug = FALSE), legend.title = element_text(size = 18)) +
  theme(legend.key.size =  unit(0.5, "in")) +
    #scale_fill_gradient(name = expression("R"[0]), low = "yellow", high = "red") +


head(rnott_plot.m)
save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/hist_R0.pdf", plot = hist.rnott, base_height = 8, base_aspect_ratio = 1.5)  


########### script to test the id of a county
county.test <- read.csv("../csvs/county_test.csv")
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_master, by = "id", all.x = TRUE)
test.plot <- merge.texas.county[order(merge.texas.county$id),]

test.sub.county <- test.plot[test.plot$Geography == "Medina County, Texas", ]


plot <- ggplot(test.plot, aes(x = long, y = lat)) + geom_polygon(data = test.plot, aes(group = group, fill = "red"), color = "black", size = .25)
layer <- plot + geom_polygon(data =test.sub.county, aes(x = long, y = lat, group = group, fill = "grey"))
save_plot(filename = paste0(fig_path, "testing_id.pdf"), plot = layer, base_height = 8, base_aspect_ratio = 1.1)


library(weathermetrics)
temperature <- county_master$temperature
county_master$temperature.c <- fahrenheit.to.celsius(county_master$temperature, round = 0)
temperature.c
unique(temperature.c)


    