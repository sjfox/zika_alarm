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
county_sensitive <- county_master
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
range(MF.expected)
MF.low <- exp(logMF.low)
MF.high.slope <- exp(logMF.high.slope)

mosquito.high <- MF.high * mosquito.abundance
mosquito.low <- MF.low * mosquito.abundance
mosquito.expected <- MF.expected * mosquito.abundance
mosquito.high.slope <- MF.high.slope * mosquito.abundance


##### Calculate R0 
mos.hum.transmission = .85
c_over_r = 9 *.77
alpha = .63
mos.mortality = 1/14



#RNOTT SCENARIOS 

rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
                                     exp(-mos.mortality*county_master$aug.eip))/mos.mortality

#### Combining to make dataframe of all R0s
county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected
write.csv(county_master, file= "Documents/zika_alarm/csvs/county_master.csv")


############### Plotting #################################################
fig_path <- "../ExploratoryFigures/"

rnott_plot.m <- melt(data = county_master, id.vars = c("id", "Geography"), measure.vars = c("rnott.expected", "rnott.low", "rnott.high", "rnott.high.slope"))
colnames(rnott_plot.m) <- c("id", "Geography", "rnott.scenario", "rnott")

texas.county <- readShapeSpatial('Documents/zika_alarm/TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]


colFunc <- colorRampPalette(c("red4","darkorange", "gold", "light yellow"))
breaks.trigger = seq(from = 0, to = max(rnott.data$rnott), by  = 0.5)
head(rnott.data)

labels <- c(rnott.expected = "Expected GDP Effect", rnott.low = "Stronger GDP Effect", 
            rnott.high = "Weaker GDP Effect", rnott.high.slope = "Heterogeneous GDP Effect")

ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))

plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color = "grey", size = .1) +
    facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    scale_fill_gradientn(name = expression("R"[0]), colors=rev(colFunc(12)),
                       na.value = "white", breaks = breaks.trigger) +
    geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
    geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right", legend.key.size =  unit(0.5, "in")) +
    theme(strip.text.x = element_text(size = 14)) 

save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/rnott_uncertainty.pdf", plot = plot.rnott,  base_height = 8, base_aspect_ratio = 1.1) 


####################### Code for R0 Sensitivity Plots ######################################################

rnott_plot.m <- melt(data = county_sensitive, id.vars = c("id", "Geography"), 
                     measure.vars = c("decrease.b", "August", "increase.b"))




colnames(rnott_plot.m) <- c("id", "Geography", "variable", "rnott")
levels(rnott_plot.m$variable) <- c("b == 0.214","b == 0.634", "b == 0.8")

merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]



plot.sensitive.b <- ggplot(rnott.data, aes(x=long, y = lat)) +
  geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color =  "grey", size = .1) +
  facet_wrap(~variable, labeller = label_parsed) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_gradientn(name = expression("R"[0]), colors=rev(colFunc(12)),
                       na.value = "white", breaks = breaks.trigger) +
  geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black")+
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank(), 
    legend.position = c(0.37, 0.05), 
    legend.title.align = 0.5) + 
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10)) +
  theme(strip.text.x = element_text(size = 16)) 

   
save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/R0sensitive_b.pdf", plot = plot.sensitive.b,  base_height = 4, base_aspect_ratio = 2.4)  



#### Code for making Histograms of R0 distributions 
#breaks <- c(0,1,6)
#mycolors <- c("yellow", "red")
labels <- c(rnott.expected = "Expected GDP Effect", rnott.low = "Stronger GDP Effect", rnott.high = "Weaker GDP Effect", rnott.high.slope = "Hetergeneous GDP Effect")
colFunc <- colorRampPalette(c("red4","darkorange", "gold", "light yellow"))(12)
breaks.trigger = seq(from = 0, to = 5.5, by  = 0.5)

rnott_plot.m$b <- cut(rnott_plot.m$rnott, breaks.trigger)

hist.rnott <- ggplot(rnott_plot.m, aes(rnott, fill = rnott_plot.m$b)) +
  geom_histogram(binwidth = .1, color = "black") +
  facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
  geom_vline(xintercept = 1, colour="red", linetype = "longdash") +
  scale_x_continuous(name = expression("R"[0]), breaks = seq(0,5,.5)) + 
  scale_y_continuous(name = "Count", breaks = seq(0,100, 10)) +
  scale_fill_manual(breaks = levels(rnott_plot.m$b), values = rev(colFunc),
                                    name = expression("R"[0]), drop = FALSE) + 
  guides(fill = FALSE) +
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + 
  theme(strip.text.x = element_text(size = 16)) 

save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/hist_R0.pdf", plot = hist.rnott, base_height = 8, base_aspect_ratio = 1.5)  



###########  Code for Making Summary Figure for Sensitivity #############
county_sensitive <- read.csv("~/Documents/zika_alarm/csvs/county_sensitive.csv")
risk_counties <- apply(county_sensitive, MARGIN = 2, function(x) length(which(x > 1)))
risk_counties <- as.data.frame(risk_counties)
risk_counties$scenario <- c("Strong GDP", "Weak GDP", "Heterogeneous GDP", "High c", "Low c", "Low b",
                            "High b", "High alpha", "Low alpha", "July",
                            "September", "October", "November", "August")


barplot <- ggplot(risk_counties, aes(x = reorder(scenario, risk_counties), y = risk_counties)) + 
  geom_bar(stat = "identity") + coord_flip() + ylab(expression(paste("Number of Counties:","R"[0], "> 1"))) +
  xlab("Scenario") + 
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) + 
  geom_text(stat='identity', aes(label = risk_counties), hjust = - .25) 

save_plot(filename = "Documents/zika_alarm/ExploratoryFigures/sensitive_barplot.pdf", plot = barplot,  base_height = 8, base_aspect_ratio = 1.5)  



########### Testing the id of a county #########################################
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



