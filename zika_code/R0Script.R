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


county_master <- read.csv("../csvs/county_master.csv")

#########################################################
# Step 1: Getting mosquito abundance data for each county 
#########################################################

# Importing the aegypti occurrence projections from Kraemer 
world.aegypti <- raster("~/Downloads/aegypti.tif.ovr")  ## TODO Will need to figure out where this is and put on repo

## setting the projection to match Texas shape files  
bb <- extent(-180, 180, -90, 90)
extent(world.aegypti) <- bb
world.aegypti <- setExtent(world.aegypti, bb, keepres=TRUE)
projection(world.aegypti) <- CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +no_defs")

# Cutting the world layer to Texas 
world.cut <- crop(world.aegypti, extent(texas.county))
texas.cut <- mask(world.cut, texas.county)

# Calculating Mosquito Occurrence for each county as the 
# Average occurrence for each cell in a county 
# TODO will need to figure out where mean.county function is 
mosquito.occurrence <- mean.county(raster = texas.cut, shapefile = texas.county, name = "occurrence") 

# Converting to proxy 
mosquito.abundance <- -log(1-mosquito.occurrence)
mosquito.abundance <- as.numeric(mosquito.abundance)

# Save to County spreadsheet 
#county_master$mosquito.abundance <- mosquito.abundance

#######################################################
# Step 2: Adjust mosquito occurrence for economic effect
# Based on Perkins 2016
########################################################

# Calculating Ranges of GDP to set Ranges of Multiplication Factor 
GDP <- log(county_master$GDP)-10
logMF.expected <- -1.79-(.07/.5)*GDP # Log Linear economic relationship 
MF.expected <- exp(logMF.expected) # Multiplication Factor 
mosquito.expected <- MF.expected * mosquito.abundance # Adjusted mosquito abundance 


#######################################################
# Step 3: Calculate R0 
#######################################################

# R0 equation is set to use the august eip from the county master 
# Baseline parameters 

mos.hum.transmission = .85
c_over_r = 9 *.77
alpha = .63
mos.mortality = 1/14
eip = county_master$aug.eip


rnott.expected <- c_over_r *(mosquito.expected * mos.hum.transmission*alpha^2*
                                     exp(-mos.mortality*eip))/mos.mortality

county_master$rnott.expected <- NA; county_master$rnott.expected <- rnott.expected
write.csv(county_master, file= "Documents/zika_alarm/csvs/county_master.csv")


#####################################
# Plotting 
####################################

fig_path <- "../ExploratoryFigures/"

rnott_plot.m <- melt(data = county_master, id.vars = c("id", "Geography"),
                     measure.vars = c("rnott.expected", "rnott.low", "rnott.high", "rnott.high.slope"))
colnames(rnott_plot.m) <- c("id", "Geography", "rnott.scenario", "rnott")

texas.county <- readShapeSpatial('~/Documents/projects/zika_alarm/TexasCountyShapeFiles/texas.county.shp',
                                 proj4string = CRS("+proj=longlat +datum=WGS84"))
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]


colorends <- c("white", "darkseagreen", "yellow", "red")
gradientends <- c(0,1,1.01,max(rnott.data$rnott))

labels <- c(rnott.expected = "Expected GDP Effect", rnott.low = "Stronger GDP Effect", 
            rnott.high = "Weaker GDP Effect", rnott.high.slope = "Heterogeneous GDP Effect")

ind <- which(map_data$Name %in% c("Austin", "Houston", "Dallas", "San Antonio"))

plot.rnott <- ggplot(rnott.data, aes(x=long, y = lat)) +
    geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color = "grey", size = .1) +
    facet_wrap(~rnott.scenario, nrow = 2, dir = "h", labeller = labeller(rnott.scenario = labels)) +
    scale_x_continuous("", breaks=NULL) + 
    scale_y_continuous("", breaks=NULL) + 
    #scale_fill_gradientn(name = expression("R"[0]), colors=rev(colFunc(12)),
    #                   na.value = "white", breaks = breaks.trigger) +
    scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
    geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
    geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black") +
    theme_cowplot() + theme(strip.background = element_blank()) + 
    theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
    theme(legend.position = "right", legend.key.size =  unit(0.5, "in")) +
    theme(strip.text.x = element_text(size = 14)) 

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/rnott_uncertainty.pdf", plot = plot.rnott,  base_height = 8, base_aspect_ratio = 1.1) 


####################### Code for R0 Sensitivity R0 Calculations ######################################################

# Sensitivity Analysis 

##### Additional Economic Relationship
logMF.high <- -.9-(.07/.5)*GDP #High
logMF.low <- -2.6-(.07/.5)*GDP #Low 
logMF.high.slope <- -1.35-(.9/.5)*GDP



rnott_plot.m <- melt(data = county_sensitive, id.vars = c("id", "Geography"), 
                     measure.vars = c("decrease.cr", "August", "increase.cr"))

colnames(rnott_plot.m) <- c("id", "Geography", "variable", "rnott")
levels(rnott_plot.m$variable) <- c("c == 0.6", "c == 0.77", "c == 0.95")
#levels(rnott_plot.m$variable) <- c("July", "August", "September", "October", "November")

merge.texas.county <- merge(texas.county.f, rnott_plot.m, by = "id", all.x = TRUE)
rnott.data <- merge.texas.county[order(merge.texas.county$id),]

colorends <- c("white", "darkseagreen", "yellow", "red")
gradientends <- c(0,1,1.01,max(rnott.data$rnott))
head(rnott.data)


plot.sensitive.c <- ggplot(rnott.data, aes(x=long, y = lat)) +
  geom_polygon(data = rnott.data, aes(group = group, fill = rnott), color =  "grey", size = .1) +
  facet_wrap(~variable, labeller = label_parsed) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
 # scale_fill_gradientn(name = expression("R"[0]), colors=rev(colFunc(12)),
  #                     na.value = "white", breaks = breaks.trigger) +
  scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
  geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black")+
  theme_cowplot() + theme(strip.background = element_blank()) + 
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank(), 
    legend.position = c(0.37, 0.05), 
    legend.title.align = 0.5) + 
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10)) +
  theme(strip.text.x = element_text(size = 16)) 


save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/R0sensitive_c.pdf", plot = plot.sensitive.c,  base_height = 4, base_aspect_ratio = 2.4)  


r0.plot <- ggplot(final.plot[which(final.plot$scenario=="importation.projected"),], aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = rnott.expected), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
  # scale_fill_gradient2(name = expression("R"[0]), low = "lightgreen", mid =  "yellow", high = "red", midpoint = 1, na.value = "white") +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10))


 
save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/R0sensitive_b.pdf", plot = plot.sensitive.b,  base_height = 4, base_aspect_ratio = 2.4)  



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
county_sensitive <- read.csv("~/Documents/projects/zika_alarm/csvs/county_sensitive.csv")
county_sensitive2 <- county_sensitive
county_sensitive <- county_sensitive[,-1]
county_sensitive <- round(county_sensitive[,2:15], 1)

risk_counties <- apply(county_sensitive, MARGIN = 2, function(x) length(which(x > 1)))
risk_counties <- as.data.frame(risk_counties)
risk_counties$scenario <- c("Strong GDP", "Weak GDP", "Heterogeneous GDP", "High c", "Low c", "Low b",
                            "High b", "High alpha", "Low alpha", "July",
                            "September", "October", "November", "August (baseline)")

library(ggplot2)
library(cowplot)

barplot <- ggplot(risk_counties, aes(x = reorder(scenario, risk_counties), y = risk_counties)) + 
  geom_bar(stat = "identity") + coord_flip() + ylab(expression(paste("Number of Counties:","R"[0], "> 1"))) +
  theme_cowplot() + theme(strip.background = element_blank()) + xlab("") +
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14),
        axis.title.x = element_text(size = 16)) + theme(legend.title=element_blank()) +
  theme(legend.position="none") +
  geom_text(stat='identity', aes(label = risk_counties, size = 14), hjust = - .25) 

barplot

save_plot(filename = "~/Documents/projects/zika_alarm/ExploratoryFigures/sensitive_barplot.pdf", plot = barplot,  base_height = 8, base_aspect_ratio = 1.5)  






colorends <- c("white", "darkseagreen", "yellow", "red")
gradientends <- c(0,1,1.01,max(final.plot$rnott.expected))


r0.plot <- ggplot(final.plot[which(final.plot$scenario=="importation.projected"),], aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = rnott.expected), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_gradientn(name = expression("R"[0]), colours = colorends, values = rescale(gradientends)) +
  # scale_fill_gradient2(name = expression("R"[0]), low = "lightgreen", mid =  "yellow", high = "red", midpoint = 1, na.value = "white") +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10))

# print(r0.plot)

import_r0_plot <- plot_grid(plot.importation.log, r0.plot, nrow = 1, labels = "AUTO", label_size = 20)
save_plot(filename = "../ExploratoryFigures/fig3_import_r0.pdf", plot = import_r0_plot, base_height = 4, base_aspect_ratio = 2.2)



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



