rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


sapply(c('branch.functions.R','plot.functions.R', 'analyze_saved_sims.R'), source)
library(plyr)
library(maptools)
library(rgeos)
library(rgdal)
library(raster)
library(plyr)
library(cowplot)


#### Read in file with county R0 and Texas shape file
county_plot <- read.csv("../csvs/county_master.csv")
texas.county <- readShapeSpatial('../TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))

## Melt R0
county_plot.m <- melt(data = county_plot, id.vars = c("id", "Geography", "rnott.expected.round", "importation_probability"), 
                      measure.vars = c("importation.current", "importation.projected", "importation.worse.projected"))
colnames(county_plot.m) <- c("id", "geography", "rnott.expected", "importation_probability", "scenario", "import.rate")

### Calculate Surveillance Triggers For R0s and Import Rates 
rnots <- sort(unique(county_plot.m$rnott.expected))
import.rates <- sort(unique(county_plot.m$import.rate))

triggers <- get_trigger_data(rnots, intro = import.rates,
                             disc = 0.0224, threshold = 20, confidence = .8, num_necessary = 10)


#### Match R0/Import Rate with Trigger
county_plot.m <- merge(x = county_plot.m, y = triggers[,c("r_not", "intro_rate", "prev_trigger")], 
                       by.x=c("rnott.expected", "import.rate"), by.y=c("r_not", "intro_rate"), all.x=TRUE, sort=FALSE)


#### Fortifying Data to ShapeFile for ggplot 
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_plot.m, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]
final.plot$import.type <- factor(final.plot$scenario, 
                                 levels = c("importation.current", "importation.projected", 
                                            "importation.worse.projected"))




###### Trigger Maps, Faceted by Scenario 
breaks.trigger = seq(from = 0, to = max(final.plot$prev_trigger,na.rm=TRUE), by = 20)
plot.trial <- ggplot(final.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = prev_trigger), color = "black", size = .25) +
  facet_wrap(~scenario, nrow = 2, dir = "h") +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Surveillance \nTrigger \n", low = "red", high = "yellow", 
                        na.value = "grey", breaks = breaks.trigger) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.key.size =  unit(0.5, "in")) 
## Plot just trigger maps
save_plot(filename = "../ExploratoryFigures/figure3_triggers.pdf", plot = plot.trial, base_height = 8, base_aspect_ratio = 1.3)


#### Importation Probability Map 
#merge.texas.county.import <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
#import.plot <- merge.texas.county[order(merge.texas.county$id),]

#max.import = max(import.plot$importation_probability,na.rm=TRUE)
#plot.importation <- ggplot(import.plot, aes(x = long, y = lat)) +
#  geom_polygon(data = import.plot, aes(group = group, fill = importation_probability), color = "black", size = .25) +
#  scale_x_continuous("", breaks=NULL) + 
#  scale_y_continuous("", breaks=NULL) + 
#  scale_fill_continuous(name = "Import \nProbability", low = "light yellow", high = "red", 
#                        na.value = "white", breaks= c(0,0.03, 0.06,0.09, 0.12, 0.25)) +
#  theme_cowplot() %+replace% theme(strip.background=element_blank(),
#                                   strip.text.x = element_blank(),
#                                   legend.key.width =  unit(0.5, "in"),
#                                   legend.position = c(0.15, 0.8)) 
                             


merge.texas.county.import <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
import.plot <- merge.texas.county[order(merge.texas.county$id),]
import.plot$importation_probability.log <- log(import.plot$importation_probability)

map_data <- read.csv("../csvs/metro_geo.csv")
map_data2 <- transform(map_data, lat2 = lat + .3)

plot.importation.log <- ggplot(import.plot, aes(x = long, y = lat)) +
  geom_polygon(data = import.plot, aes(group = group, fill = importation_probability.log), color = "grey", size = .05) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Import \nProbability", low ="light blue", high = "blue", 
                        na.value = "white") +
  geom_point(data = map_data2, aes(x = lon, y = lat), color = "black", show.legend = FALSE) +
  geom_text(data = map_data2, aes(x=lon, y = lat2, label = Name), size = 3)+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.key.width =  unit(0.3, "in"),
                                   legend.position = c(0.15, 0.8)) 

save_plot(filename = "../ExploratoryFigures/figure3_importationlog.pdf", plot = plot.importation.log, base_height = 4, base_aspect_ratio = 1.2)            

g <- ggplotGrob(plot.trial)
gl <- g$layout
## Don't need to mess with axis, because not applicable for maps, but code for that from stackoverflow is below
gl[2, 1:4] <- c(4,7,4,7)
gl[3, 1:4] <- c(8,4,8,4)
gl[4, 1:4] <- c(8,7,8,7)
g$layout <- gl
fig3_all <- ggdraw() + draw_plot(g)+ 
  draw_plot(plot = plot.importation.log, x = 0, y=0.48, width=0.46, height=0.52)+
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.45, 0, 0.45), c(1, 1, 0.5, 0.5), size = 20)

save_plot(filename = "../ExploratoryFigures/figure3_combined.pdf", plot = fig3_all, base_height = 8, base_aspect_ratio = 1.3)

#Geting the lat long of the metro areas

summary.worse.projected <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.worse.projected"])
summary.projected <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.projected"])
summary.current <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.current"])





################### SO code, helps move axis and titles with plots
# head(diamonds)
# plot <- ggplot(diamonds, aes(carat, price)) + facet_wrap(~cut) + geom_point() 
# g <- ggplotGrob(plot)
# 
# gl <- g$layout
# g$layout
# idcol <- gl$r == (ncol(g)-2)
# g$layout[idcol & gl$b < 5, c("t", "b")] <- gl[idcol & gl$b < 5, c("t", "b")] + 4
# g$layout
# plot(g)
