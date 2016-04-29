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
county_plot <- read.csv("../csvs/county_plot.csv")
texas.county <- readShapeSpatial('../TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))

## Melt R0
county_plot.m <- melt(data = county_plot, id.vars = c("id", "Geography", "metro_round", "importation_probability"), 
                      measure.vars = c("importation.current", "importation.projected", "importation.worse.projected"))
colnames(county_plot.m) <- c("id", "geography", "metro_round", "importation_probability", "scenario", "import.rate")

### Calculate Surveillance Triggers For R0s and Import Rates 
rnots <- sort(unique(county_plot.m$metro_round))
import.rates <- sort(unique(county_plot.m$import.rate))

triggers <- get_trigger_data(rnots, intro = import.rates,
                             disc = 0.0224, threshold = 20, confidence = .8)

#### Match R0/Import Rate with Trigger
county_plot.m <- merge(x = county_plot.m, y = triggers[,c("r_not", "intro_rate", "prev_trigger")], 
                       by.x=c("metro_round", "import.rate"), by.y=c("r_not", "intro_rate"), all.x=TRUE, sort=FALSE)


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
merge.texas.county.import <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
import.plot <- merge.texas.county[order(merge.texas.county$id),]
head(import.plot)


plot.importation <- ggplot(import.plot, aes(x = long, y = lat)) +
  geom_polygon(data = import.plot, aes(group = group, fill = importation_probability), color = "black", size = .25) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Import \nProbability", low = "light yellow", high = "red", 
                        na.value = "grey", breaks= c(0,0.03, 0.06,0.09, 0.12)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.key.width =  unit(0.5, "in"),
                                   legend.position = c(0.15, 0.8)) 

## Plot just importations
save_plot(filename = "../ExploratoryFigures/figure3_importation.pdf", plot = plot.importation, base_height = 4, base_aspect_ratio = 1.3)


# fig3_all <- ggdraw() +
#   draw_plot(plot.trial, x =  0,y =  0, width =  1,height =  1) +
#   draw_plot_label(c("A", "B", "C", "D"), c(0, 0.45, 0, 0.45), c(1, 1, 0.5, 0.5), size = 20)
# 
# save_plot(filename = "../ExploratoryFigures/figure3_combined.pdf", plot = fig3_all, base_height = 8, base_aspect_ratio = 1.3)


