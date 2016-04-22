rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


sapply(c('branch.functions.R','plot.functions.R', 'analyze_saved_sims.R'), source)
library(plyr)
library(ggplot2)
library(maptools)
library(ggplot2)
library(rgeos)
library(rgdal)
library(raster)
library(plyr)


county_plot <- read.csv("~/Documents/zika_alarm/county_plot.csv")
setwd('..'); setwd('TexasCountyShapeFiles')
texas.county <- readShapeSpatial('texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
setwd('../zika_code/')


#### How many times do you want to compare to current levels 
#### What is your current Rate
#### Calculate Importation Rate based on probability 
testing_rate = 1
testing_rate_high = 2
current_rate = 30/112
imports <- round(county_plot$importation_probability*testing_rate*current_rate, digits = 3)
county_plot <- cbind(county_plot, imports)
high.imports <- round(county_plot$importation_probability*testing_rate_high*current_rate, digits = 3) 
county_plot <- cbind(county_plot, high.imports)
county_plot$constant.import <- .10


county_plot.m <- melt(data = county_plot, id.vars = c("id", "Geography", "metro_round"), measure.vars = c("imports", "high.imports", "constant.import"))
colnames(county_plot.m) <- c("id", "geography", "metro_round", "import.type", "import_rate")

###### Calculate Surveillance Triggers For R0s 
r_nots <- c(0.8, 0.9)
intro_rate <- 0.1
dirPaths <- get_vec_of_files(dir_path, r_nots, intro_rate, disc_prob)

r_nots <- c(unique(sort(county_plot.m$metro_round)))
intro_rate <- unique(sort(county_plot.m$import_rate))
disc_prob <- c( 0.068, 0.011)


dir_path <- "~/Documents/zika_alarm/data/introductions/"
r0.triggers <- calculate_all_triggers(dir_path, r_nots = r_nots, intro_rate = intro_rate, 
                                      disc_prob = disc_prob, threshold = 20, confidence = .8)


# Cleaning UP r0.triggers for easier analysis 
colnames(r0.triggers) <- c("run", "r0", "detect", "intro", "threshold", "confidence", "trigger")
trigger.maxdetect <- which(r0.triggers[,"trigger"] > 100)
r0.triggers[trigger.maxdetect, "trigger"] <- NA


########Combining with Texas Maps 
#Calculate Only triggers interested in
# Add to data frame column for holding the triggers 
county_plot.m$trigger <- NA


#BaseLine First Plot 
desired_dect = 0.0110

triggers <- ddply(.data = r0.triggers, .variables = "r0", function (x) {
  row = which(x[,"detect"] == desired_dect)
  return(x[row,])
})


# Merges The Data With the County Data for Plotting 
## Matches by R0 and intro 
# Matching by two 
for (i in 1:nrow(triggers)) {
  indices = which((county_plot.m$metro_round == triggers[i,"r0"]) & 
                   (county_plot.m$import_rate== triggers[i, "intro"])) 
  county_plot.m$trigger[indices] = triggers[i,"trigger"]
}

county_plot.m$trigger <- as.numeric(county_plot.m$trigger)



#### TRYING TO PLOT 
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_plot.m, by = "id", all.x = TRUE)


final.plot <- merge.texas.county[order(merge.texas.county$id),]
final.plot$import.type <- factor(final.plot$import.type, 
                                 levels = c("constant.import", "imports", "high.imports"))


#Actual plotting 
library(cowplot)

breaks = seq(from = min(final.plot$trigger, na.rm = TRUE),to = max(final.plot$trigger,na.rm=TRUE), by = 5)

plot.trial <- ggplot(final.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = trigger), color = "black", size = .25) +
  facet_wrap(~import.type, nrow = 2, dir = "h") +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Detected \n Cases", low = "red", high = "yellow", 
                        na.value = "grey", breaks) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
  theme(legend.key.size =  unit(0.3, "in")) 

ggsave(plot.trial, filename = "low.vs.high.heteroimport.pdf", height = 7, width = 8, units = "in")




