require(maptools)
require(ggplot2)
require(rgeos)
require(rgdal)
require(raster)
require(plyr)

data <- read.csv("~/Documents/zika_alarm/county.importation.r0.csv")
setwd('..'); setwd('TexasCountyShapeFiles')
texas.county <- readShapeSpatial('texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
setwd('../zika_code/')

county_plot <- read.csv(file = "~/Documents/zika_alarm/county_ids.csv")
colnames(county_plot) <- c("id", "metro_round", "intro")


carol.rates.round.4 <- round(data$Carol_Rates, digits = 4)
carol.rates.round.3 <- round(data$Carol_Rates, digits = 3)

length(unique(carol.rates.round.4))
length(unique(carol.rates.round.3))


###### Calculate Surveillance For R0s 
r_nots <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 1.6, 1.9) 
r_nots <- c(0.2, 0.5, 1.1)
disc_prob <- c( 0.068, 0.011)
intro_rate <- c(0.01, 0.1)


#Calculate the triggers  
r0.triggers <- calculate_all_triggers(dir_path, r_nots = r_nots, intro_rate = intro_rate, 
                                      disc_prob = disc_prob, threshold = 20, confidence = .8)

colnames(r0.triggers) <- c("run", "r0", "detect", "intro", "threshold", "confidence", "trigger")

trigger.maxdetect <- which(r0.triggers[,"trigger"] > 100)
r0.triggers[trigger.maxdetect, "trigger"] <- NA


########Combining with Texas Maps 
#Calculate Only triggers interested in
# Base Case 
county_plot$med.total.intro <- NA
county_plot$high.total.intro <- NA


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
  indices = which((triggers[i,"r0"] == county_plot$metro_round) & 
                   (triggers[i, "intro"] == county_plot$intro) 
  county_plot$med.total.intro[indices] = triggers[i,"triggers"]
}

county_plot$med.total.intro <- as.numeric(county_plot$med.total.intro)


### Different Detection 
desired_dect = 0.0110

triggers <- ddply(.data = r0.triggers, .variables = "r0", function (x) {
  row = which(x[,"detect"] == desired_dect)
  return(x[row,])
})

# Merges The Data With the County Data for Plotting 
## Matches by R0 and intro 
# Matching by two 
for (i in 1:nrow(triggers)) {
  indices = which((triggers[i,"r0"] == county_plot$metro_round) & 
                    (triggers[i, "intro"] == county_plot$intro) 
                  county_plot$high.total.intro[indices] = triggers[i,"triggers"]
}

county_plot$high.total.intro <- as.numeric(county_plot$high.total.intro)


#### SOME WAY TO PLOT THIS 
county.m <- melt(data = county_reduced, id.vars = "id", measure.vars = c("Prev.Cases.1", "Prev.Cases.2"))
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_ids, by = "id", all.x = TRUE)
merge.texas.county <- merge(texas.county.f, county.m, by = "id", all.x = TRUE)

final.plot <- merge.texas.county[order(merge.texas.county$id),]
head(final.plot)


#Actual plotting 
max(r0.triggers[,"trigger"], na.rm = TRUE)
legend_breaks <- round(seq(0, 100, 10))

plot.trial <- ggplot(final.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = value)) +
  facet_wrap(~variable)+
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL)

ggsave

plotworst <- ggplot()+geom_polygon(data = final.plot, aes_string(x="long", y = "lat", group = "group", fill = "Prev.Cases"),
                                   color = "black", size = .25) + coord_map() +
  scale_fill_continuous(name = "Detected \n Cases", low = "red", high = "yellow", 
                        na.value = "grey") + #pretty_breaks(n = length(legend_breaks))) +
  #theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
  theme(legend.key.size =  unit(0.3, "in")) 


plot.r0 <- ggplot()+geom_line(data = r0.triggers, aes(r0, trigger))+facet_grid(detect ~ intro)

plotworst
plot.badintro_gooddetect
plot.goodintro_baddetect
plotbaseline


plot_grid(plotbaseline, plot.badintro_gooddetect, 
          plot.goodintro_baddetect, plotworst, labels = c("A", "B", "C", "D"))




