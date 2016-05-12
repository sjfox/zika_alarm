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
library(ggrepel)

#### Read in file with county R0 and Texas shape file
county_plot <- read.csv("../csvs/county_master.csv")
texas.county <- readShapeSpatial('../TexasCountyShapeFiles/texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))

## Melt R0
county_plot.m <- melt(data = county_plot, id.vars = c("id", "Geography", "rnott.expected.round", "importation_probability"), 
                      measure.vars = c("importation.projected", "importation.worse.projected"))
colnames(county_plot.m) <- c("id", "geography", "rnott.expected", "importation_probability", "scenario", "import.rate")

### Calculate Surveillance Triggers For R0s and Import Rates 
rnots <- sort(unique(county_plot.m$rnott.expected))
import.rates <- sort(unique(county_plot.m$import.rate))

triggers <- get_trigger_data(rnots, intro = import.rates, disc = 0.0224, confidence = .5, num_necessary = 10)


#### Match R0/Import Rate with Trigger
county_plot.m <- merge(x = county_plot.m, y = triggers[,c("r_not", "intro_rate", "prev_trigger")], 
                       by.x=c("rnott.expected", "import.rate"), by.y=c("r_not", "intro_rate"), all.x=TRUE, sort=FALSE)


#### Fortifying Data to ShapeFile for ggplot 
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_plot.m, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]
final.plot$import.type <- factor(final.plot$scenario, 
                                 levels = c("importation.projected", "importation.worse.projected"))

## Get metropolitan city data
map_data <- read.csv("../csvs/metro_geo.csv")
map_data <- map_data[1:10,]

#############################################
## R0 and importation figure first
############################################
#### Importation Probability Map 
merge.texas.county.import <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
import.plot <- merge.texas.county[order(merge.texas.county$id),]
import.plot$importation_probability.log <- log(import.plot$importation_probability)


plot.importation.log <- ggplot(import.plot, aes(x = long, y = lat)) +
  geom_polygon(data = import.plot, aes(group = group, fill = importation_probability.log), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Import Probability", low ="white", high = "blue", breaks = c(-6, -4, -2), labels = round(exp(c(-6, -4, -2)), digits = 3),
                        na.value = "white") +
  geom_point(data = map_data, aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data, aes(x=lon, y = lat, label = Name), size = 5,force=0.75, segment.color = "black")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10))
# print(plot.importation.log)

r0.plot <- ggplot(final.plot[which(final.plot$scenario=="importation.projected"),], aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = rnott.expected), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = expression("R"[0]), low = "white", high = "darkgreen", na.value = "white") +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10))

import_r0_plot <- plot_grid(plot.importation.log, r0.plot, nrow = 1, labels = "AUTO", label_size = 20)
save_plot(filename = "../ExploratoryFigures/fig3_import_r0.pdf", plot = import_r0_plot, base_height = 4, base_aspect_ratio = 2.2)

###########################################
## TEXAS RISK ASSESSMENT
###########################################
ind <- which(map_data$Name %in% c("Austin", "Houston"))
###### Trigger Maps, Faceted by Scenario 
breaks.trigger = seq(from = 0, to = max(final.plot$prev_trigger,na.rm=TRUE), by = 40)
plot.trial <- ggplot(final.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = prev_trigger), color = "grey", size = .1) +
  facet_wrap(~scenario, nrow = 1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  scale_fill_continuous(name = "Surveillance Trigger", low = "red", high = "white", 
                        na.value = "white", breaks = breaks.trigger) +
  geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.title.align = 0.5,
                                   legend.position=c(0.55,0.15)) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 12))
## Plot just trigger maps
# save_plot(filename = "../ExploratoryFigures/figure3_triggers.pdf", plot = plot.trial, base_height = 4, base_aspect_ratio = 2.2)
# plot(plot.trial)
# trigger.legend <- get_legend(plot.trial)
# plot.trial <- plot.trial+theme(legend.position="none")



# save_plot(filename = "../ExploratoryFigures/figure3_importationlog.pdf", plot = plot.importation.log, base_height = 4, base_aspect_ratio = 1.2)            

# fig3_all <- ggdraw() + draw_plot(plot.trial, x = 0, y=0, width=1, height=0.5)+ 
#   draw_plot(plot = plot.importation.log, x = 0, y=0.5, width=.5, height=0.5)+
#   draw_plot(plot = r0.plot, x = 0.5, y=0.48, width=0.5, height=0.5)+
#   draw_plot(plot = trigger.legend, x = 0.42, y=0.05, width=0.3, height=0.1)+
#   draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 20)
# 
# save_plot(filename = "../ExploratoryFigures/figure3_combined.pdf", plot = fig3_all, base_height = 8, base_aspect_ratio = 1.1)

# #Geting the lat long of the metro areas
# 
# summary.worse.projected <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.worse.projected"])
# summary.projected <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.projected"])
# summary.current <- summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.current"])

load("../data/rand_trigger.Rdata")
rand_triggers$disc <- paste0(calculate.discover(rand_triggers$disc), "%")

rand_trigger_plot <- ggplot(rand_triggers, aes(confidence, value, color=variable, linetype=as.factor(disc))) + geom_line(size=1)+
  coord_cartesian(xlim = c(0,0.75), expand=FALSE)+
  scale_color_manual(values=c("Grey", "Black"))+
  theme_cowplot()%+replace% theme(legend.position=c(0.3, 0.6),
                                  legend.box.just = "left")+
  labs(x = "Risk Tolerance", 
       y = "Trigger (Reported Cases)", 
       linetype = "Detection\nProbability",
       color="Trigger Type")
print(rand_trigger_plot)

fig4 <- ggdraw() + draw_plot(plot.trial, x = 0, y=0, width=0.66, height=1)+
  draw_plot(plot = rand_trigger_plot, x = 0.66, y=0, width=0.33, height=1)+
  draw_plot_label(c("A", "B", "C"), c(0, 0.33, 0.66), c(1, 1, 1), size = 20)
save_plot(filename = "../ExploratoryFigures/fig4_texas_risk.pdf", plot = fig4, base_height = 4, base_aspect_ratio = 3)

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
# g <- ggplotGrob(plot.trial)
# gl <- g$layout
# ## Don't need to mess with axis, because not applicable for maps, but code for that from stackoverflow is below
# gl[2, 1:4] <- c(4,7,4,7)
# gl[3, 1:4] <- c(8,4,8,4)
# gl[4, 1:4] <- c(8,7,8,7)
# g$layout <- gl


##########################################
#######################
# RANDOM SELECT VALUES
# imports <- read.csv("../csvs/county_master.csv")
# dir_path <-"~/projects/zika_alarm/data/sep_intros/"
# projected_combos <- as.character(unique(interaction(imports$rnott.expected.round, imports$importation.projected,sep = "_")))
# 
# values <- unlist(strsplit(projected_combos, "_"))
# ## R0s are the first item in each combo
# r_nots <- as.numeric(values[seq(1,length(values), by=2)])
# intro_rates <- as.numeric(values[seq(2,length(values), by=2)])
# indices <- 1:length(r_nots)
# 
# num_samps <- 10000
# samps <- sample(indices, size = num_samps, replace = T)
# samps <- table(samps)
# rand_trials <- vector("list", length = num_samps)
# list_ind <- 1
# 
# for(i in 1:length(samps)){
#   index <- as.numeric(names(samps[i]))
#   num_trials <- samps[i]
#   #  print(i)
#   #  print(get_vec_of_files(dir_path, r_nots = r_nots[index], disc_probs = 0.011, intro_rates = intro_rates[index]))
#   load(get_vec_of_files(dir_path, r_nots = r_nots[index], disc_probs = 0.0224, intro_rates = intro_rates[index]))
#   rand_trials[list_ind:(list_ind+num_trials-1)] <- trials[sample(1:10000, size=num_trials)]
#   list_ind <- list_ind+num_trials
# }
# conf <- seq(0.05, 0.95, by=0.05)
# epi_triggers <- c()
# prev_triggers <- c()
# 
# for(ci in conf){
#   temp <- get_epidemic_trigger(trials = rand_trials, threshold = 50, confidence = ci, max_detect = 300, num_necessary = 10)  
#   temp2 <- get_surveillance_trigger(trials=rand_trials, threshold = 20, confidence=ci, max_detect = 300, num_necessary=10)
#   epi_triggers <- c(epi_triggers, temp)
#   prev_triggers <- c(prev_triggers, temp2)
# }
# df <- data.frame(confidence=conf, epi_trigger = epi_triggers, prev_trigger = prev_triggers[length(prev_triggers):1])
# df <- melt(df, id.vars=c("confidence"))
# df$variable <- ifelse(df$variable=="epi_trigger", "Nowcasting", "Forecasting")
# rand_trigger_plot <- ggplot(df, aes(confidence, value, color=variable)) + geom_line(size=1)+
#   coord_cartesian(xlim = c(0,1), expand=FALSE)+
#   scale_color_manual(values=c("Grey", "Black"))+
#   theme_cowplot()%+replace% theme(legend.position=c(0.3, 0.8))+
#   labs(x = "Risk Tolerance", 
#        y = "Trigger (Reported Cases)", 
#        color = "")
# 
# print(rand_trigger_plot)
# df$disc <- 0.0224
# df10disc$disc <- 0.011
# rand_triggers <- rbind(df, df10disc)
# save(list = c("rand_triggers"), file = "../data/rand_trigger.Rdata")
