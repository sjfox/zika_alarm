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

triggers <- get_trigger_data(rnot = rnots, intro = import.rates, disc = 0.0224, confidence = .3, num_necessary = 100)
# triggers$prev_trigger <- ifelse(triggers$prev_trigger>200, 200, triggers$prev_trigger)

#### Match R0/Import Rate with Trigger
county_plot.m <- merge(x = county_plot.m, y = triggers[,c("r_not", "intro_rate", "epi_trigger")], 
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
  scale_fill_continuous(name = "Import Probability", low ="white", high = "blue", breaks = log(c(0.002, 0.02, 0.2)), labels = c(0.002,0.02,0.2),
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

colFunc <- colorRampPalette(c("red4","darkorange", "gold", "yellow"))
breaks.trigger = seq(from = 0, to = max(final.plot$epi_trigger,na.rm=TRUE), by = 1)
plot.trial <- ggplot(final.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = epi_trigger), color = "grey", size = .1) +
  facet_wrap(~scenario, nrow = 1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  # scale_fill_continuous(name = "Trigger (Nowcasting)", low = "red", high = "white",
  scale_fill_gradientn(name = "Trigger (Reported Cases)", colors=colFunc(10),
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
save_plot(filename = "../ExploratoryFigures/trigger_maps_test.pdf", plot = plot.trial, base_height = 4, base_aspect_ratio = 2.2)
# plot(plot.trial)
# trigger.legend <- get_legend(plot.trial)
# plot.trial <- plot.trial+theme(legend.position="none")



temp <- get_trigger_data(c(0.8, 0.85, seq(0.9, 1.2, by=0.05)), intro=c(0.0, 0.01, 0.05, 0.1, 0.3), disc = c(0.0224), confidence = 0.3, num_necessary = 100)

ggplot(temp, aes(intro_rate, epi_trigger, color=r_not, group=r_not)) + geom_line() +scale_y_log10()


# save_plot(filename = "../ExploratoryFigures/figure3_importationlog.pdf", plot = plot.importation.log, base_height = 4, base_aspect_ratio = 1.2)            

# fig3_all <- ggdraw() + draw_plot(plot.trial, x = 0, y=0, width=1, height=0.5)+ 
#   draw_plot(plot = plot.importation.log, x = 0, y=0.5, width=.5, height=0.5)+
#   draw_plot(plot = r0.plot, x = 0.5, y=0.48, width=0.5, height=0.5)+
#   draw_plot(plot = trigger.legend, x = 0.42, y=0.05, width=0.3, height=0.1)+
#   draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 20)
# 
# save_plot(filename = "../ExploratoryFigures/figure3_combined.pdf", plot = fig3_all, base_height = 8, base_aspect_ratio = 1.1)

load("../data/rand_trigger.Rdata")
rand_triggers1 <- rand_triggers[which(rand_triggers$disc_prob==0.0224),]
rand_triggers1$disc_prob<-NULL
rand_triggers1$r_not <- "All"
rand_triggers1$variable <- ifelse(rand_triggers1$variable=="Forecasting", "Epidemic", "Prevalence")

# test <- get_trigger_data(1.05, 0.1, 0.0224, confidence=seq(0.05, .95, by=0.05), num_necessary=10)
# test <- melt(test, id.vars=c("confidence"), measure.vars = c("epi_trigger", "prev_trigger"))
# test$confidence <- ifelse(test$variable == "epi_trigger", test$confidence, 1-test$confidence)
# test$variable <- ifelse(test$variable=="epi_trigger", "Epidemic", "Prevalence")
# test$r_not <- "1.05"
# rand_trigger_dat <- rbind(test, rand_triggers)
load("../data/rand_trigger_extra_high_risk.Rdata")
rand_triggers2 <- rand_triggers[which(rand_triggers$disc_prob==0.0224),]
rand_triggers2$disc_prob<-NULL
rand_triggers2$r_not <- "High Risk"
rand_triggers2$variable <- ifelse(rand_triggers2$variable=="Forecasting", "Epidemic", "Prevalence")

rand_trigger_dat <- rbind(rand_triggers1, rand_triggers2)
rand_trigger_dat$variable <- factor(rand_trigger_dat$variable, levels = c("Prevalence", "Epidemic"))
# rand_trigger_dat <- rand_trigger_dat[which(rand_trigger_dat$variable=="Epidemic"), ]
rand_trigger_plot <- ggplot(rand_trigger_dat, aes(confidence, value, color=variable, linetype=r_not)) + geom_line(size=1)+
  coord_cartesian(xlim = c(0,1.025), ylim=c(0,50),expand=FALSE)+
  scale_color_manual(values=c("Grey", "Black"))+
  scale_linetype_manual(values=c(10,1))+
  theme_cowplot()%+replace% theme(legend.position=c(0.2, 0.75),
                                  legend.box.just = "left")+
  labs(x = "Risk Tolerance", 
       y = "Trigger (Reported Cases)", 
       linetype = expression("R"[0]),
       color="Trigger Type")
print(rand_trigger_plot)

fig4 <- ggdraw() + draw_plot(plot.trial, x = 0, y=0, width=0.66, height=1)+
  draw_plot(plot = rand_trigger_plot, x = 0.66, y=0, width=0.33, height=1)+
  draw_plot_label(c("A", "B", "C"), c(0, 0.33, 0.66), c(1, 1, 1), size = 20)
save_plot(filename = "../ExploratoryFigures/fig4_texas_risk.pdf", plot = fig4, base_height = 4, base_aspect_ratio = 3)




summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.projected"])
summary(county_plot.m$prev_trigger[county_plot.m$scenario == "importation.worse.projected"])
triggers_both_scenarios <- dcast(county_plot.m, geography ~ scenario, value.var = "prev_trigger")

sum(!is.na(triggers_both_scenarios$importation.projected - triggers_both_scenarios$importation.worse.projected))
mean(triggers_both_scenarios$importation.projected - triggers_both_scenarios$importation.worse.projected, na.rm=T)



#####################################
## FIGURE FOR ALISON
#####################################
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

triggers <- get_trigger_data(rnot = rnots, intro = import.rates, disc = 0.0224, confidence = .7, num_necessary = 100)
# triggers$prev_trigger <- ifelse(triggers$prev_trigger>200, 200, triggers$prev_trigger)

#### Match R0/Import Rate with Trigger
county_plot.m <- merge(x = county_plot.m, y = triggers[,c("r_not", "intro_rate", "epi_trigger")], 
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
  scale_fill_continuous(name = "Import Probability", low ="white", high = "blue", breaks = log(c(0.002, 0.02, 0.2)), labels = c(0.002,0.02,0.2),
                        na.value = "white") +
  geom_point(data = map_data, aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data, aes(x=lon, y = lat, label = Name), size = 5,force=1.1, segment.color = "black")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 8))


texas.plot <- final.plot[which(final.plot$scenario=="importation.worse.projected"), ]
ind <- which(map_data$Name %in% c("Austin", "Houston"))
###### Trigger Maps, Faceted by Scenario 
colFunc <- colorRampPalette(c("red4","darkorange", "gold", "yellow"))
breaks.trigger = seq(from = 0, to = max(final.plot$epi_trigger,na.rm=TRUE), by = 10)
texas.plot.trial <- ggplot(texas.plot, aes(x=long, y = lat)) +
  geom_polygon(data = final.plot, aes(group = group, fill = epi_trigger), color = "grey", size = .1) +
  scale_x_continuous("", breaks=NULL) + 
  scale_y_continuous("", breaks=NULL) + 
  # scale_fill_continuous(name = "Trigger (Nowcasting)", low = "red", high = "white",
  scale_fill_gradientn(name = "Trigger (Reported Cases)", colors=colFunc(10),
                       na.value = "white", breaks = breaks.trigger) +
  geom_point(data = map_data[ind,], aes(x = lon, y = lat), color = "black", size=1, show.legend = FALSE) +
  geom_text_repel(data = map_data[ind,], aes(x=lon, y = lat, label = Name), force=100, nudge_x = ifelse(map_data[ind,]$Name=="Austin",-1,1),
                  size = 5, point.padding = unit(0.25, "lines"), segment.color = "black")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),
                                   strip.text.x = element_blank(),
                                   legend.position = c(0.2, 0.17), 
                                   legend.title.align = 0.5) +
  guides(fill = guide_colorbar(label.position = "bottom", title.position="top",direction = "horizontal", barwidth = 10))


###############################
load("../data/rand_trigger.Rdata")
rand_triggers1 <- rand_triggers[which(rand_triggers$disc_prob==0.0224),]
rand_triggers1$disc_prob<-NULL
rand_triggers1$r_not <- "All"
rand_triggers1$variable <- ifelse(rand_triggers1$variable=="Forecasting", "Epidemic", "Prevalence")

load("../data/rand_trigger_extra_high_risk.Rdata")
rand_triggers2 <- rand_triggers[which(rand_triggers$disc_prob==0.0224),]
rand_triggers2$disc_prob<-NULL
rand_triggers2$r_not <- "High Risk"
rand_triggers2$variable <- ifelse(rand_triggers2$variable=="Forecasting", "Epidemic", "Prevalence")

rand_trigger_dat <- rbind(rand_triggers1, rand_triggers2)
rand_trigger_dat <- rand_trigger_dat[which(rand_trigger_dat$variable=="Epidemic"), ]
rand_trigger_plot <- ggplot(rand_trigger_dat, aes(confidence, value, color=r_not)) + geom_line(size=1)+
  coord_cartesian(xlim = c(0,1.025), ylim=c(0,45),expand=FALSE)+
  scale_color_manual(values=c("Grey", "Black"))+
  scale_linetype_manual(values=c(10,1))+
  theme_cowplot()%+replace% theme(legend.position=c(0.25, 0.75),
                                  legend.box.just = "left")+
  labs(x = "Probability of an Epidemic", 
       y = "Trigger (Reported Cases)", 
       linetype = expression("R"[0]),
       color="Trigger Type")
print(rand_trigger_plot)

alison_fig <- plot_grid(plot.importation.log, texas.plot.trial, rand_trigger_plot, nrow = 1, labels = "AUTO")
save_plot(filename = "../ExploratoryFigures/alison_fig.pdf", plot = alison_fig, base_height = 4, base_aspect_ratio = 3)


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

get_rand_trials <- function(dir_path, import_loc, disc_prob, num_trials, high_risk=FALSE){
  imports <- read.csv(import_loc)
  dir_path <-"~/projects/zika_alarm/data/sep_intros/"
  
  projected_combos <- as.character(unique(interaction(imports$rnott.expected.round, imports$importation.projected,sep = "_")))
  
  values <- unlist(strsplit(projected_combos, "_"))
  ## R0s are the first item in each combo
  r_nots <- as.numeric(values[seq(1,length(values), by=2)])
  intro_rates <- as.numeric(values[seq(2,length(values), by=2)])
  
  ## If you want to do only high risk include this line:
  if(high_risk){
    indices <- which(r_nots>=1)  
    r_nots <- r_nots[indices]
    intro_rates <- intro_rates[indices]
  } 
  indices <- 1:length(r_nots)
  num_samps <- num_trials
  samps <- sample(indices, size = num_samps, replace = T)
  samps <- table(samps)
  rand_trials <- vector("list", length = num_samps)
  list_ind <- 1
  
  for(i in 1:length(samps)){
    index <- as.numeric(names(samps[i]))
    num_trials <- samps[i]
    
    load(get_vec_of_files(dir_path, r_nots = r_nots[index], disc_probs = disc_prob, intro_rates = intro_rates[index]))
    rand_trials[list_ind:(list_ind+num_trials-1)] <- trials[sample(1:10000, size=num_trials)]
    list_ind <- list_ind+num_trials
  }
  rand_trials  
}

#############################################
## Creating Random county trigger data
#############################################
## Generate random trial data
set.seed(808)
rand_trials <- get_rand_trials(dir_path, "../csvs/county_master.csv", disc_prob=0.0224, 10000)
save(list = c("rand_trials"), file = "../data/rand_trials/all_risk_county_trials_0.0224.Rdata")
rand_trials <- get_rand_trials(dir_path, "../csvs/county_master.csv", disc_prob=0.011, 10000)
save(list = c("rand_trials"), file = "../data/rand_trials/all_risk_county_trials_0.011.Rdata")
rand_trials <- get_rand_trials(dir_path, "../csvs/county_master.csv", disc_prob=0.0224, 10000, high_risk = TRUE)
save(list = c("rand_trials"), file = "../data/rand_trials/high_risk_county_trials_0.0224.Rdata")
rand_trials <- get_rand_trials(dir_path, "../csvs/county_master.csv", disc_prob=0.011, 10000, high_risk = TRUE)
save(list = c("rand_trials"), file = "../data/rand_trials/high_risk_county_trials_0.011.Rdata")

## Analyze saved trials and create data frame of epidemic and prevalence probabilities


data.files <- list.files(path="../data/rand_trials", pattern="*.Rdata", full.names=T, recursive=FALSE)
prob_data <- ldply(data.files, function(x) {
  load(x)
  ## Get the R0
  disc_prob <- get_disc_prob_rand(x) 
  risk_level <- get_risk_level_rand(x)
  prob_below <- get_prob_below_threshold(trials = rand_trials, f=totalprev_by_totaldetects, threshold=20, max_detect = 200)
  prob_below$prob_below <- 1 - prob_below$prob_below
  prob_epidemic <- get_epidemic_prob_by_d(trials = rand_trials, prev_threshold = 50, cum_threshold = 2000, max_detect = 200, num_necessary = 100)
  cbind(risk_level=risk_level, disc_prob=disc_prob, prob_below, data.frame(prob_epidemic=prob_epidemic$prob_epidemic))
})  
prob_data <- melt(prob_data, measure.vars = c("prob_below", "prob_epidemic"))
save(list = c("prob_data"), file = "../data/rand_county_prob_data.Rdata")

#############################################################################


# 
# load("../data/rand_trials/all_risk_county_trials_0.011.Rdata")
# test <- get_prob_below_threshold(rand_trials,totalprev_by_totaldetects, threshold = 20, max_detect = 200)
# 
# rand_trials
# 
# disc_probs <- c(0.011, 0.0224)
# indices <- 1:length(r_nots)
# df_all <- data.frame()
# for(disc_prob in disc_probs){
#   num_samps <- 10000
#   samps <- sample(indices, size = num_samps, replace = T)
#   samps <- table(samps)
#   rand_trials <- vector("list", length = num_samps)
#   list_ind <- 1
# 
#   for(i in 1:length(samps)){
#     index <- as.numeric(names(samps[i]))
#     num_trials <- samps[i]
#     #  print(i)
#     #  print(get_vec_of_files(dir_path, r_nots = r_nots[index], disc_probs = 0.011, intro_rates = intro_rates[index]))
#     load(get_vec_of_files(dir_path, r_nots = r_nots[index], disc_probs = disc_prob, intro_rates = intro_rates[index]))
#     rand_trials[list_ind:(list_ind+num_trials-1)] <- trials[sample(1:10000, size=num_trials)]
#     list_ind <- list_ind+num_trials
#   }
# 
#   conf <- seq(0.05, 0.95, by=0.05)
#   epi_triggers <- c()
#   prev_triggers <- c()
# 
#   for(ci in conf){
#     temp <- get_epidemic_trigger(trials = rand_trials, threshold = 50, confidence = ci, max_detect = 300, num_necessary = 100)
#     temp2 <- get_surveillance_trigger(trials=rand_trials, threshold = 20, confidence=ci, max_detect = 300, num_necessary=100)
#     epi_triggers <- c(epi_triggers, temp)
#     prev_triggers <- c(prev_triggers, temp2)
#   }
#   df <- data.frame(confidence=conf, epi_trigger = epi_triggers, prev_trigger = prev_triggers[length(prev_triggers):1])
#   df <- melt(df, id.vars=c("confidence"))
#   df$variable <- ifelse(df$variable=="epi_trigger",  "Forecasting", "Nowcasting")
#   df$disc_prob <- disc_prob
#   df_all <- rbind(df_all, df)
# }
# rand_triggers <- df_all
# 
# # # ggplot(df_all, aes(confidence, value, color=variable, linetype=as.factor(disc_prob)))+geom_line()
# save(list = c("rand_triggers"), file = "../data/rand_trigger_extra_high_risk.Rdata")
