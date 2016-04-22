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
# Code I've written just testing things out-not worth saving in the main files 
#Parameters 




prev.map <- plotheatmaps(thres.matrix.prev, type = "Prevalence", names = as.factor(dect.cases.range),
                         R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.prev, 1))
cum.map <- plotheatmaps(thres.matrix.cum, type = "Cumulative", names = as.factor(dect.cases.range),
                        R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.cumulative, 1))



 

############################################
dir_path <- "~/Documents/zika_alarm/data/introductions/"
save_path <- "~/Doucments/zika_alarm/data"

r_nots <- c(1.9, 1.1, 0.8)
disc_prob <- c( 0.068, 0.011)
intro_rate <- c(0.01, 0.1, 0.3)


#expected.cases.local <- calculate_expect_vs_detect_local(dir_path, r_nots, intro_rate, disc_prob)
#colnames(expected.cases.local) <- c("run", "r0","dect", "intro", "cases_dect", "avg", "sd" )

expected.cases.total <- calculate_expect_vs_detect_total(dir_path, r_nots, intro_rate, disc_prob)
colnames(expected.cases.total) <- c("run", "r0","dect", "intro", "cases_dect", "avg", "sd" )

# If want to split the Results to certain detection values/intro values/R0
indices.total = which((expected.cases.total$intro == 0.01 | expected.cases.total$intro == 0.3) & expected.cases.total$r0 != 1.1)
detection.total = expected.cases.total[indices.total,]

#breaks_y = seq(from = 0, to = 100, by = 10)
breaks_y_log = c(1,2,3,4, 5,10,15,20,25,50,75,100) 
breaks_y = seq(0, 200, 20)
breaks_x = seq(from = 0, to = 50, by = 5)
head(detection.total)
tail(detection.total)

plot_log <- ggplot(detection.total, aes(cases_dect, avg, fill = as.factor(r0), 
                                              color = as.factor(r0), group=interaction(as.factor(dect), r0)))  + 
  geom_line(size=1.5, aes(linetype = as.factor(dect))) + facet_wrap(~intro, nrow = 1) +
  geom_ribbon(aes(ymin = avg-sd, ymax=avg+sd), alpha=.2, color = NA) + 
  scale_y_log10(breaks = breaks_y_log, limits = c(1,150))  +
  scale_color_brewer(palette = "Set1", guide = FALSE, direction = -1) +
  scale_fill_brewer(palette="Set1", guide_legend(title = "R0"), direction = -1) +
  scale_x_continuous(name = "Cumulative Number of Detected Cases", breaks = breaks_x, limits = c(0,40)) +
  labs(y = "Expected Total Current Cases", linetype = "Detection \n Rate") #+
  #theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) #+
  #theme(axis.title.y = element_text(size=28), axis.text.y = element_text(size = 22)) + 
  #theme(axis.title.x = element_text(size=28), axis.text.x= element_text(size=22)) +
  
  #theme(legend.text=element_text(size=22, margin = margin(), debug = FALSE), legend.title = element_text(size = 28)) #+ 
#
plot_log

plot_continuous <- ggplot(detection.total, aes(cases_dect, avg, fill = as.factor(r0), 
                                              color = as.factor(r0), group=interaction(as.factor(dect), r0)))  + 
  geom_line(size=1.5, aes(linetype = as.factor(dect))) + facet_wrap(~intro, nrow = 1) +
  geom_ribbon(aes(ymin = avg-sd, ymax=avg+sd), alpha=.2, color = NA) + 
  scale_y_continuous(breaks = breaks_y, limits = c(0,150))  +
  scale_color_brewer(palette = "Set1", guide = FALSE, direction = -1) +
  scale_fill_brewer(palette="Set1", guide_legend(title = "R0"), direction = -1) +
  scale_x_continuous(name = "Cumulative Number of Detected Cases", breaks = breaks_x, limits = c(0,40)) +
  labs(y = "Expected Total Current Cases", linetype = "Detection \n Rate") #+
#theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) #+
#theme(axis.title.y = element_text(size=28), axis.text.y = element_text(size = 22)) + 
#theme(axis.title.x = element_text(size=28), axis.text.x= element_text(size=22)) +
#theme(legend.text=element_text(size=22, margin = margin(), debug = FALSE), legend.title = element_text(size = 28)) #+ 
plot_continuous
plot_grid(plot.1sd.local, plot.1sd.total, labels = c("A", "B")) 


##### TRYING TO REWRITE TEXAS TRIGGER THRESHOLD


###### Calculate Surveillance For R0s 
r_nots <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0) 
disc_prob <- c(0.011, 0.0052)
intro_rate <- c(0.1)


#Calculate the triggers  
r0.triggers <- calculate_all_triggers(dir_path, r_nots = r_nots, intro_rate = intro_rate, 
                                      disc_prob = disc_prob, threshold = 20, confidence = .8)

colnames(r0.triggers) <- c("run", "r0", "detect", "intro", "threshold", "confidence", "trigger")
trigger.maxdetect <- which(r0.triggers[,"trigger"] > 100)
r0.triggers[trigger.maxdetect, "trigger"] <- NA

#Quick plot 
plot(r0.triggers$r0[r0.triggers$detect == .0052], r0.triggers$trigger[r0.triggers$detect == .0052], col = "blue", pch = 19)
points(r0.triggers$r0[r0.triggers$detect == .011], r0.triggers$trigger[r0.triggers$detect == .011], col = "red", pch = 19)



#Combining with Texas Maps 
setwd('..'); setwd('TexasCountyShapeFiles')
texas.county <- readShapeSpatial('texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
setwd('../zika_code/')
county_ids <- read.csv(file = "~/Documents/zika_alarm/county_ids.csv")

#Calculate Only triggers interested in
# Base Case 
county_plot$constant.import <- NA

#BaseLine
desired_dect = 0.011; desired_intro = .1
#Good Detect Bad Intro
desired_dect = 0.011; desired_intro = .3
#Bad Detect Good Intro
desired_dect = 0.0052; desired_intro = .1
#Worst
desired_dect = 0.0052; desired_intro = .3


triggers <- ddply(.data = r0.triggers, .variables = "r0", function (x) {
  row = which(x[,"detect"] == desired_dect & x[,"intro"] == desired_intro)
  return(x[row,])
})

# Merges The Data With the County Data for Plotting By R0 
for (i in 1:nrow(triggers)) {
  indices = which(county_plot$metro_round == triggers[i,"r0"]) 
  county_plot$constant.import[indices] = triggers[i,"trigger" ]
}

county_plot$constant.import <- as.numeric(county_plot$constant.import)

texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_plot, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]


#Actual plotting 
max(r0.triggers[,"trigger"], na.rm = TRUE)
legend_breaks <- round(seq(0, 55, 10))

plot.constant <- ggplot()+geom_polygon(data = final.plot, aes_string(x="long", y = "lat", group = "group", fill = "constant.import"),
                                   color = "black", size = .25) + coord_map() +
  scale_fill_continuous(name = "Detected \n Cases", low = "red", high = "yellow", 
                      na.value = "grey") + #+ pretty_breaks(n = length(legend_breaks)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(legend.position = "right") +
  theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
  theme(legend.key.size =  unit(0.3, "in")) 

ggsave(filename = "plot.constant.pdf", plot.constant, width = 4, height = 4, unit = "in")

plotworst
plot.baddetect_goodintro
plot.gooddetect_badintro
plotbaseline


combined <- plot_grid(plotbaseline, plot.baddetect_goodintro, 
          plot.gooddetect_badintro, plotworst, labels = c("A", "B", "C", "D"))

r0.triggers












