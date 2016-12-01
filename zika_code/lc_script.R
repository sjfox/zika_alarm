rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/projects/zika_alarm/zika_code/')
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



#############################################
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








##################### Working on Time Distributions

dir_path = "~/Documents/projects/zika_alarm/data/zika_sims/"
fig_path = "~/Documents/projects/zika_alarm/ExploratoryFigures/"

time_by_detect <- function(df, max_detect=100, epi_thres = 50){
  ## Takes in a data frame trials, and for each
  ## returns time lapse between detections
  
  max.prev <- max_nonintro_prevalence(df)
  
  if (max.prev > epi_thres) {
  
  local.detections <- df$New_Detections-df$New_Intro_Detections 
  indicies <- which(local.detections > 0 )
  difference <- diff(indicies)
  #difference <- c(indicies[1], difference) # Adding in the time spance of the first 
  return(data.frame(difference))
  }
}

get_time_between_detects_all <- function(x, max_detect=100){
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, time_by_detect, max_detect)
}


dotplot_tdiff_data <- function(dir_path, r_nots, disc_probs, intro_rates, local){
  # Function to get all the differences between two local detection
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(dirPaths, function(x) {
    load(x)
    if(local){
      differences <- get_time_between_detects_all(trials) 
      if (nrow(differences) > 0) {
        parms <- get_parms(x)
        cbind(as.data.frame(parms), t_diff=differences$difference[sample(seq(0,10000), 1000)])
      }
    }else{
      differences <- get_time_between_detects_all(trials)  
    }
    #parms <- get_parms(x)
    #cbind(as.data.frame(parms), t_diff=differences$difference[sample(seq(0,10000), 1000)])
  })
}



time_by_detect_summary <- function(df, max_detect = 100) {
  # takes in a data frame of trials and for each
  # returns the summary statistics of the time lapses
  local.detections <- df$New_Detections-df$New_Intro_Detections 
  indicies <- which(local.detections > 0 )
  #indicies <- which(df[, "New_Detections"- "New_Intro_Detections"] > 0 )
  
  difference <- diff(indicies)
  summary <- unname(summary(difference))
  return(data.frame(median = summary[3], third = summary[5]))
}

time_by_detect_summary(df)


get_tdiff_summary_all <- function(x, max_detect = 100) {
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, time_by_detect_summary, max_detect)
}

trial$ <- get_tdiff_summary_all(trials)



dotplot_tdiff_summary_data <- function(dir_path, r_nots, disc_probs, intro_rates, type){
  # Function to get all the differences between two local detection
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(dirPaths, function(x) {
    load(x)
    differences <- get_tdiff_summary_all(trials)  
    if(type == 'median'){
      differences <- differences$median
    }else{
      differences <- differences$third 
    }
    parms <- get_parms(x)
    cbind(as.data.frame(parms), t_diff=differences)
  })
}



r_nots <- c(seq(0.7, 1.2, by=0.05))
disc_probs <- c(0.0224)
intros <- c(0.1)

#Harris County

r_nots <- c(0.8)
intros <- c(0.269)


tdiff.dots <- dotplot_tdiff_data(dir_path, r_nots, disc_probs, intro_rates = intros, local=T)

tdiff.dots.cleaned <- tdiff.dots[!is.na(tdiff.dots$t_diff),]
head(tdiff.dots.cleaned)

fake.row <- data.frame(cbind(0.8, 0.0224, 0.1, NA))
colnames(fake.row) <- colnames(tdiff.dots.cleaned)
plot.data.trail <- rbind(fake.row, tdiff.dots.cleaned)


tdiff_plot <- ggplot(plot.data.trail, aes(r_not, t_diff)) + #+ facet_wrap(~intro_rate, nrow = 1) + 
  geom_point(position="jitter", shape=20,alpha=0.5) + 
  geom_hline(yintercept=c(15), color=c("blue")) +
 # geom_hline(yintercept = c(30), color = c("red")) + 
  theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.01,0.01), breaks = r_nots)+
  panel_border(size=0.5, colour="black") +
  labs(x = expression("R"[0]), y= "Days Between Local Detection")


tdiff_plot_box <- ggplot(plot.data.trail, aes(factor(r_not), t_diff)) + 
  geom_boxplot() + 
  geom_hline(yintercept=c(15), color=c("blue")) +
  #geom_hline(yintercept = c(30), color = c("red")) + 
  theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
  scale_y_continuous(expand=c(0.01,0.01))+
  #scale_x_continuous(expand=c(0.01,0.01), breaks = r_nots)+
  panel_border(size=0.5, colour="black") +
  labs(x = expression("R"[0]), y= "Days between local detection")


tdiff_plot

save_plot(paste0(fig_path, "local_tdiff_points.pdf"), tdiff_plot, base_height = 4, base_aspect_ratio = 1.8)


## median plot 
median.dots <- dotplot_tdiff_data(dir_path, r_nots, disc_probs, intro_rates = intros, type = 'median')
head(median.dots)
head(median.dots.cleaned)
mean(median.dots.cleaned$t_diff)

median.dots.cleaned <- median.dots[!is.na(median.dots$t_diff),]

tdiff_median_plot <- ggplot(median.dots.cleaned, aes(r_not, t_diff)) + facet_wrap(~intro_rate, nrow = 1) + 
  geom_point(position="jitter", shape=20,alpha=0.5) + 
  geom_hline(yintercept=c(14), color=c("blue")) +
  theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.01,0.01), breaks = r_nots)+
  panel_border(size=0.5, colour="black") +
  labs(x = expression("R"[0]), y= "Median days between local detection")

save_plot(paste0(fig_path, "local_median_tdiff_plot.pdf"), tdiff_median_plot, base_height = 4, base_aspect_ratio = 1.8)



# third quarter plot 
#third.quarter.dots <- dotplot_tdiff_data(dir_path, r_nots, disc_probs, intro_rates = intros, type = 'third')
#third.quarter.cleaned <- third.quarter.dots[!is.na(third.quarter.dots$t_diff),] 

#tdiff_third_plot <- ggplot(third.quarter.cleaned, aes(r_not, t_diff)) + facet_wrap(~intro_rate, nrow = 1) + 
#  geom_point(position="jitter", shape=20,alpha=0.5) + 
#  geom_hline(yintercept=c(14), color=c("blue")) +
#  theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
#  scale_y_continuous(expand=c(0.01,0.01))+
#  scale_x_continuous(expand=c(0.01,0.01), breaks = r_nots)+
#  panel_border(size=0.5, colour="black") +
#  labs(x = expression("R"[0]), y= "Third Quarter Summary of days between local detection")

#save_plot(paste0(fig_path, "local_third_tdiff_plot.pdf"), tdiff_third_plot , base_height = 4, base_aspect_ratio = 1.8)



# Normal Plot 
dots <- dotplot_tdiff_data(dir_path, r_nots, disc_probs, intros, local = T)
hist(dots$t_diff)
summary(dots)

baseline = dots[dots$intro_rate == 0,]
head(baseline)

baseline.plot <- ggplot(baseline, aes(t_diff)) + geom_histogram() + 
 labs(x = 'Time Between Local Detections', y = "Count") +
  geom_vline(xintercept = c(15), color = c("blue")) +
  scale_x_continuous(breaks = seq(0,160, 10))

save_plot(paste0(fig_path, "baseline_tdiff.pdf"), baseline.plot, base_height = 4, base_aspect_ratio = 1.8)





plot_tdiff_dots <- function(dir_path, r_nots, disc_probs, intros, local=FALSE){
  
  dots <- dotplot_data(dir_path, r_nots, disc_probs, intros, local)
  
  ggplot(dots, aes(r_not, max_prev)) + facet_wrap(~intro_rate, nrow=1)+ 
    geom_point(position="jitter", shape=20,alpha=0.5) + 
    geom_hline(yintercept=c(50), color=c("blue")) +
    theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
    scale_y_continuous(expand=c(0.01,0.01))+
    scale_x_continuous(expand=c(0.01,0.01), breaks = r_nots)+
    panel_border(size=0.5, colour="black") +
    labs(x = expression("R"[0]), y= ifelse(local,"Maximum Local Daily Prevalence", "Maximum Daily Prevalences"))
}










# print(threshold_plot)
save_plot(paste0(fig_path, "local_threshold_plot.pdf"), threshold_plot, base_height = 4, base_aspect_ratio = 1.8)

projected_worst_combos <- data.frame(sort(projected_worst_combos$projected_worst_combos, decreasing = FALSE))
load(get_vec_of_files(dir_path))

load(get_vec_of_files(dir_path, 0.9, 0.0224, 0.432))

finaldetection <- all_last_cumdetect_local_values(trials)
max_local_prev <- all_max_nonintro_prevalence(trials)
#durations <- all_sim_duration(trials)

cum.cases.by.detect <- get_cumcases_by_detects_all(trials)

plot(durations, finaldetection)
plot(finaldetection, max_local_prev)


differences.trials <- get_time_between_detects_all(trials)
summary(differences.trials$difference)
hist(differences.trials$difference)



# ## Supplementary choosing threshold plot -- should be 20
r_nots <- c(seq(0.9, 1.1, by=0.05))
disc_probs <- c(0.011)
intros <- c(0.01, 0.1)
threshold_plot <- plot_dots(dir_path, r_nots, disc_probs, intros, local=T)
# print(threshold_plot)
save_plot(paste0(fig_path, "local_threshold_plot.pdf"), threshold_plot, base_height = 4, base_aspect_ratio = 1.8)



epi.data <- get_epi_data(trials, 3000)

load(get_vec_of_files(dir_path, 1.01, 0.0224, 0))

#epi.data_0.7 <- get_epi_data(trials, 3000)

outbreak_plot_1.1 

breaks <- c(1,2,3,4,5,10,50, 100)

outbreak_plot_0.7 <- ggplot(epi.data_0.6, aes(time, Cum_Detections)) + geom_line(alpha = 0.15, size = 1) +xlim (0,100) + scale_y_log10()
outbreak_plot_0.7
outbreak_plot_1.1

                                      
                                      #, group=interaction(.id, risk_level), color=risk_level)) + 
  geom_line(alpha=0.15,size=1) + 
 # scale_color_brewer(palette="Set1", direction = 1)+
  coord_cartesian(xlim=c(0,150), ylim=c(0,50), expand=FALSE) +
  guides(color=guide_legend(override.aes=list(alpha=1))) +
  theme(legend.position=c(0.3,0.8))+
  labs(x = "Time (days)", 
       y = "Reported Autochthonous Cases", 
       color = expression("R"[0]))
print(outbreak_plot)

# Looking at the distributions of time



