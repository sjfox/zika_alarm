rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

# county <- read.csv("../county_data.csv")
# 
# hist(county$HabitatSuitability, breaks=100)
# county[which.max(county$HabitatSuitability),]

sapply(c('branch.functions.R','plot.functions.R', 'incubation_branch.R', 'analyze_saved_sims.R'), source)
library(plyr)
library(cowplot)
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
branch_params <- function(r_not = 0.9,
                          infBoxes = 3,
                          incBoxes = 6,
                          recov_p = 0.3040571/(3/infBoxes),
                          incub_rate = 0.583917,
                          prop_p =  r_not*recov_p/infBoxes, 
                          e_thresh = 1000,
                          prob_symp = 1,
                          dis_prob_symp = .01,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000)
  return(as.list(environment()))

dir_path <- "~/projects/zika_alarm/data/sep_intros/"
trigger_dir_path <- "~/projects/zika_alarm/data/new_triggers/"
save_path <- "~/projects/zika_alarm/data/"
fig_path <- "~/projects/zika_alarm/ExploratoryFigures/"



# plot_final_sizes(trials)
# plot_intro_final_sizes(trials)
# plot_local_final_sizes(trials)
# plot_max_nonintro_prevalences(trials)
# plot_max_prevalences(trials)


######### Combine triggers after a tacc run/download
combine_triggers(trigger_dir_path, save_path)


get_trigger_data(1, 0.3, 0.011, 20, c(0.5,0.7,0.6,0.8))

intros <- c(0.01)
det_probs <- c(0.011)
r_nots <- c(0.8, 0.85, seq(0.9, 1, by=0.01))

dotplot_data_detect <- function(dir_path, r_nots, disc_probs, intro_rates){
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(dirPaths, function(x) {
    load(x)
    detects <- all_last_cumdetect_values(trials)
    prevs <- all_max_prevalence(trials)
    thresh_vals <- which(prevs>20)
    # prevs <- all_last_cuminfect_values(trials)
    # thresh_vals <- which(prevs>=1000)
    
    parms <- get_parms(x)
    
    if(length(thresh_vals)>0){
      cbind(as.data.frame(parms), detects=detects[thresh_vals])  
    } else{
      cbind(as.data.frame(parms), detects=NA)  
    }
  })
}

test <- dotplot_data_detect(dir_path, r_nots = r_nots, disc_probs = det_probs, intro_rates = intros)
ggplot(test, aes(r_not, detects)) + facet_wrap(~intro_rate, nrow=1)+ 
  geom_point(position="jitter", shape=20,alpha=0.5) + geom_hline(yintercept=20, color="red") +
  theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.01,0.01))+
  panel_border(size=0.5, colour="black") +
  labs(x = expression("R"[0]), y= "Cumulative Detections")




################################
## Code to Make Figure 2
################################
##### Panel A
r_nots <- c(0.9, 1.3)
disc_probs <- c(0.011, 0.0224)
intros <- 0.1 ## Harris county has 0.129, but this close
prev_plot_data <- get_prev_by_detects_plot(dir_path, r_nots = r_nots, disc_probs = disc_probs, intro_rates = intros)
prev_plot <- plot_prevalences(prev_plot_data)
print(prev_plot)  

##### Panel B
r_nots <- c(seq(0.91, 1.2, by=0.01), 1.25, seq(1.3, 2, by=0.1))
intros <- c(0.01, .1)
det_probs <- c(0.011, 0.0224)

prev_triggers <- get_trigger_data(r_nots, intros, det_probs, threshold=20, confidence=0.8)


## Remove noise from low R0 values
prev_triggers$prev_trigger[which(prev_triggers$prev_trigger>100 | is.na(prev_triggers$prev_trigger))] <- Inf

prev_trigger_plot <- ggplot(prev_triggers, aes(r_not, prev_trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) + 
  geom_line(size=1) + scale_y_continuous(expand=c(0.01,0.01), limits=c(1,100))+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0.9,1.6))+
  scale_color_manual(values=c("Grey", "Black"))+
  theme_cowplot()%+replace% theme(legend.position="none")+
  labs(x = expression("R"[0]), 
       y = "Prevalence Trigger", 
       color = "Importation\nRate", 
       linetype= "Detection \nProbability")
print(prev_trigger_plot)
######## Panel C
epi_triggers <- prev_triggers
# epi_triggers$epi_trigger <- ifelse(is.na(epi_triggers$epi_trigger), Inf, epi_triggers$epi_trigger)
epi_triggers$disc_prob <- paste0(calculate.discover(epi_triggers$disc_prob), "%")

## Remove outlier from low trigger data ( few runs, and so not infinity solely due to stochasticity)
epi_triggers$epi_trigger[which(epi_triggers$epi_trigger>100| is.na(prev_triggers$epi_trigger))] <- Inf

epi_prob_plot <- ggplot(epi_triggers, aes(r_not, epi_trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) + 
  geom_line(size=1) + scale_y_continuous(expand=c(0.01,0.01), limits=c(0,100))+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0.9,1.3))+
  scale_color_manual(values=c("Grey", "Black"))+
  labs(x = expression("R"[0]), 
       y = "Future Epidemic Trigger", 
       color = "Importation\nRate", 
       linetype= "Detection \nProbability")
print(epi_prob_plot)

################
# Print out the panels all together and gridded
figure2_panels <- plot_grid(prev_plot, prev_trigger_plot, epi_prob_plot, labels=c("A", "B", "C"), nrow = 1, rel_widths = c(1.2,1,1.2))
save_plot(paste0(fig_path, "figure2_panels.pdf"), figure2_panels, base_height = 4, base_aspect_ratio = 3)





# ######### Supplemental figures
### Porbability below threshold graph
# thresholds <- c(25)
# r_nots <- c(0.9, 1.1, 1.3)
# disc_probs <- c(0.011,.0287, 0.068)
# intro_rates <- c(.05, 0.3)
# temp <- get_prob_below_plot(dir_path, thresholds, r_nots, disc_probs, intro_rates)
# 
# prob_below <- plot_prob_below(temp)
# save_plot(paste0(fig_path, "probability_below_supplemental.pdf"),prob_below, base_height = 5, base_aspect_ratio = 2)

### Epidemic probability graph
# r_nots <- c(0.9, 1, 1.3)
# disc_probs <- c(0.011, 0.0287, 0.068)
# intro_rates <- c(.05, 0.1)
# temp <- get_epidemic_prob_plot(dir_path, prev_threshold = 25, cum_threshold = 1000, r_nots, disc_probs, intro_rates)
# 
# epi_plot <- plot_epidemic_prob(temp)
# save_plot(paste0(fig_path, "epi_prob_supplemental.pdf"),epi_plot, base_height = 5, base_aspect_ratio = 1.7)
# 


# get_rows_final_sizes <- function(final_sizes, r_nots, disc_probs, intro_rates){
#   final_sizes[which(final_sizes$r_not %in% r_nots & final_sizes$disc_prob %in% disc_probs & final_sizes$intro_rate %in% intro_rates), ]
# }
# load("../data/final_sizes.Rdata")
# r_nots <- c(0.8, 1.2)
# disc_probs <- c(0.011)
# intro_rates <- c(0, .01, 0.3, 1.92)
# 
# temp <- get_rows_final_sizes(final_sizes, r_nots, disc_probs, intro_rates)
# 
# final_plot <- ggplot(temp, aes(cumI)) + geom_histogram(bins=100)+
#   facet_grid(intro_rate~r_not, scales = "free_y") + 
#   scale_y_continuous(expand=c(0,0)) + 
#   scale_x_continuous(expand=c(0,0))+
#   labs(x = "Final Epidemic Sizes")
# save_plot(filename = "../ExploratoryFigures/final_sizes.pdf", plot = final_plot, base_height = 10, base_width = 12)
# 
# 

#####################################################
#### Figure 1, understanding the overarching question!
#####################################################
get_symp <- function(n){
  sum(sample(c("symp", "asymp"), n, prob = c(0.15, 0.85), replace = T) =="symp")  
}

set.seed(808)
head(trials[[5]], n = 100)$Cumulative_Infections
zoom_data <- trials[[5]]
zoom_data <- zoom_data[zoom_data$time <=125,]
zoom_data$Symptomatic <- sapply(zoom_data$Total_Infections, get_symp)
zoom_data$Asymptomatic <- zoom_data$Total_Infections-zoom_data$Symptomatic
zoom_data$min_asymp <- zoom_data$Symptomatic
zoom_data$min_symp <- 0
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Total_Infections", "Symptomatic"))
zoom_melt$variable <- ifelse(zoom_melt$variable=="Total_Infections", "Asymptomatic", "Symptomatic")
zoom_melt$mins <- c(zoom_data$min_asymp, zoom_data$min_symp)
arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$Total_Infections[which(zoom_data$New_Intro!=0)]
det_arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Detections!=0)])
det_arrow_data$yval <- zoom_data$Total_Infections[which(zoom_data$New_Detections!=0)]

zoom_plot <- ggplot(zoom_melt, aes(x=time, ymax=value, ymin=mins, color=variable, fill=variable)) + geom_ribbon()+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+1.5, yend=yval), 
               arrow = arrow(length = unit(0.05, "npc"), angle = 35), color="red", size=1, inherit.aes=FALSE)+
  geom_vline(data=det_arrow_data, aes(xintercept=time), linetype=2, color="red")+ 
  theme_cowplot() %+replace% theme(legend.position="none")+
  scale_y_continuous(expand=c(0.0,0.0))+
  scale_x_continuous(expand=c(0.01,0.01),limits=c(0,120))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Infections", fill="")
print(zoom_plot)

ind <- which(all_last_cuminfect_values(trials) > 20 & all_last_cuminfect_values(trials) <30)
##  8 may work
zoom_data <- trials[[ind[9]]]
# zoom_data <- zoom_data[zoom_data$time <=75,]
zoom_data$Symptomatic <- sapply(zoom_data$Total_Infections, get_symp)
zoom_data$Asymptomatic <- zoom_data$Total_Infections-zoom_data$Symptomatic
zoom_data$min_asymp <- zoom_data$Symptomatic
zoom_data$min_symp <- 0
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Total_Infections", "Symptomatic"))
zoom_melt$variable <- ifelse(zoom_melt$variable=="Total_Infections", "Asymptomatic", "Symptomatic")
zoom_melt$mins <- c(zoom_data$min_asymp, zoom_data$min_symp)
arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$Total_Infections[which(zoom_data$New_Intro!=0)]
det_arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Detections!=0)])
det_arrow_data$yval <- zoom_data$Total_Infections[which(zoom_data$New_Detections!=0)]

ylims <- ggplot_build(zoom_plot)$panel$ranges[[1]]$y.range

zoom_plot2 <- ggplot(zoom_melt, aes(x=time, ymax=value, ymin=mins, color=variable, fill=variable)) + geom_ribbon()+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+1.5, yend=yval), 
               arrow = arrow(length = unit(0.05, "npc"), angle = 35), color="red", size=1, inherit.aes=FALSE)+
  geom_vline(data=det_arrow_data, aes(xintercept=time), linetype=2, color="red")+ 
  theme_cowplot() %+replace% theme(legend.position=c(0.9,.8))+ 
  scale_y_continuous(expand=c(0.001,0.001), limits = c(0,ylims[2]))+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0,120))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Infections", fill="")
print(zoom_plot2)


all_ts <- ggdraw() +
  draw_plot(zoom_plot2, x =  0,y =  0.5,width =  1,height =  0.5) +
  draw_plot(zoom_plot, x =  0,y =  0,width =  1, height = .5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = 15)

save_plot(filename = "../ExploratoryFigures/all_ts2.pdf", plot = all_ts, base_height=4, base_aspect_ratio = 1.8)

# ## Supplementary choosing threshold plot -- should be 20
# r_nots <- c(0.8, 0.9, 1, 1.1, 1.2)
# disc_probs <- c(0.011)
# intros <- c(0.01, 0.1, 0.3)
# threshold_plot <- plot_dots(dir_path, r_nots, disc_probs, intros)
# save_plot(paste0(fig_path, "threshold_plot.pdf"), threshold_plot, base_height = 4, base_aspect_ratio = 1.8)

