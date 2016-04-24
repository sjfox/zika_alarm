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
save_path <- "~/projects/zika_alarm/data/"
fig_path <- "~/projects/zika_alarm/ExploratoryFigures/"

load(get_vec_of_files(dir_path, 1.2, .068, 0.05))

get_epidemic_trigger(trials, 25, 0.8)

get_surveillance_trigger(trials, 25, 0.8)








print(prev_trigger_plot)
# plot_final_sizes(trials)
# plot_intro_final_sizes(trials)
# plot_local_final_sizes(trials)
# plot_max_nonintro_prevalences(trials)
# plot_max_prevalences(trials)


##### Panel A
r_nots <- c(0.9, 1.3)
disc_probs <- c(0.011, 0.0287)
intros <- 0.1
prev_plot_data <- get_prev_by_detects_plot(dir_path, r_nots = r_nots, disc_probs = disc_probs, intro_rates = intros)
prev_plot <- plot_prevalences(prev_plot_data)
# print(prev_plot)  

##### Panel B
r_nots <- seq(0.5, 1.9, by=0.1)
intros <- c(0.01, .1)
det_probs <- c(0.011, 0.0287)

prev_triggers <- calculate_all_triggers(dir_path, r_nots, intro_rate = intros, disc_prob = det_probs, threshold = 25, confidence = 0.8)
prev_triggers$trigger <- ifelse(is.na(prev_triggers$trigger), Inf, prev_triggers$trigger)

prev_trigger_plot <- ggplot(prev_triggers, aes(r_not, trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) + 
  geom_line(size=1) + scale_y_continuous(expand=c(0.01,0.01), limits=c(1,150))+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0.9,1.9))+
  scale_color_manual(values=c("Grey", "Black"))+
  theme_cowplot()%+replace% theme(legend.position="none")+
  labs(x = expression("R"[0]), 
       y = "Prevalence Trigger", 
       color = "Importation\nRate", 
       linetype= "Detection \nProbability")

######## Panel C
epi_triggers <- calculate_all_epidemics(dir_path, r_nots, intros, det_probs, 25, 0.80)
epi_triggers$trigger <- ifelse(is.na(epi_triggers$trigger), Inf, epi_triggers$trigger)
epi_triggers$disc_prob <- calculate.discover(epi_triggers$disc_prob)

epi_prob_plot <- ggplot(epi_triggers, aes(r_not, trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) + 
  geom_line(size=1) + scale_y_continuous(expand=c(0.01,0.01), limits=c(1,150))+
  scale_x_continuous(expand=c(0.01,0.01), limits=c(0.9,1.4))+
  scale_color_manual(values=c("Grey", "Black"))+
  labs(x = expression("R"[0]), 
       y = "Future Epidemic Trigger", 
       color = "Importation\nRate", 
       linetype= "Detection \nProbability")

figure2_panels <- plot_grid(prev_plot, prev_trigger_plot, epi_prob_plot, labels=c("A", "B", "C"), nrow = 1, rel_widths = c(1.2,1,1.2))
save_plot(paste0(fig_path, "figure2_panels.pdf"), figure2_panels, base_height = 4, base_aspect_ratio = 3)





# ######### Supplemental figures
### Porbability below threshold graph
thresholds <- c(25)
r_nots <- c(0.9, 1.1, 1.3)
disc_probs <- c(0.011,.0287, 0.068)
intro_rates <- c(.05, 0.3)
temp <- get_prob_below_plot(dir_path, thresholds, r_nots, disc_probs, intro_rates)

prob_below <- plot_prob_below(temp)
save_plot(paste0(fig_path, "probability_below_supplemental.pdf"),prob_below, base_height = 5, base_aspect_ratio = 2)

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


example <- get_vec_of_files(dir_path, 1.1, 0.011, .1)
load(example)

## Data for growing epidemic
plotDat <- trials[[7]]
plotDat <- melt(plotDat, id.vars = c("time"), measure.vars = c("New_Detections", "Total_Infections"))
plotDat <- plotDat[plotDat$time<100, ]
plotDat$variable <- ifelse(plotDat$variable=="Total_Infections", "Prevalence", "New Detections")
## Data for dying outbreak
plotDat2 <- trials[[13]]
plotDat2 <- melt(plotDat2, id.vars = c("time"), measure.vars = c("New_Detections", "Total_Infections"))
plotDat2 <- plotDat2[plotDat2$time<100, ]
plotDat2$variable <- ifelse(plotDat2$variable=="Total_Infections", "Prevalence", "New Detections")

plotDat <- rbind(cbind(plotDat, outbreak="Epidemic"), cbind(plotDat2, outbreak="Outbreak"))
plotDat <- plotDat[plotDat$time<75, ]

ind <- which(plotDat$variable=="Prevalence")
detects <- plotDat[-ind,]
ts_plot <- ggplot(plotDat[ind,], aes(time, value, color=outbreak, fill=outbreak)) + geom_line(size=1) +
  geom_bar(data=detects, aes(time, value), position="dodge",stat="identity", width = 2)+
  theme_cowplot() %+replace% theme(legend.position=c(0.2,0.8))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.03,0.05))+
  scale_color_manual(values=c("red", "black"), guide=FALSE) +
  scale_fill_manual(values=c("red", "black")) +
  labs(x = "Time (days)", y="Individuals", fill="", color = "", linetype="")
print(ts_plot)
# save_plot(filename = "../ExploratoryFigures/new_ts.pdf", plot = ts_plot, base_height=3, base_aspect_ratio = 2)


get_symp <- function(n){
  sum(sample(c("symp", "asymp"), n, prob = c(0.15, 0.85), replace = T) =="symp")  
}

set.seed(808)
zoom_data <- trials[[7]]
zoom_data <- zoom_data[zoom_data$time <75,]
zoom_data$Symptomatic <- cumsum(sapply(zoom_data$New_Infection, get_symp))
zoom_data$Asymptomatic <- zoom_data$Cumulative_Infections-zoom_data$Symptomatic
zoom_data$min_asymp <- zoom_data$Symptomatic
zoom_data$min_symp <- 0
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Cumulative_Infections", "Symptomatic"))
zoom_melt$variable <- ifelse(zoom_melt$variable=="Cumulative_Infections", "Asymptomatic", "Symptomatic")
zoom_melt$mins <- c(zoom_data$min_asymp, zoom_data$min_symp)
arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$Cumulative_Infections[which(zoom_data$New_Intro!=0)]

zoom_plot <- ggplot(zoom_melt, aes(x=time, ymax=value, ymin=mins, color=variable, fill=variable)) + geom_ribbon()+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+8, yend=yval), 
               arrow = arrow(length = unit(0.03, "npc")), color="green", size=1.5,inherit.aes=FALSE)+
  theme_cowplot() %+replace% theme(legend.position=c(0.2,0.6))+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.03,0.05))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Individuals", fill="")
print(zoom_plot)

zoom_data <- trials[[13]]
zoom_data <- zoom_data[zoom_data$time <75,]
zoom_data$Symptomatic <- cumsum(sapply(zoom_data$New_Infection, get_symp))
zoom_data$Asymptomatic <- zoom_data$Cumulative_Infections-zoom_data$Symptomatic
zoom_data$min_asymp <- zoom_data$Symptomatic
zoom_data$min_symp <- 0
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Cumulative_Infections", "Symptomatic"))
zoom_melt$variable <- ifelse(zoom_melt$variable=="Cumulative_Infections", "Asymptomatic", "Symptomatic")
zoom_melt$mins <- c(zoom_data$min_asymp, zoom_data$min_symp)
arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$Cumulative_Infections[which(zoom_data$New_Intro!=0)]
zoom_plot2 <- ggplot(zoom_melt, aes(x=time, ymax=value, ymin=mins, color=variable, fill=variable)) + geom_ribbon()+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+3, yend=yval), 
               arrow = arrow(length = unit(0.05, "npc")), color="green", size=1.5,inherit.aes=FALSE)+
  theme_cowplot() %+replace% theme(legend.position="none")+
  scale_y_continuous(expand=c(0.01,0.01))+
  scale_x_continuous(expand=c(0.03,0.05), limits=c(0,74))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Individuals", fill="", color = "", linetype="")
print(zoom_plot2)


all_ts <- ggdraw() +
  draw_plot(ts_plot, 0, 0, .5, 1) +
  draw_plot(zoom_plot, 0.5, 0.5, .5, .5) +
  draw_plot(zoom_plot2, 0.5, 0, .5, .5) +
  draw_plot_label(c("A", "B", "C",  "Outbreak", "Epidemic"), c(0, 0.5, 0.5, 0.65, 0.65), c(1, 1, 0.5,0.5,1), size = 15, color=c(rep("black", 4),"red") )

save_plot(filename = "../ExploratoryFigures/all_ts.pdf", plot = all_ts, base_height=4, base_aspect_ratio = 2)
