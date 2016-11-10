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
library(gridExtra)

dir_path <- "~/projects/zika_alarm/data/all_trials/"
trigger_dir_path <- "~/projects/zika_alarm/data/triggers_local/"
save_path <- "~/projects/zika_alarm/data/"
fig_path <- "~/projects/zika_alarm/ExploratoryFigures/"

######### Combine triggers after a tacc run/download #####################
# combine_triggers(trigger_dir_path, save_path)
# get_trigger_data(0.7, intro = 2, disc = 0.0224, confidence=0.5, num_necessary=100)
load(get_vec_of_files(dir_path, 0.9, 0.0224, 0.3))
sum(all_last_cuminfect_local_values(trials)>2000 & all_max_nonintro_prevalence(trials)>50)

################################
## Code to Make Figure 3
################################
get_epi_data <- function(trials, n){
  ## Returns first n trials in data frame form
  names(trials) <- seq_along(trials)
  ldply(trials[1:n], data.frame)
}

data.files <- list.files(path="../data/rand_trials", pattern="*.Rdata", full.names=T, recursive=FALSE)
load(data.files[4])
# plot_final_sizes(rand_trials)
temp <- ldply(data.files, function(x) {
  load(x)
  disc_prob <- get_disc_prob_rand(x)
  risk_level <- get_risk_level_rand(x)
  data <- get_epi_data(rand_trials, n= 3000)
  cbind(risk_level, disc_prob, data)
})  
# temp <- temp[which(temp$disc_prob==0.0224),]
temp <- temp[which(temp$risk_level=="high_risk" & temp$disc_prob==0.0224),]
temp$disc_prob <- paste0(calculate.discover(temp$disc_prob), "%")
# temp$disc_prob <- factor(temp$disc_prob, levels = c("20%", "10%"))

load(get_vec_of_files(dir_path, 1.1, 0.0224, 0.01))
known_rnot <- cbind(data.frame(risk_level="1.1", disc_prob=0.0224), get_epi_data(trials, 3000))
known_rnot$disc_prob <- paste0(calculate.discover(known_rnot$disc_prob), "%")
both_rnots <- rbind(temp, known_rnot)
both_rnots$risk_level <- ifelse(both_rnots$risk_level=="high_risk", "High Risk", "1.1")
both_rnots$risk_level <- factor(both_rnots$risk_level, levels = c("High Risk", "1.1"))

outbreak_plot <- ggplot(both_rnots, aes(time, Cum_Detections, group=interaction(.id, risk_level), color=risk_level)) + 
  geom_line(alpha=0.15,size=1) + 
  scale_color_brewer(palette="Set1", direction = 1)+
  coord_cartesian(xlim=c(0,150), ylim=c(0,50), expand=FALSE) +
  guides(color=guide_legend(override.aes=list(alpha=1))) +
  theme(legend.position=c(0.3,0.8))+
  labs(x = "Time (days)", 
       y = "Reported Autochthonous Cases", 
       color = expression("R"[0]))
print(outbreak_plot)



unknown_prev <- ldply(data.files, function(x) {
  load(x)
  disc_prob <- get_disc_prob_rand(x)
  risk_level <- get_risk_level_rand(x)
  prevalences <- get_prev_by_detects_all(rand_trials, f=totalprev_by_totaldetects)  
  
  prevalences <- ddply(prevalences, .(detected), .fun = function(x){ 
    quants <-  quantile(x = x$prevalence, probs = c(0.5, 0.25, 0.75), names=FALSE) 
    data.frame(median=quants[1], min = quants[2], max = quants[3])
  })
  cbind(risk_level, disc_prob, prevalences)
})
unknown_prev <- unknown_prev[which(unknown_prev$risk_level=="high_risk" ),]

known_prev <- get_prev_by_detects_plot(dir_path, r_nots = 1.1, disc_probs = c(0.011,0.0224), intro_rates = 0.01)
known_prev$intro_rate <- NULL
names(known_prev)[1] <- "risk_level"
known_prev$risk_level <- "1.1"

prev_plot_data <- rbind(unknown_prev, known_prev)

prev_plot_data$disc_prob <- paste0(calculate.discover(prev_plot_data$disc_prob), "%")
# prev_plot_data$disc_prob <- factor(prev_plot_data$disc_prob, levels = c("20%", "10%"))
prev_plot_data$risk_level <- ifelse(prev_plot_data$risk_level=="high_risk", "High Risk", "1.1")
prev_plot_data$risk_level <- factor(prev_plot_data$risk_level, levels = c("High Risk", "1.1"))


prev_plot <- ggplot(prev_plot_data, aes(detected, median, color=risk_level, fill=risk_level, linetype=as.factor(disc_prob), group = interaction(risk_level, disc_prob))) + 
  geom_line(size=1)+
  #geom_hline(yintercept=20)+
  geom_ribbon(aes(ymax=max, ymin=min), alpha=0.1, color=NA)+
  scale_y_log10(expand=c(0,0),limits=c(1,50), breaks = c(5,10,25,50))+
  coord_cartesian(xlim = c(0,15))+
  scale_x_continuous(expand=c(0.01,0.01))+
  theme(legend.position = c(0.3,0.79),
        #legend.direction = "horizontal",
        legend.box="horizontal")+
  scale_color_brewer(palette="Set1", direction = 1)+
  scale_fill_brewer(palette="Set1", direction=1)+
  guides(linetype=guide_legend(title.hjust = 0, override.aes=list("fill"=NA), title="Reporting Rate"),
         color=FALSE,fill=FALSE)+
  labs(x = "Reported Cases", 
       y = "Cases (log scale)", 
       color = expression("R"[0]), 
       fill = expression("R"[0]))
print(prev_plot)  
# prev_plot_data[prev_plot_data$detected==10,]

load("../data/rand_county_prob_data.Rdata")
# prob_data <- prob_data[which(prob_data$disc_prob==0.0224),]
prob_data <- prob_data[which(prob_data$risk_level=="high_risk"), ]
prob_data <- prob_data[which(prob_data$variable=="prob_epidemic"), ]
prob_data$variable <- NULL

load(get_vec_of_files(dir_path = dir_path, r_nots = 1.1, disc_probs = c(0.0224), intro_rates = 0.01))
temp_prev_below <- cbind("1.1", 0.0224, get_epidemic_prob_by_d(trials,20, 2000, 100, 100))
load(get_vec_of_files(dir_path = dir_path, r_nots = 1.1, disc_probs = c(0.011), intro_rates = 0.01))
temp_prev_below2 <- cbind("1.1", 0.011, get_epidemic_prob_by_d(trials,20, 2000, 100, 100))
names(temp_prev_below2) <- names(temp_prev_below)
temp_prev_below <- rbind(temp_prev_below, temp_prev_below2)
rm(temp_prev_below2)
names(temp_prev_below)[c(1:2, 4)] <- c("risk_level", "disc_prob", "value")
# temp_prev_below$prob_below <- 1 - temp_prev_below$prob_below
# temp_prev_below <- melt(temp_prev_below, measure.vars = c("prob_epidemic"))

prob_data <- rbind(prob_data, temp_prev_below)
prob_data$risk_level <- ifelse(prob_data$risk_level=="1.1", "1.1", "High Risk")
# prob_data$variable <- ifelse(prob_data$variable=="prob_below", "Prevalence", "Epidemic")
prob_data$risk_level <- factor(prob_data$risk_level, levels = c("High Risk", "1.1"))

prob_data$disc_prob <- paste0(calculate.discover(prob_data$disc_prob), "%")
# prob_data <- prob_data[seq(1,nrow(prob_data),by=3),]

prob_plot <- ggplot(prob_data, aes(detected, value, linetype=as.factor(disc_prob), color=risk_level)) + 
  geom_line(size=1) + 
  coord_cartesian(xlim=c(0,15), ylim=c(0,1), expand=FALSE)+
  geom_vline(xintercept=2, size=0.5, linetype=2)+
  scale_color_brewer(palette="Set1", direction = 1) +
  background_grid(major = "xy", minor = "none")+
  theme(legend.position="none",
        legend.box.just="left")+
  labs(x = "Reported Cases", 
       y = "Epidemic Probability", 
       color = "Trigger Type",
       linetype= "County Risk")
# print(prob_plot)
# fig2 <- plot_grid(outbreak_plot, prev_plot, prob_plot, nrow = 1, labels="AUTO")

fig2 <- ggdraw() + draw_plot(outbreak_plot, x = 0, y=0, width=.33, height=1)+
  draw_plot(plot = prev_plot, x = 0.33, y=0.0, width=0.33, height=1)+
  draw_plot(plot = prob_plot, x = 0.66, y=0.0, width=0.33, height=1)+
  draw_plot_label(c("A", "B", "C"), c(0, 0.33, 0.66), c(1, 1, 1), size = 20)

save_plot(paste0(fig_path, "figure3_local.pdf"), fig2, base_height = 4, base_aspect_ratio = 3)


################################################################################################
## End code for Figure 3.
################################################################################################

#######################################################
# 0.25, 0.5, 0.75, and 1.0 (in addition to the 2 you already included)
## Fig for Michael Johannson
# load("../data/rand_county_prob_data.Rdata")
# prob_data <- prob_data[which(prob_data$risk_level=="high_risk"), ]
# prob_data <- prob_data[which(prob_data$variable=="prob_epidemic"), ]
# prob_data$variable <- NULL
# 
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.2, disc_probs = c(0.0224), intro_rates = 0.001))
# temp_prev_below <- cbind("0.2", 0.0224, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.2, disc_probs = c(0.011), intro_rates = 0.001))
# temp_prev_below2 <- cbind("0.2", 0.011, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.2, disc_probs = c(0.0505), intro_rates = 0.001))
# temp_prev_below3 <- cbind("0.2", 0.0505, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# 
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.5, disc_probs = c(0.0224), intro_rates = 0.001))
# temp_prev_below4 <- cbind("0.5", 0.0224, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.5, disc_probs = c(0.011), intro_rates = 0.001))
# temp_prev_below5 <- cbind("0.5", 0.011, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.5, disc_probs = c(0.0505), intro_rates = 0.001))
# temp_prev_below6 <- cbind("0.5", 0.0505, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# 
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.7, disc_probs = c(0.0224), intro_rates = 0.001))
# temp_prev_below7 <- cbind("0.7", 0.0224, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.7, disc_probs = c(0.011), intro_rates = 0.001))
# temp_prev_below8 <- cbind("0.7", 0.011, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 0.7, disc_probs = c(0.0505), intro_rates = 0.001))
# temp_prev_below9 <- cbind("0.7", 0.0505, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# 
# load(get_vec_of_files(dir_path = dir_path, r_nots = 1, disc_probs = c(0.0224), intro_rates = 0.001))
# temp_prev_below10 <- cbind("1.0", 0.0224, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 1, disc_probs = c(0.011), intro_rates = 0.001))
# temp_prev_below11 <- cbind("1.0", 0.011, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# load(get_vec_of_files(dir_path = dir_path, r_nots = 1, disc_probs = c(0.0505), intro_rates = 0.001))
# temp_prev_below12 <- cbind("1.0", 0.0505, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
# 
# 
# names(temp_prev_below2) <-names(temp_prev_below3)<-names(temp_prev_below4)<-names(temp_prev_below5)<-names(temp_prev_below6)<-
#   names(temp_prev_below7)<-names(temp_prev_below8)<-names(temp_prev_below9)<-names(temp_prev_below10)<-names(temp_prev_below11)<- names(temp_prev_below12)<-names(temp_prev_below)
# 
# temp_prev_below <- rbind(temp_prev_below, temp_prev_below2, temp_prev_below3,temp_prev_below4,
#                          temp_prev_below5,temp_prev_below6,temp_prev_below7,temp_prev_below8,
#                          temp_prev_below9,temp_prev_below10,temp_prev_below11,temp_prev_below12)
# rm(temp_prev_below2)
# names(temp_prev_below)[c(1:2, 4)] <- c("risk_level", "disc_prob", "value")
# # temp_prev_below$prob_below <- 1 - temp_prev_below$prob_below
# # temp_prev_below <- melt(temp_prev_below, measure.vars = c("prob_epidemic"))
# 
# prob_data <- rbind(prob_data, temp_prev_below)
# 
# prob_data$risk_level <- ifelse(prob_data$risk_level=="high_risk", "High Risk", as.character(prob_data$risk_level))
# # prob_data$variable <- ifelse(prob_data$variable=="prob_below", "Prevalence", "Epidemic")
# # prob_data$risk_level <- factor(prob_data$risk_level, levels = c("High Risk", "1.1"))
# 
# # prob_data <- prob_data[seq(1,nrow(prob_data),by=3),]
# prob_data <- prob_data[which(prob_data$disc_prob!=0.0505),]
# prob_data$disc_prob <- paste0(calculate.discover(prob_data$disc_prob), "%")
# 
# mj_prob_below <- ggplot(prob_data, aes(detected, value, linetype=as.factor(disc_prob), color=risk_level)) + 
#   geom_line(size=1) + 
#   coord_cartesian(xlim=c(0,30), ylim=c(0,1), expand=FALSE)+
#   # geom_hline(yintercept=0.5, size=0.5, linetype=2)+
#   scale_color_brewer(palette="Set1", direction = 1) +
#   background_grid(major = "xy", minor = "none")+
#   theme(#legend.position="none",
#     legend.box.just="left")+
#   labs(x = "Cumulative Reported Cases", 
#        y = "Probability of Sustained Transmission", 
#        color = expression("R"[0]),
#        linetype= "Reporting\nRate")
# save_plot(filename = "../ExploratoryFigures/prob_below.pdf", mj_prob_below, base_height = 5, base_aspect_ratio = 1.3)

#################################################################




temp <- get_rows_final_sizes(final_sizes, r_nots, disc_probs, intro_rates)

final_plot <- ggplot(temp, aes(cumI)) + geom_histogram(bins=100)+
  facet_grid(intro_rate~r_not, scales = "free_y") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(expand=c(0,0))+
  labs(x = "Final Epidemic Sizes")
save_plot(filename = "../ExploratoryFigures/final_sizes.pdf", plot = final_plot, base_height = 10, base_width = 12)


#####################################################
#### Figure 1, understanding the overarching question!
#####################################################
get_symp <- function(new_intros, diedout=FALSE){
  # sum(sample(c("symp", "asymp"), n, prob = c(0.15, 0.85), replace = T) =="symp")
  df <- data.frame(new_intros=new_intros)
  df$symp <- 0
  df$asymp <- 0
  for(row in 1:nrow(df)){

    if(df$new_intros[row] !=0 ){
      # browser()
      samps <- sample(c("symp", "asymp"), df$new_intros[row], prob = c(0.15, 0.85), replace = T)
      new_symp <- sum(samps=="symp")
      new_asymp <- df$new_intros[row] - new_symp
      symp_time <- round(rnorm(1, mean = 8, sd = 2))
      asymp_time <- round(rnorm(1, mean = 8, sd = 2))
      while(row + symp_time > nrow(df)){
        symp_time <- symp_time-1
      }
      while(row + asymp_time > nrow(df)){
        asymp_time <- asymp_time-1
      }
      df$symp[row:(row+symp_time)] <- df$symp[row:(row+symp_time)] + new_symp
      df$asymp[row:(row+asymp_time)] <- df$asymp[row:(row+asymp_time)] + new_asymp


    }
  }

  if(diedout){
    df[nrow(df),] <- 0
  }

  df
}
load(get_vec_of_files(dir_path, 1.1, 0.011, 0.1))
set.seed(808)
ind <- which(all_last_cuminfect_values(trials) > 2000)
zoom_data <- trials[[ind[27]]]
zoom_data <- zoom_data[zoom_data$time <=110,]

symp_data <- get_symp(zoom_data$New_Infection, diedout = TRUE)

zoom_data$Symptomatic <- symp_data$symp
zoom_data$Asymptomatic <- symp_data$asymp
zoom_data$tot <- zoom_data$Symptomatic + zoom_data$Asymptomatic
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Asymptomatic", "Symptomatic"))
zoom_melt <- zoom_melt[order(zoom_melt$variable, decreasing = T),]

arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$tot[which(zoom_data$New_Intro!=0)]
det_arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Detections!=0)])
det_arrow_data$yval <- zoom_data$tot[which(zoom_data$New_Detections!=0)]

zoom_plot <- ggplot(zoom_melt, aes(x=time, y=value, color=NA, fill= variable)) + geom_bar(stat="identity", width=1.05)+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+1.5, yend=yval),
               arrow = arrow(length = unit(0.05, "npc"), angle = 35), color="red", size=1, inherit.aes=FALSE)+
  geom_vline(data=det_arrow_data, aes(xintercept=time), linetype=2, color="red")+
  theme_cowplot() %+replace% theme(legend.position="none")+
  scale_y_continuous(expand=c(0.0,0.0))+
  scale_x_continuous(expand=c(0.01,0.01),limits=c(0,101))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Prevalence", fill="")
print(zoom_plot)

# ex_growing_epidemic_data <- zoom_melt

ind <- which(all_last_cuminfect_values(trials) > 20 & all_last_cuminfect_values(trials) <30)
##  8 may work
zoom_data <- trials[[ind[1]]]
# zoom_data <- zoom_data[zoom_data$time <=75,]
symp_data <- get_symp(zoom_data$New_Infection, diedout = TRUE)

zoom_data$Symptomatic <- symp_data$symp
zoom_data$Asymptomatic <- symp_data$asymp
zoom_data$tot <- zoom_data$Symptomatic + zoom_data$Asymptomatic
zoom_melt <- melt(zoom_data, id.vars = c("time"), measure.vars = c("Asymptomatic", "Symptomatic"))
zoom_melt <- zoom_melt[order(zoom_melt$variable, decreasing = T),]

arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Intro!=0)])
arrow_data$yval <- zoom_data$tot[which(zoom_data$New_Intro!=0)]
det_arrow_data <- data.frame(time=zoom_data$time[which(zoom_data$New_Detections!=0)])
det_arrow_data$yval <- zoom_data$tot[which(zoom_data$New_Detections!=0)]

ylims <- ggplot_build(zoom_plot)$panel$ranges[[1]]$y.range
zoom_plot2 <- ggplot(zoom_melt, aes(x=time, y=value, color=NA, fill= variable)) + geom_bar(stat="identity", width=1.05)+
  geom_segment(data=arrow_data, aes(x=time, xend=time, y=yval+1.5, yend=yval),
               arrow = arrow(length = unit(0.05, "npc"), angle = 35), color="red", size=1, inherit.aes=FALSE)+
  geom_vline(data=det_arrow_data, aes(xintercept=time), linetype=2, color="red")+
  theme_cowplot() %+replace% theme(legend.position=c(0.9,0.85))+
  scale_y_continuous(expand=c(0.0,0.0), limits = c(0,ylims[2]))+
  scale_x_continuous(expand=c(0.01,0.01),limits=c(0,101))+
  scale_color_manual(values=c("black", "grey"), guide=FALSE) +
  scale_fill_manual(values=c("black", "grey")) +
  labs(x = "Time (days)", y="Prevalence", fill="")
print(zoom_plot2)

# ex_dying_epidemic_data <- zoom_melt

all_ts <- ggdraw() +
  draw_plot(zoom_plot2, x =  0,y =  0.5,width =  1,height =  0.5) +
  draw_plot(zoom_plot, x =  0,y =  0,width =  1, height = .5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = 15)

save_plot(filename = "../ExploratoryFigures/all_ts2.pdf", plot = all_ts, base_height=4, base_aspect_ratio = 1.8)

#########################################################################

# ## Supplementary choosing threshold plot -- should be 20
r_nots <- c(seq(0.9, 1.1, by=0.05))
disc_probs <- c(0.011)
intros <- c(0.01, 0.1)
threshold_plot <- plot_dots(dir_path, r_nots, disc_probs, intros, local=T)
# print(threshold_plot)
save_plot(paste0(fig_path, "local_threshold_plot.pdf"), threshold_plot, base_height = 4, base_aspect_ratio = 1.8)






#################################################################
## TESTING EPIDEMIC THRESHOLD IDEAS
#################################################################

is_epidemic <- function(trial){
  ## Function should only run if the run hit cumulative threshold and
  ## the prevalence threshold
  if(last_cuminfect_value(trial) >= 1000 & max_prevalence(trial)>20){
    times <- 1:nrow(trial)

    rows <- split(times, cut(times, breaks=3))

    mid_chunk <- rows[[2]]
    final_chunk <- rows[[3]]

    mid_data <-trial[mid_chunk, c("Total_Infections")]
    final_data <-trial[final_chunk, c("Total_Infections")]

    mod <- t.test(mid_data, final_data)
    data.frame(statistic=mod$statistic, p_val = mod$p.value)
  } else{
    return(data.frame(statistic = NA, p_val= NA))
  }
}
all_is_epidemic <- function(trials, ...){
  ## Returns all the betas for the 10000 simulations
  ldply(trials, is_epidemic, ...)
}
thirds <- all_is_epidemic(trials)
sum(thirds$p_val<0.05/10000, na.rm=T)/10000
thirds_ind <- which(thirds$p_val<0.05/10000)

plot_trial_data(trials[thirds_ind], n = 5)


get_beta <- function(trial, num_rows=50){
  ## Function should only run if the run hit cumulative threshold and
  ## the prevalence threshold
  if(last_cuminfect_value(trial) >= 1000 & max_prevalence(trial)>20){
    n <- nrow(trial)
    rows_used <- (n-num_rows):n
    rows_used <- rows_used[which(rows_used>0)]
    data <-trial[rows_used, c("time","Total_Infections")]
    data$logged <- ifelse(data$Total_Infections==0, NA, log(data$Total_Infections))
    model <- lm(Total_Infections~ time, data = data)

    summary(model)$coefficients["time",]
  } else{
    return(c(Estimate = NA, "Std. Error"= NA,  "t value"= NA, "Pr(>|t|)"=NA))
  }
}

all_get_beta <- function(trials, ...){
  ## Returns all the betas for the 10000 simulations
  ldply(trials, get_beta, ...)
}
load(get_vec_of_files(dir_path, 1.03, 0.011, 0.1))

ind <- which(all_last_cuminfect_values(trials)>=2000 & all_max_prevalence(trials) >=50)

plot_trial_data(trials[-ind], n = 5)



test <- all_get_beta(trials, num_rows=200)
no_epis <- sum(is.na(test$Estimate))
ind <- which(test$`Pr(>|t|)`< (0.05/(10000-no_epis)) & test$Estimate>0)
length(ind)/10000
ind <- which(all_last_cuminfect_values(trials)>=1000 & all_max_prevalence(trials)>=20)
length(ind)/10000


plot_trial_data(trials[ind], n = 5)


plot(all_max_prevalence(trials), test$Estimate)

plot(all_last_cuminfect_intro_values(trials),all_last_cuminfect_local_values(trials))

temp <- all_last_cuminfect_intro_values(trials)/all_last_cuminfect_local_values(trials)
plot(temp, all_max_prevalence(trials), pch=20)


pos_beta <- which(test$Estimate>0)
blacks <- rep("red", 10000)
blacks[pos_beta] <- "black"
plot(test$Estimate, all_max_prevalence(trials),  col=blacks,pch=20)

## Plot 3 non epidemics with 3 epidemics
load(get_vec_of_files(dir_path, 1.05, 0.011, 0.1))
samples <- sample(x = 1:length(trials), size = 10, replace = F)  
df <- do.call("rbind", trials[samples]) ## Combine into data frame
df$r_not <- "1.1" ## Add index for each simulation
df$index <- rep(seq_along(trials[samples]), sapply(trials[samples], nrow))
load(get_vec_of_files(dir_path, .95, 0.011, 0.1))
samples <- sample(x = 1:length(trials), size = 10, replace = F)  
df <- rbind(df, cbind(do.call("rbind", trials[samples]),r_not="0.9", index=rep(seq_along(trials[samples])+10, sapply(trials[samples], nrow))))

df <- melt(df, id.vars=c("time", "index", "r_not"), measure.vars = c("Total_Infections"))  
ggplot(df, aes(time, value, color=as.factor(r_not), group=as.factor(index))) + geom_line() + guides(color=FALSE)

####################################################




###############################
## Risk tolerance plot
r_nots <- c(1.1)#, 1.3)
intros <- c(0.01, .1)
det_probs <- c(0.011, 0.0224)
conf <-seq(0.05, .95, by=0.05)
risk_tol <- get_trigger_data(r_nots, intros, det_probs, threshold=20, confidence=conf, num_necessary = 10)
risk_tol$risk_tolerance <- 1-risk_tol$confidence
risk_tol$disc_prob <- paste0(calculate.discover(risk_tol$disc_prob), "%")
risk_tol_plot <- ggplot(risk_tol, aes(risk_tolerance, prev_trigger, color=as.factor(intro_rate), linetype=as.factor(disc_prob))) + 
  geom_line(size=1) + 
  scale_color_manual(values=c("Grey", "Black"))+
  coord_cartesian(ylim=c(0, 150), xlim = c(0,1), expand = FALSE)+
  labs(x = "Risk Tolerance", 
       y = "Threshold (Nowcasting)", 
       color = "Importation\nRate", 
       linetype= "Detection \nProbability")

print(risk_tol_plot)

