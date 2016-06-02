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
trigger_dir_path <- "~/projects/zika_alarm/data/triggers50/"
save_path <- "~/projects/zika_alarm/data/"
fig_path <- "~/projects/zika_alarm/ExploratoryFigures/"

######### Combine triggers after a tacc run/download #####################
# combine_triggers(trigger_dir_path, save_path)


################################
## Code to Make Figure 2
################################
##### Panel A


get_epi_data <- function(trials, n){
  ## Returns first n trials in data frame form
  names(trials) <- seq_along(trials)
  ldply(trials[1:n], data.frame)
}

data.files <- list.files(path="../data/rand_trials", pattern="*.Rdata", full.names=T, recursive=FALSE)
load(data.files[4])
plot_final_sizes(rand_trials)
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
known_rnot <- cbind(data.frame(risk_level="1.1", disc_prob=0.0224), get_epi_data(trials, 2000))
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
       y = "Cumulative Reported Cases", 
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

prev_plot_data <- prev_plot_data[seq(1,nrow(prev_plot_data), by=2),]
prev_plot <- ggplot(prev_plot_data, aes(detected, median, color=risk_level, fill=risk_level, linetype=as.factor(disc_prob), group = interaction(risk_level, disc_prob))) + 
  geom_line(size=1)+
  #geom_hline(yintercept=20)+
  geom_ribbon(aes(ymax=max, ymin=min), alpha=0.1, color=NA)+
  scale_y_log10(expand=c(0,0),limits=c(1,200), breaks = c(5,10,25,50,100))+
  coord_cartesian(xlim = c(0,30))+
  scale_x_continuous(expand=c(0.01,0.01))+
  theme(legend.position = c(0.3,0.79),
        #legend.direction = "horizontal",
        legend.box="horizontal")+
  scale_color_brewer(palette="Set1", direction = 1)+
  scale_fill_brewer(palette="Set1", direction=1)+
  guides(linetype=guide_legend(title.hjust = 0, override.aes=list("fill"=NA), title="Reporting Rate"),
         color=FALSE,fill=FALSE)+
  labs(x = "Cumulative Reported Cases", 
       y = "Prevalence (log scale)", 
       color = expression("R"[0]), 
       fill = expression("R"[0]))
print(prev_plot)  


load("../data/rand_county_prob_data.Rdata")
# prob_data <- prob_data[which(prob_data$disc_prob==0.0224),]
prob_data <- prob_data[which(prob_data$risk_level=="high_risk"), ]
prob_data <- prob_data[which(prob_data$variable=="prob_epidemic"), ]
prob_data$variable <- NULL

load(get_vec_of_files(dir_path = dir_path, r_nots = 1.1, disc_probs = c(0.0224), intro_rates = 0.01))
temp_prev_below <- cbind("1.1", 0.0224, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
load(get_vec_of_files(dir_path = dir_path, r_nots = 1.1, disc_probs = c(0.011), intro_rates = 0.01))
temp_prev_below2 <- cbind("1.1", 0.011, get_epidemic_prob_by_d(trials,50, 2000, 200, 100))
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
  coord_cartesian(xlim=c(0,30), ylim=c(0,1), expand=FALSE)+
  geom_hline(yintercept=0.5, size=0.5)+
  scale_color_brewer(palette="Set1", direction = 1) +
  background_grid(major = "xy", minor = "none")+
  theme(legend.position="none",
        legend.box.just="left")+
  labs(x = "Cumulative Reported Cases", 
       y = "Threshold Probability", 
       color = "Trigger Type",
       linetype= "County Risk")
# print(prob_plot)
# fig2 <- plot_grid(outbreak_plot, prev_plot, prob_plot, nrow = 1, labels="AUTO")

fig2 <- ggdraw() + draw_plot(outbreak_plot, x = 0, y=0, width=.33, height=1)+
  draw_plot(plot = prev_plot, x = 0.33, y=0.0, width=0.33, height=1)+
  draw_plot(plot = prob_plot, x = 0.66, y=0.0, width=0.33, height=1)+
  draw_plot_label(c("A", "B", "C"), c(0, 0.33, 0.66), c(1, 1, 1), size = 20)

save_plot(paste0(fig_path, "figure3_new.pdf"), fig2, base_height = 4, base_aspect_ratio = 3)





##### Panel B
thresholds <- c(20)
r_nots <- c(1.1)
disc_probs <- c(.0224)
intro_rates <- c(.01, .1)
prob_below <- get_prob_below_plot(dir_path, thresholds, r_nots, disc_probs, intro_rates)
prob_below$prob_below <- 1-prob_below$prob_below
epi_prob <- get_epidemic_prob_plot(dir_path, prev_threshold = 50, cum_threshold = 2000, r_nots, disc_probs, intro_rates)
both_combined <- merge(x = prob_below, y= epi_prob, by.x= c("r_not", "intro_rate", "disc_prob", "detected"), by.y=c("r_not", "intro_rate", "disc_prob", "detected"))
both_combined$threshold <- NULL
df <- melt(both_combined, measure.vars = c("prob_below", "prob_epidemic"))
df$variable <- ifelse(df$variable=="prob_below", "Prevalence", "Epidemic")
df$intro_rate <- factor(df$intro_rate, levels=c(0.1,0.01))
df <- df[seq(1,nrow(df),by=3),]
epi_prev_threshold <- ggplot(df, aes(detected, value, linetype=intro_rate, color = variable)) + 
  geom_line(size=1) +
  geom_hline(yintercept=0.5, size=0.5)+
  coord_cartesian(ylim=c(0,1), xlim = c(0,100), expand = FALSE)+
  scale_color_manual(values=c("Black", "Grey")) +
  background_grid(major = "xy", minor = "none")+
  theme(legend.position="none")+
  labs(x = "Cumulative Reported Cases", 
       y = "Trigger Probability", 
       color = "Trigger Type",
       linetype= "Importation\nRate")
# print(epi_prev_threshold)

# save_plot(paste0(fig_path, "trigger_probability.pdf"),epi_prev_threshold, base_height = 5, base_aspect_ratio = 1.3)

r_nots <- c(0.8, 0.85, seq(0.9, 1.2, by=0.01), 1.25, seq(1.3, 2, by=0.1))
intros <- c(0.01, .1)
det_probs <- c(0.0224)

triggers <- get_trigger_data(r_nots, intros, det_probs, confidence=0.5, num_necessary = 100)


## Remove NA to make plots drift off to top of figure
triggers$prev_trigger[which(is.na(triggers$prev_trigger))] <- 200
triggers$epi_trigger[which(is.na(triggers$epi_trigger))] <- 200
temp <- melt(triggers, measure.vars = c("epi_trigger", "prev_trigger"))
temp$variable <- ifelse(temp$variable=="prev_trigger", "Prevalence", "Epidemic")
trigger_plot <- ggplot(temp, aes(r_not, value, linetype=as.factor(intro_rate), color=variable)) + 
  geom_line(size=1) + 
  scale_color_manual(values=c("Grey", "Black"))+
  coord_cartesian(ylim=c(0, 150), xlim = c(0.95,1.4), expand = FALSE)+
  scale_x_continuous(expand=c(0.01,0.01))+
  theme(legend.position=c(0.75,0.7), legend.box.just="left")+
  labs(x = expression("R"[0]), 
       y = "Trigger (Reported Cases)", 
       linetype = "Importation\nRate", 
       color = "Trigger Type")
print(trigger_plot)
# save_plot(paste0(fig_path, "triggers.pdf"), trigger_plot, base_height = 5, base_aspect_ratio = 1.3)

## Align figures
p1 <-ggplot_gtable(ggplot_build(prev_plot))
p2 <- ggplot_gtable(ggplot_build(trigger_plot))
p2$heights <- p1$heights
figure2 <- plot_grid(p1, epi_prev_threshold, p2, labels=c("A", "B", "C"), nrow = 1, rel_widths = c(1, 1, 1))
save_plot(paste0(fig_path, "figure2_panels.pdf"), figure2, base_height = 4, base_aspect_ratio = 3)

# 
# prev_trigger_plot <- ggplot(prev_triggers, aes(r_not, prev_trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) + 
#   geom_line(size=1) + 
#   scale_color_manual(values=c("Grey", "Black"))+
#   coord_cartesian(ylim=c(0, 150), xlim = c(0.9,1.4), expand = FALSE)+
#   theme_cowplot()%+replace% theme(legend.position="none")+
#   labs(x = expression("R"[0]), 
#        y = "Trigger (Nowcasting)", 
#        color = "Importation\nRate", 
#        linetype= "Reporting \nRate")
# print(prev_trigger_plot)
# ######## Panel C
# ### Epidemic probability plot
# r_nots <- c(0.8, 0.85, seq(0.9, 1.2, by=0.01), 1.25, seq(1.3, 2, by=0.1))
# intros <- c(0.01, .1)
# det_probs <- c(0.011, 0.0224)
# 
# epi_triggers  <- get_trigger_data(r_nots, intros, det_probs, confidence=0.5, num_necessary = 100)
# 
# # epi_triggers$epi_trigger <- ifelse(is.na(epi_triggers$epi_trigger), Inf, epi_triggers$epi_trigger)
# epi_triggers$disc_prob <- paste0(calculate.discover(epi_triggers$disc_prob), "%")
# 
# ## Remove outlier from low trigger data ( few runs, and so not infinity solely due to stochasticity)
# epi_triggers$epi_trigger[which(is.na(epi_triggers$epi_trigger))] <- 200
# 
# epi_prob_plot <- ggplot(epi_triggers, aes(r_not, epi_trigger, linetype=as.factor(disc_prob), color=as.factor(intro_rate))) +
#   geom_line(size=1) + scale_color_manual(values=c("Grey", "Black"))+
#   coord_cartesian(ylim=c(0, 150), xlim = c(1,1.25), expand = FALSE)+
#   labs(x = expression("R"[0]),
#        y = "Trigger (Forecasting)",
#        color = "Importation\nRate",
#        linetype= "Reporting\nRate")
# print(epi_prob_plot)

# ################
# # Print out the panels all together and gridded
# epi_legend <- get_legend(epi_prob_plot)
# epi_prob_plot <- epi_prob_plot + theme(legend.position="none")
# prev_legend <- get_legend(prev_plot)
# prev_plot <- prev_plot + theme(legend.position="none")
# 
# p1 <-ggplot_gtable(ggplot_build(prev_plot))
# p2 <- ggplot_gtable(ggplot_build(epi_prob_plot))
# max_height <- grid::unit.pmax(p1$heights, p2$heights)
# p1$heights <- max_height
# p2$hights <- max_height
# plot_grid(p1, p2)
# 
# figure2_panels <- plot_grid(p1, prev_trigger_plot, p2, labels=c("A", "B", "C"), nrow = 1, rel_widths = c(1, 1, 1))
# figure2_panels <- ggdraw()+ draw_plot(figure2_panels, x =  0,y =  0,width =  .95,height =  1) +
#   draw_plot(epi_legend, x =  0.81, y =  0.2,width =  .2, height = 1) +
#   draw_plot(prev_legend, x =  0.79, y =  -0.125, width =  .2, height = 1)
# 
# save_plot(paste0(fig_path, "figure2_panels.pdf"), figure2_panels, base_height = 4, base_aspect_ratio = 3)





# ######### Supplemental figures
### Porbability below threshold graph
thresholds <- c(20)
r_nots <- c(0.9, 1.1, 1.3)
disc_probs <- c(0.011,.0224)
intro_rates <- c(.05, 0.3)
temp <- get_prob_below_plot(dir_path, thresholds, r_nots, disc_probs, intro_rates)

prob_below <- plot_prob_below(temp)
save_plot(paste0(fig_path, "probability_below_supplemental.pdf"),prob_below, base_height = 5, base_aspect_ratio = 2)

### Epidemic probability graph
r_nots <- c(1.1, 1.3)
disc_probs <- c(0.011, 0.0224)
intro_rates <- c(.05, 0.1)
temp <- get_epidemic_prob_plot(dir_path, prev_threshold = 25, cum_threshold = 1000, r_nots, disc_probs, intro_rates)

epi_plot <- plot_epidemic_prob(temp)
save_plot(paste0(fig_path, "epi_prob_supplemental.pdf"),epi_plot, base_height = 5, base_aspect_ratio = 1.7)



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

zoom_data <- trials[[18]]
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


ind <- which(all_last_cuminfect_values(trials) > 20 & all_last_cuminfect_values(trials) <30)
##  8 may work
zoom_data <- trials[[ind[7]]]
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


all_ts <- ggdraw() +
  draw_plot(zoom_plot2, x =  0,y =  0.5,width =  1,height =  0.5) +
  draw_plot(zoom_plot, x =  0,y =  0,width =  1, height = .5) +
  draw_plot_label(c("A", "B"), c(0, 0), c(1, 0.5), size = 15)

save_plot(filename = "../ExploratoryFigures/all_ts2.pdf", plot = all_ts, base_height=4, base_aspect_ratio = 1.8)

# ## Supplementary choosing threshold plot -- should be 20
r_nots <- c(seq(0.95, 1.1, by=0.01))
disc_probs <- c(0.011)
intros <- c(0.01, 0.1, 0.3)
threshold_plot <- plot_dots(dir_path, r_nots, disc_probs, intros)
threshold_plot <- threshold_plot + geom_hline(yintercept=20, color="blue")
# print(threshold_plot)
save_plot(paste0(fig_path, "threshold_plot.pdf"), threshold_plot, base_height = 4, base_aspect_ratio = 1.8)






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

