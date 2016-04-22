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

load(get_vec_of_files(dir_path, 1.2, .068, 0.1))


# plot_final_sizes(trials)
# plot_intro_final_sizes(trials)
# plot_local_final_sizes(trials)
# plot_max_nonintro_prevalences(trials)
# plot_max_prevalences(trials)

get_prev_by_detects_plot <- function(dir_path, r_nots, disc_probs, intro_rates){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    prevalences <- get_prev_by_detects_all(trials, f=totalprev_by_totaldetects)  
    
    prevalences <- ddply(prevalences, .(detected), .fun = function(x){ 
                      quants <-  quantile(x = x$prevalence, probs = c(0.5, 0.25, 0.75), names=FALSE) 
                      data.frame(median=quants[1], min = quants[2], max = quants[3])
                    })
    cbind(as.data.frame(parms), prevalences)
  })  
}

r_nots <- c(0.9, 1.3)
disc_probs <- c(0.011, 0.0287)
intros <- 0.1
prev_plot_data <- get_prev_by_detects_plot(dir_path, r_nots = r_nots, disc_probs = disc_probs, intro_rates = intros)
prev_plot <- plot_prevalences(prev_plot_data)
# print(prev_plot)  















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


# example <- get_vec_of_files(dir_path, 1.1, 0.011, .1)
# load(example)
# plotDat <- trials[[3]]
# plotDat <- melt(plotDat, id.vars = c("time"), measure.vars = c("New_Detections","New_Infectious"))
# plotDat <- plotDat[plotDat$time<100, ]
# plotDat$variable <- ifelse(plotDat$variable=="New_Detections", "Newly Detected", "Newly Infectious")
# ts_plot <- ggplot(plotDat, aes(time, y=value, fill=variable)) + geom_bar(stat="identity", position="dodge") +
#   theme_cowplot() %+replace% theme(legend.position=c(0.2,0.8))+
#   scale_y_continuous(expand=c(0.0,0.0))+
#   scale_x_continuous(expand=c(0.03,0.05))+
#   scale_fill_manual(values=c("black", "grey")) + 
#   scale_color_manual(values=c("black", "grey")) + 
#   labs(x = "Time (days)", y="Individuals", color = "", fill="") 
# print(ts_plot)
# save_plot(filename = "../ExploratoryFigures/ex_ts_bars.pdf", plot = ts_plot, base_height=2, base_aspect_ratio = 2.8)
# 
# plotDat <- trials[[3]]
# plotDat <- melt(plotDat, id.vars = c("time"), measure.vars = c("Cum_Detects", "Cumulative_Infections"))
# plotDat <- plotDat[plotDat$time<100, ]
# plotDat$min_val <- ifelse(plotDat$variable=="Cum_Detects", 0, NA)
# plotDat$min_val[which(is.na(plotDat$min_val))] <- plotDat$value[which(is.na(plotDat$min_val))-nrow(plotDat)/2]
# plotDat$variable <- ifelse(plotDat$variable=="Cum_Detects", "Detected Cases", "Undetected Cases")
# ts_plot <- ggplot(plotDat, aes(time, ymax=value,ymin=min_val,  fill=variable)) + geom_ribbon() +
#   theme_cowplot() %+replace% theme(legend.position=c(0.2,0.8))+
#   scale_y_continuous(expand=c(0.01,0.01))+
#   scale_x_continuous(expand=c(0.03,0.05))+
#   scale_fill_manual(values=c("black", "grey")) + 
#   labs(x = "Time (days)", y="Individuals", color = "", fill="") 
# print(ts_plot)
# 
# save_plot(filename = "../ExploratoryFigures/ex_ts.pdf", plot = ts_plot, base_height=2, base_aspect_ratio = 2.8)
