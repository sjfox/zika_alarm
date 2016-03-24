rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

sapply(c('branch.functions.R','plot.functions.R'), source)
library(plyr)
library(cowplot)
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
branch_params <- function(prop_p = 1.7/7 , 
                          recov_p = 1.0/7,
                          d_thres = 5,
                          e_thresh = 300,
                          prob_symp = 1,
                          incub_p = 1,
                          dis_prob_symp = 1,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000)
  return(as.list(environment()))

# Average Run 
trials <- run_branches(num_reps = 1000, branch_params())
lastvalues <- all_last_cuminfect_values(trials)
hist(lastvalues, breaks = 100)
hist(lastvalues, breaks = 100, ylim = c(0,100))
sort(lastvalues)

# R0 = 1 Run 
trials2 <- run_branches(num_reps = 1000, branch_params(prop_p=.5/7))
lastvalues2 <- all_last_cuminfect_values(trials2)
hist(lastvalues2, breaks = 100)
hist(lastvalues2, breaks = 100, ylim = c(0,100))
sort(lastvalues2)


count_escapes(trials, 1, branch_params()$e_thresh)

plot_final_sizes(trials)


escape_prob <- prob_ext(prop_p = prop_p,recov_p = recov_p, incub_p = incub_p, prob_symp = prob_symp, 
         d_thres=d_thres,e_thresh = e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100)
escape_prob

escape_underDoverE <- prob_underD.overE(prop_p = prop_p,recov_p = recov_p, incub_p = incub_p, prob_symp = prob_symp, 
                        d_thres,e_thresh = e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100)


# Given we have detected X cases- what is the average number of Total Infected Cases and Cumulative Cases 
#at.detect <- all_detect_rows(trials)


## Plot grid of histograms, not really of interest anymore
# dis_prob_symps <- seq(0.01, 1, length.out = 5)
# prop_ps <- seq(1.1/7, 4/7, length.out = 5)
# all_runs <- data.frame()
# for(dis_prob_symp in dis_prob_symps){
#   for(prop_p in prop_ps){
#     trials <- run_branches(num_reps = 1000, prop_p, recov_p, incub_p, prob_symp, 
#                            d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate)
#     all_runs <- rbind(all_runs, data.frame(dis_prob_symp=dis_prob_symp, prop_p=prop_p, final_size = all_last_cuminfect_values(trials)))
#   }
# }
# 
# grid_hist <- ggplot(all_runs, aes(final_size)) + geom_histogram(binwidth = 20) + 
#   facet_grid(dis_prob_symp~prop_p, labeller=label_both, scales = "free_y")+
#   scale_y_continuous(expand=c(0,0))+
#   scale_x_continuous(expand=c(0.01,0.01))
# 
# ggsave(filename = "../ExploratoryFigures/Final_size_hist.pdf",plot = grid_hist, width=16, height=10)
