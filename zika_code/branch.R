rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

source("branch.functions.R")
library(plyr)
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
prop_p <- 1.7/7  
recov_p <- 1.0/7
d_thres <- 25
e_thresh <- 300
prob_symp <- 1
incub_p <- 1
dis_prob_symp <- 0.01
dis_prob_asymp <- 0.00 
intro_rate <- .000


trials <- run_branches(num_reps = 1000, prop_p, recov_p, incub_p, prob_symp, 
                       d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate)


escape_prob <- prob_ext(prop_p = prop_p,recov_p = recov_p, incub_p = incub_p, prob_symp = prob_symp, 
         d_thres,e_thresh = e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=1000)

escape_underDoverE <- prob_underD.overE(prop_p = prop_p,recov_p = recov_p, incub_p = incub_p, prob_symp = prob_symp, 
                        d_thres,e_thresh = e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100)


# Given we have detected X cases- what is the average number of Total Infected Cases and Cumulative Cases 
at.detect <- all_detect_rows(trials)

