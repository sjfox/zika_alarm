rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

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


run_simple_branches <- function(num_reps, ...){
  rlply(.n = num_reps, .expr = run_branch_simple(...) )
}

getMaxCumI <- function(x){
  return(x$cumI[nrow(x)])
}
all_getMaxCumI <- function(x) {
  laply(x, getMaxCumI)
}

trials <- run_simple_branches(1000, prop_p=1.7/7, recov_p=1/7, disc_p=.01, d_thresh=5, e_thresh=200)
final.sizes <- all_getMaxCumI(trials)
hist(final.sizes, breaks=100)

trials <- run_branches(num_reps = 1000, branch_params(prop_p=.5/7, dis_prob_symp = .01))

p <- plot_final_sizes(trials)
print(p)

ggsave(filename = "../ExploratoryFigures/r0_0.2_disc_0.01_hist.pdf", plot = p, width=6, height=5)
