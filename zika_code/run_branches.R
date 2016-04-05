rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

library(plyr)
args <- (commandArgs(TRUE)) ## load arguments from R CMD BATCH

## Needs 
if(length(args)>0)  { ## Then cycle through each element of the list and evaluate the expressions.
  print(paste0('loading in ', args, ' from R CMD BATCH'))
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

source("incubation_branch.R")

branch_params <- function(r_not = 1.1,
                          infBoxes = 6,
                          incBoxes = 6,
                          recov_p = 0.6071267,
                          incub_rate = 0.5777811,
                          prop_p =  r_not*recov_p/infBoxes, 
                          d_thres = 5,
                          e_thresh = 500,
                          prob_symp = 1,
                          dis_prob_symp = .01,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000)
return(as.list(environment()))

params <- branch_params(r_not=desired_Rnot, dis_prob_symp = disc_prob, intro_rate = intro_rate )

trials <- run_branches_inc(num_reps = 10000, params)

save(list = c("trials", "params"), file = paste0("../../workfolder/data/zika/", "zika_sims_", desired_Rnot, "_", disc_prob, "_", intro_rate, ".Rdata"))
