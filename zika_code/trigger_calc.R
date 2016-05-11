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

sapply(c('branch.functions.R','plot.functions.R', 'incubation_branch.R', 'analyze_saved_sims.R'), source)

prevalence <- 20
epi_prevalence <- 50
confidences <- seq(0.05, 0.95, by=0.05)
num_necessary <- c(10, 20, 100)

load(data.file)
params <- get_parms(data.file)
triggers <- data.frame()
for(num in num_necessary){
  for(conf in confidences){
    epi <- get_epidemic_trigger(trials = trials, threshold = epi_prevalence, confidence = conf, max_detect=300, num_necessary=num)  
    prev <- get_surveillance_trigger(trials = trials,threshold =  prevalence, confidence = conf, max_detect=300, num_necessary=num) 
    triggers <- rbind(triggers, cbind(params, data.frame(prev_threshold= prevalence, epi_threshold=epi_prevalence, confidence=conf, num_necessary=num, epi_trigger=epi, prev_trigger=prev)))
  }
}

save(list = c("triggers"), file = paste0("../../workfolder/data/zika_triggers_all50/", "zika_triggers_", params$r_not, "_", params$disc_prob, "_", params$intro_rate, ".Rdata"))

