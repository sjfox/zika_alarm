rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

run_type <- "importations"

## Want total discovery rates of 5%, 10%, 15%, 20%, 40% 
## Calculated by total discovery probability = (1-(1-daily_prob)^9.88)
disc_probs <- c(0.0052, 0.011, 0.01635, 0.0224, 0.0505) 
if(run_type== "importations"){
  imports <- read.csv("../csvs/county_plot.csv")
  current_combos <- as.character(unique(interaction(imports$R0_round, imports$importation.current,sep = "_")))
  projected_combos <- as.character(unique(interaction(imports$R0_round, imports$importation.projected,sep = "_")))
  projected_worst_combos <- as.character(unique(interaction(imports$R0_round, imports$importation.worse.projected,sep = "_")))
  
  all_combos <- unique(c(current_combos, projected_combos, projected_worst_combos))
  
  values <- unlist(strsplit(all_combos, "_"))
  ## R0s are the first item in each combo
  r_nots <- values[seq(1,length(values), by=2)]
  intro_rates <- values[seq(2,length(values), by=2)]
  
  sink('../launcher/run_import_branches.txt')
  for(disc_prob in disc_probs){
    for(ii in 1: length(r_nots)){
      startCmd <- "R CMD BATCH '--no-restore --no-save --args"
      paramCmd <- paste0(' desired_Rnot=', r_nots[ii], ' disc_prob=', disc_prob, ' intro_rate=', intro_rates[ii], " '")
      endCmd <- " ../zika_code/run_branches.R"
      full_cmd <- paste0(startCmd, paramCmd, endCmd)
      cat(full_cmd)              
      cat('\n')              
    }
  }
  sink()  
  
}else{
  r_nots <- c(0.8, 0.85, seq(0.9, 1.2, by=0.01), 1.25, seq(1.3, 2, by=0.1))
  intro_rates <- c(0.0, 0.01, 0.05, 0.1, 0.3)
  
  sink('../launcher/run_branches.txt')
  for(r_not in r_nots){
    for(disc_prob in disc_probs){
      for(intro_rate in intro_rates){
        startCmd <- "R CMD BATCH '--no-restore --no-save --args"
        paramCmd <- paste0(' desired_Rnot=', r_not, ' disc_prob=', disc_prob, ' intro_rate=', intro_rate, " '")
        endCmd <- " ../zika_code/run_branches.R"
        full_cmd <- paste0(startCmd, paramCmd, endCmd)
        cat(full_cmd)              
        cat('\n')              
      }
    }
  }
  sink()  
}
