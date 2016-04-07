rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

r_nots <- c(seq(0.8, 2, by=0.05), seq(2.2, 3, by=.5))

## Want total discovery rates of 5%, 10%, 50%, 90%
disc_probs <- c(0.0052, 0.011, 0.068, 0.21) ## Change

intro_rates <- c(0.0, .01, 0.3, 1, 1.92)

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