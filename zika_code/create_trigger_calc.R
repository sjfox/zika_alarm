rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


dirPath <- "../../workfolder/data/zika_sensitivity/"

data.files <- list.files(path=dirPath, pattern="*.Rdata", full.names=T, recursive=FALSE)
sink('../launcher/trigger_calc.txt')
for(file in data.files){
  startCmd <- "R CMD BATCH '--no-restore --no-save --args"
  fileCmd <- paste0(' data.file="', file, '"')
  endCmd <- "' ../zika_code/trigger_calc.R"
  full_cmd <- paste0(startCmd, fileCmd, endCmd)
  # print(full_cmd)
  cat(full_cmd)               # add command
  cat('\n')              # add new line
}
sink()
