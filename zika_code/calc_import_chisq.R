rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

library(tidyverse)
library(cowplot)
library(stringr)

cty_preds <- read_csv("../csvs/county_master.csv")

cty_data <- read_csv("../csvs/texas_county_imports.csv")
cty_data <- cty_data %>% mutate(imports = ifelse(is.na(local), cases, cases-local)) %>%
              rename(Geography=county)

county_data <- cty_preds %>% mutate(Geography = str_replace_all(Geography, pattern = " County, Texas", "")) %>%
  left_join(cty_data) %>% 
  select(Geography, importation_probability, imports) %>%
  mutate(imports = if_else(is.na(imports), as.integer(0), imports),
         pred_imports = importation_probability * sum(imports)) 

county_data %>% ggplot(aes(imports, pred_imports)) + geom_point()+
  geom_abline(intercept= 0 , slope=1)

lm(imports ~ pred_imports, data=county_data) %>% summary()

tes <- county_data %>% filter(pred_imports>=1)
chisq.test(tes$imports, p = tes$importation_probability/sum(tes$importation_probability))

