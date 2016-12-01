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

rnots <- read_csv(file = "../csvs/county_sensitive.csv")
rnots <- rnots%>%mutate(Geography = str_replace_all(Geography, pattern = " County, Texas", ""))

rnots_ranked <- rnots %>% mutate_at(.funs = funs(rank), .cols = vars(rnott.expected:August), ties.method = "average")

order_x_axis <- c("rnott.expected", "rnott.low", "rnott.high", "rnott.high.slope", "increase.cr", 
                  "decrease.cr", "decrease.b", "increase.b", "increase.alpha", "decrease.alpha", 
                  "July",  "August", "September",  "October",  "November")


rank_plot <- rnots_ranked %>% arrange(desc(rnott.expected)) %>%
  mutate(Geography = factor(Geography, levels = Geography)) %>%
  slice(1:15) %>%
  gather(key, value, rnott.expected:August) %>%
  mutate(key = factor(key, levels = order_x_axis), value = 251-value) %>%
  ggplot(aes(Geography, key, fill = value)) + geom_tile() +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
    scale_fill_gradient(low="purple", high="yellow") +
    labs(x = "County", y = "R0 Sensitivity Analysis", fill="rank")

save_plot("../ExploratoryFigures/county_rank.pdf", rank_plot, base_height = 5, base_aspect_ratio = 1.2)


cty_corr_plot <- rnots %>% select(-(1:2)) %>%
  cor(method = "spearman") %>%
  as_data_frame() %>%
  mutate(key1 = colnames(.)) %>%
  gather(key2, value, rnott.expected:August) %>%
  mutate(key1 = factor(key1, levels = order_x_axis), key2 = factor(key2, levels = order_x_axis)) %>%
  ggplot(aes(key2, key1, fill=value)) + geom_tile() +
    coord_cartesian(expand=FALSE) + 
    theme(axis.text.x = element_text(angle=45, vjust = 1, hjust = 1)) +
    scale_fill_gradient(low="purple", high="yellow") +
    labs(x = "R0 Sensitivity Analysis", y = "R0 Sensitivity Analysis", fill="Correlation")

save_plot("../ExploratoryFigures/county_corr.pdf", cty_corr_plot, base_height = 5, base_aspect_ratio = 1.3)

