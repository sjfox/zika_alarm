rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

sapply(c('branch.functions.R','plot.functions.R', 'incubation_branch.R'), source)
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
                          e_thresh = 200,
                          prob_symp = 1,
                          incub_rate = 1/16.5,
                          dis_prob_symp = 1,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000,
                          zeroInc_prob = 6.825603e-08)
  return(as.list(environment()))



getEscapebyD <- function(trials, e_thresh){
  d <- 0:50
  esc_data <- count_escapes_vec(trials, d, e_thresh)
  return(data.frame(d_thresh=d, probEsc=esc_data))
}
  
discoveries <- seq(0.01, .95, length.out = 4)
prop_ps <- seq(0.9, 1.3, by=.1)/7
escape_data <- data.frame()
for(disc_p in discoveries){
  for(prop in prop_ps){
    trials <- run_branches_inc(num_reps = 1000, branch_params(prop_p=prop, dis_prob_symp = disc_p))
    escape_data <- rbind(escape_data, cbind(disc_p=disc_p, prop_p=prop, getEscapebyD(trials, branch_params()$e_thresh)))   
  }
}
escape_data$r_not <- escape_data$prop_p*7
plot1 <- ggplot(escape_data, aes(d_thresh, probEsc, color = as.factor(r_not))) + geom_line(size=2) + facet_wrap(~disc_p)+
  scale_y_continuous(expand=c(0,0)) +
  scale_color_brewer(palette="Set1")
print(plot1)

ggsave(filename = "../ExploratoryFigures/d_thresh_plot.pdf", plot = plot1, width=10, height=8)


trials <- run_branches_inc(num_reps = 1000, branch_params(prop_p=1.1/7, dis_prob_symp = .8))

p <- plot_final_sizes(trials)
print(p)
count_escapes(trials, 5, 300)

ggsave(filename = "../ExploratoryFigures/r0_0.2_disc_0.01_hist.pdf", plot = p, width=6, height=5)



# run_simple_branches <- function(num_reps, ...){
#   rlply(.n = num_reps, .expr = run_branch_simple(...) )
# }
# 
# getMaxCumI <- function(x){
#   return(x$cumI[nrow(x)])
# }
# all_getMaxCumI <- function(x) {
#   laply(x, getMaxCumI)
# }
# 
# trials <- run_simple_branches(1000, prop_p=.85/7, recov_p=1/7, disc_p=.01, d_thresh=5, e_thresh=200)
# final.sizes <- all_getMaxCumI(trials)
# hist(final.sizes, breaks=100)
