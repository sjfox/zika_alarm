rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

# county <- read.csv("../county_data.csv")
# 
# hist(county$HabitatSuitability, breaks=100)
# county[which.max(county$HabitatSuitability),]

sapply(c('branch.functions.R','plot.functions.R', 'incubation_branch.R', 'analyze_saved_sims.R'), source)
library(plyr)
library(cowplot)
# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

#Parameters 
branch_params <- function(r_not = 1.1,
                          infBoxes = 3,
                          incBoxes = 6,
                          recov_p = 0.3040571/(3/infBoxes),
                          incub_rate = 0.583917,
                          prop_p =  r_not*recov_p/infBoxes, 
                          d_thres = 5,
                          e_thresh = 500,
                          prob_symp = 1,
                          dis_prob_symp = .01,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000)
  return(as.list(environment()))

dir_path <- "~/projects/zika_alarm/data/first_runs/"
save_path <- "~/projects/zika_alarm/data/"

save_final_sizes(dir_path, save_path)


r_nots <- c(0.9, 1.2, 1.8)
intro_rate <- c(0.3)
disc_prob <- c(0.011)

files <- paste("zika_sims", r_nots, intro_rate, disc_prob, sep="_")
list.files(path=dir_path, pattern="*0.9_0.011_0.3.Rdata", full.names=T, recursive=FALSE)

trials <- run_branches_inc(num_reps = 500, branch_params(r_not=.8))
plot_final_sizes(trials)

temp <- all_last_cuminfect_values(trials)

load("../zika_sims_0.85_0.1_0.Rdata")
plot_final_sizes(trials)


getEscapebyD <- function(trials, e_thresh){
  d <- 0:50
  esc_data <- count_escapes_vec(trials, d, e_thresh)
  return(data.frame(d_thresh=d, probEsc=esc_data))
}
  
discoveries <- seq(0.01, .1, length.out = 2)
prop_ps <- c(1.2)/7
escape_data <- data.frame()
for(disc_p in discoveries){
  for(prop in prop_ps){
    trials <- run_branches_inc(num_reps = 1000, branch_params(prop_p=prop, dis_prob_symp = disc_p))
    escape_data <- rbind(escape_data, cbind(disc_p=disc_p, prop_p=prop, getEscapebyD(trials, branch_params()$e_thresh)))   
  }
}
escape_data$r_not <- escape_data$prop_p*7

plot1 <- ggplot(escape_data, aes(d_thresh, probEsc, color = as.factor(r_not))) + 
  geom_line(size=1.5, aes(linetype=as.factor(disc_p))) + #facet_wrap(~disc_p)+
  scale_y_continuous(expand=c(0.01,0.01)) +
  scale_color_brewer(palette="Set1")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank())+
  labs(x = "Cumulative Number of Detected Cases", y = "Probability of an Epidemic", color = expression("R"[0]))+
  guides(linetype= FALSE )

print(plot1)

save_plot(filename = "../ExploratoryFigures/d_thresh_plot.pdf", plot = plot1, base_aspect_ratio = 1.5)

temp <- escape_data[which(escape_data$disc_p==.1), ]


trials <- run_branches_inc(num_reps = 1000, branch_params(prop_p=1.1/7, dis_prob_symp = .8))

p <- plot_final_sizes(trials)
print(p)
count_escapes(trials, 5, 300)

ggsave(filename = "../ExploratoryFigures/r0_0.2_disc_0.01_hist.pdf", plot = p, width=6, height=5)


###################################################
## Determining incubation period boxes and infectious period boxes
## simulate distributions
# incubation_time <- function(nboxes, prob){
#   box = 1
#   t = 0
#   while(box <= nboxes){
#     if(runif(1) < prob){
#       box = box+1
#     }
#     t = t + 1
#   }
#   t
# }
# rincubation <- function(num_reps, ...) {
#   raply(.n = num_reps, .expr = incubation_time(...) )
# }

fitDist <- function(prob, size, mean) {
  draws <- rnbinom(n=100000, size = size, prob = prob) + size
  mean - mean(draws)
}

## Look at weibull distribution from lessler 
qweibull(c(0.025, 0.5, 0.975), shape=1.97, scale=10.87)
hist(rweibull(100000, shape=1.97, scale=10.87), breaks=100)

## First fit the weibll viral clearance distribution
## Lessler et al biorxiv key times paper
## mean = 9.88, cI : 
n_boxes <- 3
infectious_best_fit <- uniroot(f = fitDist, interval = c(0.001,.999), size= n_boxes, mean=9.88)
infectious_times <- rnbinom(100000, size = n_boxes, prob = infectious_best_fit$root) + n_boxes
hist(infectious_times, breaks=50)
quantile(infectious_times, c(.025, .5, .975))

## Now determine the incubation period boxes
## Use the serial interval from brownstein paper
## 10-23 days 
fitDistBoth <- function(prob, size, min, max, infectious_times) {
  draws <- rnbinom(n=100000, size = size, prob = prob) + size
  qs <- quantile((draws+infectious_times/2), probs = c(.025, .975))
  sum(qs - c(min,max))^2
}
inc_boxes = 6
incubation_best_fit <- optimize(f = fitDistBoth, interval = c(0.001,.999), size= inc_boxes, min=10, max=23, infectious_times=infectious_times)
incubation_times <- rnbinom(n = 100000, size = inc_boxes, prob = incubation_best_fit$minimum) + inc_boxes
quantile((incubation_times+infectious_times/2), probs = c(.025, .5, .975))
hist(incubation_times+infectious_times/2, breaks=100, xlim=c(0,45))

infectious_best_fit$root
incubation_best_fit$minimum

mean(infectious_times)







###############################################################
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
