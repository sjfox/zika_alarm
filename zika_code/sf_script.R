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

# save_final_sizes(dir_path, save_path)

r_nots <- c(0.9, 1.1, 1.5)
disc_probs <- c(0.011, 0.068)
intro_rates <- c(.01, 0.05, 0.1, .3)

test <- get_escape_prob_by_d(dir_path, r_nots, disc_probs, intro_rates)

plot1 <- ggplot(test, aes(d_thresh, prob_esc, color = as.factor(r_not))) + 
  geom_line(size=1, aes(linetype=as.factor(disc_prob))) + facet_wrap(~intro_rate)+
  scale_y_continuous(expand=c(0.01,0.01)) +
  scale_color_brewer(palette="Set1")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank())+
  labs(x = "Cumulative Number of Detected Cases", y = "Probability of an Epidemic", color = expression("R"[0]))+
  guides(linetype= FALSE)

print(plot1)
save_plot(filename = "../ExploratoryFigures/epi_prob_by_detected.pdf", plot = plot1, base_aspect_ratio = 1.5)




######################################
## Plot maps of Texas







# get_rows_final_sizes <- function(final_sizes, r_nots, disc_probs, intro_rates){
#   final_sizes[which(final_sizes$r_not %in% r_nots & final_sizes$disc_prob %in% disc_probs & final_sizes$intro_rate %in% intro_rates), ]
# }
# load("../data/final_sizes.Rdata")
# r_nots <- c(0.8, 1.2)
# disc_probs <- c(0.011)
# intro_rates <- c(0, .01, 0.3, 1.92)
# 
# temp <- get_rows_final_sizes(final_sizes, r_nots, disc_probs, intro_rates)
# 
# final_plot <- ggplot(temp, aes(cumI)) + geom_histogram(bins=100)+
#   facet_grid(intro_rate~r_not, scales = "free_y") + 
#   scale_y_continuous(expand=c(0,0)) + 
#   scale_x_continuous(expand=c(0,0))+
#   labs(x = "Final Epidemic Sizes")
# save_plot(filename = "../ExploratoryFigures/final_sizes.pdf", plot = final_plot, base_height = 10, base_width = 12)
# 
# 
