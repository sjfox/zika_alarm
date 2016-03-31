rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


sapply(c('branch.functions.R','plot.functions.R'), source)
library(plyr)
library(ggplot2)
# Code I've written just testing things out-not worth saving in the main files 
#Parameters 
}

branch_params <- function(prop_p = 1.7/7 , 
                          recov_p = 1.0/7,
                          d_thres = 5,
                          e_thresh = 200,
                          prob_symp = 1,
                          incub_rate = 1/16.5,
                          dis_prob_symp = 0.0246,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000,
                          zeroInc_prob = 6.825603e-08)
  return(as.list(environment()))


# THIS WORKS FOR 
prevalence.long <- data.frame()
cumulative.long <- data.frame()
disc_p <- seq(from = .01, to = .1, length.out = 2)

prop_range <- seq(from = 1.5/7, to = 1.5/7, length.out = 1) 
#disc_range <- seq(from =  0.1, to =  0.1, by = .1) 
confidence.range <- seq(from = .6, to = .9, by = .1)

trials <- run_branches_inc(num_reps = 1000, branch_params(dis_prob_symp=.01, prop_p = 1.5/7, e_thresh = 500))
lastdetected <- all_last_cumdetect_values(trials)

R0 = 1.5
disc_p =.01

for (prop_p in prop_range)) {  
  for (confidence in confidence.range) {
    
    R0  = prop_p*7
    percent.discover <- calculate.discover(disc_p) 

    #Running trials
    trials <- run_branches_inc(num_reps = 1000, branch_params(dis_prob_symp=disc_p, prop_p = prop_p, e_thresh = 500))
    lastdetected <- all_last_cumdetect_values(trials)
    
    max <- set.max.bin(max(lastdetected)) 
    dect.cases.range <- seq(1:max)
    d_thres <- max
        
    #setting up bins to calculate frequencies 
    bins.prev <- set.prev.bins(d_thres, trials)
    bins.cumulative <- set.cum.bins(d_thres, trials)
    
    #Setting Up data frames and vectors to store 
    thres.matrix.prev <- data.frame()
    thres.matrix.cum <- data.frame()
    
    median.vec <- numeric()
    average.vec <- numeric()
  
    #Writing the values 
    # Could split this up into two analysis: average/median and frequencies 
    for (dect.case in dect.cases.range) {
      d_thres <- dect.case
      dataframe <- all_detect_rows(trials) # takes already whatever the current thresholds 
      
      median.vec <- c(median.vec, median(dataframe$Cumulative_Infections))
      average.vec <- c(average.vec, mean(dataframe$Cumulative_Infections))
      
      frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
      thres.matrix.prev <- cbind.all(thres.matrix.prev, frequencies.prev)
      
      frequencies.cum <- bin.frequency(dataframe[,8], bins.cumulative)
      thres.matrix.cum <- data.frame(cbind.all(thres.matrix.cum, frequencies.cum))
    }      
    
    #dect.case.median <- function(dect.case, d_thres) {
    #  d_thres = dect.case
    #  dataframe = all_detect_rows(trials)
    #  median.vec <- median(dataframe$Cumulative_Infections)
    #  return(median.vec)
    #}
    #dect.case.median(13)
 
     #dect.range.median.analysis(dect.cases.range)
     
    
     
     #dect.range.median.analysis <- function(dect.cases.range) {
    #  return(ddply(.data = dect.cases.range,.variables = dect.cases.range, .fun =  dect.case.median))
    #}
    
    #  average.vec <- c(average.vec, mean(dataframe$Cumulative_Infections))
      
    #  frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
    #  frequencies.cum <- bin.frequency(dataframe[,8], bins.cumulative)
      
     # summary <- list(median = median.vec, average = average.vec, frequencies.prev = frequencies.prev, frequencies.cum=frequencies.cum)
    #  return(summary)
    #}

  
    
    # ANALYSIS FOR TRIGGER THRESHOLD
    integer.prev <- find_thres_cases(bins.prev, thres = 10, df=thres.matrix.prev, threshold_value=confidence)
    integer.cumulative <- find_thres_cases(bins = bins.cumulative, thres = 100, df = thres.matrix.cum, threshold_value = confidence)
     
    integer.cumulative <- find_thres_cases(bins = bins.cumulative, thres = 100, df = thres.matrix.cum, threshold_value = .1)
    prevalence.long <- rbind(prevalence.long, cbind(disc_p=disc_p, prop_p=prop_p, confidence = confidence, cases = unname(integer.prev)))
    cumulative.long <- rbind(cumulative.long, cbind(disc_p=disc_p, prop_p=prop_p, confidence = confidence, cases = unname(integer.cumulative)))
  

    # IF WANT TO WRITE AND SAVE             
    #Writing and saving the files - Not necessary to save all of these for now to check 
    #filename.prev <- paste(name.generator(R0, percent.discover, "Prev"), "csv", sep = ".")
    #write.csv(x = thres.matrix.prev, file = filename.prev, row.names = TRUE)
    
    #filename.cum <- paste(name.generator(R0, percent.discover, "Cumulative"), "csv", sep = ".")
    #write.csv(x = thres.matrix.cum, file = filename.cum, row.names = TRUE)
    
    # FUNCTIONS FOR PLOTTING A PARTICULAR PARAMETER SET 
    prev.map <- plotheatmaps(thres.matrix.prev, type = "Prevalence", names = as.factor(dect.cases.range), R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.cumulative, 1))
    filename.prev.map <- paste(name.generator(R0, percent.discover, "Prev_Incubation"), "pdf", sep = ".")
    ggsave(filename = filename.prev.map, plot = prev.map, width=14, height=9)
    
    cum.map <- plotheatmaps(thres.matrix.cum, type = "Cumulative", names = as.factor(dect.cases.range), R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.cumulative, 1))
    filename.cum.map <- paste(name.generator(R0, percent.discover, "Cumulative_Incubation"), "pdf", sep = ".")
    ggsave(filename = filename.cum.map, plot = cum.map, width=14, height=9)
    
    middle.stats.average <- data.frame(cbind(dect.cases, average.vec))
    plot.average <- ggplot(middle.stats.average, aes(dect.cases, average.vec))
    plot.average = plot.average + geom_line(size = 2) + labs(x = "Cumulative Detected Cases", y = "Cumulative Total Cases") + theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size = 20)) + 
      scale_x_continuous(breaks = x.breaks.seq) +
      scale_y_log10(breaks = c(1,5,10, 25,50,75, 100)) +  expand_limits(y=c(0, max = tail(bins.cumulative, 1)))
    
    plot.grid = plot_grid(cum.map, plot.average, labels = c("A", "B"), ncol = 1)
    
   } 
}





plot1 <- ggplot(prevalence.long, aes(prop_p, cases, color = as.factor(confidence))) + geom_line(size=2) + facet_wrap(~disc_p)+
  scale_y_continuous(expand=c(0,0)) +
  scale_color_brewer(palette="Set1") +
  labs(x = "Daily Propogation Rate", y = "Detected Cases") +
  theme(legend.text = element_text(size = 16)) +
  guides(fill=guide_legend(title="Confidence Probability")) + 
  ggtitle("Prevalence Threshold of 10 current cases")

plot2 <- ggplot(cumulative.long, aes(prop_p, cases, color = as.factor(confidence))) + geom_line(size=2) + facet_wrap(~disc_p)+
  scale_y_continuous(expand=c(0,0)) +
  scale_color_brewer(palette="Set1") + 
  labs(x = "Daily Propogation Rate", y = "Detected Cases") +
  theme(legend.text = element_text(size = 16)) +
  guides(fill=guide_legend(title="Confidence Probability")) + 
  ggtitle("Cumulative threshold of 100 cases")

middle.stats.m <- melt(middle.stats)

x.breaks.seq <- seq(0,25,1)
y.breaks.seq <- seq(0,60,10)



plot.average

middle.stats.median <- data.frame(cbind(dect.cases, median.vec))
plot.median <- ggplot(middle.stats.median, aes(dect.cases, median.vec))
plot.median + geom_line(size = 2)



####### Comparing histograms of outbreak size with and without incubation period
trials.simple<- run_simple_branches(10000, prop_p=1.7/7, recov_p=1/7, disc_p=.01, d_thresh=5, e_thresh=200)
#trials.no.inc <- run_branches_noinc(num_reps = 1000, branch_params())
trials.incubation <-  run_branches(num_reps = 10000, branch_params(prop_p = 1.7/7, e_thresh = 200, incub_p = 1/15.5))
#trials.with <- run_branches(num_reps = 1000, branch_params(incub_p = 1/16.5))

par(mfrow = c(1,2))
final.sizes <- all_getMaxCumI(trials.simple)
hist(final.sizes, breaks=100,  main = "Simple Model", xlab = "Final Size")
hist(final.sizes, breaks=100, ylim = c(0,2000), main = "Simple Model-Close Up", xlab = "Final Size")
median(final.sizes)
mean(final.sizes)
sort(final.sizes) # Largest outbreak is 52


outbreak.detections.incubation <- all_last_cumdetect_values(trials.incubation) # Highest small case 31 
hist(outbreak.detections.incubation, breaks = 100, main = "Complex Model", xlab = "Final Detection")
hist(outbreak.detections.incubation, breaks = 100, ylim = c(0,200), main = "Complex Model - Close Up", xlab = "Outbreak Detections")
mean(outbreak.size.inc1)
median(outbreak.size.inc1)
sorted.outbreak.size <- sort(outbreak.size.inc1)
sorted.outbreak.size # Highest small case 50


