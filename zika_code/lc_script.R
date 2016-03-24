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
branch_params <- function(prop_p = 1.7/7 , 
                          recov_p = 1.0/7,
                          d_thres = 5,
                          e_thresh = 300,
                          prob_symp = 1,
                          incub_p = 1,
                          dis_prob_symp = 1,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000)
  return(as.list(environment()))

last_infected <- all_last_instantInf_values(trials)
end_cum_infected <- all_last_cuminfect_values(trials)
last_cum_detected <- all_last_cumdetect_values(trials)

sort(last_cum_detected)
hist(end_cum_infected)

detect.vs.cumulative <- cbind(last_cum_detected, end_cum_infected)


#Plotting Preliminary Results
lastvalues <- all_last_cuminfect_values(trials)
hist(lastvalues, main = paste( "1000 Runs, R0 =", prop_p/recov_p, sep = ""), breaks = 10, ylim = c(1,1000))

plot(last_infected, lastvalues)
last_detect_values <- all_last_cumdetect_values(trials)
hist(last_detect_values)
plot(last_detect_values, lastvalues,  xlab = "Cumulative Detected", ylab =  "Cumulative Infected", main = "1000 Runs-Running until Instantaneous Infected > e_thres")



## Writing Functions for Plotting Heat Maps




prop_range <- seq(from = 1.5/7, to = 1.6/7, by = 0.1/7) #Changed to =  2.0/7 for testing
disc_range <- seq(from = 0.01, to = 0.02, by = .01) # Changed discover range to - 0.10



 
for (m in 1:length(prop_range)) {
  print(paste("Starting Prop Range", m, "-"))
  
  for (j in 1:length(disc_range)) {
    
    R0  = prop_range[m]*7
    
    #Running trials
    trials <- run_branches(num_reps = 1000, branch_params(dis_prob_symp=disc_range[j], prop_p = prop_range[m]))
    lastdetected <- all_last_cumdetect_values(trials)
    max <- max(lastdetected)*.25
    dect.cases <- seq(1:max)
        
    print(paste("Finished Trials", j, sep = "-"))
    
    #setting up bins to calculate frequencies 
    d_thres <- dect.cases[length(dect.cases)] #Highest Number of Cases to Consider 
    
    bins.prev <- set.prev.bins(d_thres, trials)
    bins.cumulative <- set.cum.bins(d_thres, trials)
    
    #resetting d_thres for trials

    
    
    #Setting Up matrices 
    thres.matrix.prev <- data.frame(matrix(nrow = length(dect.cases), ncol = length(bins.prev)-1))
    colnames(thres.matrix.prev) <- paste("< ", bins.prev[2:length(bins.prev)]); rownames(thres.matrix.prev) <- paste("Detected Cases = ", dect.cases)
    
    thres.matrix.cum <- data.frame(matrix(nrow = length(dect.cases), ncol = length(bins.cumulative)-1))
    colnames(thres.matrix.cum) <- paste("<", bins.cumulative[2:length(bins.cumulative)]); rownames(thres.matrix.cum) <- paste("Detected Cases = ", dect.cases)
    
    
    #Writing the values 
    for (i in 1:length(dect.cases)) {
      d_thres <- dect.cases[i]
      dataframe <- all_detect_rows(trials) # takes already whatever the current thresholds 
      frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
      frequencies.cum <- bin.frequency(dataframe[,8], bins.cumulative)
      print(paste("Sum of Prevalence Frequences =", sum(frequencies.prev), sep = ""))
      print(paste("Sum of Cumulative Frequences =", sum(frequencies.cum), sep = ""))
      
      thres.matrix.prev[i,] <- frequencies.prev
      thres.matrix.cum[i,] <- frequencies.cum
    }
    print("Calculated Frequencies")
    
    #Writing and saving the files - Not necessary to save all of these for now to check 
    filename.prev <- paste(name.generator(R0, disc_range[j], "Prev"), "csv", sep = ".")
    write.csv(x = thres.matrix.prev, file = filename.prev, row.names = TRUE)
    
    filename.cum <- paste(name.generator(R0, disc_range[j], "Cumulative"), "csv", sep = ".")
    write.csv(x = thres.matrix.cum, file = filename.cum, row.names = TRUE)
    
    # Plotting Function to go here
    prev.map <- plotheatmaps(thres.matrix.prev, type = "Prevalence", names = as.factor(dect.cases), R0 = R0, disc_value = disc_range[j])
    filename.prev.map <- paste(name.generator(R0, disc_range[j], "Prev"), "pdf", sep = ".")
    ggsave(filename = filename.prev.map, plot = prev.map, width=14, height=9)
    
    
    cum.map <- plotheatmaps(thres.matrix.cum, type = "Cumulative", names = as.factor(dect.cases), R0 = R0, disc_value = disc_range[j])
    filename.cum.map <- paste(name.generator(R0, disc_range[j], "Cumulative"), "pdf", sep = ".")
    ggsave(filename = filename.cum.map, plot = cum.map, width=14, height=9)
  } 
}

