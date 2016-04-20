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




prev.map <- plotheatmaps(thres.matrix.prev, type = "Prevalence", names = as.factor(dect.cases.range),
                         R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.prev, 1))
cum.map <- plotheatmaps(thres.matrix.cum, type = "Cumulative", names = as.factor(dect.cases.range),
                        R0 = R0, percent.discover = percent.discover, max.infect = tail(bins.cumulative, 1))



 

############################################
dir_path <- "~/Documents/zika_alarm/data/introductions/"
save_path <- "~/Doucments/zika_alarm/data"

r_nots <- c(1.5, 1.1, 0.9)
disc_prob <- c( 0.068, 0.011)
intro_rate <- c(0.01, 0.1, 0.3)


expected.cases.local <- calculate_expect_vs_detect_local(dir_path, r_nots, intro_rate, disc_prob)
expected.cases.total <- calculate_expect_vs_detect_total(dir_path, r_nots, intro_rate, disc_prob)

colnames(expected.cases.total) <- c("run", "r0","dect", "intro", "cases_dect", "avg", "sd" )
colnames(expected.cases.local) <- c("run", "r0","dect", "intro", "cases_dect", "avg", "sd" )

# If want to split the Results to certain detection values/intro values/R0
indices.local = which(expected.cases.local$intro == 0.3 | expected.cases.local$intro == 0.1)
indices.total = which(expected.cases.total$intro == 0.3 | expected.cases.total$intro == 0.1)

detection.local = expected.cases.local[indices.local, ]
detection.total = expected.cases.total[indices.total,]

breaks_y = seq(from = 0, to = 200, by = 25 )
breaks_y_log = c(0, 5,10,15,20,25,50,75,100,150,200) #, 200, 300,400,500) #, 200, 300, 400, 500) # Set according to max of value + se
breaks_x = seq(from = 0, to = 100, by = 10)


plot.1sd.local <- ggplot(detection.local, aes(cases_dect, avg, fill = as.factor(r0), color = as.factor(r0), group=interaction(as.factor(dect), r0)))  + 
  geom_line(size=1.5, aes(linetype = as.factor(dect))) + facet_wrap(~intro, ncol = 1) +
  geom_ribbon(aes(ymin = avg-sd, ymax=avg+sd), alpha=.2, color = NA) + 
  scale_y_continuous(breaks = breaks_y)  +
  scale_color_brewer(palette = "Set1", guide = FALSE, direction = -1) +
  scale_fill_brewer(palette="Set1", guide_legend(title = "R0"), direction = -1) +
  scale_x_continuous(name = "Cumulative Number of Detected Cases", breaks = breaks_x, limits = c(0,100)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) #+
  #theme(axis.title.y = element_text(size=28), axis.text.y = element_text(size = 22)) + 
  #theme(axis.title.x = element_text(size=28), axis.text.x= element_text(size=22)) +
  #labs(y = "Expected Total Current Cases", linetype = "Detection \n Rate") +
  #theme(legend.text=element_text(size=22, margin = margin(), debug = FALSE), legend.title = element_text(size = 28)) #+ 
#

plot.1sd.total <- ggplot(detection.total, aes(cases_dect, avg, fill = as.factor(r0), color = as.factor(r0), group=interaction(as.factor(dect), r0)))  + 
  geom_line(size=1.5, aes(linetype = as.factor(dect))) + facet_wrap(~intro, ncol = 1) +
  geom_ribbon(aes(ymin = avg-sd, ymax=avg+sd), alpha=.2, color = NA) + 
  scale_y_continuous(breaks = breaks_y)  +
  scale_color_brewer(palette = "Set1", guide = FALSE, direction = -1) +
  scale_fill_brewer(palette="Set1", guide_legend(title = "R0"), direction = -1) +
  scale_x_continuous(name = "Cumulative Number of Detected Cases", breaks = breaks_x, limits = c(0,100)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank())


plot_grid(plot.1sd.local, plot.1sd.total, labels = c("A", "B")) 


##### TRYING TO REWRITE TEXAS TRIGGER THRESHOLD

r_nots <- 1.1
intro_rate = .01
disc_prob = .011
load(dirPaths)

calculate_expect_vs_detect_local <- function(dir_path, r_nots, intro_rate, disc_prob) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  calculate_average_sd <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    detection.df <- get_prev_by_detects_all(trials, localprev_by_localdetects)
    average.prevalence <- ddply(.data = detection.df, .variables = "detected", summarize, mean(prevalence), sd(prevalence)) 
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate) 
    result <- cbind(as.data.frame(matrix(parms,ncol=3)), average.prevalence)
    return(result)
  })
}



###### Calculate Surveillance For R0s 
freq_at_detect <- function(df, detected) {
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of 
  ## prevalence for a specific detection 
  rows <- which(df[,"detected"] == detected)
  if(length(rows)==0){
    return(NA)
  }
  table.rows <- df[rows,]  
  frequency <- count(table.rows, vars = c("detected", "prevalence"))
  frequency.dis <- frequency$freq/sum(frequency$freq)
  frequency$freq <- frequency.dis
  return(frequency)
}


#freq_at_detect_vec <- Vectorize(freq_at_detect, vectorize.args = "detected")


get_freq_at_detect <- function(trials, f, type = "all", max_detect=50) {
  ## Returns the frequency distribution of actual cases in a given set of trials for each detection
 
   detected <- seq(0,max_detect)
  
   data <- get_prev_by_detects_all(trials, f)
   freq.dis <- freq_at_detect_vec(data, detected)
  return(data.frame(detected=detected, prob_below=probs))
}





calculate_surveillance_prevalence <- function(dir_path,r_nots, intro_rate, disc_prob, threshold.prevalence) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
 
   frequency_at_threshold <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    max.cumulative <- max(all_last_cuminfect_values(trials))
    max.prev <- max(all_max_prevalence(trials))
    lastdetected <- all_last_cumdetect_values(trials)
    max <- set.max.bin(max(lastdetected))
    dect.cases.range <- seq(1:max)
    
    #splits up trials into detection 
    trials_by_detection <- alply(.data = dect.cases.range, .margins = 1, function (x) {
      d_thres = as.numeric(x)
      dataframe <- all_detect_rows(trials, threshold = d_thres ) 
      dataframe = na.omit(dataframe)
      return(dataframe)
    })    
    
    #setting up bins to calculate frequencies
    bins.prev <- set.prev.bins(max.prev)
    bins.cumulative <- set.cum.bins(max.cumulative)    
    
    #Frequency Calculations 
    frequency.cumulative <- ldply(.data = trials_by_detection, function (x) {
      frequency = bin.frequency(x[,"Cumulative_Infections"], bins.cumulative)
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    frequency.prevalence <- ldply(.data = trials_by_detection, function (x) {
      frequency = unname(bin.frequency(x[,"Total_Infections"], bins.prev))
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    # Clean UP 
    colnames(frequency.prevalence) <- bins.prev; rownames(frequency.prevalence) <- dect.cases.range
    colnames(frequency.cumulative) <- bins.cumulative; rownames(frequency.cumulative) <- dect.cases.range
    frequency.prevalence <- frequency.prevalence[,-1]
    frequency.cumulative <- frequency.cumulative[,-1]
    
    #return(list(cumulative = frequency.cumulative, current = frequency.prevalence))
    # calculate frequency at specified threshold 
    threshold.frequency.prev <- frequency_threshold(bins = bins.prev,threshold.cases = threshold.prevalence, df = frequency.prevalence)
    threshold.frequency.cum <- frequency_threshold(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative)
    
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, threshold.prevalence, threshold.cumulative)
    cbind(as.data.frame(matrix(parms,ncol=5)),dect.cases.range, threshold.frequency.prev, threshold.frequency.cum)
  })
  return(frequency_at_threshold)
}

# functin to calculate the trigger based on a specificed number of cases and risk tolerance level 
calculate_surveillance_triggers <- function(dir_path,r_nots, intro_rate, disc_prob, threshold.cumulative, threshold.prevalence) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  triggers_threshold <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    max.cumulative <- max(all_last_cuminfect_values(trials))
    max.prev <- max(all_max_prevalence(trials))
    lastdetected <- all_last_cumdetect_values(trials)
    max <- set.max.bin(max(lastdetected))
    dect.cases.range <- seq(1:max)
    
    #splits up trials into detection 
    trials_by_detection <- alply(.data = dect.cases.range, .margins = 1, function (x) {
      d_thres = as.numeric(x)
      dataframe <- all_detect_rows(trials, threshold = d_thres ) 
      dataframe = na.omit(dataframe)
      return(dataframe)
    })    
    
    #setting up bins to calculate frequencies
    bins.prev <- set.prev.bins(max.prev)
    bins.cumulative <- set.cum.bins(max.cumulative)    
    
    #Frequency Calculations 
    frequency.cumulative <- ldply(.data = trials_by_detection, function (x) {
      frequency = bin.frequency(x[,"Cumulative_Infections"], bins.cumulative)
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    frequency.prevalence <- ldply(.data = trials_by_detection, function (x) {
      frequency = unname(bin.frequency(x[,"Total_Infections"], bins.prev))
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    # Clean UP 
    colnames(frequency.prevalence) <- bins.prev; rownames(frequency.prevalence) <- dect.cases.range
    colnames(frequency.cumulative) <- bins.cumulative; rownames(frequency.cumulative) <- dect.cases.range
    frequency.prevalence <- frequency.prevalence[,-1]
    frequency.cumulative <- frequency.cumulative[,-1]
    
    
    integer.prev <- unname(find_thres_cases(bins.prev,  threshold.cases = threshold.prevalence, df=frequency.prevalence, 
                                            confidence.value = confidence))
    integer.cumulative <- unname(find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative,
                                                  confidence.value = confidence))
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, confidence, threshold.prevalence, threshold.cumulative) 
    cbind(as.data.frame(matrix(parms,ncol=6)), dect.cases.range, integer.prev, integer.cumulative)
  })
  return(triggers_threshold)
}



