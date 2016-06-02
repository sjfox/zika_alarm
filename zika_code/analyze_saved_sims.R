
get_vec_of_files <- function(dir_path, r_nots, disc_probs, intro_rates){
  data_files <- c()
  for(r_not in r_nots){
    for(intro_rate in intro_rates){
      for(disc_prob in disc_probs){
        pattern <- paste0("*_", paste(r_not,  disc_prob, intro_rate, sep="_"), ".Rdata")  
        data_files <- c(data_files, list.files(path=dir_path, pattern=pattern, full.names=T, recursive=FALSE))
      }
    }
  }
  data_files
}

save_final_sizes <- function(dirPath, saveLoc, saveResults=TRUE){
  data.files <- list.files(path=dirPath, pattern="*.Rdata", full.names=T, recursive=FALSE)
  final_sizes <- ldply(data.files, function(x) {
    load(x)
    # Setup the fixed parameters from the run being analyzed
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate)
    # Calculate the final infected and detected
    cumI <- all_last_cuminfect_values(trials)
    cumD <- all_last_cumdetect_values(trials)
    
    # Return parms and cumI/cumD in data.frame
    cbind(as.data.frame(matrix(parms,ncol=3)), cumI=cumI, cumD=cumD)
  })
  colnames(final_sizes) <- c("r_not", "disc_prob", "intro_rate", "cumI", "cumD")
  # Save the results or simply return them
  if(saveResults){
    save( list = c('final_sizes'), file = file.path(saveLoc, paste0("final_sizes.Rdata")))  
  } else {
    final_sizes
  }
}


get_parms <- function(path){
  ## Get parms of a run from the pathway to the run
  temp <- strsplit(path, split="/")[[1]]
  temp <- temp[length(temp)]
  temp <- strsplit(temp, split="_")[[1]]
  
  r_not <- as.numeric(temp[3])
  disc_prob <- as.numeric(temp[4])
  intro_rate_data <- strsplit(temp[5], "\\.")[[1]]
  
  intro_rate <- as.numeric(paste0(intro_rate_data[1:length(intro_rate_data)-1],collapse ="."))
  
  list(r_not=r_not, disc_prob=disc_prob, intro_rate=intro_rate) 
}

####################################################
## Functions for getting escape probability by detection threshold
cumcases_by_detects <- function(df, max_detect=100){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_total(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=max_detect]
  
  data.frame(detected = unique_detects, cum_infections = last_cuminfect_value(df), max_prevalence = max_prevalence(df))
}

get_cumcases_by_detects_all <- function(x, max_detect=100){
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, cumcases_by_detects, max_detect)
}


freq_above_thresh <- function(df, detected, cum_threshold, prev_threshold, num_necessary=10){
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of times that
  ## prevalence for a specific detection criteria is below a threshold
  rows <- which(df[,"detected"] == detected)
  if(length(rows)<=num_necessary){
    return(NA)
  }else{
    ## Return number of rows that excede both thresholds divided by the total rows
    sum(df[rows, "cum_infections"] >= cum_threshold & df[rows,"max_prevalence"] >= prev_threshold) / length(rows)
  }
}
freq_above_thresh_vec <- Vectorize(freq_above_thresh, vectorize.args = "detected")

get_epidemic_prob_by_d <- function(trials, prev_threshold, cum_threshold, max_detect=50, num_necessary){
  ## Returns the probability in a given set of trials that the prevalence is below
  ## a specified threshold when X number of cases have been detected
  
  detected <- seq(0, max_detect)
  
  data <- get_cumcases_by_detects_all(trials, max_detect = max_detect)

  probs <- freq_above_thresh_vec(data, detected, cum_threshold, prev_threshold, num_necessary)
  return(data.frame(detected=detected, prob_epidemic=probs))
}


get_epidemic_prob_plot <- function(dir_path, prev_threshold, cum_threshold, r_nots, disc_probs, intro_rates,max_detect=150, num_necessary=100){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    prob_belows <- get_epidemic_prob_by_d(trials = trials, prev_threshold=prev_threshold, cum_threshold=cum_threshold, num_necessary=num_necessary, max_detect=max_detect)  
    cbind(as.data.frame(parms), prob_belows)
  })  
}


#############################################
## Get porbability below thresholds by detection functions


totalprev_by_totaldetects <- function(df, max_detect){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_total(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=max_detect]
  
  matches <- match(unique_detects, all_detects)
  data.frame(detected = unique_detects, prevalence=prevalence_total(df)[matches])
}

get_prev_by_detects_all <- function(x, f, max_detect=100){
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, f, max_detect)
}

freq_below_thresh <- function(df, detected, threshold, num_necessary){
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of times that
  ## prevalence for a specific detection criteria is below a threshold
  rows <- which(df[,"detected"] == detected)
  if(length(rows)<=num_necessary){
    return(NA)
  }
  sum(df[rows, "prevalence"] < threshold)/ length(rows)
}

freq_below_thresh_vec <- Vectorize(freq_below_thresh, vectorize.args = "detected")

get_prob_below_threshold <- function(trials, f, threshold, max_detect=50){
  ## Returns the probability in a given set of trials that the prevalence is below
  ## a specified threshold when X number of cases have been detected 
  detected <- seq(0,max_detect) 
  data <- get_prev_by_detects_all(trials, f) 
  probs <- freq_below_thresh_vec(data, detected, threshold, num_necessary=100)
  return(data.frame(detected=detected, prob_below=probs))
}


get_prob_below_plot <- function(dir_path, thresholds, r_nots, disc_probs, intro_rates){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    data <- data.frame()
    for(threshold in thresholds){
      prob_belows <- get_prob_below_threshold(trials = trials, f=totalprev_by_totaldetects, threshold=threshold, max_detect = 150)  
      data <- rbind(data, cbind(as.data.frame(parms), data.frame(threshold=threshold), prob_belows))
    }
    data
  })  
}


get_prev_by_detects_plot <- function(dir_path, r_nots, disc_probs, intro_rates){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    prevalences <- get_prev_by_detects_all(trials, f=totalprev_by_totaldetects)  
    
    prevalences <- ddply(prevalences, .(detected), .fun = function(x){ 
      quants <-  quantile(x = x$prevalence, probs = c(0.5, 0.25, 0.75), names=FALSE) 
      data.frame(median=quants[1], min = quants[2], max = quants[3])
    })
    cbind(as.data.frame(parms), prevalences)
  })  
}


#########################################################
## Lauren functions
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

calculate_expect_vs_detect_total <- function(dir_path, r_nots, intro_rate, disc_prob) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  calculate_average_sd <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    detection.df <- get_prev_by_detects_all(trials, totalprev_by_totaldetects)
    average.prevalence <- ddply(.data = detection.df, .variables = "detected", summarize, mean(prevalence), sd(prevalence)) 
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate) 
    result <- cbind(as.data.frame(matrix(parms,ncol=3)), average.prevalence)
    return(result)
  })
}



#### Calculate Frequency Distributions For Each Number of Detected Cases
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
  frequency
}

freq_at_detect_vec <- Vectorize(freq_at_detect, vectorize.args = "detected", SIMPLIFY = FALSE)

get_freq_at_detect <- function(trials, f, type = "all", max_detect=50) {
  ## Returns the frequency distribution of actual cases in a given set of trials for each detection
  detected <- seq(0,max_detect)  
  data <- get_prev_by_detects_all(trials, f)
  freq.ds <- freq_at_detect_vec(data, detected)
  return(freq.ds)
}


## Function to get all surveillance for all R0s 
calculate_all_triggers <- function(dir_path, r_nots, intro_rate, disc_prob,threshold, confidence) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  ldply(dirPaths, function(x) {
    load(x)
    trigger <- get_surveillance_trigger(trials, threshold, confidence, max_detect = 200)
    parms <- get_parms(x)
    cbind(parms, data.frame(threshold=threshold, confidence=confidence, trigger=trigger))
  })
}

calculate_all_triggers_vec <- Vectorize(calculate_all_triggers, vectorize.args="confidence")


get_surveillance_trigger <- function(trials, threshold, confidence, max_detect=200, num_necessary=10){
  ## Returns the max number of detected cases based on
  ## a specified threshold when X number of cases have been detected 
  ## and a tolerance for being X sure
  detected <- seq(0, max_detect) 
  data <- get_prev_by_detects_all(x = trials, f=totalprev_by_totaldetects, max_detect = max_detect) 
  probs <- freq_below_thresh_vec(data, detected, threshold=threshold, num_necessary)
  # threshold.probs <- probs - confidence
  # threshold.positive <- threshold.probs[threshold.probs > 0]
  temp <- which(probs < confidence)
  
  if (is.na(temp[1])) { 
    # All were negatives-already took off
    threshold = NA
  } else if (temp[1]==1)  { 
    # Never Hit Threshold 
    threshold = 0
  } else {
    ## subtract 2, because 1 due to detecteds starting at 0, second due to finding
    ## the first match for being less than 0.8, and we want match just before 0.8
    threshold <- temp[1] - 2
  }
  return(threshold)
}


## Function to get all surveillance for all R0s 
calculate_all_epidemics <- function(dir_path, r_nots, intro_rate, disc_prob, threshold, confidence) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  ldply(dirPaths, function(x) {
    load(x)
    trigger <- get_epidemic_trigger(trials, threshold, confidence, max_detect = 200)
    parms <- get_parms(x)
    cbind(parms, data.frame(threshold=threshold, confidence=confidence, trigger=trigger))
  })
}

calculate_all_epidemics_vec <- Vectorize(calculate_all_epidemics, vectorize.args="confidence")


get_epidemic_trigger <- function(trials, threshold, confidence, max_detect=200, num_necessary=10){
  ## Returns the max number of detected cases based on
  ## a specified threshold when X number of cases have been detected 
  ## and a tolerance for being X sure
  detected <- seq(0, max_detect) 
  data <- get_epidemic_prob_by_d(trials = trials, prev_threshold = threshold, cum_threshold = 2000, max_detect = max_detect, num_necessary) 
  temp <- which(data$prob_epidemic >= confidence)
  if (is.na(temp[1])) { 
    # Never hit the threshold
    trigger = NA
  } else if (temp[1]==1)  { 
    # hit threshold right away
    trigger = 1
  } else {
    ## subtract 2, because 1 due to detecteds starting at 0
    trigger <- temp[1] - 1
  }
  return(trigger)
}

#############################
## Get saved trigger data
#############################
get_trigger_data <- function(rnot, intro, disc, prev_threshold=c(20), epi_threshold=c(50), confidence, num_necessary){
  ## type should equal "prevalence" or "epidemic"
  ## Gets saved trigger data, and returns data from requested runs
  load("../data/all_triggers.Rdata")
  
  df <- all_triggers
  df[which(as.character(df$r_not)%in%rnot & 
             as.character(df$intro_rate) %in% intro & 
             as.character(df$disc_prob)%in%disc  & 
             as.character(df$prev_threshold) %in% prev_threshold & 
             as.character(df$epi_threshold) %in% epi_threshold & 
             as.character(df$confidence) %in% confidence & 
             as.character(df$num_necessary) %in% num_necessary), ]  
  
}

combine_triggers <- function(dir_path, save_path) {
  ## Function that combines TACC trigger data output into one file
  data.files <- list.files(path=dir_path, pattern="*.Rdata", full.names=T, recursive=FALSE)
  all_triggers <- ldply(data.files, function(x) {
    load(x)
    triggers
  })
  save(list = c("all_triggers"), file = paste0(save_path, "all_triggers.Rdata"))
}

















# convert_trigger_data <- function(df){
#   ## Converts output of running all trigger runs to a dataframe of
#   ## specific format
#   
#   data.frame(r_not = unlist(df$r_not),
#            disc_prob=unlist(df$disc_prob),
#            intro_rate = unlist(df$intro_rate),
#            threshold=unlist(df$threshold),
#            confidence= unlist(df$confidence),
#            trigger= unlist(df$trigger))
# }

# save_calc_epidemic_triggers <- function(dir_path, save_path, threshold, confidences){
#   ## CAlculates and saves all epidemic probability triggers
#   r_nots <- seq(0.1, 2, by=0.1)
#   disc_probs <- c(0.0052, 0.011, 0.01635, .0287, 0.068) 
#   intro_rates <- c(0.0, 0.01, 0.05, 0.1, 0.2, 0.3)
#   epidemic_triggers <- as.data.frame(t(calculate_all_epidemics_vec(dir_path = dir_path, r_nots = r_nots, 
#                                                                    intro_rate = intro_rates, disc_prob = disc_probs, 
#                                                                    threshold = threshold, confidence=confidences)))
#   epidemic_triggers <- convert_trigger_data(epidemic_triggers)
#   save( list = c('epidemic_triggers'), file = file.path(save_path, paste0("epidemic_triggers.Rdata")))  
# }
# 
# save_calc_prev_triggers <- function(dir_path, save_path, threshold, confidences){
#   ## Calculates and saves all prevalence triggers
#   r_nots <- seq(0.1, 2, by=0.1)
#   disc_probs <- c(0.0052, 0.011, 0.01635, .0287, 0.068) 
#   intro_rates <- c(0.0, 0.01, 0.05, 0.1, 0.2, 0.3)
#   prevalence_triggers <- as.data.frame(t(calculate_all_triggers_vec(dir_path = dir_path, r_nots = r_nots, 
#                                                                     intro_rate = intro_rates, disc_prob = disc_probs, 
#                                                                     threshold = threshold, confidence=confidences)))
#   prevalence_triggers <- convert_trigger_data(prevalence_triggers)
#   save( list = c('prevalence_triggers'), file = file.path(save_path, paste0("prevalence_triggers.Rdata")))  
# }
