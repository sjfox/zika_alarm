
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
cumcases_by_detects <- function(df){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_total(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=100]
  
  data.frame(detected = unique_detects, cum_infections = last_cuminfect_value(df), max_prevalence = max_prevalence(df))
}

get_cumcases_by_detects_all <- function(x){
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, cumcases_by_detects)
}


freq_above_thresh <- function(df, detected, cum_threshold, prev_threshold){
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of times that
  ## prevalence for a specific detection criteria is below a threshold
  rows <- which(df[,"detected"] == detected)
  if(length(rows)==0){
    return(NA)
  }else{
    ## Return number of rows that excede both thresholds divided by the total rows
    sum(df[rows, "cum_infections"] >= cum_threshold & df[rows,"max_prevalence"] >= prev_threshold) / length(rows)
  }
}
freq_above_thresh_vec <- Vectorize(freq_above_thresh, vectorize.args = "detected")

get_epidemic_prob_by_d <- function(trials, prev_threshold, cum_threshold, max_detect=50){
  ## Returns the probability in a given set of trials that the prevalence is below
  ## a specified threshold when X number of cases have been detected
  
  detected <- seq(0, max_detect)
  
  data <- get_cumcases_by_detects_all(trials)

  probs <- freq_above_thresh_vec(data, detected, cum_threshold, prev_threshold)
  return(data.frame(detected=detected, prob_epidemic=probs))
}


get_epidemic_prob_plot <- function(dir_path, prev_threshold, cum_threshold, r_nots, disc_probs, intro_rates){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    prob_belows <- get_epidemic_prob_by_d(trials = trials, prev_threshold=prev_threshold, cum_threshold=cum_threshold)  
    cbind(as.data.frame(parms), prob_belows)
  })  
}


#############################################
## Get porbability below thresholds by detection functions

# localprev_by_localdetects <- function(df){
#   ## Takes in a data frame trials, and for each
#   ## First instance of a new local detection, returns the local prevalence
#   
#   all_detects <- cum_detect_local(df)
#   unique_detects <- unique(all_detects)
#   ## Only  interested in maximum of 100 detections
#   unique_detects <- unique_detects[unique_detects<=100]
#   
#   matches <- match(unique_detects, all_detects)
#   data.frame(detected = unique_detects, prevalence = prevalence_local(df)[matches])
# }
# 
# totalprev_by_localdetects <- function(df){
#   ## Takes in a data frame trials, and for each
#   ## First instance of a new local detection, returns the total prevalence
#   
#   all_detects <- cum_detect_local(df)
#   unique_detects <- unique(all_detects)
#   ## Only  interested in maximum of 100 detections
#   unique_detects <- unique_detects[unique_detects<=100]
#   
#   matches <- match(unique_detects, all_detects)
#   data.frame(detected = unique_detects, prevalence=prevalence_total(df)[matches])
# }  


totalprev_by_totaldetects <- function(df){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_total(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=100]
  
  matches <- match(unique_detects, all_detects)
  data.frame(detected = unique_detects, prevalence=prevalence_total(df)[matches])
}

# localprev_by_totaldetects <- function(df){
#   ## Takes in a data frame trials, and for each
#   ## First instance of a new local detection, returns the total prevalence
#   
#   all_detects <- cum_detect_total(df)
#   unique_detects <- unique(all_detects)
#   ## Only  interested in maximum of 100 detections
#   unique_detects <- unique_detects[unique_detects<=100]
#   
#   matches <- match(unique_detects, all_detects)
#   data.frame(detected = unique_detects, prevalence=prevalence_local(df)[matches])
# }

get_prev_by_detects_all <- function(x, f){
  ## Returns data frame of all prevalence by detections for all trials
  ldply(x, f)
}

freq_below_thresh <- function(df, detected, threshold){
  ## Takes in dataframe of all prevalence by detect
  ## Returns a single frequency of times that
  ## prevalence for a specific detection criteria is below a threshold
  rows <- which(df[,"detected"] == detected)
  if(length(rows)==0){
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
  probs <- freq_below_thresh_vec(data, detected, threshold)
  return(data.frame(detected=detected, prob_below=probs))
}


get_prob_below_plot <- function(dir_path, thresholds, r_nots, disc_probs, intro_rates){
  data.files <- get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    data <- data.frame()
    for(threshold in thresholds){
      prob_belows <- get_prob_below_threshold(trials = trials, f=totalprev_by_totaldetects, threshold=threshold)  
      data <- rbind(data, cbind(as.data.frame(parms), data.frame(threshold=threshold), prob_belows))
    }
    data
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
  calculate_triggers <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    trigger <- get_surveillance_trigger(trials, f = totalprev_by_totaldetects, threshold, confidence, max_detect = 100)
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate) 
    result <- cbind(as.data.frame(matrix(parms,ncol=3)), threshold, confidence, trigger)
    return(result)
  })
}
  

get_surveillance_trigger <- function(trials, f, threshold, confidence, type="all", max_detect=100){
  ## Returns the max number of detected cases based on
  ## a specified threshold when X number of cases have been detected 
  ## and a tolerance for being X sure
  detected <- seq(0,max_detect) 
  data <- get_prev_by_detects_all(trials, f=totalprev_by_totaldetects) 
  probs <- freq_below_thresh_vec(data, detected, threshold=20)
  threshold.probs <- probs - confidence
  threshold.positive <- threshold.probs[which(threshold.probs > 0)]
  if (length(threshold.positive) == 0) { 
    # All were negatives-already took off
    threshold = 0
  } else if (threshold.positive[length(threshold.positive)] == (1-confidence) & is.na(threshold.probs[length(threshold.positive) + 1]))  { 
    # Never Hit Threshold 
    threshold = NA
  } else 
    threshold <- length(threshold.positive)
  return(threshold)
}



