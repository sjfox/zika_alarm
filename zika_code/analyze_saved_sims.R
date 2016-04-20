
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
test_escape <- function(df, d_thres, e_thresh, prev_thresh=20){
  if(last_cumdetect_local_value(df) < d_thres) return(NA)
  if(last_cuminfect_value(df)>=e_thresh) {
    if(max_nonintro_prevalence(df) > prev_thresh){
      TRUE  
    } else{
      FALSE
    }
  } else{
    FALSE
  }
}

count_escapes <- function(x, d_thres, e_thresh){
  ## Function to get probability of escape, if detecteds
  ## are greater than the d_thres
  ## x must be list of runs
  
  escapes <- laply(x, test_escape, d_thres, e_thresh)
  numEscape <- sum(escapes, na.rm=T)
  numPossible <- sum(!is.na(escapes), na.rm=T)
  if(numPossible ==0) {
    NA
  }else {
    numEscape/numPossible
  }
}

count_escapes_vec <- Vectorize(count_escapes, vectorize.args = "d_thres")

calc_escape_prob_by_d <- function(trials, e_thresh){
  d <- 0:50
  esc_data <- count_escapes_vec(trials, d, e_thresh)
  return(data.frame(d_thresh=d, prob_esc=esc_data))
}

get_escape_prob_by_d <- function(dir_path, r_nots, disc_probs, intro_rates, e_thresh=1000){
  data.files <- get_vec_of_files(dir_path, r_nots,  disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    escape_probs <- calc_escape_prob_by_d(trials, e_thresh)
    
    cbind(as.data.frame(parms), escape_probs)
  })
}


#############################################
## Get porbability below thresholds by detection functions

localprev_by_localdetects <- function(df){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the local prevalence
  
  all_detects <- cum_detect_local(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=100]
  
  matches <- match(unique_detects, all_detects)
  data.frame(detected = unique_detects, prevalence = prevalence_local(df)[matches])
}

totalprev_by_localdetects <- function(df){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_local(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=100]
  
  matches <- match(unique_detects, all_detects)
  data.frame(detected = unique_detects, prevalence=prevalence_total(df)[matches])
}  


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

localprev_by_totaldetects <- function(df){
  ## Takes in a data frame trials, and for each
  ## First instance of a new local detection, returns the total prevalence
  
  all_detects <- cum_detect_total(df)
  unique_detects <- unique(all_detects)
  ## Only  interested in maximum of 100 detections
  unique_detects <- unique_detects[unique_detects<=100]
  
  matches <- match(unique_detects, all_detects)
  data.frame(detected = unique_detects, prevalence=prevalence_local(df)[matches])
}

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

get_prob_below_threshold <- function(trials, f, threshold, type="all", max_detect=50){
  ## Returns the probability in a given set of trials that the prevalence is below
  ## a specified threshold when X number of cases have been detected
  
  detected <- seq(0,max_detect)
  
  data <- get_prev_by_detects_all(trials, f)
  
  probs <- freq_below_thresh_vec(data, detected, threshold)
  return(data.frame(detected=detected, prob_below=probs))
}



#########################################################
## Lauren functions

calculate_expect_vs_detect <- function(dir_path, r_nots, intro_rate, disc_prob) {
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)
  calculate_average_sd <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
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
    average.cumulative <- ldply(.data = trials_by_detection, function (x) {
      mean.cumulative <- mean(x[,"Cumulative_Infections"])
      sd.cumulative <- sd(x[,"Cumulative_Infections"])
      return(c(mean.cumulative, sd.cumulative))
    })
    average.prevalence <- ldply(.data = trials_by_detection, function (x) {
      mean.prevalence <- mean(x[,"Total_Infections"])
      sd.prevalence<- sd(x[,"Total_Infections"])
      return(c(mean.prevalence, sd.prevalence))   
    })
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate) 
    result <- cbind(as.data.frame(matrix(parms,ncol=3)), dect.cases.range, 
                    average.cumulative[,2], average.cumulative[, 3], average.prevalence[,2], average.prevalence[,3])
    return(result)    
  })
}



# Calculate the Probability of Exceeding a specified number of cases 
calculate_frequency_threshold <- function(dir_path,r_nots, intro_rate, disc_prob, threshold.cumulative, threshold.prevalence) {
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
    
    #return(list(cumulative = frequency.cumulative, current = frequency.prevalence))
   
    integer.prev <- unname(find_thres_cases(bins.prev,  threshold.cases = threshold.prevalence, df=frequency.prevalence, 
                                            confidence.value = confidence))
    integer.cumulative <- unname(find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative,
                                                  confidence.value = confidence))
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, confidence, threshold.prevalence, threshold.cumulative) 
    cbind(as.data.frame(matrix(parms,ncol=6)), dect.cases.range, integer.prev, integer.cumulative)
  })
  return(triggers_threshold)
}

