
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
    # All were negatives
    threshold = 0
  } else if (threshold.positive[length(threshold.positive)] == (1-confidence)) { 
    # Were always above 1 
    threshold = NA
  } else 
    threshold <- length(threshold.positive)
  return(threshold)
}



