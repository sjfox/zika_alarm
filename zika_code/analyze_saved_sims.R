
get_vec_of_files <- function(dir_path, r_nots, disc_probs, intro_rates){
  data_files <- c()
  for(r_not in r_nots){
    for(intro_rate in intro_rates){
      for(disc_prob in disc_probs){
        pattern <- paste0("*", paste(r_not,  disc_prob, intro_rate, sep="_"), ".Rdata")
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
test_escape <- function(df, d_thres, e_thresh){
  if(last_cumdetect_value(df) < d_thres) return(NA)
  if(last_cuminfect_value(df)>=e_thresh) {
    TRUE
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

get_escape_prob_by_d <- function(dir_path, r_nots, disc_probs, intro_rates, e_thresh=500){
  data.files <- get_vec_of_files(dir_path, r_nots,  disc_probs, intro_rates)
  ldply(data.files, function(x) {
    load(x)
    parms <- get_parms(x)
    escape_probs <- calc_escape_prob_by_d(trials, e_thresh)
    
    cbind(as.data.frame(parms), escape_probs)
  })
}







