
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
  
  # Save the results or simply return them
  if(saveResults){
    save( list = c('final_sizes'), file = file.path(saveLoc, paste0("final_sizes.Rdata")))  
  } else {
    final_sizes
  }
}










