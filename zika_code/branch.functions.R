run_branch_simple <- function(prop_p, recov_p, disc_p, d_thresh, e_thresh) {
  UI = 1; DI = 0; D = 0
  I = UI + DI
  cumI = UI
  time_record <- data.frame(matrix(data = c(I,D, cumI), ncol = 3))
  colnames(time_record) <- c("I", "D", "cumI")
  
  while  (((I < e_thresh) & (I > 0)) | ((I > 0) & (D < d_thresh))) {
    #while Number of infected is below epidemic threshold hold and more than 0 infected
    # or while number of infecteds is above 0 and the number of detected is below threshold
    
    newDI_draws = runif(DI) #New discovered infected possibilities 
    newDI_count = sum(newDI_draws < prop_p) #Number of new infected by detected individuals 
    
    decDI_draws = runif(DI) #Detected indivduals-will they recover?
    decDI_count = sum(decDI_draws < recov_p)
    
    newUI_draws = runif(UI) #Number of possible new Undiscovred infecteds 
    newUI_count =  sum(newUI_draws < prop_p)  #Number of new Infecteds 
    
    decUI_draws = runif(UI) #Undiscovered infecteds that stay infected
    decUI_count = sum(decUI_draws < recov_p)
    
    disUI_draws = runif(UI)
    disUI_count = sum(disUI_draws < disc_p)  #probability that undetected individuals are discovered 
    
    removeUI = sum((decUI_draws < recov_p) | (disUI_draws < disc_p)) 
    
    #Undiscovered infecteds can remove by being recovered or by being discovoered 
    
    removeUDI =  sum((decUI_draws < recov_p) & ( disUI_draws < disc_p))
    
    #remove undiscovered infecteds that are both no longer infected and were eventually discovered 
    
    UI = UI + newUI_count + newDI_count - removeUI 
    #Undiscovered = previously undiscovered + newly infected by UI + newly infected by DI - either recovered or discovered 
    DI = DI - decDI_count + disUI_count - removeUDI
    #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered 
    D = D + disUI_count #Cumulatve Detected = those detected + newly detected
    I = UI + DI #Undiscovered and Discovered Infected 
    cumI = cumI + newUI_count + newDI_count
    time_record <- rbind(time_record,c(I,D,cumI))
  }
  return(time_record)
}

# Kinds of counters:
# UI - Undiscovered Infecteds, including UI_Symp and UI_Asymp
# DI - Discovered Infecteds, including DI_Symp and DI_Asymp
# I - Number of Infectious 
# IncI - Incubation Infecteds
# D - Cumulative number of discoveries
#
# Data to Store 
# Run Number 
# New_Exposed
# New_Detected 
# New_Infectious
# Total_Infected 
# Introductions 
# Comulative Infections
# Simulate until we have either:
# - we run out of infecteds
# - e_thresh number of infecteds, and d_thresh number of discoveries

run_branch <- function(params) {
  with(params,{
    UI_Symp = 1; UI_Asymp = 0; DI_Symp = 0; DI_Asymp = 0 
    UI = UI_Symp + UI_Asymp; DI = DI_Symp + DI_Asymp
    D = 0  
    incubationInfecteds = 0
    I = UI + DI
    cumI = UI
    time_record <- data.frame(matrix(data = 0, ncol = 7))
    colnames(time_record) <- c("New_Exposed", "New_Infectious", "Intros", "New_Detections", "Cum_Detects", "Total_Infected", "Cumulative_Infections") 
    time_record$Total_Infected <- 1
    time_record$Cumulative_Infections <- 1
    while  (((I < e_thresh) & (I > 0)) | ((I > 0) & (D < d_thres))) {
      #while Number of infected is below epidemic threshold hold and more than 0 infected
      # or while number of infecteds is above 0 and the number of detected is below threshold
            
      ########################### First Introudction Rate
      intro_draws = runif(1) #some possible number of events
      intro_count = sum(intro_draws < intro_rate) #Number of Introductions 
      intro_type_draw = runif(intro_count)
      intro_Symp = sum(intro_type_draw < prob_symp) # Proportion of Symptomatic Introductions
      intro_Asymp = intro_count-intro_Symp # Asym Into = Total Intro-symptomatic 
      
      
      ############################## INCUBATION
      
      leavingInc_draws = runif(incubationInfecteds) # Probabiity of leaving Incubation
      leavingInc_count = sum(leavingInc_draws < incub_p) # Cout of those leaving incubation
      
      ############################## DETERMING ASYMP/SYMP
      newUI_Symp_draws = runif(leavingInc_count) #New Symptomatic Infections based on leaving Inc
      newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
      newUI_Asymp_count = leavingInc_count - newUI_Symp_count # New Asymptomatic Infectiosn leaving Inc
      
      
      ############################## INFECTION
      
      #Infection counts don't matter if they're done by symptomatic or asymptomatic individuals at this point
      # becaue the rates are the same
      
      #Infection for Detected Individuals 
      
      newDI_draws = runif(DI) #New exposed possibilities 
      newDI_count = sum(newDI_draws < prop_p) #Number of new exposed by detected individuals 
      
      #newDI_Symp_draws = runif(newDI_count)
      #newDI_Symp_count = sum(newDI_Symp_draws < prob_symp)
      #newDI_Asymp_count = newDI_count - newDI_Symp_count
      
      # Recovering for Detected Symptomatic 
      recoveredDI_Symp_draws = runif(DI_Symp) # Recovered possibilities for detected symptomatic
      recoveredDI_Symp_count = sum(recoveredDI_Symp_draws < recov_p) # number of recovered DI_symp
      
      # Recovering for Detected Asymptomatic 
      recoveredDI_Asymp_draws = runif(DI_Asymp) #  Recovered possibilities for detected asymptomatic
      recoveredDI_Asymp_count = sum(recoveredDI_Asymp_draws < recov_p)
      
      
      #Innfection for Undetected Individuals 
      newUI_draws = runif(UI) # New exposed probabilities by Undiscovred infecteds 
      newUI_count =  sum(newUI_draws < prop_p)  # Number of new Infecteds 
      
      #newUI_Symp_draws = runif(newUI_count)
      #newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
      #newUI_Asymp_count = newUI_count - newUI_Symp_count
      
      # Recovering for Undetected Symptomatic 
      recoveredUI_Symp_draws = runif(UI_Symp) #Detected indivduals-will they recover?
      recoveredUI_Symp_count = sum(recoveredDI_Symp_draws < recov_p)
      
      # Recovering for Undetected Asymptomatic 
      recoveredUI_Asymp_draws = runif(UI_Asymp) #Detected indivduals-will they recover?
      recoveredUI_Asymp_count = sum(recoveredUI_Asymp_draws < recov_p)
      
      ####################################################################
      # DETECTION
      
      # Detection - Symptomatic 
      discoveredUI_Symp_draws = runif(UI_Symp)  #probability that undetected symptomatic individuals are discovered 
      discoveredUI_Symp_count = sum(discoveredUI_Symp_draws < dis_prob_symp) 
     
      # Detection - Asymptomatic 
      discoveredUI_Asymp_draws = runif(UI_Asymp)  #probability that undetected asymptomatic  individuals are discovered 
      discoveredUI_Asymp_count = sum(discoveredUI_Asymp_draws < dis_prob_asymp) 
      
      # Updating - UI- Undiscovered infecteds can be removed by being recovered or by being discovoered 
      removeUI_Symp = sum((recoveredUI_Symp_draws < recov_p) | (discoveredUI_Symp_draws < dis_prob_symp)) 
      removeUI_Asymp = sum((recoveredUI_Asymp_draws < recov_p) | (discoveredUI_Asymp_draws < dis_prob_asymp)) 
      
      # remove undiscovered infecteds that are both no longer infected and were eventually discovered 
      removeUDI_Symp =  sum((recoveredUI_Symp_draws < recov_p) & ( discoveredUI_Symp_draws < dis_prob_symp))
      removeUDI_Asymp =  sum((recoveredUI_Asymp_draws < recov_p) & ( discoveredUI_Asymp_draws < dis_prob_asymp))
      
      #########################################################
      # UPDATING 
          
      # Incubating
      incubationInfecteds = incubationInfecteds +  newDI_count + newUI_count - leavingInc_count
      # previously incubation + new infections by DI and UI - those that are now infectious
      
      # Newly Exposed 
      newEx = newUI_count + newDI_count
      # newly exposed are those infected by UI and DI 
      
      # Newly Infectious
      #newInf = newUI_Symp_count + newUI_Asymp_count + intro_Asymp + intro_Symp  # or leaving inc+ intro 
      newInf = leavingInc_count + intro_count
      # newly infectious = those that have left incubation + intro counts 
      cumI = cumI + newInf 
          
      # Newly Detected = symptomatic and asymptomatic cases discovered
      newlyDisc = discoveredUI_Asymp_count + discoveredUI_Symp_count
      
      #Cumulatve Detected = those already detected + newly detected
      D = D + newlyDisc 
      
      #Undetected = previously + those introduced + those leaving incuation - recovered or detected 
      UI_Symp = UI_Symp + intro_Symp + newUI_Symp_count - removeUI_Symp
      UI_Asymp = UI_Asymp + intro_Asymp  + newUI_Asymp_count - removeUI_Asymp
      
      UI = UI_Symp + UI_Asymp 
      
      # Detected = previously + discovered - recovered discovered - in one time were discovered/recovered
      DI_Symp = DI_Symp + discoveredUI_Symp_count - recoveredDI_Symp_count - removeUDI_Symp
      DI_Asymp = DI_Asymp + discoveredUI_Asymp_count - recoveredDI_Asymp_count - removeUDI_Asymp
      
      DI = DI_Symp + DI_Asymp
      
      #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered       
      I = UI + DI #Undiscovered and Discovered Infected 
         
    
    #adding time step data 
    time_record <- rbind(time_record, c(newEx, newInf, intro_count, newlyDisc, D , I, cumI))     
  }
  time <- data.frame(seq(1:nrow(time_record)))
  colnames(time) <- "time"
  time_record <- cbind(time, time_record)
  
  return(time_record)
  })
}


# Call Run Branch and Save multiple runs 
# Takes in all the parameters and replicates 
run_branches <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch(...) ) 
  
}

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
#####################################################################

# Analysis functions for extracting various elements from the trials 
last_cuminfect_value <- function(x) {
  row <- tail(x, 1) 
  value <- row[8]
  return(value)
}
all_last_cuminfect_values <- function(x) {
  return(unlist(laply(x, last_cuminfect_value)))
}


last_cumdetect_value <- function(x) {
  x[nrow(x), "Cum_Detects"]
}
all_last_cumdetect_values <- function(x) {
  laply(x, last_cumdetect_value)
}

last_instantInf_value <- function(x) {
  row <- tail(x, 1) 
  value <- row[7]
  return(value)
}
all_last_instantInf_values <- function(x) {
  return(unlist(laply(x, last_instantInf_value)))
}



## Set of functions to calculate given I have X cases, what is the distribution of cases I see 
library(ggplot2)
library(reshape2)
library(scales)

## Set of functions to calculate given I have X cases, what is the distribution of cases I see 
detection_rows <- function(x, detect_thres=d_thres) {
  detections <- x[,6]
  rows <- which(detections == detect_thres)
  reduced <- x[rows,]
  return(reduced)
}

all_detect_rows <- function(x, d_thres) {   # Function to return all rows in trials that match detection threshold
  return(ldply(x, detection_rows))
}



# Additional Post Processing Functions 
prob_ext <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100) {
  escapes = 0
  escapes_prob <- c(escapes)
  i = 1 

  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) #Run the simulation 
    Final.I = record[nrow(record),7]
    Final.D = record[nrow(record),6] ## should 6 be changed to 5? # no it's cumulative detection, instantaneous infectious

    if (Final.D < d_thres)  {
      next
    }
    if (Final.I > e_thresh) {
      escapes = escapes + 1   
    }
    escapes_prob <- c(escapes_prob, escapes/i)
    i = i+1    
  }
  return(escapes_prob)
}



# If you have not crossed the detection threshold, what is the probability that you have crossed your epidemic threshold? Seems to never happen
prob_underD.overE <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=10) {
  escapes = 0
  escapes_prob <- c(escapes)
  i = 1 
  
  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) #Run the simulation 
    Final.I = record[nrow(record),7]
    Final.D = record[nrow(record),6]
    
    if (Final.D < d_thres)  { # Haven't Crossed Detection 
      if (Final.I > e_thresh) { # Is my Instantaneous Number of Infections more than my epidemic threshold
        escapes = escapes + 1   
        escapes_prob <- c(escapes_prob, escapes/i)
        
      }
    }
    i = i+1   
  }
  return(escapes_prob)
}




##### Functions for generating the heat maps 
name.generator <- function(R0, disc_value, type) {
  filename.base <- paste(R0, disc_value, sep = "_")
  filename.type <- paste(type, filename.base, sep = "_")
  return(filename.type)
} 


set.cum.bins <- function(highest.thres, trials) {
  highest.dect <- dect.cases[length(dect.cases)]
  d_thres <- highest.thres
  at.high.dect <- all_detect_rows(trials)
  max.cum <- max(at.high.dect[,8])
  max.bin <- max.cum + 25
  bins <- seq(from = 0, to = max.bin, by = 25)
  return(bins)
}


set.prev.bins <- function(highest.thres, trials) {
  d_thres <- highest.thres
  at.high.dect <- all_detect_rows(trials)
  max.prev <- max(at.high.dect[,7])
  max.bin <- max.prev + 5
  bins <- seq(from = 0, to = max.bin, by = 5)
  return(bins)
}



#Function to take the row and calculate proabilities for bins
bin.frequency <- function(df, bins) {
  df.cut <- cut(df, breaks = bins, right = FALSE)
  df.freq = table(df.cut)
  df.freq <- cbind(df.freq)
  df.freq <- df.freq/sum(df.freq)
  return(df.freq)
}

# Keep Things Ordered Here ording by a certain variable

plotheatmaps <- function(df, type, names, disc_value, R0) {
  df <- cbind(names, df)
  df$names <- factor(df$names, levels = df$names[order(df$names)])
  df.m <- melt(df)
  title <- c("R0 = ", R0, " ; Detect Prob = ", disc_value)
  title <- paste(title, collapse = "")
  
  p <- ggplot(df.m, aes(variable, names)) 
  p <- p +  geom_tile(aes(fill = value), colour = "white") + theme_bw()+ scale_fill_gradient(low = "lightyellow",high = "red", name = "Frequency") +labs(x = type, y = "Detected Cases")
  p <- p + theme(axis.title.x = element_text(size=15), axis.text.x  = element_text(size=9), axis.title.y = element_text(size=15), axis.text.y  = element_text(size=9) ) + ggtitle(title) + theme(plot.title = element_text(lineheight=.8, face="bold"))
  return(p)
  p
}


