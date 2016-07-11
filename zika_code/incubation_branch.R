run_branches_inc <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch_inc(...) ) 
}

run_branches_inc_old <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch_inc_old(...) ) 
}

update_state <- function(n, prob){
  ## Draws n uniform random numbers
  ## Returns a sum of all numbers less than the probability
  sum(runif(n = n) < prob)
}

update_states <- Vectorize(update_state, vectorize.args = c("n"))

track_infectious <- function(UI, recov_prob, disc_prob){
  ## Issues with UI and DI movement and boxes
  ## Takes in UI vectors
  ## Each vector is same length, and contains sum individuals
  ## in each box
  ## Outputs a list:
  ## UI leaving each box
  ## DI to add from UI recoverin
  ##    DI to add is vector length DI+1
  ##    The last item is the number that should 
  ##    recover but count as detected
  
  UI_progressing = vector("numeric", length(UI))  ## The number of UI moving to next UI box
  UI_detecting = vector("numeric", length(UI))    ## The number of UI that were detected
  DI_adding = vector("numeric", length(UI)+1)     ## The number of UI that are detected by stage of recovery
  
  for(box in 1:length(UI)){
    UI_recovery_draws = runif(UI[box])
    UI_detected_draws = runif(UI[box])
    
    ## Get the total numbers recovering, detected, and (detected and recovered)
    UI_recovering = sum(UI_recovery_draws < recov_prob)
    UI_detected = sum(UI_detected_draws < disc_prob)
    UI_detected_and_recovered = sum( (UI_detected_draws < disc_prob) & (UI_recovery_draws < recov_prob) )
    
    ## Keep track of total leaving the current UI box into another UI box
    UI_progressing[box] = UI_recovering - UI_detected_and_recovered
    
    ## Neeed to keep track of those leaving due to detection
    UI_detecting[box] = UI_detected
    
    ## Now update those entering the DI class
    ## Those that didn't recover, but were detected, are added to current DI box
    DI_adding[box] = DI_adding[box] + UI_detected - UI_detected_and_recovered    
    
    ## Those that recovered and detected are added to next DI class
    DI_adding[box+1] = DI_adding[box+1] + UI_detected_and_recovered    
  }
  
  return(list(UI_progressing=UI_progressing, UI_detecting=UI_detecting, DI_adding=DI_adding))
}



run_branch_inc <- function(params) {
  with(params,{
    ## The 4 infected classes all have that number of boxes
    UI_Symp = vector("numeric", length=infBoxes) 
    UI_Asymp = vector("numeric", length=infBoxes)
    DI_Symp = vector("numeric", length=infBoxes) 
    DI_Asymp = vector("numeric", length=infBoxes) 
    
    ## Classes for storing the introduced infectious
    UI_Intro_Symp = vector("numeric", length=infBoxes) 
    UI_Intro_Asymp = vector("numeric", length=infBoxes)
    DI_Intro_Symp = vector("numeric", length=infBoxes) 
    DI_Intro_Asymp = vector("numeric", length=infBoxes) 
    
    ## Initialize the epidemic with one undetected infectious
    UI_Intro_Symp[1] = 1
    
    ## UI and DI keep track of the total number of undetected and detected
    ## cases at any given time
    UI = sum(UI_Symp, UI_Asymp, UI_Intro_Symp, UI_Intro_Asymp) 
    DI = sum(DI_Symp, DI_Asymp, DI_Intro_Symp, DI_Intro_Asymp)
    
    ## Only those that are introduced
    UI_Intro = sum(UI_Intro_Symp, UI_Intro_Asymp) 
    DI_Intro = sum(DI_Intro_Symp, DI_Intro_Asymp) 
    
    ## I keeps track of total infectious
    I = UI + DI 
    
    ## Only introduced ones
    I_Intro = UI_Intro + DI_Intro
    
    cumD = 0  ## Cum detected
    cumD_Intro = 0  ## Cum detected of introduced
    
    ## Keeps track of all exposed classes
    incubationInfecteds = vector("numeric", length = incBoxes)
    
    ## Cumulative incidence
    cumI = I
    cumI_Intro = I_Intro
    
    ## Data.frame for holding time series
    time_record = data.frame(matrix(data = 0, ncol = 11))
    colnames(time_record) = c("New_Exposed", "New_Infection", "Total_Infections", "Cumulative_Infections", 
                              "New_Intro", "Total_Intro_Infections", "Cumulative_Intro_Infections", 
                              "New_Detections", "Cum_Detections", "New_Intro_Detections", "Cum_Intro_Detections") 
    
    ## Initialize appropriately
    time_record$New_Infection = time_record$Total_Infections = time_record$Cumulative_Infections = I
    
    while  (( (cumI-cumI_Intro) < e_thresh) & ((I+sum(incubationInfecteds)) > 0)) {
      #while Number of infected is below epidemic threshold and more than 0 infected
      # or while number of infecteds is above 0 and the number of detected is below threshold
      
      ########################### First Introduction Rate
      ## intro_rate is rate of introductions, drawn from poisson distribution daily
      intro_count = rpois(n = 1, lambda = intro_rate) 
      
      ## get total introductions that are symptomatic
      intro_symp = update_state(n=intro_count, prob = prob_symp)
      intro_asymp = intro_count-intro_symp 
      
      ############################## INCUBATION 
      ## Keep track of transitions out of each incubation box
      leavingInc = update_states(n = incubationInfecteds, prob = incub_rate)
    
      ## DETERMING ASYMP/SYMP of those leaving last incubation box (becoming infectious)
      newUI_Symp = update_state(n = leavingInc[incBoxes], prob = prob_symp)
      newUI_Asymp = leavingInc[incBoxes] - newUI_Symp
      
      ############################## INFECTION
      #Calculate the number of new exposed from all infectious
      ## Add in second line if decide that symp/asymp should have different transmission rates
      ## potentially detected/undetected as well.
      newExposed = sum(rpois(n = I, lambda = prop_p))
      
      ############################## INFECTIOUS CLASSES MOVEMENT
      # Infectious classes movement Undetected infectious
      ## Can progress in UI boxes, or get transferred
      ## to DI boxes
      infectious_update_Symp = track_infectious(UI = UI_Symp, recov_prob = recov_p, disc_prob = dis_prob_symp)
      infectious_update_Asymp = track_infectious(UI = UI_Asymp, recov_prob = recov_p, disc_prob = dis_prob_asymp)
      
      ## Update introduction infectious classes
      infectious_intro_update_Symp = track_infectious(UI = UI_Intro_Symp, recov_prob = recov_p, disc_prob = dis_prob_symp)
      infectious_intro_update_Asymp = track_infectious(UI = UI_Intro_Asymp, recov_prob = recov_p, disc_prob = dis_prob_asymp)
    
      ## DI movement (can only progress towards recovery)
      DI_symp_leaving <- update_states(n=DI_Symp, prob=recov_p)
      DI_asymp_leaving <- update_states(n=DI_Asymp, prob=recov_p)

      DI_intro_symp_leaving <- update_states(n=DI_Intro_Symp, prob=recov_p)
      DI_intro_asymp_leaving <- update_states(n=DI_Intro_Asymp, prob=recov_p)
      
      #########################################################
      # UPDATING 
      # Incubation classes
      for(i in 1:incBoxes){
        if(i==1){
          incubationInfecteds[i] = incubationInfecteds[i] + newExposed - leavingInc[i]
        } else {
          incubationInfecteds[i] = incubationInfecteds[i] + leavingInc[i-1] - leavingInc[i]
        }
      }
      
      # UI_leaving=UI_leaving, DI_adding=DI_adding
      # Infectious classes updating
      for(i in 1:infBoxes){
        if(i==1){
          ## Update the UI first box
          UI_Symp[i] = UI_Symp[i] +  newUI_Symp - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + newUI_Asymp - infectious_update_Asymp$UI_progressing[i] -  infectious_update_Asymp$UI_detecting[i]

        } else {
          ## Update the other UI boxes
          UI_Symp[i] = UI_Symp[i] + infectious_update_Symp$UI_progressing[i-1] - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + infectious_update_Asymp$UI_progressing[i-1] - infectious_update_Asymp$UI_progressing[i]-  infectious_update_Asymp$UI_detecting[i]
          
        }
        
        ## Update the DI boxes (all the same process, so not in conditionals)
        DI_Symp[i] = DI_Symp[i] + infectious_update_Symp$DI_adding[i] - DI_symp_leaving[i]
        DI_Asymp[i] = DI_Asymp[i] + infectious_update_Asymp$DI_adding[i] - DI_asymp_leaving[i]
      }
      
      for(i in 1:infBoxes){
        if(i==1){
          ## Update the UI first box
          UI_Intro_Symp[i] = UI_Intro_Symp[i] + intro_symp - infectious_intro_update_Symp$UI_progressing[i] -  infectious_intro_update_Symp$UI_detecting[i]
          UI_Intro_Asymp[i] = UI_Intro_Asymp[i] + intro_asymp +  infectious_intro_update_Asymp$UI_progressing[i] -  infectious_intro_update_Asymp$UI_detecting[i]
          
        } else {
          ## Update the other UI boxes
          UI_Intro_Symp[i] = UI_Intro_Symp[i] + infectious_intro_update_Symp$UI_progressing[i-1] - infectious_intro_update_Symp$UI_progressing[i] -  infectious_intro_update_Symp$UI_detecting[i]
          UI_Intro_Asymp[i] = UI_Intro_Asymp[i] + infectious_intro_update_Asymp$UI_progressing[i-1] - infectious_intro_update_Asymp$UI_progressing[i]-  infectious_intro_update_Asymp$UI_detecting[i]
          
        }
        
        ## Update the DI boxes (all the same process, so not in conditionals)
        DI_Intro_Symp[i] = DI_Intro_Symp[i] + infectious_intro_update_Symp$DI_adding[i] - DI_intro_symp_leaving[i]
        DI_Intro_Asymp[i] = DI_Intro_Asymp[i] + infectious_intro_update_Asymp$DI_adding[i] - DI_intro_asymp_leaving[i]
      }
      
      
      # Newly Infectious
      newInf = leavingInc[incBoxes] + intro_symp + intro_asymp 
      newInf_Intro = intro_symp + intro_asymp 
      
      # Update cumulative infectious
      cumI = cumI + newInf 
      cumI_Intro = cumI_Intro + newInf_Intro 
      
      # Newly Detected = symptomatic and asymptomatic cases discovered
      newlyDisc = sum(infectious_update_Symp$DI_adding, infectious_update_Asymp$DI_adding, infectious_intro_update_Symp$DI_adding, infectious_intro_update_Asymp$DI_adding)
      newlyDisc_Intro = sum(infectious_intro_update_Symp$DI_adding, infectious_intro_update_Asymp$DI_adding)
      
      #Cumulatve Detected = those already detected + newly detected
      cumD = cumD + newlyDisc 
      cumD_Intro = cumD_Intro + newlyDisc_Intro
      
      UI = sum(UI_Symp, UI_Asymp, UI_Intro_Symp, UI_Intro_Asymp) 
      DI = sum(DI_Symp, DI_Asymp, DI_Intro_Symp, DI_Intro_Asymp)
      
      ## Only those that are introduced
      UI_Intro = sum(UI_Intro_Symp, UI_Intro_Asymp) 
      DI_Intro = sum(DI_Intro_Symp, DI_Intro_Asymp) 
      
      #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered       
      I = UI + DI #Undiscovered and Discovered Infected 
      I_Intro = UI_Intro + DI_Intro

      #adding time step data 
      time_record <- rbind(time_record, c(newExposed, newInf, I, cumI, newInf_Intro, I_Intro, cumI_Intro, newlyDisc, cumD, newlyDisc_Intro, cumD_Intro))     
    }
    time <- data.frame(time=seq(1:nrow(time_record)))
    time_record <- cbind(time, time_record)
    
    return(time_record)
  })
}

run_branch_inc_old <- function(params) {
  with(params,{
    ## The 4 infected classes all have that number of boxes
    UI_Symp = vector("numeric", length=infBoxes) 
    UI_Asymp = vector("numeric", length=infBoxes)
    DI_Symp = vector("numeric", length=infBoxes) 
    DI_Asymp = vector("numeric", length=infBoxes) 
    
    ## Initialize the epidemic with one undetected infectious
    UI_Symp[1] = 1
    
    ## UI and DI keep track of the total number of undetected and detected
    ## cases at any given time
    UI = sum(UI_Symp, UI_Asymp) 
    DI = sum(DI_Symp, DI_Asymp)
    
    ## I keeps track of total infectious
    I = UI + DI
    
    cumD = 0  ## Cum detected
    
    ## Keeps track of all exposed classes
    incubationInfecteds = vector("numeric", length = incBoxes)
    
    ## Cumulative incidence
    cumI = UI + DI
    
    ## Data.frame for holding time series
    time_record = data.frame(matrix(data = 0, ncol = 7))
    colnames(time_record) = c("New_Exposed", "New_Infectious", "Intros", "New_Detections", "Cum_Detects", "Total_Infected", "Cumulative_Infections") 
    
    ## Initialize appropriately
    time_record$Total_Infected = I
    time_record$New_Infectious = I
    time_record$Cumulative_Infections = cumI
    
    while  ((cumI < e_thresh) & ((I+sum(incubationInfecteds)) > 0)) {
      #while Number of infected is below epidemic threshold and more than 0 infected
      # or while number of infecteds is above 0 and the number of detected is below threshold
      
      ########################### First Introduction Rate
      ## intro_rate is rate of introductions, drawn from poisson distribution daily
      intro_count = rpois(n = 1, lambda = intro_rate) 
      
      ## get total introductions that are symptomatic
      intro_symp = update_state(n=intro_count, prob = prob_symp)
      intro_asymp = intro_count-intro_symp 
      
      ############################## INCUBATION 
      ## Keep track of transitions out of each incubation box
      leavingInc = update_states(n = incubationInfecteds, prob = incub_rate)
      
      ## DETERMING ASYMP/SYMP of those leaving last incubation box (becoming infectious)
      newUI_Symp = update_state(n = leavingInc[incBoxes], prob = prob_symp)
      newUI_Asymp = leavingInc[incBoxes] - newUI_Symp
      
      ############################## INFECTION
      #Calculate the number of new exposed from all infectious
      ## Add in second line if decide that symp/asymp should have different transmission rates
      ## potentially detected/undetected as well.
      newExposed = sum(rpois(n = I, lambda = prop_p))
      
      ############################## INFECTIOUS CLASSES MOVEMENT
      # Infectious classes movement Undetected infectious
      ## Can progress in UI boxes, or get transferred
      ## to DI boxes
      infectious_update_Symp = track_infectious(UI = UI_Symp, recov_prob = recov_p, disc_prob = dis_prob_symp)
      infectious_update_Asymp = track_infectious(UI = UI_Asymp, recov_prob = recov_p, disc_prob = dis_prob_asymp)
      
      ## DI movement (can only progress towards recovery)
      DI_symp_leaving <- update_states(n=DI_Symp, prob=recov_p)
      DI_asymp_leaving <- update_states(n=DI_Asymp, prob=recov_p)
      
      #########################################################
      # UPDATING 
      # Incubation classes
      for(i in 1:incBoxes){
        if(i==1){
          incubationInfecteds[i] = incubationInfecteds[i] + newExposed - leavingInc[i]
        } else {
          incubationInfecteds[i] = incubationInfecteds[i] + leavingInc[i-1] - leavingInc[i]
        }
      }
      
      # UI_leaving=UI_leaving, DI_adding=DI_adding
      # Infectious classes updating
      for(i in 1:infBoxes){
        if(i==1){
          ## Update the UI first box
          UI_Symp[i] = UI_Symp[i] + intro_symp + newUI_Symp - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + intro_asymp + newUI_Asymp - infectious_update_Asymp$UI_progressing[i] -  infectious_update_Asymp$UI_detecting[i]
          
        } else {
          ## Update the other UI boxes
          UI_Symp[i] = UI_Symp[i] + infectious_update_Symp$UI_progressing[i-1] - infectious_update_Symp$UI_progressing[i] -  infectious_update_Symp$UI_detecting[i]
          UI_Asymp[i] = UI_Asymp[i] + infectious_update_Asymp$UI_progressing[i-1] - infectious_update_Asymp$UI_progressing[i]-  infectious_update_Asymp$UI_detecting[i]
          
        }
        
        
        ## Update the DI boxes (all the same process, so not in conditionals)
        DI_Symp[i] = DI_Symp[i] + infectious_update_Symp$DI_adding[i] - DI_symp_leaving[i]
        DI_Asymp[i] = DI_Asymp[i] + infectious_update_Asymp$DI_adding[i] - DI_asymp_leaving[i]
      }
      
      
      # Newly Infectious
      newInf = leavingInc[incBoxes] + intro_symp + intro_asymp 
      
      # Update cumulative infectious
      cumI = cumI + newInf 
      
      # Newly Detected = symptomatic and asymptomatic cases discovered
      newlyDisc = sum(infectious_update_Symp$DI_adding, infectious_update_Asymp$DI_adding)
      
      #Cumulatve Detected = those already detected + newly detected
      cumD = cumD + newlyDisc 
      
      UI = sum(UI_Symp, UI_Asymp)
      
      # Detected = previously + discovered - recovered discovered - in one time were discovered/recovered
      DI = sum(DI_Symp, DI_Asymp)
      
      #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered       
      I = UI + DI #Undiscovered and Discovered Infected 
      
      
      #adding time step data 
      time_record <- rbind(time_record, c(newExposed, newInf, intro_count, newlyDisc, cumD , I, cumI))     
    }
    time <- data.frame(seq(1:nrow(time_record)))
    colnames(time) <- "time"
    time_record <- cbind(time, time_record)
    
    return(time_record)
  })
}
