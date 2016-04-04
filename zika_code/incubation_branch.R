run_branches_inc <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch_inc(...) ) 
}

run_branches_modified_inc <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch_modified_inc(...) ) 
}

# branch_params <- function(prop_p = 1.7/7 , 
#                           recov_p = 1.0/7,
#                           d_thres = 5,
#                           e_thresh = 200,
#                           prob_symp = 1,
#                           incub_rate = 1/16.5,
#                           dis_prob_symp = 0.0246,
#                           dis_prob_asymp = 0.00 ,
#                           intro_rate = 0.000,
#                           zeroInc_prob = 6.825603e-08)
#   return(as.list(environment()))

# 
# run_branch_modified_inc <- function(params) {
#   with(params,{
#     UI_Symp = 1; UI_Asymp = 0; DI_Symp = 0; DI_Asymp = 0 
#     UI = UI_Symp + UI_Asymp; DI = DI_Symp + DI_Asymp
#     D = 0  
#     incubationInfecteds = 0
#     I = UI + DI
#     cumI = UI
#     time_record <- data.frame(matrix(data = 0, ncol = 7))
#     colnames(time_record) <- c("New_Exposed", "New_Infectious", "Intros", "New_Detections", "Cum_Detects", "Total_Infected", "Cumulative_Infections") 
#     time_record$Total_Infected <- 1
#     time_record$Cumulative_Infections <- 1
#     while  (((cumI < e_thresh) & ((I+incubationInfecteds) > 0)) | (((I+incubationInfecteds) > 0) & (D < d_thres))) {
#       #while Number of infected is below epidemic threshold hold and more than 0 infected
#       # or while number of infecteds is above 0 and the number of detected is below threshold
#       
#       ########################### First Introudction Rate
#       intro_draws = runif(1) #some possible number of events
#       intro_count = sum(intro_draws < intro_rate) #Number of Introductions 
#       intro_type_draw = runif(intro_count)
#       intro_Symp = sum(intro_type_draw < prob_symp) # Proportion of Symptomatic Introductions
#       intro_Asymp = intro_count-intro_Symp # Asym Into = Total Intro-symptomatic 
#       
#       
#       ############################## INCUBATION
#       
#       leavingInc_draws = runif(incubationInfecteds) # Probabiity of leaving Incubation
#       leavingInc_count = sum(leavingInc_draws < incub_rate) # Cout of those leaving incubation
#       
#       ############################## DETERMING ASYMP/SYMP
#       newUI_Symp_draws = runif(leavingInc_count) #New Symptomatic Infections based on leaving Inc
#       newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
#       newUI_Asymp_count = leavingInc_count - newUI_Symp_count # New Asymptomatic Infectiosn leaving Inc
#       
#       
#       ############################## INFECTION
#       
#       #Infection counts don't matter if they're done by symptomatic or asymptomatic individuals at this point
#       # becaue the rates are the same
#       
#       #Infection for Detected Individuals 
#       
#       newDI_draws = runif(DI) #New exposed possibilities 
#       newDI_count = sum(newDI_draws < prop_p) #Number of new exposed by detected individuals 
#       
#       zeroInc_DI_draws = runif(newDI_count)
#       zeroInc_DI_count = sum(zeroInc_DI_draws < zeroInc_prob)
#       
#       
#       #newDI_Symp_draws = runif(newDI_count)
#       #newDI_Symp_count = sum(newDI_Symp_draws < prob_symp)
#       #newDI_Asymp_count = newDI_count - newDI_Symp_count
#       
#       # Recovering for Detected Symptomatic 
#       recoveredDI_Symp_draws = runif(DI_Symp) # Recovered possibilities for detected symptomatic
#       recoveredDI_Symp_count = sum(recoveredDI_Symp_draws < recov_p) # number of recovered DI_symp
#       
#       # Recovering for Detected Asymptomatic 
#       recoveredDI_Asymp_draws = runif(DI_Asymp) #  Recovered possibilities for detected asymptomatic
#       recoveredDI_Asymp_count = sum(recoveredDI_Asymp_draws < recov_p)
#       
#       
#       #Innfection for Undetected Individuals 
#       newUI_draws = runif(UI) # New exposed probabilities by Undiscovred infecteds 
#       newUI_count =  sum(newUI_draws < prop_p)  # Number of new Infecteds 
#       
#       zeroInc_UI_draws = runif(newUI_count)
#       zeroInc_UI_count = sum(zeroInc_UI_draws < zeroInc_prob)
#       
#       #newUI_Symp_draws = runif(newUI_count)
#       #newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
#       #newUI_Asymp_count = newUI_count - newUI_Symp_count
#       
#       # Recovering for Undetected Symptomatic 
#       recoveredUI_Symp_draws = runif(UI_Symp) 
#       recoveredUI_Symp_count = sum(recoveredUI_Symp_draws < recov_p)
#       
#       # Recovering for Undetected Asymptomatic 
#       recoveredUI_Asymp_draws = runif(UI_Asymp) #Detected indivduals-will they recover?
#       recoveredUI_Asymp_count = sum(recoveredUI_Asymp_draws < recov_p)
#       
#       
#       # Determination of zero incubation strains
#       zeroInc_count = zeroInc_UI_count + zeroInc_DI_count
#       zero_Symp_draws = runif(zeroInc_count) #New Symptomatic Infections based on leaving Inc
#       zero_Symp_count = sum(zero_Symp_draws < prob_symp)
#       zero_Asymp_count = zeroInc_count - zero_Symp_count
#       
#       ####################################################################
#       # DETECTION
#       
#       # Detection - Symptomatic 
#       discoveredUI_Symp_draws = runif(UI_Symp)  #probability that undetected symptomatic individuals are discovered 
#       discoveredUI_Symp_count = sum(discoveredUI_Symp_draws < dis_prob_symp) 
#       
#       # Detection - Asymptomatic 
#       discoveredUI_Asymp_draws = runif(UI_Asymp)  #probability that undetected asymptomatic  individuals are discovered 
#       discoveredUI_Asymp_count = sum(discoveredUI_Asymp_draws < dis_prob_asymp) 
#       
#       # Updating - UI- Undiscovered infecteds can be removed by being recovered or by being discovoered 
#       removeUI_Symp = sum((recoveredUI_Symp_draws < recov_p) | (discoveredUI_Symp_draws < dis_prob_symp)) 
#       removeUI_Asymp = sum((recoveredUI_Asymp_draws < recov_p) | (discoveredUI_Asymp_draws < dis_prob_asymp)) 
#       
#       # remove undiscovered infecteds that are both no longer infected and were eventually discovered 
#       removeUDI_Symp =  sum((recoveredUI_Symp_draws < recov_p) & ( discoveredUI_Symp_draws < dis_prob_symp))
#       removeUDI_Asymp =  sum((recoveredUI_Asymp_draws < recov_p) & ( discoveredUI_Asymp_draws < dis_prob_asymp))
#       
#       #########################################################
#       # UPDATING 
#       # Some of newly exposed are those infected by UI and DI will go straight 
#       # Newly Exposed 
#       newEx = newUI_count + newDI_count
#       
#       # Incubating
#       incubationInfecteds = incubationInfecteds +  newDI_count + newUI_count - leavingInc_count - zeroInc_count
#       # previously incubation + new infections by DI and UI - those that are now infectious
#       
#       
#       
#       
#       # Newly Infectious
#       #newInf = newUI_Symp_count + newUI_Asymp_count + intro_Asymp + intro_Symp  # or leaving inc+ intro 
#       newInf = leavingInc_count + intro_count + zeroInc_count
#       # newly infectious = those that have left incubation + intro counts 
#       cumI = cumI + newInf 
#       
#       # Newly Detected = symptomatic and asymptomatic cases discovered
#       newlyDisc = discoveredUI_Asymp_count + discoveredUI_Symp_count
#       
#       #Cumulatve Detected = those already detected + newly detected
#       D = D + newlyDisc 
#       
#       #Undetected = previously + those introduced + those leaving incuation - recovered or detected 
#       UI_Symp = UI_Symp + intro_Symp + newUI_Symp_count - removeUI_Symp + zero_Symp_count
#       UI_Asymp = UI_Asymp + intro_Asymp  + newUI_Asymp_count - removeUI_Asymp + zero_Asymp_count
#       
#       UI = UI_Symp + UI_Asymp 
#       
#       # Detected = previously + discovered - recovered discovered - in one time were discovered/recovered
#       DI_Symp = DI_Symp + discoveredUI_Symp_count - recoveredDI_Symp_count - removeUDI_Symp
#       DI_Asymp = DI_Asymp + discoveredUI_Asymp_count - recoveredDI_Asymp_count - removeUDI_Asymp
#       
#       DI = DI_Symp + DI_Asymp
#       
#       #Detected = Detected + Recovered Individuals + Newly Detected - Previously UI that have been detected and Recovered       
#       I = UI + DI #Undiscovered and Discovered Infected 
#       
#       
#       #adding time step data 
#       time_record <- rbind(time_record, c(newEx, newInf, intro_count, newlyDisc, D , I, cumI))     
#     }
#     time <- data.frame(seq(1:nrow(time_record)))
#     colnames(time) <- "time"
#     time_record <- cbind(time, time_record)
#     
#     return(time_record)
#   })
# }


run_branch_modified_inc <- function(params) {
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
    while  (((cumI < e_thresh) & ((I+incubationInfecteds) > 0)) | (((I+incubationInfecteds) > 0) & (D < d_thres))) {
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
      leavingInc_count = sum(leavingInc_draws < incub_rate) # Cout of those leaving incubation
      
      ############################## DETERMING ASYMP/SYMP
      newUI_Symp_draws = runif(leavingInc_count) #New Symptomatic Infections based on leaving Inc
      newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
      newUI_Asymp_count = leavingInc_count - newUI_Symp_count # New Asymptomatic Infectiosn leaving Inc
      
      
      ############################## INFECTION
      
      #Infection counts don't matter if they're done by symptomatic or asymptomatic individuals at this point
      # becaue the rates are the same
      
      #Infection for Detected Individuals 
      
      newDI_count = sum(rpois(DI, lambda = prop_p)) #New exposed possibilities 
      
      zeroInc_DI_draws = runif(newDI_count)
      zeroInc_DI_count = sum(zeroInc_DI_draws < zeroInc_prob)
      
      
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
      newUI_count = sum(rpois(UI, lambda = prop_p)) # New exposed probabilities by Undiscovred infecteds 
      
      
      zeroInc_UI_draws = runif(newUI_count)
      zeroInc_UI_count = sum(zeroInc_UI_draws < zeroInc_prob)
      
      #newUI_Symp_draws = runif(newUI_count)
      #newUI_Symp_count = sum(newUI_Symp_draws < prob_symp)
      #newUI_Asymp_count = newUI_count - newUI_Symp_count
      
      # Recovering for Undetected Symptomatic 
      recoveredUI_Symp_draws = runif(UI_Symp) 
      recoveredUI_Symp_count = sum(recoveredUI_Symp_draws < recov_p)
      
      # Recovering for Undetected Asymptomatic 
      recoveredUI_Asymp_draws = runif(UI_Asymp) #Detected indivduals-will they recover?
      recoveredUI_Asymp_count = sum(recoveredUI_Asymp_draws < recov_p)
      
      
      # Determination of zero incubation strains
      zeroInc_count = zeroInc_UI_count + zeroInc_DI_count
      zero_Symp_draws = runif(zeroInc_count) #New Symptomatic Infections based on leaving Inc
      zero_Symp_count = sum(zero_Symp_draws < prob_symp)
      zero_Asymp_count = zeroInc_count - zero_Symp_count
      
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
      # Some of newly exposed are those infected by UI and DI will go straight 
      # Newly Exposed 
      newEx = newUI_count + newDI_count
      
      # Incubating
      incubationInfecteds = incubationInfecteds +  newDI_count + newUI_count - leavingInc_count - zeroInc_count
      # previously incubation + new infections by DI and UI - those that are now infectious
      
      
      
      
      # Newly Infectious
      #newInf = newUI_Symp_count + newUI_Asymp_count + intro_Asymp + intro_Symp  # or leaving inc+ intro 
      newInf = leavingInc_count + intro_count + zeroInc_count
      # newly infectious = those that have left incubation + intro counts 
      cumI = cumI + newInf 
      
      # Newly Detected = symptomatic and asymptomatic cases discovered
      newlyDisc = discoveredUI_Asymp_count + discoveredUI_Symp_count
      
      #Cumulatve Detected = those already detected + newly detected
      D = D + newlyDisc 
      
      #Undetected = previously + those introduced + those leaving incuation - recovered or detected 
      UI_Symp = UI_Symp + intro_Symp + newUI_Symp_count - removeUI_Symp + zero_Symp_count
      UI_Asymp = UI_Asymp + intro_Asymp  + newUI_Asymp_count - removeUI_Asymp + zero_Asymp_count
      
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
    
    while  (((cumI < e_thresh) & ((I+sum(incubationInfecteds)) > 0)) | 
            (((I+sum(incubationInfecteds)) > 0) & (cumD < d_thres))) {
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

