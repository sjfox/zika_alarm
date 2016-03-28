run_branches_inc <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch_modified_inc(...) ) 
}

branch_params <- function(prop_p = 1.7/7 , 
                          recov_p = 1.0/7,
                          d_thres = 5,
                          e_thresh = 200,
                          prob_symp = 1,
                          incub_rate = 1/16.5,
                          dis_prob_symp = 1,
                          dis_prob_asymp = 0.00 ,
                          intro_rate = 0.000,
                          zeroInc_prob = 6.825603e-08)
  return(as.list(environment()))


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
    while  (((I < e_thresh) & ((I+incubationInfecteds) > 0)) | (((I+incubationInfecteds) > 0) & (D < d_thres))) {
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
      
      newDI_draws = runif(DI) #New exposed possibilities 
      newDI_count = sum(newDI_draws < prop_p) #Number of new exposed by detected individuals 
      
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
      newUI_draws = runif(UI) # New exposed probabilities by Undiscovred infecteds 
      newUI_count =  sum(newUI_draws < prop_p)  # Number of new Infecteds 
      
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

#trials.incubation <-  run_branches_inc(num_reps = 1000, branch_params(prop_p = 1.7/7, e_thresh = 300, incub_rate = 1, zeroInc_prob = 1))
