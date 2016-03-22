
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

run_branch <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) {
  UI = 1; DI = 0;  D = 0; newI = 0; IncI = 0 
  UI_Symp = 1; UI_Asymp = 0; DI_Symp = 0; DI_Asymp = 0 
  CurrentInfecteds = UI + DI
  I = CurrentInfecteds
  
  time_record <- data.frame(matrix(data = 0, ncol = 7))
  colnames(time_record) <- c("New_Exposed", "New_Infectious", "Intros", "New_Detections", "Cum_Detects", "Total_Infected", "Cumulative_Infections") 
  time_record$Total_Infected <- 1
  time_record$Cumulative_Infections <- 1
  
  while  ( ((CurrentInfecteds > 0) & (D < d_thres))   |   ((CurrentInfecteds < e_thresh) & (CurrentInfecteds > 0)) ) {
    
    #i = 1
    #for(i in 1:100) { 
    #while Number of infected is below epidemic threshold hold and more than 0 infected
    # or while number of infecteds is above 0 and the number of detected is below threshold
    
    ######## Step 1 - Infection ##################
    #Could potentially have 4 different rates here , but assuming that transmissibility determined by detected or not 
    newDI_draws = runif(DI) #New infected by detected prop_p 
    newDI_count = sum(newDI_draws < prop_p) # For Right now prop_p is the same re
    
    newUI_draws = runif(UI) #Number infected by Undetected  
    newUI_count =  sum(newUI_draws < prop_p)  
    
    newI = newUI_count + newDI_count #New Infections (Exposed)
    
    ####### Step 2 - Exposed Period to Undetected Infectious ##########   
    IncI_draws = runif(IncI) #Probability of Exiting Infectious Period
    ExitI = sum(IncI_draws < incub_p) #Num Infectious 
    
    
    # Symptomatic vs Asymptomatic
    newSymp_draws = runif(ExitI)
    newSymp_count = sum(newSymp_draws < prob_symp)
    newAsymp_count = ExitI-newSymp_count
    
    #Potential Introduction
    intro_prob = runif(1) #some possible number of events
    intro_num = sum(intro_prob < intro_rate) #Number of Introductions 
    intro_type_draw = runif(intro_num)
    intro_type_Symp = sum(intro_type_draw < prob_symp)
    intro_type_Asymp = intro_num - intro_type_Symp
    
    
    #### step 3 Detection  ######################### 
    # Detection of new Cases according to Symptomatic / Asymptomatic Probabilities  
    dectSymp_draws = runif(UI_Symp) 
    dectSymp_count = sum(dectSymp_draws < dis_prob_symp) 
    
    dectAsymp_draws = runif(UI_Asymp)
    dectAsymp_count = sum(dectAsymp_draws < dis_prob_asymp)
    
    NewlyDisc = dectSymp_count + dectAsymp_count # Number of New Detecteds 
    
    
    ###### step 4 Recovery  ##################
    #Detected 
    recDI_Symp_draws = runif(DI_Symp) #Detected Symptomatic indivduals-will they recover?
    recDI_Symp_count = sum(recDI_Symp_draws < recov_p)
    
    recDI_Asymp_draws = runif(DI_Asymp) #Detected Asymptomatic indivduals-will they recover?
    recDI_Asymp_count = sum(recDI_Asymp_draws < recov_p)
    
    
    #Undetected 
    recUI_Symp_draws = runif(UI_Symp) #Undetected Symptomatic indivduals-will they recover?
    recUI_Symp_count = sum(recUI_Symp_draws < recov_p)
    
    recUI_Asymp_draws = runif(UI_Asymp) #Undetected Asymptomatic indivduals-will they recover?
    recUI_Asymp_count = sum(recUI_Asymp_draws < recov_p)
    
    ##### step 5  Updating final counters and time record ############      
    IncI = IncI + newI - ExitI # Incubated are the number there, plus new, minus those that left 
    
    
    removeUI_Symp = sum((recUI_Symp_draws  < recov_p) | ( dectSymp_draws < dis_prob_symp) )
    removeUDI_Symp = sum((recUI_Symp_draws  < recov_p) & ( dectSymp_draws < dis_prob_symp) )
    
    UI_Symp = UI_Symp + newSymp_count + intro_type_Symp - removeUI_Symp
    DI_Symp = DI_Symp + dectSymp_count - removeUDI_Symp
    
    removeUI_Asymp = sum(( recUI_Asymp_draws  < recov_p) | ( dectAsymp_draws < dis_prob_asymp) )
    removeUDI_Asymp = sum(( recUI_Asymp_draws  < recov_p) & ( dectAsymp_draws < dis_prob_asymp) )
     
    UI_Asymp = UI_Asymp + newAsymp_count + intro_type_Asymp - removeUI_Asymp
    DI_Asymp = DI_Asymp + dectAsymp_count - removeUDI_Asymp

    
    UI = UI_Symp + UI_Asymp 
    DI = DI_Symp + DI_Asymp
   
    D = D + NewlyDisc # Cumulative Detected 
    
    CurrentInfecteds = UI + DI #Undiscovered and Discovered Infecteds 
    if (CurrentInfecteds < 0) CurrentInfecteds = 0
    I = I + ExitI
    
    #adding time step data 
    time_record <- rbind(time_record, c(newI, ExitI, intro_num, NewlyDisc, D , CurrentInfecteds, I))     
    # i = i + 1
  }
  time <- data.frame(seq(1:nrow(time_record)))
  colnames(time) <- "time"
  time_record <- cbind(time, time_record)
  return(time_record)
}



# Call Run Branch and Save multiple runs 
# Takes in all the parameters and replicates 
run_branches <- function(num_reps, ...) {
  rlply(.n = num_reps, .expr = run_branch(...) ) 
}








# Post Processing Functions 


# analysis function









prob_ext <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100) {
  escapes = 0
  i = 1 

  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) #Run the simulation 
    Final.I = record[nrow(record),7]
    Final.D = record[nrow(record),6]

    if (Final.D < d_thres)  {
      next
    }
    if (Final.I > e_thresh) {
      escapes = escapes + 1   
    }
    print(paste('Estimate', escapes/i , sep = ": "))
    i = i+1    
  }
}
