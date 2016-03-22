
# Three kinds of counters:
# UI - Undiscovered Infecteds
# DI - Discovered Infecteds
# IncI - Incubation Infecteds
# D - Cumulative number of discoveries
#
# Simulate until we have either:
# - we run out of infecteds
# - e_thresh number of infecteds, and d_thresh number of discoveries

run_branch <- function(prop_p, recov_p, disc_p, incub_p, d_thres, e_thresh) {
  UI = 1; DI = 0; D = 0; IncI = 0
  I = UI + DI
  time_record <- data.frame(matrix(data = c(I,D), ncol = 2))
  colnames(time_record) <- c("I", "D")
  
  while  (((I < e_thresh) & (I > 0)) | ((I > 0) & (D < d_thresh))) {
    #while Number of infected is below epidemic threshold hold and more than 0 infected
    # or while number of infecteds is above 0 and the number of detected is below threshold
    
    ######## Step 1 - Infection ##################
    newDI_draws = runif(DI) #New infected by Detected prop_p 
    newDI_count = sum(newDI_draws < prop_p) 
    
    newUI_draws = runif(UI) #Number infected by Undetected  
    newUI_count =  sum(newUI_draws < prop_p)  
    
    IncI = newUI_count + newDI_count #New Infections
    
    
    ####### Step 2 - Exposed Period to Undetected Infectious ##########
    IncI_draws = reunif(InI) #Probability of Exiting Infectious Period
    ExitI = sum(IncI_draws < incub_p)
    
    
    # Symptomatic vs Asymptomatic
    newSymp_draws = runif(ExitI)
    newSymp_count = sum(newSymp_draws < prob_symp)
    newAsymp_count = ExitI-newSympCount
   
    #Potential Introduction
    intro_prob = reunif(1) #some possible number of events
    intro_num = sum(intro_prob < intro_rate)
    intro_type_draw = reunif(intro_num)
    intro_type_Symp = sum(intro_type_draw < prob_symp)
    intro_type_Asymp = intro_num - intro_type_Symp
    
    #Add new Undetected Infectious Cases and Intro Cases Together 
    UI_Symp = UI_Symp + newSymp_count + intro_type_Symp
    UI_Asymp = UI_Aysmp + newAsymp_count + intro_type_Asymp
  
    
    #### step 3 Detection  ######################### 
    # Detection of new Cases according to Symptomatic / Asymptomatic Probabilities  
    dectSymp_draws = runif(UI_Symp) 
    dectSymp_count = sum(dectSymp_draws < dis_prob_symp) 
    
    dectAsymp_draws = reunif(UI_Asymp)
    dectAsymp_count = sum(dectAsymp_draws < dis_prob_asymp)
    
    DI_Symp = DI_Symp + dectSymp_count
    DI_Asymp = DI_Aymp + dectAsymp_count
    
    NewlyDisc = DI_Symp + DI_Asymp 
    
    #update UI 
    UI_Symp = UI_Symp - detectSymp_count
    UI_Asymp = UI_Asymp - detectAsymp_count
  
    ###### step 4 Recovery and Update ##################
    #Detected 
    recDI_Symp_draws = runif(DI_Symp) #Detected Symptomatic indivduals-will they recover?
    recDI_Symp_count = sum(recDI_Symp_draws < recov_p)
    
    recDI_Asymp_draws = runif(DI_Asymp) #Detected Asymptomatic indivduals-will they recover?
    recDI_Asymp_count = sum(recDI_Asymp_draws < recov_p)
    
    DI_Symp = DI_Symp - recDI_Symp_count
    DI_Asymp = DI_Asymp - recDI_Asymp_count
    
    #Undetected 
    recUI_Symp_draws = runif(UI_Symp) #Undetected Symptomatic indivduals-will they recover?
    recUI_Symp_count = sum(recUI_Symp_draws < recov_p)
    
    recUI_Asymp_draws = runif(UI_Asymp) #Undetected Asymptomatic indivduals-will they recover?
    recUI_Asymp_count = sum(recUI_Asymp_draws < recov_p)
    
    UI_Symp = UI_Symp - recUI_Symp_count
    UI_Asymp = UI_Asymp - recUI_Asymp_count
    
    ##### step 5  Updating final counters ############
    removeUDIR_Asymp = sum((recUI_Asymp_draws < recov_p) & (dectAsymp_draws < dis_prob_asymp))
    removeUDIR_Symp = sum((recUI_Symp_draws < recov_p) & (dectSsymp_draws < dis_prob_symp)) 
    
    #Not sure if I need these....or if it's better to do the double 
   
    UI = UI_Symp  + UI_Asymp 
    DI = DI_Symp + DI_Asymp
    D = D + NewlyDisc
    I = UI + DI #Undiscovered and Discovered Infecteds 
    time_record <- rbind(time_record,c(I,D))
  }
  return(time_record)
}





prob_ext <- function(prop_p, recov_p, disc_p, d_thresh, e_thresh, nsamples=10000) {
  escapes = 0
  i = 1 
  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, disc_p, d_thresh, e_thresh) #Run the simulation 
    Final.I = record[nrow(record),1]
    # print(paste("Final.I", Final.I, sep = "-"))
    Final.D = record[nrow(record),2]
    # print(paste("Final.D", Final.D, sep = "-"))
    if (Final.D < d_thresh)  
      next
    if (Final.I > e_thresh) {
      escapes = escapes + 1
    }
    i = i+1
    print(paste('Estimate', escapes / i, sep = ": "))
  }
}
