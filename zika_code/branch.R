
# Three kinds of counters:
# UI - Undiscovered Infecteds
# DI - Discovered Infecteds
# D - Cumulative number of discoveries
#
# Simulate until we have either:
# - we run out of infecteds
# - e_thresh number of infecteds, and d_thresh number of discoveries
run_branch <- function(prop_p, recov_p, disc_p, d_thres, e_thresh) {
  UI = 1; DI = 0; D = 0
  I = UI + DI
  time_record <- data.frame(matrix(data = c(I,D), ncol = 2))
  colnames(time_record) <- c("I", "D")
  
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
    time_record <- rbind(time_record,c(I,D))
  }
  return(time_record)
}


# Estimate the probability we will reach e_thresh number of I's (as opposed to
# the infection dying out) given we have discovered at least d_thresh cases
#
# prop_p -- probability an I infects a new individual in a time period
# recov_p -- probability an I recovers in a time period
# disc_p -- probability we discover an I
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative)

prob_ext <- function(prop_p, recov_p, disc_p, d_thresh, e_thresh, nsamples=10000) {
  escapes = 0
  i = 1 
  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, disc_p, d_thresh, e_thresh) #Run the simulation 
    Final.I = record[nrow(record),1]
    # print(paste("Final.I", Final.I, sep = "-"))
    Final.D = record[nrow(record),2]
    # print(paste("Final.D", Final.D, sep = "-"))
    if (Final.D < d_thresh) {
      next 
    }
    if (Final.I > e_thresh){  
      escapes = escapes + 1
    }
    i = i+1
    print(paste('Estimate', escapes / i, sep = ": "))
  }
}

# prop_p -- probability an I infects a new individual in a time period: how to determine this? 
# recov_p -- probability an I recovers in a time period: human recovery 
# disc_p -- probability we discover an I: symptom ratio 
# d_thresh -- unused... number of discoveries
# e_thresh -- total instantaneous I's that count as epidemic escape (not cumulative):

prop_p <- 2/7
recov_p <- 1.0/7
disc_p <- .01
d_thresh <- 5
e_thresh <- 150

prob_ext(prop_p, recov_p, disc_p, d_thresh, e_thresh)
  
