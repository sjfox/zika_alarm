# Code I've written just testing things out-not worth saving in the main files 
)
last_infected <- all_last_instantInf_values(trials)
end_cum_infected <- all_last_cuminfect_values(trials)
last_cum_detected <- all_last_cumdetect_values(trials)

sort(last_cum_detected)
hist(end_cum_infected)

detect.vs.cumulative <- cbind(last_cum_detected, end_cum_infected)


#Plotting Preliminary Results
lastvalues <- all_last_cuminfect_values(trials)
hist(lastvalues, main = paste( "1000 Runs, R0 =", prop_p/recov_p, sep = ""), breaks = 10, ylim = c(1,1000))

plot(last_infected, lastvalues)
last_detect_values <- all_last_cumdetect_values(trials)
hist(last_detect_values)
plot(last_detect_values, lastvalues,  xlab = "Cumulative Detected", ylab =  "Cumulative Infected", main = "1000 Runs-Running until Instantaneous Infected > e_thres")



## Writing Functions for Plotting Heat Maps

dect.cases <- c(0:20)
prop_range <- seq(from = 0.8/7, to = .9/7, by = 0.1/7) #Changed to =  2.0/7 for testing
disc_range <- seq(from = 0.01, to = 0.03, by = .01) # Changed discover range to - 0.10



 
for (m in 1:length(prop_range)) {
  print(paste("Starting Prop Range", m, "-"))
  
  for (j in 1:length(disc_range)) {
    
    R0  = prop_range[m]*7
    
    #Running trials
    trials <- run_branches(num_reps = 1000, branch_params(dis_prob_symp=disc_range[j], prop_p = prop_range[m]))
    print(paste("Finished Trials", j, sep = "-"))
    
    #setting up bins to calculate frequencies 
    d_thres <- dect.cases[length(dect.cases)] #Highest Number of Cases to Consider 
    bins.prev <- set.prev.bins(highest.thres, trials)
    bins.cumulative <- set.cum.bins(highest.thres, trials)
    
    #resetting d_thres for trials

    
    
    #Setting Up matrices 
    thres.matrix.prev <- data.frame(matrix(nrow = length(dect.cases), ncol = length(bins.prev)-1))
    colnames(thres.matrix.prev) <- paste("< ", bins.prev[2:length(bins.prev)]); rownames(thres.matrix.prev) <- paste("Detected Cases = ", dect.cases)
    
    thres.matrix.cum <- data.frame(matrix(nrow = length(dect.cases), ncol = length(bins.cumulative)-1))
    colnames(thres.matrix.cum) <- paste("<", bins.cumulative[2:length(bins.cumulative)]); rownames(thres.matrix.cum) <- paste("Detected Cases = ", dect.cases)
    
    
    #Writing the values 
    for (i in 1:length(dect.cases)) {
      d_thres <- dect.cases[i]
      print(d_thres)
      dataframe <- all_detect_rows(trials) # takes already whatever the current thresholds 
      frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
      frequencies.cum <- bin.frequency(dataframe[,8], bins.cumulative)
      
      thres.matrix.prev[i,] <- frequencies.prev
      thres.matrix.cum[i,] <- frequencies.cum
    }
    print("Calculated Frequencies")
    
    #Writing and saving the files - Not necessary to save all of these for now to check 
    filename.prev <- paste(name.generator(R0, disc_range[j], "Prev"), "csv", sep = ".")
    write.csv(x = thres.matrix.prev, file = filename.prev, row.names = TRUE)
    
    filename.cum <- paste(name.generator(R0, disc_range[j], "Cumulative"), "csv", sep = ".")
    write.csv(x = thres.matrix.cum, file = filename.cum, row.names = TRUE)
    
    # Plotting Function to go here
    prev.map <- plotheatmaps(thres.matrix.prev, type = "Prevalence", names = as.factor(dect.cases), R0 = R0, disc_value = disc_range[j])
    filename.prev.map <- paste(name.generator(R0, disc_range[j], "Prev"), "pdf", sep = ".")
    ggsave(filename = filename.prev.map, plot = prev.map, width=14, height=9)
    
    
    cum.map <- plotheatmaps(thres.matrix.cum, type = "Cumulative", names = as.factor(dect.cases), R0 = R0, disc_value = disc_range[j])
    filename.cum.map <- paste(name.generator(R0, disc_range[j], "Cumulative"), "pdf", sep = ".")
    ggsave(filename = filename.cum.map, plot = prev.map, width=14, height=9)
  } 
}


 
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
  
  paste(sdata, collapse = '')
  p <- ggplot(df.m, aes(variable, names)) 
  p <- p +  geom_tile(aes(fill = value), colour = "white") + theme_bw()+ scale_fill_gradient(low = "lightyellow",high = "red", name = "Frequency") +labs(x = type, y = "Detected Cases")
  p <- p + theme(axis.title.x = element_text(size=15), axis.text.x  = element_text(size=9), axis.title.y = element_text(size=15), axis.text.y  = element_text(size=9) ) + ggtitle(title) + theme(plot.title = element_text(lineheight=.8, face="bold"))
  return(p)
  p
}

