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

dect.thresholds <- c(0:20)
bins.prev <- seq(from = 0, to = 80, by = 10)
bins.cumulative <- seq(from = 0, to = 250, by = 25)


prev.names <- paste("< ", bins.prev[2:length(bins.prev)])
cumulative.names <- paste("<", bins.cumulative[2:length(bins.cumulative)])

thres.matrix.prev <- data.frame(matrix(nrow = length(dect.thresholds), ncol = length(bins.prev)-1))
colnames(thres.matrix.prev) <- prev.names ; rownames(thres.matrix.prev) <- paste("Detected Cases = ", dect.thresholds)
thres.matrix.cum <- data.frame(matrix(nrow = length(dect.thresholds), ncol = length(bins.cumulative)-1))
colnames(thres.matrix.cum) <- cumulative.names; rownames(thres.matrix.cum) <- paste("Detected Cases = ", dect.thresholds)


R0 = prop_p/recov_p
for (i in 1:length(dect.thresholds)) {
  d_thres <- dect.thresholds[i]
  print(d_thres)
  dataframe <- all_detect_rows(trials) # takes already whatever the current thresholds 
  frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
  frequencies.cum <- bin.frequency(dataframe[,8], bins.cumulative)
  
  thres.matrix.prev[i,] <- frequencies.prev
  thres.matrix.cum[i,] <- frequencies.cum
}



filename <- paste(R0, dis_prob_symp, sep = "_")

filename.prev <- paste("Prev", filename, sep = ".")
filename.prev <- paste(filename.prev, "csv", sep = ".")
write.csv(x = thres.matrix.prev, file = filename.prev, row.names = TRUE)



filename.cum <- paste("Cumulative", filename, sep = ".")
filename.cum <- paste(filename.cum, "csv", sep = ".")

write.csv(x = thres.matrix.cum, file = filename.cum, row.names = TRUE)



#Function to take the row and calculate proabilities for bins
bin.frequency <- function(df, bins) {
  df.cut <- cut(df, breaks = bins, right = FALSE)
  df.freq = table(df.cut)
  df.freq <- cbind(df.freq)
  df.freq <- df.freq/sum(df.freq)
  return(df.freq)
}


