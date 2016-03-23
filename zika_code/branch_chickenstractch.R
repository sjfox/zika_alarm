# Code I've written just testing things out-not worth saving in the main files 




last_infected <- all_last_instantInf_values(trials)
hist(last_infected)



hist(escape_prob[200:1000], breaks = 10, main = "d_thres = 10")  
plot(last_infected, lastvalues)


#Plotting Preliminary Results
lastvalues <- all_last_cuminfect_values(trials)
hist(lastvalues, main = paste( "1000 Runs, R0 =", prop_p/recov_p, sep = ""), breaks = 10, ylim = c(1,1000))

plot(last_infected, lastvalues)
last_detect_values <- all_last_cumdetect_values(trials)
hist(last_detect_values)
plot(last_detect_values, lastvalues,  xlab = "Cumulative Detected", ylab =  "Cumulative Infected", main = "1000 Runs-Running until Instantaneous Infected > e_thres")


par(mfrow = c(2,4))

for(i in 1:length(d_thres_seq)) {
  d_thres = d_thres_seq[i]
  at.detect <- all_detect_rows(trials)
  hist(at.detect$Cumulative_Infections, main = paste("Detection Threshold = ", d_thres, sep = ""), xlab = "Cumulative Infections")
  mean.currentinfected[i] <- mean(at.detect$Cumulative_Infections)
}



mean.totalinfected[i] <- mean(at.detect$Total_Infected)
plot(d_thres_seq, mean.totalinfected, xlab = "Detection Threshold", ylab = "Mean Current Infecteds", main = "Given X cases at Threhsold, How Many Current Infecteds")
plot(d_thres_seq, mean.currentinfected, xlab = "Detection Threshold", ylab = "Mean Cumulative Infections", main = "Given X cases at Threhsold, How Many Cumulative Infections")
