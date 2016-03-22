# Code I've written just testing things out-not worth saving in the main files 


thousand.total_Infected_20 <- mean(at.detect$Total_Infected)
thousand.cum_Infected_20 <- mean(at.detect$Cumulative_Infections)


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
