
#####################################################################

# Analysis functions for extracting various elements from the trials 
last_cuminfect_value <- function(x) {
  x[nrow(x), "Cumulative_Infections"]
}
all_last_cuminfect_values <- function(x) {
  return(unlist(laply(x, last_cuminfect_value)))
}

last_cuminfect_local_value <- function(x) {
  x[nrow(x), "Cumulative_Infections"] - x[nrow(x), "Cumulative_Intro_Infections"]
}
all_last_cuminfect_local_values <- function(x) {
  return(unlist(laply(x, last_cuminfect_local_value)))
}

last_cuminfect_intro_value <- function(x) {
  x[nrow(x), "Cumulative_Intro_Infections"]
}
all_last_cuminfect_intro_values <- function(x) {
  return(unlist(laply(x, last_cuminfect_intro_value)))
}

## Get last prevalence
last_infectious_value <- function(x) {
  x[nrow(x), "Total_Infections"]
}
all_last_infectious_values <- function(x) {
  return(unlist(laply(x, last_infectious_value)))
}

## Get last cumulative detections
last_cumdetect_value <- function(x) {
  x[nrow(x), "Cum_Detections"]
}
all_last_cumdetect_values <- function(x) {
  laply(x, last_cumdetect_value)
}


## Get last cumulative introduced detections
last_cumdetect_intro_value <- function(x) {
  x[nrow(x), "Cum_Intro_Detections"]
}
all_last_cumdetect_intro_values <- function(x) {
  laply(x, last_cumdetect_intro_value)
}

## Get last cumulative local detections
last_cumdetect_local_value <- function(x) {
  x[nrow(x), "Cum_Detections"] - x[nrow(x), "Cum_Intro_Detections"]  
}
all_last_cumdetect_local_values <- function(x) {
  laply(x, last_cumdetect_local_value)
}


## Find the max prevalence
max_prevalence <- function(x){
  return(max(x[,"Total_Infections"]))
}
all_max_prevalence <- function(x) {
  return(unlist(laply(x, max_prevalence)))
}

## Find the max non introduced prevalence
max_nonintro_prevalence <- function(x){
  return(max(x[,"Total_Infections"]-x[,"Total_Intro_Infections"]))
}
all_max_nonintro_prevalence <- function(x) {
  return(unlist(laply(x, max_nonintro_prevalence)))
}

cum_infect_local <- function(x){
  ## Returns column of cumulative local infections
  x[, "Cumulative_Infections"] - x[, "Cumulative_Intro_Infections"]
}

cum_detect_local <- function(x){
  ## Returns column of cumulative local detections
  x[, "Cum_Detections"] - x[, "Cum_Intro_Detections"]
}

cum_detect_total <- function(x){
  ## Returns column of cumulative total detections
  x[, "Cum_Detections"]
}

prevalence_local <- function(x){
  ## Returns column of local prevalence
  x[, "Total_Infections"] - x[, "Total_Intro_Infections"]
}

prevalence_total <- function(x){
  ## Returns column of local prevalence
  x[, "Total_Infections"]
}

## Get Time of max prevalence
time_max_prevalence <- function(x){
  return(x[which.max(x[,"Total_Infections"]), "time"])
}
all_time_max_prevalence <- function(x) {
  return(unlist(laply(x, time_max_prevalence)))
}

## Get duration of simulation
sim_duration <- function(x){
  return(nrow(x))
}
all_sim_duration <- function(x) {
  return(unlist(laply(x, sim_duration)))
}


### Calculate Secondary Transmitted Detections for whole time series
calculate_second_detections <- function(x) {
  second_detections <- x[,"Cum_Detections"] - x[,"Cum_Intro_Detections"]
  x <- cbind(x, second_detections)
  return(x)
}

all_second_detections <- function(x) {
  return(llply(x, calculate_second_detections))
}

### Calculate Secondary Transmitted Cases for whole time series 
calculate_second_prevalence <- function(x) {
  second_infections <- x[,"Total_Infections"] - x[,"Total_Intro_Infections"]
  x <- cbind(x, second_infections)
  return(x)
}

all_second_current_infections <- function(x) {
  return(llply(x, calculate_second_prevalence))
}




## Calculate Secondary Cumulative Cases for whole time series 




## Set of functions to calculate given I have X cases, what is the distribution of cases I see 
library(ggplot2)
library(reshape2)
library(scales)

## Set of functions to calculate given I have X cases, what is the distribution of cases I see 
detection_rows <- function(trial, threshold) {
  detections <- trial[,"Cum_Detections"]
  rows <- which(detections == threshold)
  reduced <- trial[rows[1],]
  return(reduced)
}

# Function to return all rows in trials that match detection threshold
all_detect_rows <- function(trials, threshold) {
  all_rows_detection = ldply(.data = trials, .fun = detection_rows, threshold = threshold) 
  return(all_rows_detection)
}


metro_areas <- list(dallas =c("Collin", "Dallas", "Denton", "Ellis", "Hood", "Hunt", "Johnson", "Kaufman", "Parker", "Rockwall", "Somervell", "Tarrant", "Wise"),
                    houston = c("Harris", "Fort Bend", "Montgomery", "Brazoria", "Galveston", "Liberty", "Waller", "Chambers", "Austin"),
                    san_antonio = c("Atascosa", "Bandera", "Bexar", "Comal", "Guadalupe", "Kendall", "Medina", "Wilson"),
                    austin = c("Bastrop", "Caldwell", "Hays", "Travis", "Williamson"),
                    el_paso = c("El Paso", "Hudspeth"),
                    mcallen = c("Hidalgo"),
                    corpus_christi = c("Aransas", "Nueces", "San Patricio"),
                    brownsville = c("Cameron"),
                    killeen = c("Bell", "Coryell", "Lampasas"),
                    beaumont = c("Hardin", "Jefferson", "Newton", "Orange"),
                    lubbock = c("Crosby", "Lubbock", "Lynn"),
                    laredo = c("Webb"),
                    amarillo = c("Armstrong", "Carson", "Potter", "Randall"),
                    waco = c("McLennan", "Falls"),
                    college_station = c("Brazos", "Burleson", "Robertson"))

# Additional Post Processing Functions 
prob_ext <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=100) {
  escapes = 0
  escapes_prob <- c(escapes)
  i = 1 

  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) #Run the simulation 
    Final.I = record[nrow(record),7]
    Final.D = record[nrow(record),6] ## should 6 be changed to 5? # no it's cumulative detection, instantaneous infectious

    if (Final.D < d_thres)  {
      next
    }
    if (Final.I > e_thresh) {
      escapes = escapes + 1   
    }
    escapes_prob <- c(escapes_prob, escapes/i)
    i = i+1    
  }
  return(escapes_prob)
}



# If you have not crossed the detection threshold, what is the probability that you have crossed your epidemic threshold? Seems to never happen
prob_underD.overE <- function(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate, nsamples=10) {
  escapes = 0
  escapes_prob <- c(escapes)
  i = 1 
  
  while (i < nsamples) {
    record = run_branch(prop_p, recov_p, incub_p, prob_symp, d_thres, e_thresh, dis_prob_asymp, dis_prob_symp, intro_rate) #Run the simulation 
    Final.I = record[nrow(record),7]
    Final.D = record[nrow(record),6]
    
    if (Final.D < d_thres)  { # Haven't Crossed Detection 
      if (Final.I > e_thresh) { # Is my Instantaneous Number of Infections more than my epidemic threshold
        escapes = escapes + 1   
        escapes_prob <- c(escapes_prob, escapes/i)
        
      }
    }
    i = i+1   
  }
  return(escapes_prob)
}




##### Functions for generating the heat maps 
name.generator <- function(R0, disc_value, type) {
  filename.base <- paste(R0, disc_value, sep = "_")
  filename.type <- paste(type, filename.base, sep = "_")
  return(filename.type)
} 


set.cum.bins <- function(max.cum) {
  #at.high.dect <- all_detect_rows(trials)
  #max.cum <- max(at.high.dect[,8])
  max.bin <- max.cum + 10 # normally 25
  #beginning.bins <- seq(from = 0, to = 18, by = 2) When I had higher resolution 
  #bins.1 <- seq(from = 20, to = max.bin, by = 20)
  #bins <- c(beginning.bins, bins.1)
  bins <- seq(from = 0, to = max.bin, by = 10)
  return(bins)
}


set.prev.bins <- function(max.prev) {
  #at.high.dect <- all_detect_rows(trials)
  #max.prev <- max(at.high.dect[,7])
  max.bin <- max.prev + 1
  bins <- seq(from = 0, to = max.bin, by = 1)
  #if (max.bin > 10) {
  #bins.1 <- seq(from = 10, to = max.bin, by = 5)
  #bins <- c(bins, bins.1)
  #}
  return(bins)
}

#Function to take the row and calculate proabilities for bins
bin.frequency <- function(df, bins) {
  df.cut <- cut(df, breaks = bins, right = TRUE)
  df.freq = table(df.cut)
  df.freq <- cbind(df.freq)
  df.freq <- df.freq/sum(df.freq)
  return(df.freq)
}

# Keep Things Ordered Here ording by a certain variable

plotheatmaps <- function(df, type, names, percent.discover, R0, max.infect) {
  df <- cbind(names, df)
  df$names <- factor(df$names, levels = df$names[order(df$names)])
  df.m <- melt(df)
  
  x.breaks.seq = seq(0, length(names), by = 2)
  
  if (type == "Cumulative") {
    y.breaks.seq <- seq(0, max.infect+20, by = 20)
    title = "Cumulative Total Cases"
  } else {
    y.breaks.seq <- seq(0, max.infect+5, by = 5)
    title = "Current Infected Cases"
  } 
  
  p <- ggplot(df.m, aes(as.factor(names), variable))
  p <- p + geom_tile(aes(fill = value), colour = "white") + 
    theme_bw()+ scale_fill_gradient(low = "lightyellow",high = "red", name = "Probability") + 
    labs(x = "Cumulative Detected Cases", y = title) + 
    theme(axis.text.x = element_text(vjust = 1, hjust=1))    
  p <- p + theme(axis.title.x = element_text(size=20), axis.text.x= element_text(size=14))
  p <- p + theme(axis.title.y = element_text(size=20), axis.text.y = element_text(size = 14)) + 
    theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
    scale_y_discrete(breaks = y.breaks.seq)+
    scale_x_discrete(breaks = x.breaks.seq)
  return(p)
}

## Function to find number of detected cases for X% sure is X or less
find_thres_cases <- function(bins, threshold.cases, df, confidence.value) {
  df.t <- t(df)
  if(max(bins) < threshold.cases) {
    detection.thres <- NA
  }
  else {
    marker <- which(bins == threshold.cases)-1
    df.t <- df.t[1:marker,] # Cut 
    col_sum <- colSums(x = df.t)  
    detection.thres.candidates <- which(col_sum > confidence.value)
    
    if (length(detection.thres.candidates) == 0) {
      prob = max(col_sum) 
      #detection.thres <- which(col_sum == prob)
      #detection.thres <- min(detection.thres) # If there's multiple that have the same 
      detection.thres <- 0
    } else {
      prob <- min(col_sum[detection.thres.candidates])
      if (round(prob, digits = 2) == 1) {
        detection.thres <- NA
      } else {
        detection.thres.multiple <- which(col_sum == prob)
        detection.thres <- min(detection.thres.multiple) #If there's multiple that have the same 
      }
    }
  }
  return(case.trigger = detection.thres)
  #return(list(thres.int = detection.thres, prob = prob))
}



frequency_threshold <- function(bins, threshold.cases, df) {
  df.t <- t(df)
  if(max(bins) < threshold.cases) {
    detection.thres <- NA
  }
  else {
    marker <- which(bins == threshold.cases)-1
    df.t <- df.t[1:marker,] # Cut 
    col_sum <- colSums(x = df.t)  
    return(col_sum) 
  }
}


set.max.bin <- function(max) {
  if (max >= 200) {
    max = round(max * .2)
  } else if (max > 100) {
    max = round(max*.3)
  } else if (max > 80 & max <= 100) {
    max = round(max * 0.40)
  }  else if (max > 60 & max <= 80) {
    max = round(max* 0.5)
  } else if (max > 40 & max <= 60) {
    max = round(max * 0.75)
  } else max = max
  return(max)
}



### Function to take in the relative sustained and scale relative to 1.7
scale_rnott <- function(relative, max) {
  relative.rnott <- relative/max * 1.7
  return(relative.rnott) 
}

# Function to calculate discovery percentage 
calculate.discover <- function(disc_p){
  round(1-((1-disc_p)^9.88), digits = 2)*100
} 



  