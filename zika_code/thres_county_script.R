rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')


# Texas Script
# 1) Takes in a list of relative R0 from each metropolitan area and returns values for R0 to use in the simulation
# 2) Uses output of step 1 to run the analysis, stores the number of cases until detection 
# 3) Maps the number of cases onto the 
require(maptools)
require(ggplot2)
require(rgeos)
require(rgdal)
require(raster)
require(plyr)
require(ggthemes)

# 1 County Data
# However we want to scale R0 right now using a log 
# Done in R0 script-will stay the same 
# If wanted to do it by metro area 
R0_range = ddply(test.R0s,.variables = 'MetroArea', summarise, mean(scaledrnott))
colnames(R0_range) = c("MetroArea", "meanR0")


#2 Take in R0 values and calculate the trigger number based on an epidemic threshold and confidence level
# Will separate this into running the trials for each R0 and then for each list calculate the value....
#

# Will want to read in here the different files 
#disc_p = .0246
confidence = .9
threshold.prevalence = 25
threshold.cumulative = 100


#Chosen for analysis 
r_nots <- c(0.9, 1.2, 1.4, 1.6, 1.8, 2, 2.2)
intro_rate <- c(0.3)
disc_prob <- c(0.011, 0.21, 0.0052, 0.068)
dir_path <- "~/Documents/zika_alarm/data/first_runs/"
dirPaths = get_vec_of_files(dir_path, r_nots,  disc_prob, intro_rate)
saveLoc  <- "~/Documents/zika_alarm/data/"

debug(threshold_R0)

calculate_threshold <- threshold_R0(dir_path = dir_path, saveLoc = saveLoc, saveResults = FALSE, r_nots = r_nots,
             intro_rate = intro_rate, disc_prob = disc_prob, confidence = confidence,
             threshold.prevalence = threshold.prevalence, threshold.cumulative = threshold.cumulative)


threshold_R0 <- function(dir_path, saveLoc, saveResults=TRUE, r_nots, intro_rate, disc_prob, 
                         confidence, threshold.prevalence, threshold.cumulative) {
  
  dirPaths = get_vec_of_files(dir_path, r_nots,  disc_prob, intro_rate)
  
  
  calculate_threshold <- adply(.data = dirPaths, .margins = 1, .fun = function (x) {
    load(x)
    
    max.cumulative <- max(all_last_cuminfect_values(trials))
    max.prev <- max(all_last_instantInf_values(trials))
    
    
    lastdetected <- all_last_cumdetect_values(trials)
    max <- set.max.bin(max(lastdetected))
    dect.cases.range <- seq(1:max)
    
    
    #splits up trials into detection 
    trials_by_detection <- alply(.data = dect.cases.range, .margins = 1, function (x) {
      d_thres = as.numeric(x)
      dataframe <- all_detect_rows(trials, threshold = d_thres )      
      return(dataframe)
    })
    
    #setting up bins to calculate frequencies
    bins.prev <- set.prev.bins(max.prev)
    bins.cumulative <- set.cum.bins(max.cumulative)
    
    average.cumulative <- laply(.data = trials_by_detection, function (x) {
      mean.cumulative <- mean(x[,8])
      return(mean.cumulative)
    })
    
    median.cumulative <- laply(.data = trials_by_detection, function (x) {
      median.cum <- median(x[,8])
      return(median.cum)
    })
    
    
    average.prevalence <- laply(.data = trials_by_detection, function (x) {
      mean.prevalence <- mean(x[,7])
      return(mean.prevalence)
    })
    
    median.prevalence <- laply(.data = trials_by_detection, function (x) {
      median.prevalence <- median(x[,7])
    })
    
    
    #Frequency Calculations 
    frequency.cumulative <- ldply(.data = trials_by_detection, function (x) {
      frequency = bin.frequency(x[,8], bins.cumulative)
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    frequency.prevalence <- ldply(.data = trials_by_detection, function (x) {
      frequency = unname(bin.frequency(x[,7], bins.prev))
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    #Clean UP 
    colnames(frequency.prevalence) <- bins.prev; rownames(frequency.prevalence) <- dect.cases.range
    colnames(frequency.cumulative) <- bins.cumulative; rownames(frequency.cumulative) <- dect.cases.range
    
    frequency.prevalence <- frequency.prevalence[,-1]
    frequency.cumulative <- frequency.cumulative[,-1]
    
    
    # ANALYSIS FOR TRIGGER THRESHOLD
    integer.prev <- unname(find_thres_cases(bins.prev,  threshold.cases = threshold.prevalence, df=frequency.prevalence, confidence.value = confidence))
    integer.cumulative <- unname(find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative,
                                                  confidence.value = confidence))
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, confidence, threshold.prevalence, threshold.cumulative)  
    cbind(as.data.frame(matrix(parms,ncol=6)), integer.prev, integer.cumulative)
  })
  
  # Save the results or simply return them
  if(saveResults){
    save( list = c('calculate_threshold'), file = file.path(saveLoc, paste0("calculate_threshold.Rdata")))  
  } else {
    calculate_threshold
  }
}

saveLoc <- "~/Documents/zika_alarm/data/"
save( list = c('calculate_threshold'), file = paste(saveLoc, "calculate_threshold.6.Rdata"))


saveLoc
#for (i in 1:nrow(R0_range)) {
#    prop_p  = R0_range$meanR0[i] / 7
 #   MetroArea = R0_range$MetroArea[i]
  #  percent.discover <- calculate.discover(disc_p)
    
    #Running trials
 #   trials <-
#      run_branches_inc(num_reps = 1000, branch_params(
#        dis_prob_symp = disc_p, prop_p = prop_p, e_thresh = 500
#      ))
 #   lastdetected <- all_last_cumdetect_values(trials)
    
  #  max <- set.max.bin(max(lastdetected))
   # dect.cases.range <- seq(1:max)
   # d_thres <- max
    
    #setting up bins to calculate frequencies
  #  bins.prev <- set.prev.bins(d_thres, trials)
  #  bins.cumulative <- set.cum.bins(d_thres, trials)
    
    #Setting Up data frames and vectors to store
    #thres.matrix.prev <-data.frame(matrix(nrow = length(dect.cases.range), ncol = length(bins.prev) - 1))
    #colnames(thres.matrix.prev) <-bins.prev[2:length(bins.prev)]; rownames(thres.matrix.prev) <-paste("Detected Cases =", dect.cases.range)
    
    #thres.matrix.cum <-data.frame(matrix(nrow = length(dect.cases.range), ncol = length(bins.cumulative) - 1))
    #colnames(thres.matrix.cum) <- bins.cumulative[2:length(bins.cumulative)]; rownames(thres.matrix.cum) <-paste("Detected Cases =", dect.cases.range)
    
    
    #Writing the values
    # Could split this up into two analysis: average/median and frequencies
    #for (dect.case in dect.cases.range) {
     # d_thres <- dect.case
      #dataframe <-
       # all_detect_rows(trials) # takes already whatever the current thresholds
      
      #frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
      #thres.matrix.prev[dect.case,] <- frequencies.prev
      
      #frequencies.cum <-
       # bin.frequency(dataframe[,8], bins.cumulative)
      #thres.matrix.cum[dect.case,] <- frequencies.cum
  #  }
    
    # ANALYSIS FOR TRIGGER THRESHOLD
    #integer.prev <- find_thres_cases(bins.prev,  threshold.cases= threshold.prevalence, df=thres.matrix.prev, confidence.value = confidence)
    #integer.cumulative <- find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = thres.matrix.cum, confidence.value = confidence)
    
   
    # When only want to return the first
    #prevalence.long <- rbind(prevalence.long, cbind(id = MetroArea, disc_p=disc_p, R0=prop_p*7, confidence = confidence, cases = unname(integer.prev)))
    #cumulative.long <- rbind(cumulative.long, cbind(id = MetroArea, disc_p=disc_p, R0=prop_p*7, confidence = confidence, cases = unname(integer.cumulative)))
#} 









# 3  Mapping Output onto Texas
setwd('..'); setwd('TexasCountyShapeFiles')
texas.county <- readShapeSpatial('texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
setwd('../zika_code/')

# Need to read in county ids 
county.ids$Prev.Cases <- NA
county.ids$Cum.Cases <- NA

#Merging the R0 with the final shapefile 
for (i in 1:nrow(prevalence.long)) {
 indices = which(prevalence.long$id[i] == county.ids$Metro) 
  county.ids$Prev.Cases[indices] = prevalence.long$cases[i]
  county.ids$Cum.Cases[indices] = cumulative.long$cases[i]
}
     
texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county.ids, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]

# Decide which type you want to lot

plot_county_threshold(final.plot, type = "Cumulative", confidence = .7, case.threshold = 100)

plot_county_threshold <- function(shp, type, confidence, case.threshold) {
  ## Requires the results to be already fotified with shape file
  if (type == "Prevalence") {
    grey_county <- shp[which(shp$Prev_Cases==0),]
    actual_county <- shp[which(shp$Prev_Cases!=0 | is.na(shp$Prev_Cases)),]
    case.type = "Prev_Cases"
    case.min <- min(shp$Prev_Cases[shp$Prev_Cases != 0], na.rm = TRUE)
    case.max <- max(shp$Prev_Cases, na.rm = TRUE)
    unique.cases <- unique(actual_county_prev$Prev_Cases, na.rm = TRUE)
    
    title = c("Max Cumulative Detected Cases to be ", confidence*100, "% certain \n the outbreak currently is fewer than", case.threshold, "infections")
    title = paste(title, sep="", collapse=" ")
  } else {
    grey_county<- shp[which(shp$Cum_Cases==0),]
    actual_county <- shp[which(shp$Cum_Cases!=0 | is.na(shp$Cum_Cases)),]
    case.type = "Cum_Cases"
    case.min <- min(shp$Cum_Cases[shp$Cum_Cases != 0], na.rm = TRUE)
    case.max <- max(shp$Cum_Cases, na.rm = TRUE)
    unique.cases <- unique(actual_county_cum$Cum_Cases, na.rm = TRUE)
    
    title = c("Max Cumulative Detected Cases to be ", confidence*100, "% certain \n the outbreak currently is fewer than", case.threshold, "infections")
    title = paste(title, sep="", collapse=" ")
  }

    plot <- ggplot()+geom_polygon(data = actual_county_prev, aes_string(x="long", y = "lat", group = "group", fill = case.type), color = "black", size = .25)+coord_map() +
      scale_fill_gradient(name = "Cases", limits = c((case.min-1) ,(case.max+1)), low = "red", high = "yellow",  
                          na.value = "white", guide = "colorbar", breaks = pretty_breaks(n = length(unique.cases))) + 
      labs(title = title) + theme(plot.title = element_text(size = 22)) +
      geom_polygon(data=grey_county_prev, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)
  return(plot)
}


### Save this to make sure everything is working 
grey_county_cum <- final.plot[which(final.plot$Cum_Cases==0),]
actual_county_cum <- final.plot[which(final.plot$Cum_Cases!=0 | is.na(final.plot$Cum_Cases)),]

cumulative.plot <- ggplot()+geom_polygon(data = actual_county_cum, aes(x=long, y = lat, group = group, fill = Cum_Cases), color = "black", size = .25)+coord_map() +
   scale_fill_gradient(name = "Cases", limits = c(min(cumulative.long$cases[cumulative.long$cases != 0])-1 , max(cumulative.long$cases)+1), low = "red", high = "yellow",  
                        na.value = "white", guide = "colorbar", breaks = pretty_breaks(length(unique(cumulative.long$cases)))) + 
  labs(title = "Max Cumulative Detected Cases to be 70% certain \n the outbreak is smaller than 100 cumulative infections") + 
  theme(plot.title = element_text(size = 22)) +
  geom_polygon(data=grey_county_cum, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)



