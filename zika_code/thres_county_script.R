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

# 1 County Data
# However we want to scale R0 right now using a log 
# Found mean of the habitat suitability and then calculated scaled R0 based on that 
habitat_range = ddply(county.ids, .variables = 'Metro', summarise, mean(HabitatSuitability))
habitat_range = habitat_range[-16,]
habitat_max = max(habitat_range$..1)

habitat_range$scaledrnott = habitat_range$..1/habitat_max * 2.5

test.R0s$logrelative <- log(test.R0s$RelativeR0+1)
max.relativeR0 = max(test.R0s$logrelative)
test.R0s$scaledrnott = test.R0s$logrelative/max.relativeR0 * 1.5



R0_range = ddply(test.R0s,.variables = 'MetroArea', summarise, mean(scaledrnott))
colnames(R0_range) = c("MetroArea", "meanR0")


#2 Take in R0 values and calculate the trigger number based on an epidemic threshold and confidence level
# Will separate this into running the trials for each R0 and then for each list calculate the value....
#

disc_p = .0246
confidence = .8
threshold.prevalence = 10
threshold.cumulative = 50
prevalence.long <- data.frame()
cumulative.long <- data.frame()


#trials_R0 <- function(df) {
#  prop_p = df[,2]/7
#  MetroArea = df[,1]/7
#  trials <-run_branches_inc(num_reps = 1000, branch_params(prop_p = prop_p))
#  return(trials)
#}

# Try this for step 1 
#trials.list <- dlply(R0_range, .variables = 'MetroArea', trials_R0)
# This may work but suppppppper slow 

for (i in 1:nrow(R0_range)) {
    prop_p  = R0_range$meanR0[i] / 7
    MetroArea = R0_range$MetroArea[i]
    percent.discover <- calculate.discover(disc_p)
    
    #Running trials
    trials <-
      run_branches_inc(num_reps = 1000, branch_params(
        dis_prob_symp = disc_p, prop_p = prop_p, e_thresh = 500
      ))
    lastdetected <- all_last_cumdetect_values(trials)
    
    max <- set.max.bin(max(lastdetected))
    dect.cases.range <- seq(1:max)
    d_thres <- max
    
    #setting up bins to calculate frequencies
    bins.prev <- set.prev.bins(d_thres, trials)
    bins.cumulative <- set.cum.bins(d_thres, trials)
    
    #Setting Up data frames and vectors to store
    thres.matrix.prev <-data.frame(matrix(nrow = length(dect.cases.range), ncol = length(bins.prev) - 1))
    colnames(thres.matrix.prev) <-bins.prev[2:length(bins.prev)]; rownames(thres.matrix.prev) <-paste("Detected Cases =", dect.cases.range)
    
    thres.matrix.cum <-data.frame(matrix(nrow = length(dect.cases.range), ncol = length(bins.cumulative) - 1))
    colnames(thres.matrix.cum) <- bins.cumulative[2:length(bins.cumulative)]; rownames(thres.matrix.cum) <-paste("Detected Cases =", dect.cases.range)
    
    
    #Writing the values
    # Could split this up into two analysis: average/median and frequencies
    for (dect.case in dect.cases.range) {
      d_thres <- dect.case
      dataframe <-
        all_detect_rows(trials) # takes already whatever the current thresholds
      
      frequencies.prev <- bin.frequency(dataframe[,7], bins.prev)
      thres.matrix.prev[dect.case,] <- frequencies.prev
      
      frequencies.cum <-
        bin.frequency(dataframe[,8], bins.cumulative)
      thres.matrix.cum[dect.case,] <- frequencies.cum
    }
    
    # ANALYSIS FOR TRIGGER THRESHOLD
    integer.prev <- find_thres_cases(bins.prev,  threshold.cases= threshold.prevalence, df=thres.matrix.prev, confidence.value = confidence)
    integer.cumulative <- find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = thres.matrix.cum, confidence.value = confidence)
    
   
    # When only want to return the first
    prevalence.long <- rbind(prevalence.long, cbind(id = MetroArea, disc_p=disc_p, R0=prop_p*7, confidence = confidence, cases = unname(integer.prev)))
    cumulative.long <- rbind(cumulative.long, cbind(id = MetroArea, disc_p=disc_p, R0=prop_p*7, confidence = confidence, cases = unname(integer.cumulative)))
} 




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
plot_county_threshold(final.plot, type = "Cumulative", confidence = .8, case.threshold = 50)
plot_county_threshold(final.plot, type = "Prevalence", confidence = .8, case.threshold = 10)

plot_county_threshold <- function(shp, type, confidence, case.threshold) {
  ## Requires the results to be already fotified with shape file
  if (type == "Prevalence") {
    grey_county <- shp[which(shp$Prev.Cases==0),]
    actual_county <- shp[which(shp$Prev.Cases!=0 | is.na(shp$Prev.Cases)),]
    case.type = "Prev.Cases"
    case.min <- min(shp$Prev.Cases[shp$Prev.Cases != 0], na.rm = TRUE)
    case.max <- max(shp$Prev.Cases, na.rm = TRUE)
    unique.cases <- unique(actual_county$Prev.Cases, na.rm = TRUE)
    
    title = c("Max Cumulative Detected Cases to be ", confidence*100, "% certain \n the outbreak currently is fewer than", case.threshold, "infections")
    title = paste(title, sep="", collapse=" ")
  } else {
    grey_county<- shp[which(shp$Cum.Cases==0),]
    actual_county <- shp[which(shp$Cum.Cases!=0 | is.na(shp$Cum.Cases)),]
    case.type = "Cum.Cases"
    case.min <- min(shp$Cum.Cases[shp$Cum.Cases != 0], na.rm = TRUE)
    case.max <- max(shp$Cum.Cases, na.rm = TRUE)
    unique.cases <- unique(actual_county$Cum.Cases, na.rm = TRUE)
    
    title = c("Max Cumulative Detected Cases to be ", confidence*100, "% certain \n the outbreak currently is fewer than", case.threshold, "total infections")
    title = paste(title, sep="", collapse=" ")
  }

    plot <- ggplot()+geom_polygon(data = actual_county, aes_string(x="long", y = "lat", group = "group", fill = case.type), color = "black", size = .25)+coord_map() +
      scale_fill_gradient(name = "Cases", limits = c((case.min-1) ,(case.max+1)), low = "red", high = "yellow",  
                          na.value = "white", guide = "colorbar", breaks = pretty_breaks(n = length(unique.cases))) + 
      labs(title = title) + theme(plot.title = element_text(size = 22)) +
      geom_polygon(data=grey_county, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)
  return(plot)
}


### Save this to make sure everything is working 
grey_county_prev <- final.plot[which(final.plot$Prev.Cases==0),]
actual_county_prev <- final.plot[which(final.plot$Prev.Cases!=0 | is.na(final.plot$Prev.Cases)),]

prevalence.plot <- ggplot()+geom_polygon(data = actual_county_prev, aes(x=long, y = lat, group = group, fill = Prev.Cases), color = "black", size = .25)+coord_map() +
   scale_fill_gradient(name = "Cases", limits = c((case.min-1) , (case.max+1)), low = "red", high = "yellow",  
                        na.value = "white", guide = "colorbar", breaks = pretty_breaks(length(unique(cumulative.long$cases)))) + 
  labs(title = "Max Cumulative Detected Cases to be 80% certain \n the number of current infections is les than 20") + 
  theme(plot.title = element_text(size = 22)) +
  geom_polygon(data=grey_county_prev, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)



