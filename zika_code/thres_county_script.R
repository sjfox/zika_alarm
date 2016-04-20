rm(list=ls())
if(grepl('spencerfox', Sys.info()['login'])) setwd('~/projects/zika_alarm/zika_code/')
if(grepl('vagrant', Sys.info()['user'])) setwd('/vagrant/zika_alarm/zika_code/')
if(grepl('sjf826', Sys.info()['login'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('tacc', Sys.info()['nodename'])) setwd('/home1/02958/sjf826/zika_alarm/zika_code/')
if(grepl('meyerslab', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')
if(grepl('laurencastro', Sys.info()['login'])) setwd('~/Documents/zika_alarm/zika_code/')

sapply(c('branch.functions.R','plot.functions.R', 'incubation_branch.R', 'analyze_saved_sims.R'), source)



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
#require(cowplot)

# 1 County Data
# However we want to scale R0 right now using a log 
# Done in R0 script-will stay the same 
# If wanted to do it by metro area 
R0_range = ddply(test.R0s,.variables = 'MetroArea', summarise, mean(scaledrnott))
colnames(R0_range) = c("MetroArea", "meanR0")

R0_metro = ddply(.data = county_ids, .variables = 'Metro', summarise, mean(R0))
R0_metro = R0_metro[-16,]


#2 Take in R0 values and calculate the trigger number based on an epidemic threshold and confidence level


######## Case Scenarios

#Best Case
# High Detect - .068 = 50%
# Low Intro - .01 = .9 cases over 90 days 

# Average Case
# Detection around Symptomatic - .011 = 10%
# Medium Intro (Harris County) - .1 (11 cases over 90 days)

# Worse Cast
# Low Detection - 0.0052 = 5%
# High Intro (State-Wide ) .3 = 27 cases over 90 days

confidence = .8
threshold.prevalence = 20
threshold.cumulative = 100


#Chosen for analysis 
# Trigger Analysis 
r_nots <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.3, 1.5, 1.6, 1.9) 

r_nots <- c(1.5, 1.1, 0.9)
disc_prob <- c( 0.068, 0.011)
intro_rate <- c(0.01, 0.1, 0.3)


#Average Analysis
#r_nots <- c(0.9, 1.5)
#intro_rate <- c(.01, .1, 0.3)
#disc_prob <- c( 0.011, 0.068)


dir_path <- "~/Documents/zika_alarm/data/introductions/"
dirPaths = get_vec_of_files(dir_path, r_nots,  disc_prob, intro_rate)
saveLoc  <- "~/Documents/zika_alarm/data/"


#calculate_threshold <- threshold_R0(dir_path = dir_path, saveLoc = saveLoc, saveResults = FALSE, r_nots = r_nots,
#                                    type = "average", intro_rate = intro_rate, disc_prob = disc_prob, confidence = confidence,
#                                    threshold.prevalence = threshold.prevalence, threshold.cumulative = threshold.cumulative)


threshold_R0 <- function(dir_path, saveLoc, saveResults=TRUE, r_nots, intro_rate, disc_prob, type, 
                         confidence, threshold.prevalence, threshold.cumulative) {
  
 
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_prob, intro_rate)


  calculate_average <- adply(.data = dirPaths, .margins = 1, .expand = TRUE, .fun = function (x) {
    load(x)
    
    max.cumulative <- max(all_last_cuminfect_values(trials))
    max.prev <- max(all_max_prevalence(trials))
    
    
    lastdetected <- all_last_cumdetect_values(trials)
    max <- set.max.bin(max(lastdetected))
    dect.cases.range <- seq(1:max)
    
    #splits up trials into detection 
    trials_by_detection <- alply(.data = dect.cases.range, .margins = 1, function (x) {
      d_thres = as.numeric(x)
      dataframe <- all_detect_rows(trials, threshold = d_thres ) 
      dataframe = na.omit(dataframe)
      return(dataframe)
    })    
    
    average.cumulative <- ldply(.data = trials_by_detection, function (x) {
      mean.cumulative <- mean(x[,"Cumulative_Infections"])
      sd.cumulative <- sd(x[,"Cumulative_Infections"])
    # return(mean.cumulative)
     return(c(mean.cumulative, sd.cumulative))
    })
    # print(average.cumulative)
    
    #print(average.cumulative)
     #median.cumulative <- ldply(.data = trials_by_detection, function (x) {
     # median.cumulative <- median(x[,"Cumulative_Infections])
    # return(median.cumulative)
    #})
    
    average.prevalence <- ldply(.data = trials_by_detection, function (x) {
     mean.prevalence <- mean(x[,"Total_Infections"])
    sd.prevalence<- sd(x[,"Total_Infections"])
    #return(mean.prevalence)
    return(c(mean.prevalence, sd.prevalence))   
    })
    #print(average.prevalence)
    
    #}
    
    #median.prevalence <- ldply(.data = trials_by_detection, function (x) {
    # median.prevalence <- median(x[,"Total_Infections])
    #return(median.prevalence)
    #})   
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate) 
    #print(parms)
    result <- cbind(as.data.frame(matrix(parms,ncol=3)), dect.cases.range, 
                    average.cumulative[,2], average.cumulative[, 3], average.prevalence[,2], average.prevalence[,3])
    # print(result)
    return(result)    
    })
    
    #, median.cumulative[,2], median.prevalence[,2])
    
    
    #setting up bins to calculate frequencies
    bins.prev <- set.prev.bins(max.prev)
    bins.cumulative <- set.cum.bins(max.cumulative)    
    
    #Frequency Calculations 
    frequency.cumulative <- ldply(.data = trials_by_detection, function (x) {
      frequency = bin.frequency(x[,"Cumulative_Infections"], bins.cumulative)
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    frequency.prevalence <- ldply(.data = trials_by_detection, function (x) {
      frequency = unname(bin.frequency(x[,"Total_Infections"], bins.prev))
      frequency = as.data.frame(matrix(frequency, nrow=1))
      return(frequency)
    })
    
    # Clean UP 
    colnames(frequency.prevalence) <- bins.prev; rownames(frequency.prevalence) <- dect.cases.range
    colnames(frequency.cumulative) <- bins.cumulative; rownames(frequency.cumulative) <- dect.cases.range
    
    frequency.prevalence <- frequency.prevalence[,-1]
    frequency.cumulative <- frequency.cumulative[,-1]
    
    # ANALYSIS FOR TRIGGER THRESHOLD
    #integer.prev <- unname(find_thres_cases(bins.prev,  threshold.cases = threshold.prevalence, df=frequency.prevalence, confidence.value = confidence))
    #integer.cumulative <- unname(find_thres_cases(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative,
    #confidence.value = confidence))
    threshold.frequency.prev <- frequency_threshold(bins = bins.prev,threshold.cases = threshold.prevalence, df = frequency.prevalence)
    threshold.frequency.cum <- frequency_threshold(bins = bins.cumulative, threshold.cases = threshold.cumulative, df = frequency.cumulative)
    #parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, confidence, threshold.prevalence, threshold.cumulative) 
    parms <- c(params$r_not, params$dis_prob_symp, params$intro_rate, threshold.prevalence, threshold.cumulative)
    
    cbind(as.data.frame(matrix(parms,ncol=5)),dect.cases.range, threshold.frequency.prev, threshold.frequency.cum)
  })


    # Save the results or simply return them
    if(saveResults){
      save( list = c('calculate_threshold'), file = file.path(saveLoc, paste0("calculate_threshold.Rdata")))  
    } else {
  
    calculate_threshold
  }
}



saveLoc <- "~/Documents/zika_alarm/data/"
save(list = c('calculate_threshold'), file = paste(saveLoc, "calculate_threshold_1000.Rdata"))


sd <- calculate_threshold_avg_sd[,6:7]
full_threshold <- cbind(calculate_threshold, sd)
head(full_threshold)


colnames(calculate_threshold) <- c("Run", "R0", "Dect", "Intro", "DectCases",
                                   "Avg.Cum",  "avg.cum.sd", "Avg.Prev","avg.prev.sd")


head(full_threshold)

# Going from wide format to long format : combining the average value with error bars 
detection.m <- melt(data = full_threshold, id.vars = c("R0", "Dect", "Intro", "DectCases"), measure.vars = c("Avg.Cum", "Avg.Prev" ))
error.m <- melt(data = full_threshold, id.vars = c("R0", "Dect", "Intro", "DectCases"), measure.vars = c("avg.cum.sd", "avg.prev.sd"))

colnames(error.m) <- c("R0", "Dect", "Intro", "DectCases", "Error", "SE")
detection.m <- cbind(detection.m, error.m[,5:6])
head(detection.m)

# If want to split the Results to certain detection values and if plotting Prevalence versus Cumulative 
indices = which(detection.m$variable == "Avg.Prev"  & (detection.m$R0 == 0.9 | detection.m$R0 == 1.5) & (detection.m$Intro == 0.3))
indices = which(detection.m$variable == "Avg.Prev") 
detection.avg = detection.m[indices, ]

max(detection.avg$value)
max(detection.avg$DectCases)
breaks = c(1,5,10,20,30,40,50, 100, 500) #, 200, 300,400,500) #, 200, 300, 400, 500) # Set according to max of value + se
breaks_x = seq(from = 0, to = 30, by = 10)



plot_trial <- ggplot(detection.avg, aes(DectCases, value, fill = as.factor(R0), color = as.factor(R0), group=interaction(as.factor(Dect), R0)))  + 
  geom_line(size=1.5, aes(linetype = as.factor(Dect))) + #facet_grid(~Intro) +
  geom_ribbon(aes(ymin = value-SE, ymax=value+SE), alpha=.2, color = NA) + 
  scale_y_log10(breaks = breaks)  +
  scale_color_brewer(palette = "Set1", guide = FALSE, direction = -1) +
  scale_fill_brewer(palette="Set1", guide_legend(title = "R0"), direction = -1) +
  scale_x_continuous(name = "Cumulative Number of Detected Cases", breaks = breaks_x, limits = c(0,30)) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.title.y = element_text(size=28), axis.text.y = element_text(size = 22)) + 
  theme(axis.title.x = element_text(size=28), axis.text.x= element_text(size=22)) +
  labs(y = "Expected Total Current Cases", linetype = "Detection \n Rate") +
  theme(legend.text=element_text(size=22, margin = margin(), debug = FALSE), legend.title = element_text(size = 28)) #+ 
#
plot_trial


 
  


# 3  Mapping Output onto Texas-Need Different Outcomes 
setwd('..'); setwd('TexasCountyShapeFiles')
texas.county <- readShapeSpatial('texas.county.shp', proj4string = CRS("+proj=longlat +datum=WGS84"))
setwd('../zika_code/')

# Need to read in county ids 
county_ids$Prev.Cases.worst <- NA
county_ids$Cum.Cases.worst <- NA

#Merging the R0 with the final shapefile 
# Put desired Target Threshold 
unique(calculate_threshold$V1)


#Collects Only the Triggers You're Interested In 
desired_dect = 0.0110; desired_intro = .01

triggers <- ddply(.data = calculate_threshold, .(V1), function (x) {
  row = which(x[,3] == desired_dect & x[,4] == desired_intro)
  return(x[row,])
})

# Merges The Data With the County Data for Plotting By R0 


for (i in 1:nrow(trigger_worst)) {
  indices = which(trigger_worst[i,2] == county_ids$metro_round) 
  county_ids$Prev.Cases.worst[indices] = trigger_worst[i, 8 ]
  county_ids$Cum.Cases.worst[indices] = trigger_worst[i,9]
}
     
county_ids$Prev.Cases.worst <- as.numeric(county_ids$Prev.Cases.worst)
county_ids$Cum.Cases.worst <- as.numeric(county_ids$Cum.Cases.worst)

texas.county.f <- fortify(texas.county, region = "ID")
merge.texas.county <- merge(texas.county.f, county_ids, by = "id", all.x = TRUE)
final.plot <- merge.texas.county[order(merge.texas.county$id),]
save(list = c('final.plot'), file = paste(saveLoc, "Avg_Worst_Plot.Rdata"))
# Decide which type you want to plot

legend_breaks <- round(seq(0, 40, 8))

p1.leg <- ggplot(data = final.plot, aes_string("long", "lat", "group", fill = "Prev.Cases.avg")) + geom_polygon() 

plotworst <- ggplot()+geom_polygon(data = final.plot, aes_string(x="long", y = "lat", group = "group", fill = "Prev.Cases.worst" ),
                              color = "black", size = .25) + coord_map() +
  scale_fill_gradient(name = "Detected Cases", low = "red", high = "yellow", 
                      na.value = "grey", breaks = legend_breaks) + #pretty_breaks(n = length(legend_breaks))) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(legend.position = "bottom") +
  theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
  theme(legend.key.size =  unit(0.5, "in")) 


plotavg <- ggplot()+geom_polygon(data = final.plot, aes_string(x="long", y = "lat", group = "group", fill = "Prev.Cases.avg" ),
                                   color = "black", size = .25) + coord_map() +
  scale_fill_gradient(name = "Detected Cases", low = "red", high = "yellow", 
                      na.value = "grey", breaks = legend_breaks) + #pretty_breaks(n = length(legend_breaks))) +
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
  theme(axis.ticks = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), line = element_blank()) +
  labs(x=NULL, y = NULL) +
  theme(legend.position = "none")
  #theme(legend.position = "right") +
  #theme(legend.text=element_text(size=14, margin = margin(), debug = FALSE), legend.title = element_text(size = 20)) +
  #theme(legend.key.size =  unit(0.5, "in")) 


#Funtion to extract one legend 
library(gridExtra)
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plotworst)

tmp <- arrangeGrob(plotavg + theme(legend.position = "none"), plotworst + 
                     theme(legend.position = "none"), layout_matrix = matrix(c(1, 2), nrow = 1))

p3 <- grid.arrange(tmp, mylegend, nrow= 1 ) #, heights = unit.c(unit(1, "npc") - lheight, lheight))

p3 <- grid.arrange(arrangeGrob(plotavg + theme(legend.position="none"),
                               plotworst + theme(legend.position="none"),
                               nrow=1),
                   mylegend, nrow=2,heights=c(10, 1))
plot_grid(plotavg, plotworst, labels = c("A", "B"), ncol = 1)



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
lamb_county <- final.plot[which(final.plot$Geography == "Lamb County, Texas"),]
grey_county_cum <- final.plot[which(final.plot$Cum_Cases==0),]


actual_county_cum <- final.plot[which(final.plot$Cum_Cases!=0 | is.na(final.plot$Cum_Cases)),]

cumulative.plot <- ggplot()+geom_polygon(data = actual_county_cum, aes(x=long, y = lat, group = group, fill = Cum_Cases), color = "black", size = .25)+coord_map() +
   scale_fill_gradient(name = "Cases", limits = c(min(cumulative.long$cases[cumulative.long$cases != 0])-1 , max(cumulative.long$cases)+1), low = "red", high = "yellow",  
                        na.value = "white", guide = "colorbar", breaks = pretty_breaks(length(unique(cumulative.long$cases)))) + 
  labs(title = "Max Cumulative Detected Cases to be 70% certain \n the outbreak is smaller than 100 cumulative infections") + 
  theme(plot.title = element_text(size = 22)) +
  geom_polygon(data=grey_county_cum, aes(x=long, y = lat, group = group), fill="grey", color = "black", size = .25, inherit.aes = FALSE)




##### Line Graphs For When  Wnat to Compare Trends Over Multiple Confidence Levels and detection Frequencies 


calculate_threshold.m <- melt(calculate_threshold, id.vars = c("V1", "V2"), measure.vars = c("integer.prev", "integer.cumulative"))


#middle.stats.average <- data.frame(cbind(dect.cases.range, average.vec))
#plot.average <- ggplot(middle.stats.average, aes(dect.cases.range, average.vec))
#plot.average + geom_point()
#plot.average = plot.average + geom_line(size = 2) + labs(x = "Cumulative Detected Cases", y = "Average Cumulative Total Cases")

plot1 <- ggplot(calculate_threshold.m, aes(V1, value, color = variable)) + 
  geom_line(size=1.5, aes(linetype = as.factor(V2))) + facet_grid(~V2) +
  scale_color_brewer(palette="Set1") +
  labs(x = "R0", y = "Max Number of Cases")
  
theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank()) +
##### Examples
plot1 <- ggplot(escape_data, aes(d_thresh, probEsc, color = as.factor(r_not))) + 
  geom_line(size=1.5, aes(linetype=as.factor(disc_p))) + facet_grid(~disc_p)+
  scale_y_continuous(expand=c(0.01,0.01)) +
  scale_color_brewer(palette="Set1")+
  theme_cowplot() %+replace% theme(strip.background=element_blank(),strip.text.x = element_blank())+
  labs(x = "Cumulative Number of Detected Cases", y = "Probability of an Epidemic", color = expression("R"[0]))+
  guides(linetype= FALSE )
