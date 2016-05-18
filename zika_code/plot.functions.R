#####################################
##
## Plotting functions for zika data
#####################################

plot_final_sizes <- function(x){
  final.sizes <- all_last_cuminfect_values(x)
  qplot(final.sizes, geom="histogram", bins=100) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Final Epidemic Sizes")
}

plot_local_final_sizes <- function(x){
  final.sizes <- all_last_cuminfect_local_values(x)
  qplot(final.sizes, geom="histogram", bins=100) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Final Local Epidemic Sizes")
}

plot_intro_final_sizes <- function(x){
  final.sizes <- all_last_cuminfect_intro_values(x)
  qplot(final.sizes, geom="histogram", bins=100) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Final Introduced Epidemic Sizes")
}

plot_max_nonintro_prevalences <- function(x){
  final.sizes <- all_max_nonintro_prevalence(x)
  qplot(final.sizes, geom="histogram", bins=50) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Max non-introduced prevalence")
}

plot_max_prevalences <- function(x){
  final.sizes <- all_max_prevalence(x)
  qplot(final.sizes, geom="histogram", bins=50) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Max TOTAL prevalence")
}

plot_detections <- function(x){
  final.sizes <- all_last_cumdetect_values(x)
  qplot(final.sizes, geom="histogram", bins=50) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Detections")
}

plot_prob_below <- function(df){
  df$disc_prob <- paste0(calculate.discover(df$disc_prob), "%")
  ggplot(df, aes(detected, prob_below, linetype=as.factor(disc_prob), color = as.factor(r_not))) + 
    geom_line(size=1) + facet_wrap(~intro_rate)+
    scale_y_continuous(expand=c(0.01,0.01)) +
    theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=1, linetype=1))+
    scale_color_brewer(palette="Set1") +
    panel_border(size=.5, colour="black") +
    background_grid(major = "y", minor = "none")+
    labs(x = "Cumulative Number of Reported Cases", 
         y = "Probability Prevalence is Below 20", 
         color = expression("R"[0]),
         linetype= "Reporting\nRate")
  
}


plot_epidemic_prob <- function(df){
  df$disc_prob <- paste0(calculate.discover(df$disc_prob), "%")
  ggplot(df, aes(detected, prob_epidemic, color = as.factor(r_not))) +
    geom_line(size=1, aes(linetype=as.factor(disc_prob))) + facet_wrap(~intro_rate)+
    scale_y_continuous(expand=c(0.01,0.01)) +
    scale_color_brewer(palette="Set1")+
    theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=1, linetype=1))+
    panel_border(size=.5, colour="black") +
    background_grid(major = "y", minor = "none")+
    labs(x = "Cumulative Number of Reported Cases", 
         y = "Probability of an Epidemic", 
         color = expression("R"[0]), 
         linetype= "Reporting\nRate")
}

plot_prevalences <- function(df){
  
  df$disc_prob <- paste0(calculate.discover(df$disc_prob), "%")
  
  ggplot(df, aes(detected, median, color=as.factor(r_not), linetype=as.factor(disc_prob), fill=as.factor(r_not))) + 
    geom_line(size=1)+
    geom_ribbon(aes(ymax=max, ymin=min), alpha=0.1)+
    scale_y_log10(expand=c(0.01,0.01), limits=c(1,100), breaks = c(5,10,25,50,100))+
    scale_x_continuous(expand=c(0.01,0.01), limits=c(0,30))+
    scale_color_brewer(palette="Set1", direction = -1)+
    scale_fill_brewer(palette="Set1", direction=-1)+
    guides(linetype=FALSE)+
    labs(x = "Cumulative Reported Cases", 
         y = "Prevalence (log scale)", 
         color = expression("R"[0]), 
         fill = expression("R"[0]))
}

dotplot_data <- function(dir_path, r_nots, disc_probs, intro_rates){
  dirPaths = get_vec_of_files(dir_path, r_nots, disc_probs, intro_rates)
  ldply(dirPaths, function(x) {
    load(x)
    prevs <- all_max_prevalence(trials)
    parms <- get_parms(x)
    cbind(as.data.frame(parms), max_prev=prevs[sample(seq(0,10000), 1000)])
  })
}
plot_dots <- function(dir_path, r_nots, disc_probs, intros){
  dots <- dotplot_data(dir_path, r_nots, disc_probs, intros)
  ggplot(dots, aes(r_not, max_prev)) + facet_wrap(~intro_rate, nrow=1)+ 
    geom_point(position="jitter", shape=20,alpha=0.5) + geom_hline(yintercept=c(20,50), color=c("red","blue")) +
    theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=0.5, linetype=1))+
    scale_y_continuous(expand=c(0.01,0.01))+
    scale_x_continuous(expand=c(0.01,0.01))+
    panel_border(size=0.5, colour="black") +
    labs(x = expression("R"[0]), y= "Maximum Total Infectious")
}


max_prev_time_data <- function(trials){
  ## Gives data back for plotting the maximum prevalence against
  ## The time difference between the max prevalence and duration
  time <- all_sim_duration(trials) - all_time_max_prevalence(trials)
  max_prev <- all_max_prevalence(trials)
  
  data.frame(duration_diff=time, max_prev=max_prev)
}

plot_trial_data <- function(trials, type="prev", n=10, rando=TRUE, seqs=1){
  ## plots trial data
  if(rando){
    samples <- sample(x = 1:length(trials), size = n, replace = F)  
  } else{
    samples <- seqs
  }
  df <- do.call("rbind", trials[samples]) ## Combine into data frame
  df$index <- rep(seq_along(trials[samples]), sapply(trials[samples], nrow)) ## Add index for each simulation

  
  if(type=="prev"){
    df <- melt(df, id.vars=c("time", "index"), measure.vars = c("Total_Infections"))  
  } else if(type=="inc"){
    df <- melt(df, id.vars=c("time", "index"), measure.vars = c("New_Infection"))  
  }
  
  ggplot(df, aes(time, value, color=as.factor(index), group=as.factor(index))) + geom_line() + guides(color=FALSE)
}

plot_trial_data_local <- function(trials, type="prev", n=10, rando=TRUE, seqs=1){
  ## plots trial data
  if(rando){
    samples <- sample(x = 1:length(trials), size = n, replace = F)  
  } else{
    samples <- seqs
  }
  df <- do.call("rbind", trials[samples]) ## Combine into data frame
  df$index <- rep(seq_along(trials[samples]), sapply(trials[samples], nrow)) ## Add index for each simulation
  
  
  if(type=="prev"){
    df$local <- df$Total_Infections - df$Total_Intro_Infections
    df <- melt(df, id.vars=c("time", "index"), measure.vars = c("local"))  
  } else if(type=="inc"){
    df <- melt(df, id.vars=c("time", "index"), measure.vars = c("New_Infection"))  
  }
  
  ggplot(df, aes(time, value, color=as.factor(index), group=as.factor(index))) + geom_line() + guides(color=FALSE)
}



get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

panel_border <- function (colour = "gray80", size = 0.5, linetype = 1, remove = FALSE) 
{
  if (remove) {
    return(theme(panel.border = element_blank()))
  }
  theme(panel.border = element_rect(colour = colour, fill = NA, 
                                    linetype = linetype, size = size))
}



