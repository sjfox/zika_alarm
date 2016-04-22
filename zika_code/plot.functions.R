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

plot_prob_below <- function(df){
  df$disc_prob <- calculate.discover(df$disc_prob)
  ggplot(df, aes(detected, prob_below, linetype=as.factor(disc_prob), color = as.factor(r_not))) + 
    geom_line(size=1) + facet_wrap(~intro_rate)+
    scale_y_continuous(expand=c(0.01,0.01)) +
    theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=1, linetype=1))+
    scale_color_brewer(palette="Set1") +
    panel_border(size=.5, colour="black") +
    background_grid(major = "y", minor = "none")+
    labs(x = "Cumulative Number of Detected Cases", 
         y = "Probability Prevalence is Below 25", 
         color = expression("R"[0]),
         linetype= "Detection \nProbability")
  
}


plot_prob_below <- function(df){
  df$disc_prob <- calculate.discover(df$disc_prob)
  ggplot(df, aes(detected, prob_epidemic, color = as.factor(r_not))) +
    geom_line(size=1, aes(linetype=as.factor(disc_prob))) + facet_wrap(~intro_rate)+
    scale_y_continuous(expand=c(0.01,0.01)) +
    scale_color_brewer(palette="Set1")+
    theme_cowplot() %+replace% theme(strip.background=element_rect(fill=NULL, color="black", size=1, linetype=1))+
    panel_border(size=.5, colour="black") +
    background_grid(major = "y", minor = "none")+
    labs(x = "Cumulative Number of Detected Cases", 
         y = "Probability of an Epidemic", 
         color = expression("R"[0]), 
         linetype= "Detection \nProbability")
}



panel_border <- function (colour = "gray80", size = 0.5, linetype = 1, remove = FALSE) 
{
  if (remove) {
    return(theme(panel.border = element_blank()))
  }
  theme(panel.border = element_rect(colour = colour, fill = NA, 
                                    linetype = linetype, size = size))
}