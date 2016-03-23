#####################################
##
## Plotting functions for zika data
#####################################

plot_final_sizes <- function(x){
  final.sizes <- all_last_cuminfect_values(x)
  qplot(final.sizes, ..density.., geom="histogram", binwidth=5) +
    scale_y_continuous(expand=c(0,0)) + 
    scale_x_continuous(expand=c(0,0))+
    labs(x = "Final Epidemic Sizes")
}
