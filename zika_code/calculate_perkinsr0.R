## Want a function that 
# 1 Takes in a list of months
# 2 Calculates all rnotts for those months
# 3 Calculate distribution of rnotts for a given temperature and saves to csv 

#############################################################
# Step 1 load Perkins Nature Microbiology data and Texas data 
#############################################################
load("../data/perkins_sims/functions_R0_AR_random_draws.RData")
county_parms <- read.csv("../csvs/county_r0_parameters.csv")
tx_temps <- read.csv("../csvs/tx_county_temps.csv")


#####################################################################
# Functions for calculating rnot according to a specified temperature 
#####################################################################

select_month <- function(temps, desired.months) {
  # Down select the temperatures 
  temps %>% dplyr::select(county, which(colnames(temps) %in% 
                    eval(parse(text = "desired.months"))))
}

align_data <- function(temps.selected, county_parms) {
  full.data <- left_join(county_parms, temps.selected, by = "county")
}


rnot_calc_cty_dist <- function(mosq_abundance, gdp, temperature,
                               a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function returns a matrix of rnots based on abundances/gdp/temperatures
  # Dimensions are 254 rows x 1000 columns
  args <- list(mosq_abundance = as.list(mosq_abundance),
               gdp = as.list(gdp),
               temperature = as.list(temperature))
  
  cty_rnot_dists.l <- purrr::pmap(args, rnot_calc_dist, a=a, b=b, c.r=c.r, 
                                mort.fun.list=mort.fun.list, eip.fun.list=eip.fun.list,scam.est.list=scam.est.list)
  
  cty_rnot_dists <- as.data.frame(matrix(unlist(cty_rnot_dists.l), ncol = 1000, byrow = TRUE))
  colnames(cty_rnot_dists) <- seq(1:1000)
  return(cty_rnot_dists)
}
 
# save for each month 
calculate_rnots_bymonth <- function(temps, desired.months, county_parms) {
  
  # select months 
  temps.selected <- select_month(temps, desired.months)
  test.data <- align_data(temps.selected, county_parms)
 
  # separate it into desired months now that its aligned 
  temps.selected <- dplyr::select(test.data, 
                                  which(colnames(test.data) %in% 
                                          eval(parse(text = "desired.months"))))

  # Will return a list of 254x1000s r0s for each temperature selected 
  rnot_ests_temp <- apply(X = temps.selected, MARGIN = 2, function(x) {
    rnot_ests <- rnot_calc_counties(test.data$mosquito.abundance, 
                                    test.data$gdp, 
                                    temperature = x, 
                                    a=a, b=b, c.r=c.r, 
                                    mort.fun.list=mort.fun, 
                                    eip.fun.list=eip.fun, 
                                    scam.est.list=scam.est.list)
    # Cleaning up the distribution 
    names <- county_parms$county
    rnot_ests$county <-  names
    return(rnot_ests)
  })
  # Saving each entry as a separate month 
  lapply(seq_along(rnot_ests_temp), function(i){
    file.name = paste0(csv.path, desired.months[i], "_r0.csv")
    write.csv(rnot_ests_temp[[i]], file.name, row.names = FALSE)
  })
}

############################################################ 
desired.months = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
csv.path = "../csvs/rnot_temp/"
calculate_rnots_bymonth(temps = temps, desired.months = desired.months, county_parms = county_parms)


#### Ploting relationship for SCAME functions 
econ <- seq(6, 15, by=0.1)
plot(econ, predict(scam.est.list[[1]], newdata=data.frame(econ=econ)), type="l", ylim=c(-2, 5))
for(ii in 2:1000){
  lines(econ, predict(scam.est.list[[ii]], newdata=data.frame(econ=econ)), type="l")
}
