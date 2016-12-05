###############################################################
# Computing rnott by row
# One Function takes a vector of eip.fun, mort.fun, and scam.fun and calculates 1000 values of rnott for a single county
# Second function goes through each county


rnot_calc_dist <- function(mosq_abundance, gdp, temperature,
                           a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function that returns the full distribution of R0s for single set (county) of abundances gdp and temperature
  require(scam)
  
  index = seq(1:length(scam.est.list))
  
  # Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
  # calculates distribution of mortality rates
  g <- 1 / sapply(index, FUN = function(x) {
    mort.fun[[x]](temperature)
  })
  
  # calculates distribution of eip lengths
  e <- sapply(index, FUN = function(x) {
    eip.fun[[x]](temperature)
  })
  
  # calculates distribution of gdp scaling multipliers 
  gdp_scaling <- unname(sapply(index, FUN = function(x) {
    predict(scam.est.list[[x]], newdata=data.frame(econ=log(gdp)))
  })) 
  
  # calculates and returns distributions of rnott 
  as.numeric(exp(gdp_scaling) * mosq_abundance * a ^ 2 * b * c.r * exp(-g * e) / g)
}


rnot_calc_cty_dist <- function(mosq_abundance, gdp, temperature,
                               a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function returns a matrix of rnots based on abundances/gdp/temperatures
  # Dimensions are 254 rows x 1000 columns
  args <- list(mosq_abundance = as.list(mosq_abundance),
               gdp = as.list(gdp),
               temperature = as.list(temperature))
  
  cty_rnot_dists.l <- purrr::pmap(args, rnot_calc_dist, a=a, b=b, c.r=c.r, 
                                  mort.fun.list=mort.fun.list, 
                                  eip.fun.list=eip.fun.list,
                                  scam.est.list=scam.est.list)
  
  
  cty_rnot_dists <- as.data.frame(matrix(unlist(cty_rnot_dists.l), ncol = length(scam.est.list), byrow = TRUE))
  colnames(cty_rnot_dists) <- seq(1:ncol(cty_rnot_dists))
  return(cty_rnot_dists)
}


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
    rnot_ests <- rnot_calc_cty_dist(test.data$mosquito.abundance, 
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


# Original Perksin rnott functions written by SJF 

rnot_calc <- function(mosq_abundance, gdp, temperature, 
                      a, b, c.r, mort.fun, eip.fun, scam.est){
  # Function that returns a single rnot estimate from single values
  require(scam)
  # Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
  g <- 1 / mort.fun(temperature)
  e <- eip.fun(temperature)
  gdp_scaling <- predict(scam.est, newdata=data.frame(econ=log(gdp)))
  as.numeric(exp(gdp_scaling) * mosq_abundance * a ^ 2 * b * c.r * exp(-g * e) / g)
}

# New function with same name
#rnot_calc_dist <- function(mosq_abundance, gdp, temperature,
# a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
# Function that returns the full distribution of R0s for single set of abundances gdp and temperature
#  require(scam)
# Mortality function actually gives lifespan, so take 1 / lifespan to get mortality
#  args <- list(mort.fun = mort.fun.list,
#               eip.fun = eip.fun.list,
#               scam.est = scam.est.list)
#  unlist(purrr::pmap(args, rnot_calc, a=a, b=b, c.r=c.r, 
#                     mosq_abundance=mosq_abundance, gdp=gdp, temperature=temperature))
#}


rnot_calc_counties_q <- function(mosq_abundance, gdp, temperature,
                                 a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function returns median, hi, low r0 for vector of abundances/gdp/temperatures
  args <- list(mosq_abundance = as.list(mosq_abundance),
               gdp = as.list(gdp),
               temperature = as.list(temperature))
  
  cty_rnot_dists <- purrr::pmap(args, rnot_calc_dist, a=a, b=b, c.r=c.r, 
                                mort.fun.list=mort.fun.list, eip.fun.list=eip.fun.list,scam.est.list=scam.est.list)
  quantile_df <- function(x, probs, na.rm =F, names = F){
    z <- quantile(x, probs, na.rm, names)
    df <- data.frame(low_rnot=numeric(1), median_rnot = numeric(1), high_rnot=numeric(1))
    df[1,] <- z
    df
  }
  cty_rnot_dists %>% purrr::map( ~ quantile_df(.x, probs = c(0.025, 0.5, 0.975))) %>%
    bind_rows() %>%
    mutate(low_rnot = ifelse(low_rnot<0, 0, low_rnot))
}

rnot_calc_counties <- function(mosq_abundance, gdp, temperature,
                               a, b, c.r, mort.fun.list, eip.fun.list, scam.est.list){
  # Function returns for vector of abundances/gdp/temperatures
  args <- list(mosq_abundance = as.list(mosq_abundance),
               gdp = as.list(gdp),
               temperature = as.list(temperature))
  
  cty_rnot_dists <- purrr::pmap(args, rnot_calc_dist, a=a, b=b, c.r=c.r, 
                                mort.fun.list=mort.fun.list, eip.fun.list=eip.fun.list,scam.est.list=scam.est.list)
}




