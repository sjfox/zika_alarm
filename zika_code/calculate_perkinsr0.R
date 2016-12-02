## Want a function that 
# 1; Takes in a list of months
# 2; Calculates all rnotts for those months
# 3 Calculate distribution of rnotts for a given temperature

select_month <- function(temps, desired.months) {
  # Down select the temperatures 
  temps %>% 
    gather(month, value, -Area.Name) %>%
    filter(month %in% desired.months) %>%
    spread(key = month, value = value)
}

align_data <- function(temps.selected, county_parms) {
  temps.selected$Area.Name = as.character(county_parms$Area.Name)
  full.data <- left_join(county_parms, temps.selected, by = "Area.Name")
}

test.data <- align_data(temp.selected, county_parms)

rnots_ests_bymonth <- function(full.data, ...) {
  
  
  # separate it into desired months 
  
  temps.selected <- dplyr::select(test.data, 
                                  which(colnames(test.data) %in% 
                                          eval(parse(text = "desired.months"))))
  
  rnot_ests_temp <- apply(X = temps.selected, MARGIN = 2, function(x) {
    browser()
    rnot_ests <- rnot_calc_counties(test.data$mosquito.abundance, 
                                    test.data$GDP, 
                                    temperature = x, 
                                    a=a, b=b, c.r=c.r, 
                                    mort.fun.list=mort.fun, 
                                    eip.fun.list=eip.fun, 
                                    scam.est.list=scam.est.list)
    
  })
  colnames(rnot_ests) <- colnames(temps.selected)
  return(rnot_ests)
}
 
  



